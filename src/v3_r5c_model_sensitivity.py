"""
R5c: 진동자 모델의 파라미터 민감도 — '짧은 쪽이 넓다'는 비대칭이 견고한가.

왜 필요한가
  계획서 예상결과 (4)에 "14시간 군(T=28h)이 10시간 군(T=20h)보다 강한 자극을
  요구한다" 고 썼다. 이 주장은 Arnold tongue 이 짧은 주기 쪽으로 넓게 기울어져
  있다는 모델 결과에 의존한다.

  그런데 모델 파라미터는 실측값이 아니라 '20-28h 주기에서 상대진폭이 최대'
  라는 자기정합 기준으로 고른 값이다. 그 선택 과정은 저장소에 재현 가능한
  형태로 남아 있지도 않다. 따라서 파라미터를 흔들었을 때 비대칭의 **방향**이
  유지되는지 확인하지 않으면, 계획서의 그 주장은 근거가 없다.

  방향이 뒤집히면 계획서에서 해당 주장을 빼야 한다.

방법
  기준 파라미터에서 주요 계수를 각각 ±20% 흔들고, 자극 세기 두 값에서
  동조 가능한 구동 주기 범위를 다시 구해 고유 주기 기준 좌우 폭을 비교한다.
  전체 651 격자 대신 필요한 세기만 계산해 비용을 줄인다.

산출: data/rhythm/model_sensitivity.csv, results/v3/r5c_sensitivity.txt
"""
import os
import sys
import copy

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "scripts"))
import oscillator_model as om  # noqa: E402

if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA = os.path.join(ROOT, "data", "rhythm")
RES = os.path.join(ROOT, "results", "v3")
os.makedirs(RES, exist_ok=True)

STRENGTHS = [0.05, 0.10]
N_PERIOD = 25          # 기준 실행(31)보다 성기지만 좌우 폭 비교에는 충분


class use_par:
    """파라미터를 실제로 주입하는 컨텍스트 매니저.

    주의 — `om.PAR = par` 로는 바뀌지 않는다.
    oscillator_model.rhs 의 시그니처가 `def rhs(y, u, kdeg, p=PAR)` 이라
    기본 인자 p 가 **함수 정의 시점의 PAR 객체에 묶여** 있기 때문이다.
    모듈 전역 PAR 을 재할당해도 rhs 는 옛 객체를 계속 본다.
    (초판 민감도 실행이 이 함정에 걸려 6개 파라미터 변형이 전부 무시됐고,
     고유주기가 소수점까지 동일하게 나와 발각됐다.)
    integrate() 는 호출 시점에 모듈 전역 rhs 를 찾으므로, rhs 자체를 갈아끼운다.
    """

    def __init__(self, par):
        self.par = par

    def __enter__(self):
        self.orig = om.rhs
        orig = self.orig
        par = self.par
        om.rhs = lambda y, u, kdeg, p=par: orig(y, u, kdeg, p)
        return self

    def __exit__(self, *a):
        om.rhs = self.orig
        return False


def scan(par, kdeg, S, tau, n_t=N_PERIOD):
    """세기 S 하나에서 동조 가능한 구동 주기 범위를 구한다."""
    with use_par(par):
        periods = np.linspace(tau * 0.80, tau * 1.20, n_t)
        N = len(periods)
        T = periods.copy()
        kd = np.full(N, kdeg)
        y = om.to_limit_cycle(kd, burn=600.0)

        def u_fn(t):
            return np.where((t % T) < 1.0, S, 0.0)

        T0 = 120 * tau
        y = om.integrate(y, kd, int(T0 / om.DT), u_fn=u_fn)
        _, tr = om.integrate(y, kd, int(50 * tau / om.DT), u_fn=u_fn, t0=T0, record="full")

        ent = np.zeros(N, bool)
        for j in range(N):
            b = tr[:, om.PHASE_VAR, j]
            alive = (b.max() - b.min()) / max(b.mean(), 1e-9) > 0.10
            pk = om.peaks(b, t0=T0)
            if len(pk) < 10 or not alive:
                continue
            rot = (len(pk) - 1) * T[j] / (pk[-1] - pk[0])
            phi = (pk % T[j]) / T[j]
            u = np.unwrap(phi * 2 * np.pi) / (2 * np.pi)
            cyc = (pk - pk[0]) / T[j]
            fit = np.polyfit(cyc, u, 1)
            resid = u - np.polyval(fit, cyc)
            ent[j] = (abs(rot - 1.0) < 0.02 and abs(fit[0]) < 0.004 and resid.std() < 0.02)
        if not ent.any():
            return np.nan, np.nan
        return float(periods[ent].min()), float(periods[ent].max())


def main():
    base = dict(om.PAR)
    kdeg = om.KDEG0

    variants = [("기준", base, kdeg)]
    for key in ["v1", "v2", "k3", "v4", "k5", "n"]:
        for f, lab in ((0.8, "-20%"), (1.2, "+20%")):
            p = dict(base)
            p[key] = base[key] * f
            variants.append((f"{key} {lab}", p, kdeg))
    for f, lab in ((0.8, "-20%"), (1.2, "+20%")):
        variants.append((f"kdeg {lab}", dict(base), kdeg * f))

    rows = []
    print(f"{'파라미터':12s} {'고유주기':>8s} {'세기':>6s} {'동조 최소':>9s} {'최대':>8s} "
          f"{'짧은쪽':>7s} {'긴쪽':>7s} {'짧은/긴':>8s}")
    for lab, par, kd in variants:
        with use_par(par):
            per0, _, rel0, _, _ = om.characterize(np.array([kd]))
        tau = float(per0[0])
        if not np.isfinite(tau) or float(rel0[0]) < 0.05:
            print(f"{lab:12s} {'진동 없음':>8s}")
            rows.append(dict(variant=lab, tau=np.nan, oscillating=False))
            continue
        for S in STRENGTHS:
            lo, hi = scan(par, kd, S, tau)
            if not np.isfinite(lo):
                print(f"{lab:12s} {tau:8.2f} {S:6.2f} {'동조 없음':>9s}")
                rows.append(dict(variant=lab, tau=tau, strength=S, oscillating=True,
                                 lo=np.nan, hi=np.nan, short=np.nan, long=np.nan))
                continue
            short, lng = tau - lo, hi - tau
            ratio = short / lng if lng > 1e-6 else np.inf
            print(f"{lab:12s} {tau:8.2f} {S:6.2f} {lo:9.2f} {hi:8.2f} "
                  f"{short:7.2f} {lng:7.2f} {ratio:8.2f}")
            rows.append(dict(variant=lab, tau=tau, strength=S, oscillating=True,
                             lo=lo, hi=hi, short=short, long=lng, short_over_long=ratio))

    df = pd.DataFrame(rows)
    df.to_csv(os.path.join(DATA, "model_sensitivity.csv"), index=False, encoding="utf-8-sig")

    ok = df[df.short_over_long.notna()]
    n_short_wider = int((ok.short_over_long > 1).sum())
    lines = []
    L = lines.append
    L("=" * 84)
    L("R5c  진동자 모델 파라미터 민감도 — 비대칭 방향이 견고한가")
    L("=" * 84)
    L("")
    L("검사 대상 주장: '짧은 주기 쪽으로 동조가 더 쉽다'")
    L("  -> 계획서의 'T=28h 가 T=20h 보다 강한 자극을 요구한다' 가 여기에 의존한다")
    L("")
    L(f"파라미터 변형 {len(variants)}종 x 자극 세기 {len(STRENGTHS)}종")
    L(f"동조 범위가 산출된 조합 {len(ok)}건 중 짧은 쪽이 더 넓은 경우: "
      f"{n_short_wider}건 ({n_short_wider/max(len(ok),1):.0%})")
    L(f"짧은쪽/긴쪽 비 중앙값 {ok.short_over_long.median():.2f} "
      f"(범위 {ok.short_over_long.min():.2f} ~ {ok.short_over_long.max():.2f})")
    L("")
    if n_short_wider == len(ok):
        L("판정: 검사한 모든 변형에서 방향이 유지된다. 계획서의 주장을 유지한다.")
    elif n_short_wider >= 0.8 * len(ok):
        L("판정: 대부분 유지되나 예외가 있다. 계획서에 '검사한 파라미터 범위에서' 를 명시한다.")
    else:
        L("판정: 방향이 파라미터에 의존한다. 계획서에서 해당 주장을 철회해야 한다.")
    L("")
    L(df.to_string(index=False))

    txt = "\n".join(lines)
    print("\n" + txt)
    with open(os.path.join(RES, "r5c_sensitivity.txt"), "w", encoding="utf-8") as fh:
        fh.write(txt + "\n")
    print(f"\n-> data/rhythm/model_sensitivity.csv, results/v3/r5c_sensitivity.txt")


if __name__ == "__main__":
    main()
