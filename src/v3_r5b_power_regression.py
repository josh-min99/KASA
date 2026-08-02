"""
R5b: 검정력 재계산 — 계획서의 실제 판정 기준(회귀선 비교)으로.

왜 다시 계산하는가
  R5(`v3_r5_power_wetlab.py`)는 각 구간의 **평균 위상**을 비교했다. 그런데 계획서 §2 의
  판정 기준은 "baseline 구간 회귀선과 처치 구간 회귀선의 **기울기 및 절편**이 유의하게
  달라지는지" 다. 둘은 같지 않다.

    - 평균 비교는 구간 안에서 위상이 표류하지 않는다고 가정한다.
      그런데 T != 24h 군에서는 위상이 하루 (T-24)시간씩 체계적으로 표류한다.
      표류하는 계열의 산술평균은 '위상' 이 아니다.
    - 회귀는 기울기를 함께 추정하므로 절편의 표준오차가 커진다.
      실제로 처치 시작일을 기준점으로 두면 합성 표준오차가 0.762h -> 1.482h 로
      1.95배 나빠진다. 즉 R5 의 필요 마리수는 낙관적이었다.

여기서 검정 두 가지를 계획서 문구 그대로 나눈다
  검정 A (절편) : baseline 회귀선을 처치 시작일로 연장한 값과 처치 회귀선의 절편 차이.
                  효과 크기는 새로 형성된 위상 이동량 dphi.
  검정 B (기울기): 두 구간 회귀선의 기울기 차이. 효과 크기는 (T - tau) 로 설계가 정한다.
                  주기 자체가 바뀌었는지를 보므로 T != 24h 군에서 강력하다.

개체별로 두 구간에 각각 최소자승 직선을 적합하고, 개체 단위 통계량에 1표본 t 검정을 한다.
개체 고유 위상은 구간 간 차분에서 상쇄되므로 자기대조 설계의 이점이 유지된다.

산출: data/rhythm/power_regression.csv, results/v3/r5b_power_regression.txt
"""
import os
import sys
import warnings

import numpy as np
import pandas as pd
from scipy import stats

warnings.filterwarnings("ignore")
if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA = os.path.join(ROOT, "data", "rhythm")
RES = os.path.join(ROOT, "results", "v3")
os.makedirs(RES, exist_ok=True)

RNG = np.random.default_rng(20260805)
N_SIM = 4000
ALPHA = 0.05
TAU = 23.7


def seg_stats(x, d0):
    """구간 x 에서 OLS 로 (1) 기준일 d0 에서의 적합값, (2) 기울기를 얻을 때의
    분산 계수를 돌려준다. 둘 다 관측값의 선형결합이므로 계수만 있으면 된다.

        Var[yhat(d0)] = sigma^2 * (1/n + (d0-xbar)^2/Sxx)
        Var[slope]    = sigma^2 / Sxx
    """
    n = len(x)
    xbar = x.mean()
    sxx = ((x - xbar) ** 2).sum()
    return 1.0 / n + (d0 - xbar) ** 2 / sxx, 1.0 / sxx


def t_power(effect, sd, n, alpha=ALPHA):
    """1표본 t 검정 검정력 (비중심 t 분포, 양측). 개체별 통계량이
    N(effect, sd^2) 를 따르므로 닫힌 형태로 정확히 계산된다."""
    if sd <= 0:
        return 1.0
    ncp = effect / sd * np.sqrt(n)
    tc = stats.t.ppf(1 - alpha / 2, n - 1)
    return float(1 - stats.nct.cdf(tc, n - 1, ncp) + stats.nct.cdf(-tc, n - 1, ncp))


def powers(n_animals, nb, nt, sd_b, sd_t, dphi, T):
    """검정 A(절편)·B(기울기) 검정력."""
    xb = np.arange(nb, dtype=float)
    xt = np.arange(nb, nb + nt, dtype=float)
    d0 = float(nb)
    cA_b, cB_b = seg_stats(xb, d0)
    cA_t, cB_t = seg_stats(xt, d0)
    sdA = np.sqrt(sd_b ** 2 * cA_b + sd_t ** 2 * cA_t)     # 절편 차의 개체 내 SD
    sdB = np.sqrt(sd_b ** 2 * cB_b + sd_t ** 2 * cB_t)     # 기울기 차의 개체 내 SD
    effB = T - TAU                                          # 기울기 차의 참값
    return t_power(dphi, sdA, n_animals), t_power(abs(effB), sdB, n_animals)


def _selfcheck():
    """닫힌 형태가 맞는지 몬테카를로로 한 칸만 검산한다."""
    nb, nt, sd_b, sd_t, dphi, T, n = 3, 14, 0.31, 2.77, 2.0, 20, 6
    xb = np.arange(nb, dtype=float); xt = np.arange(nb, nb + nt, dtype=float)
    d0 = float(nb); hitA = hitB = 0; N = 4000
    for _ in range(N):
        DA = np.empty(n); DB = np.empty(n)
        for i in range(n):
            # 귀무가설은 '처치 구간이 baseline 자유진행선을 그대로 잇는다' 이다.
            # 따라서 처치 구간의 기준일 값은 baseline 선의 기준일 값 + dphi 여야 한다.
            # 초판 검산은 이 오프셋 (TAU-24)*d0 를 빼먹어 효과를 2.0 이 아니라 2.9 로
            # 생성했고, 그래서 몬테카를로가 닫힌 형태보다 높게 나왔다.
            base_at_d0 = (TAU - 24) * d0
            yb = (TAU - 24) * xb + RNG.normal(0, sd_b, nb)
            yt = base_at_d0 + dphi + (T - 24) * (xt - d0) + RNG.normal(0, sd_t, nt)

            def ols(x, y):
                xm = x.mean(); sxx = ((x - xm) ** 2).sum()
                b = ((x - xm) * (y - y.mean())).sum() / sxx
                return y.mean() + b * (d0 - xm), b
            fb, bb = ols(xb, yb); ft, bt = ols(xt, yt)
            DA[i] = ft - fb; DB[i] = bt - bb
        for D, which in ((DA, "A"), (DB, "B")):
            t = D.mean() / (D.std(ddof=1) / np.sqrt(n))
            if 2 * (1 - stats.t.cdf(abs(t), n - 1)) < ALPHA:
                if which == "A":
                    hitA += 1
                else:
                    hitB += 1
    fa, fb_ = powers(n, nb, nt, sd_b, sd_t, dphi, T)
    return (hitA / N, fa), (hitB / N, fb_)


def main():
    prec = pd.read_csv(os.path.join(DATA, "phase_precision.csv"))
    s = (prec[prec.condition.isin(["baseline", "treatment"])]
         .groupby(["variable", "dataset", "condition"]).phase_sd_h.median().unstack())

    log = []
    P = log.append
    P("=" * 96)
    P("R5b  계획서 판정 기준(회귀선 비교)에 따른 검정력")
    P("=" * 96)
    P("")
    P("R5(구간 평균 비교)와의 차이 — 심부체온, baseline 3일 / 처치 14일 기준")
    sd_b, sd_t = 0.31, 2.77
    nb, nt = 3, 14
    xb = np.arange(nb, dtype=float); xt = np.arange(nb, nb + nt, dtype=float)
    Sxx_b = ((xb - xb.mean()) ** 2).sum(); Sxx_t = ((xt - xt.mean()) ** 2).sum()
    se_mean = np.hypot(sd_b / np.sqrt(nb), sd_t / np.sqrt(nt))
    se_reg = np.hypot(sd_b * np.sqrt(1 / nb + (nb - xb.mean()) ** 2 / Sxx_b),
                      sd_t * np.sqrt(1 / nt + (nb - xt.mean()) ** 2 / Sxx_t))
    P(f"  구간 평균 비교  합성 표준오차 {se_mean:.3f} h")
    P(f"  회귀선 절편 비교 합성 표준오차 {se_reg:.3f} h   ({se_reg/se_mean:.2f} 배)")
    P(f"  baseline 기울기 표준오차 {sd_b/np.sqrt(Sxx_b):.3f} h/일 — 3일로는 기울기를 못 잰다")
    P("")

    (mcA, fA), (mcB, fB) = _selfcheck()
    P("자체 검산 (T=20h, 심부체온, 6마리, dphi=2h) — 닫힌 형태 대 몬테카를로")
    P(f"  검정 A 절편  몬테카를로 {mcA:.3f}  닫힌형태 {fA:.3f}")
    P(f"  검정 B 기울기 몬테카를로 {mcB:.3f}  닫힌형태 {fB:.3f}")
    P("")

    rows = []
    for var, lab in [("tb_core", "심부체온"), ("tb_sub", "피하온도"), ("activity", "활동량")]:
        sub = s.loc[var] if var in s.index.get_level_values(0) else None
        if sub is None:
            continue
        for ds in sub.index:
            b, t = float(sub.loc[ds, "baseline"]), float(sub.loc[ds, "treatment"])
            # 처치 전 구간 길이도 격자에 넣는다. 3일이 기본 설계값이고,
            # 5·7일은 '동물 수 대신 관측 기간을 늘리는' 대안을 평가하기 위한 것이다.
            for nb_i in (3, 5, 7):
                for T in (24, 20, 28):
                    for dphi in (1.0, 2.0, 3.0):
                        for n in (3, 4, 5, 6, 7, 8, 10, 12, 16, 20, 24):
                            pa, pb = powers(n, nb_i, nt, b, t, dphi, T)
                            rows.append(dict(variable=var, label=lab, dataset=ds,
                                             n_days_base=nb_i, T_hours=T,
                                             dphi=dphi, n_animals=n,
                                             power_intercept=pa, power_slope=pb,
                                             sd_base=b, sd_treat=t))
    df = pd.DataFrame(rows)
    df.to_csv(os.path.join(DATA, "power_regression.csv"), index=False, encoding="utf-8-sig")

    P("-" * 96)
    P("[검정 A — 절편 차이]  '처치로 새 위상이 형성되었는가'")
    P("검정력 0.8 을 얻는 최소 마리수 (처치 14일). T 에 따라 달라지지 않는다(절편 검정이므로)")
    P("")
    P(f"  {'지표':10s} {'코호트':14s} {'dphi=1h':>9s} {'dphi=2h':>9s} {'dphi=3h':>9s}")
    for var, lab in [("tb_core", "심부체온"), ("tb_sub", "피하온도"), ("activity", "활동량")]:
        for ds in sorted(df[df.variable == var].dataset.unique()):
            cells = []
            for dphi in (1.0, 2.0, 3.0):
                q = df[(df.variable == var) & (df.dataset == ds) & (df.dphi == dphi)
                       & (df.T_hours == 24) & (df.n_days_base == 3)
                       & (df.power_intercept >= 0.8)]
                cells.append(str(int(q.n_animals.min())) if len(q) else ">24")
            P(f"  {lab:10s} {ds:14s} " + "".join(f"{c:>9s}" for c in cells))

    P("")
    P("-" * 96)
    P("[검정 B — 기울기 차이]  '주기 자체가 인가 주기로 바뀌었는가'")
    P("효과 크기는 설계가 정한다: (T - tau). T=24 -> +0.3, T=20 -> -3.7, T=28 -> +4.3 h/일")
    P("")
    P(f"  {'지표':10s} {'코호트':14s} {'T=24h':>8s} {'T=20h':>8s} {'T=28h':>8s}")
    for var, lab in [("tb_core", "심부체온"), ("tb_sub", "피하온도"), ("activity", "활동량")]:
        for ds in sorted(df[df.variable == var].dataset.unique()):
            cells = []
            for T in (24, 20, 28):
                q = df[(df.variable == var) & (df.dataset == ds) & (df.T_hours == T)
                       & (df.dphi == 2.0) & (df.n_days_base == 3)
                       & (df.power_slope >= 0.8)]
                cells.append(str(int(q.n_animals.min())) if len(q) else ">24")
            P(f"  {lab:10s} {ds:14s} " + "".join(f"{c:>8s}" for c in cells))

    P("")
    P("=" * 96)
    P("판독")
    P("=" * 96)
    P("")
    P("1) R5 의 필요 마리수는 낙관적이었다. 계획서 기준(회귀선)으로 다시 계산하면")
    P("   절편 검정의 합성 표준오차가 1.95배 커지고, 필요 마리수도 그만큼 늘어난다.")
    P("")
    P("2) 그러나 T != 24h 군에서는 기울기 검정이 압도적으로 강력하다.")
    P("   주기가 20h 나 28h 로 바뀌면 하루 4시간씩 표류하므로, 적은 마리수로도 검출된다.")
    P("   즉 10시간 · 14시간 군의 주 검정은 절편이 아니라 기울기여야 한다.")
    P("")
    P("3) 12시간 군(T=24h)은 기울기 차이가 0.3 h/일 로 작아 절편 검정에 의존한다.")
    P("   세 군 중 이 군이 가장 많은 마리수를 요구한다.")
    P("")
    P("4) 처치 전 구간을 늘리는 것의 효과는 검정마다 다르다. 아래는 심부체온 기준.")
    P("")
    P(f"   {'처치전':>6s} {'절편검정 dphi=2h':>16s} {'기울기검정 T=24h':>16s} {'기울기검정 T=20h':>16s}")
    for nb_i in (3, 5, 7):
        cells = []
        for col, f in (("power_intercept", dict(dphi=2.0, T_hours=24)),
                       ("power_slope", dict(dphi=2.0, T_hours=24)),
                       ("power_slope", dict(dphi=2.0, T_hours=20))):
            q = df[(df.variable == "tb_core") & (df.n_days_base == nb_i)]
            for k, v in f.items():
                q = q[q[k] == v]
            q = q[q[col] >= 0.8]
            cells.append(str(int(q.n_animals.min())) if len(q) else ">24")
        P(f"   {nb_i:>4d}일 {cells[0]:>16s} {cells[1]:>16s} {cells[2]:>16s}")
    P("")
    P("   절편 검정은 처치 구간 오차(2.77h)가 표준오차의 대부분을 차지하므로")
    P("   처치 전 구간을 늘려도 개선되지 않는다. 반면 기울기 검정은 처치 전 구간의")
    P("   기울기 추정이 직접 들어가므로 3일 -> 5일에서 뚜렷하게 개선되고 그 뒤로는 포화된다.")

    txt = "\n".join(log)
    print(txt)
    with open(os.path.join(RES, "r5b_power_regression.txt"), "w", encoding="utf-8") as fh:
        fh.write(txt + "\n")
    print(f"\n-> data/rhythm/power_regression.csv, results/v3/r5b_power_regression.txt")


if __name__ == "__main__":
    main()
