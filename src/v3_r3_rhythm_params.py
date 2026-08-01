"""
R3 + R3.5: 공통 리듬 파라미터 추출과 잡음 바닥.

측정 방식이 제각각인 자료(임플란트 활동카운트 / 심부체온 / 영상 강도)를 같은 축에
놓아야 하므로, 두 종류의 지표를 함께 뽑는다.

  모수 지표   cosinor(24h) -> MESOR, 진폭, acrophase, 리듬검출 p
              리듬검출 p 는 영진폭 검정(zero-amplitude test)의 F 검정이다.
  비모수 지표 IS(일간 안정성) / IV(일내 변동성) / RA(상대 진폭)
              분포 가정도 파형 가정도 없어서 측정 방식이 달라도 비교 가능하다.

자기 대조 정규화
  개체별로 baseline 구간 값을 기준으로 삼는다. 종·장비·단위 차이를 흡수한다.

R3.5 잡음 바닥 (반드시 비율 논증보다 먼저)
  PROGRESS.md [2026-07-31 20:15] 에 잡음 바닥을 확인하지 않고 '64%' 라는 비율을
  주장했다가 철회한 기록이 있다. 그래서 여기서는 어떤 대비를 계산하기 전에,
  '효과가 없을 때 이 지표가 보이는 변화의 크기'를 먼저 구한다.
  방법: baseline 구간 안에서만 두 조각을 무작위로 갈라 같은 대비를 계산하는 것을 반복한다.
        효과가 0 인 대비이므로 이 분포가 곧 잡음 바닥이다.

게이트 G2
  Helissen 논문(Life 13:844) 이 보고한 cosinor 수치를 우리 파이프라인이 재현하는가.
  실패하면 이후 단계를 진행하지 않는다.

산출: data/rhythm/params.csv, data/rhythm/noise_floor.csv, results/v3/gate_log.txt
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
OUT = os.path.join(ROOT, "data", "rhythm")
RES = os.path.join(ROOT, "results", "v3")
os.makedirs(RES, exist_ok=True)

RNG = np.random.default_rng(20260802)
N_BOOT = 2000
OMEGA = 2 * np.pi / 24


# ------------------------------------------------------------------ cosinor
def cosinor(t_hours, y, period=24.0):
    """단일 성분 cosinor. 반환: MESOR, 진폭, acrophase(h), 영진폭검정 p, n."""
    t = np.asarray(t_hours, float)
    y = np.asarray(y, float)
    ok = np.isfinite(t) & np.isfinite(y)
    t, y = t[ok], y[ok]
    n = len(y)
    if n < 6 or np.allclose(y, y[0]):
        return dict(mesor=np.nan, amplitude=np.nan, acrophase=np.nan, p_rhythm=np.nan, n=n)

    w = 2 * np.pi / period
    X = np.column_stack([np.ones(n), np.cos(w * t), np.sin(w * t)])
    try:
        beta, *_ = np.linalg.lstsq(X, y, rcond=None)
    except np.linalg.LinAlgError:
        return dict(mesor=np.nan, amplitude=np.nan, acrophase=np.nan, p_rhythm=np.nan, n=n)

    fit = X @ beta
    rss = float(((y - fit) ** 2).sum())
    tss = float(((y - y.mean()) ** 2).sum())
    amp = float(np.hypot(beta[1], beta[2]))
    # acrophase: y = M + A*cos(w(t - phi)) 의 phi (시간 단위, 0-24)
    acro = float((np.arctan2(beta[2], beta[1]) / w) % period)

    # 영진폭 검정: H0 beta1=beta2=0. df1=2, df2=n-3
    df2 = n - 3
    if df2 <= 0 or rss <= 0:
        p = np.nan
    else:
        f = ((tss - rss) / 2) / (rss / df2)
        p = float(1 - stats.f.cdf(f, 2, df2))
    return dict(mesor=float(beta[0]), amplitude=amp, acrophase=acro, p_rhythm=p, n=n)


# ------------------------------------------------------- 비모수 (측정방식 불문)
def nonparam(t_hours, y, bin_h):
    """IS / IV / RA. bin_h 는 표본 간격(시간)."""
    t = np.asarray(t_hours, float)
    y = np.asarray(y, float)
    ok = np.isfinite(t) & np.isfinite(y)
    t, y = t[ok], y[ok]
    if len(y) < 8:
        return dict(IS=np.nan, IV=np.nan, RA=np.nan)

    per_day = int(round(24 / bin_h))
    # IS: 하루 중 같은 시각끼리 모은 평균의 분산 / 전체 분산
    slot = (np.round(t / bin_h).astype(int)) % per_day
    gm = y.mean()
    tot = ((y - gm) ** 2).sum()
    if tot <= 0:
        return dict(IS=np.nan, IV=np.nan, RA=np.nan)
    slot_means = np.array([y[slot == s].mean() for s in range(per_day) if (slot == s).any()])
    n_slot = np.array([(slot == s).sum() for s in range(per_day) if (slot == s).any()])
    IS = float((n_slot * (slot_means - gm) ** 2).sum() / tot)

    # IV: 인접 표본 차이의 분산 / 전체 분산
    order = np.argsort(t)
    ys = y[order]
    IV = float((len(ys) * ((np.diff(ys) ** 2).sum())) / ((len(ys) - 1) * tot))

    # RA: M10/L5 를 시간 폭으로 정의하고 표본 간격에 맞춰 창 길이를 정한다
    w10 = max(1, int(round(10 / bin_h)))
    w5 = max(1, int(round(5 / bin_h)))
    if len(ys) < w10 + 1:
        return dict(IS=IS, IV=IV, RA=np.nan)
    roll10 = pd.Series(ys).rolling(w10).mean().dropna().values
    roll5 = pd.Series(ys).rolling(w5).mean().dropna().values
    M10, L5 = float(roll10.max()), float(roll5.min())
    RA = float((M10 - L5) / (M10 + L5)) if (M10 + L5) != 0 else np.nan
    return dict(IS=IS, IV=IV, RA=RA)


def bin_width(t):
    d = np.diff(np.sort(np.unique(np.asarray(t, float))))
    d = d[d > 0]
    return float(np.median(d)) if len(d) else np.nan


# --------------------------------------------------------------------- main
def per_segment_params(g):
    bw = bin_width(g.t_hours.values)
    out = cosinor(g.t_hours.values, g.value.values)
    out.update(nonparam(g.t_hours.values, g.value.values, bw if np.isfinite(bw) else 2.0))
    out["bin_h"] = bw
    out["span_h"] = float(g.t_hours.max() - g.t_hours.min())
    return out


def main():
    long = pd.read_csv(os.path.join(OUT, "long.csv"))
    print("=" * 78)
    print("R3  리듬 파라미터 추출")
    print("=" * 78)

    # osd595 는 명기 전용이라 24h 위상/진폭을 계산하면 안 된다 (R1 판정: LEVEL_ONLY).
    rhythm_part = long[long.dataset != "osd595"]
    print(f"\n리듬 축 대상 {rhythm_part.dataset.nunique()} 데이터셋 "
          f"(osd595 는 LEVEL_ONLY 라 리듬 계산에서 제외)")

    rows = []
    for (ds, subj, cond, var), g in rhythm_part.groupby(
            ["dataset", "subject_id", "condition", "variable"]):
        if len(g) < 6:
            continue
        p = per_segment_params(g)
        p.update(dataset=ds, subject_id=subj, condition=cond, variable=var,
                 species=g.species.iloc[0], gravity=g.gravity.iloc[0])
        rows.append(p)
    par = pd.DataFrame(rows)
    par.to_csv(os.path.join(OUT, "params.csv"), index=False, encoding="utf-8-sig")
    print(f"파라미터 {len(par)} 행 (데이터셋 x 개체 x 구간 x 변수)")

    # ---------------------------------------------------------- 게이트 G2
    log = []
    L = log.append
    L("=" * 78)
    L("게이트 G2 — Helissen 2023 (Life 13:844) 보고값 재현")
    L("=" * 78)
    L("논문 보고 (피하 텔레메트리 코호트, n=5):")
    L("  활동  MESOR 0.26+-0.11 -> 0.04+-0.01 | 진폭 0.18+-0.11 -> 0.01+-0.01 |"
      " 리듬검출 p 0.005 -> 0.420")
    L("  피하온도 MESOR 33.81+-0.72 -> 32.69+-0.51 | 진폭 0.78+-0.19 -> 0.23+-0.07 |"
      " 리듬검출 p 0.004 -> 0.129")
    L("")

    checks = []
    for ds in ["helissen2020", "helissen2021"]:
        for var, label in [("activity", "활동"), ("tb_sub", "피하온도")]:
            sub = par[(par.dataset == ds) & (par.variable == var)]
            if not len(sub):
                continue
            L(f"[{ds} / {label}]  개체 {sub.subject_id.nunique()}")
            for cond in ["baseline", "treatment", "recovery"]:
                s = sub[sub.condition == cond]
                if not len(s):
                    continue
                L(f"   {cond:10s} MESOR {s.mesor.mean():8.3f}  진폭 {s.amplitude.mean():7.3f}"
                  f"  acrophase {s.acrophase.mean():5.2f}h  리듬검출 p {s.p_rhythm.mean():.3f}"
                  f"  (개체별 p 중앙값 {s.p_rhythm.median():.3f})")
            b = sub[sub.condition == "baseline"]
            t = sub[sub.condition == "treatment"]
            if len(b) and len(t):
                checks.append(dict(ds=ds, var=var,
                                   mesor_b=b.mesor.mean(), mesor_t=t.mesor.mean(),
                                   amp_b=b.amplitude.mean(), amp_t=t.amplitude.mean(),
                                   p_b=b.p_rhythm.median(), p_t=t.p_rhythm.median()))
            L("")

    # 판정: 방향과 대략적 크기가 맞는가 (소수점 일치를 요구하지 않는다.
    #       논문은 군 평균 cosinor, 우리는 개체별 cosinor 의 평균이라 완전 일치는 기대할 수 없다)
    g2 = []
    for c in checks:
        if c["var"] == "activity":
            g2.append(("활동 MESOR 감소", c["mesor_t"] < c["mesor_b"] * 0.4))
            g2.append(("활동 진폭 감소", c["amp_t"] < c["amp_b"] * 0.4))
            g2.append(("활동 리듬검출 상실(p 상승)", c["p_t"] > c["p_b"]))
        else:
            g2.append(("피하온도 MESOR 감소", c["mesor_t"] < c["mesor_b"]))
            g2.append(("피하온도 진폭 감소", c["amp_t"] < c["amp_b"] * 0.6))
    npass = sum(1 for _, ok in g2 if ok)
    L("-" * 78)
    L("G2 판정 항목")
    for name, ok in g2:
        L(f"   {'OK ' if ok else '!! '}{name}")
    L(f"\nG2: {npass}/{len(g2)} 통과")
    gate_ok = npass >= int(0.8 * len(g2))
    L("판정: 통과 — 파이프라인이 원논문 결과를 재현한다" if gate_ok
      else "판정: 실패 — 파이프라인 점검 필요. 이후 단계 중단")

    # ------------------------------------------------- R3.5 잡음 바닥
    L("")
    L("=" * 78)
    L("R3.5  잡음 바닥 — 효과가 0 일 때 이 지표가 보이는 변화")
    L("=" * 78)
    L("방법: baseline 구간만 사용해, 같은 개체의 baseline 을 무작위로 두 조각으로 갈라")
    L("      우리가 쓸 대비(진폭비·MESOR차·acrophase차)를 계산하는 것을 2,000회 반복한다.")
    L("      효과가 0 인 대비이므로 이 분포가 곧 잡음 바닥이다.")
    L("      이후 모든 주장은 '잡음 대비 몇 배' 로 적는다.")
    L("")

    nf_rows = []
    for (ds, var), g in rhythm_part.groupby(["dataset", "variable"]):
        base = g[g.condition == "baseline"]
        if base.subject_id.nunique() < 3:
            continue
        d_amp, d_mesor, d_acro = [], [], []
        for _ in range(N_BOOT):
            a_amp, b_amp, a_m, b_m, a_p, b_p = [], [], [], [], [], []
            for subj, sg in base.groupby("subject_id"):
                t = sg.t_hours.values
                if len(t) < 14:
                    continue
                # 시간축을 무작위 지점에서 둘로 가른다 (각 조각 >= 하루)
                lo, hi = t.min() + 24, t.max() - 24
                if hi <= lo:
                    continue
                cut = RNG.uniform(lo, hi)
                s1, s2 = sg[sg.t_hours < cut], sg[sg.t_hours >= cut]
                if len(s1) < 6 or len(s2) < 6:
                    continue
                c1 = cosinor(s1.t_hours.values, s1.value.values)
                c2 = cosinor(s2.t_hours.values, s2.value.values)
                if not np.isfinite(c1["amplitude"]) or not np.isfinite(c2["amplitude"]):
                    continue
                a_amp.append(c1["amplitude"]); b_amp.append(c2["amplitude"])
                a_m.append(c1["mesor"]); b_m.append(c2["mesor"])
                a_p.append(c1["acrophase"]); b_p.append(c2["acrophase"])
            if len(a_amp) < 2:
                continue
            A1, A2 = np.mean(a_amp), np.mean(b_amp)
            if A1 > 0:
                d_amp.append(abs(A2 / A1 - 1))
            d_mesor.append(abs(np.mean(b_m) - np.mean(a_m)))
            dp = (np.mean(b_p) - np.mean(a_p) + 12) % 24 - 12
            d_acro.append(abs(dp))
        if not d_amp:
            continue
        nf_rows.append(dict(dataset=ds, variable=var,
                            amp_ratio_floor=float(np.median(d_amp)),
                            amp_ratio_p95=float(np.percentile(d_amp, 95)),
                            mesor_floor=float(np.median(d_mesor)),
                            mesor_p95=float(np.percentile(d_mesor, 95)),
                            acrophase_floor_h=float(np.median(d_acro)),
                            acrophase_p95_h=float(np.percentile(d_acro, 95)),
                            n_boot=len(d_amp)))
    nf = pd.DataFrame(nf_rows)
    nf.to_csv(os.path.join(OUT, "noise_floor.csv"), index=False, encoding="utf-8-sig")

    L("데이터셋별 잡음 바닥 (중앙값 / 95백분위)")
    L(f"{'데이터셋':16s} {'변수':10s} {'진폭비':>16s} {'MESOR차':>16s} {'acrophase차(h)':>18s}")
    for _, r in nf.iterrows():
        L(f"{r.dataset:16s} {r.variable:10s} "
          f"{r.amp_ratio_floor:7.3f}/{r.amp_ratio_p95:<7.3f} "
          f"{r.mesor_floor:7.3f}/{r.mesor_p95:<7.3f} "
          f"{r.acrophase_floor_h:8.2f}/{r.acrophase_p95_h:<8.2f}")

    txt = "\n".join(log)
    print("\n" + txt)
    with open(os.path.join(RES, "gate_log.txt"), "w", encoding="utf-8") as fh:
        fh.write(txt + "\n")
    print(f"\n-> data/rhythm/params.csv, data/rhythm/noise_floor.csv, results/v3/gate_log.txt")
    return 0 if gate_ok else 1


if __name__ == "__main__":
    sys.exit(main())
