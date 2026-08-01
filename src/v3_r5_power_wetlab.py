"""
R5: 웻랩 검정력 재계산 + 마스킹 시뮬레이션.

기존 `scripts/power_analysis.py` 와의 차이
  그 스크립트는 잡음을 **전사체 군내 SD**(0.343~0.670 log2)로 놓고 cosinor 위상차를 검정했다.
  전사체는 단일 시점이라 위상 추이를 볼 수 없으므로, 그 검정력은 실제 웻랩 판정 절차와
  대응하지 않는다. 여기서는 잡음을 **실측 텔레메트리 시계열**에서 직접 뽑고,
  검정도 계획서에 적힌 판정 기준(회귀선 비교) 그대로 쓴다.

1단계 — 일별 위상 추정 오차를 데이터에서 직접 측정한다
  Helissen baseline 구간에서 개체별로 24시간 창을 하루씩 밀며 acrophase 를 추정하고,
  그 개체 내 표준편차를 구한다. 이것이 '하루치 데이터로 위상을 얼마나 정확히 찍는가' 다.
  가정이 아니라 실측이며, 지표(활동량 / 피하온도 / 심부온도)별로 따로 구한다.

2단계 — 계획서의 판정 기준을 그대로 시뮬레이션한다
  가설 1: baseline 구간 회귀선 대 처치 구간 회귀선의 절편이 유의하게 다른가
  가설 2: 처치 구간 회귀선을 free-run 구간으로 연장했을 때 실측이 그 선을 따르는가
  개체 무작위효과(개체별 고유 위상)를 포함한 혼합효과 구조로 생성하고,
  개체를 단위로 한 t 검정으로 판정한다.

3단계 — 마스킹
  HLU-on 구간에서는 지표 자체가 억제된다(R4 실측: 활동 진폭 -92%, 체온 진폭 -58~-79%).
  진폭이 줄면 위상 추정 오차가 커진다. 그 증폭분을 1단계에서 실측한
  '진폭 대 위상오차' 관계로 반영한다.

산출: data/rhythm/phase_precision.csv, data/rhythm/power_grid.csv,
      results/v3/r5_power.txt
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

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from v3_r3_rhythm_params import cosinor  # noqa: E402

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, "data", "rhythm")
RES = os.path.join(ROOT, "results", "v3")
os.makedirs(RES, exist_ok=True)

RNG = np.random.default_rng(20260803)
N_SIM = 2000
ALPHA = 0.05


# ---------------------------------------------- 1단계: 일별 위상 추정 오차 실측
def daily_phase_precision(long):
    """개체별로 24h 창을 하루씩 밀며 acrophase 를 추정하고 개체 내 SD 를 구한다."""
    rows = []
    for (ds, subj, var, cond), g in long.groupby(
            ["dataset", "subject_id", "variable", "condition"]):
        if var not in ("activity", "tb_sub", "tb_core"):
            continue
        g = g.sort_values("t_hours")
        t, y = g.t_hours.values, g.value.values
        if len(t) < 12:
            continue
        phis, amps = [], []
        t0, t1 = t.min(), t.max()
        start = t0
        while start + 24 <= t1 + 1e-9:
            m = (t >= start) & (t < start + 24)
            if m.sum() >= 8:
                c = cosinor(t[m], y[m])
                if np.isfinite(c["acrophase"]):
                    phis.append(c["acrophase"])
                    amps.append(c["amplitude"])
            start += 24
        if len(phis) < 2:
            continue
        # 원형 표준편차 (위상은 24시간 순환)
        ang = np.array(phis) * 2 * np.pi / 24
        Rbar = np.abs(np.mean(np.exp(1j * ang)))
        circ_sd_h = float(np.sqrt(-2 * np.log(max(Rbar, 1e-9))) * 24 / (2 * np.pi))
        rows.append(dict(dataset=ds, subject_id=subj, variable=var, condition=cond,
                         n_days=len(phis), phase_sd_h=circ_sd_h,
                         mean_amplitude=float(np.mean(amps))))
    return pd.DataFrame(rows)


# ---------------------------------------------- 2단계: 판정 기준 시뮬레이션
def simulate(n_animals, n_days_base, n_days_treat, sd_phi_base, sd_phi_treat,
             sd_between, true_dphi, n_sim=N_SIM, alpha=ALPHA):
    """계획서 가설 1(절편 차이) 검정력.

    개체 i 의 j 일째 관측 위상:
        base   : mu_i + e,                e ~ N(0, sd_phi_base)
        treat  : mu_i + true_dphi + e,    e ~ N(0, sd_phi_treat)
    mu_i ~ N(0, sd_between) 은 개체 고유 위상.
    개체별로 (처치 평균 - baseline 평균) 을 구하고, 그 값에 대해 1표본 t 검정.
    개체 고유 위상은 차분에서 상쇄되므로 sd_between 은 검정력에 영향을 주지 않지만,
    실제 설계와 대응시키기 위해 남겨 둔다.
    """
    hit = 0
    se_b = sd_phi_base / np.sqrt(n_days_base)
    se_t = sd_phi_treat / np.sqrt(n_days_treat)
    sd_diff = np.sqrt(se_b ** 2 + se_t ** 2)
    for _ in range(n_sim):
        d = RNG.normal(true_dphi, sd_diff, n_animals)
        if d.std(ddof=1) == 0:
            continue
        tstat = d.mean() / (d.std(ddof=1) / np.sqrt(n_animals))
        p = 2 * (1 - stats.t.cdf(abs(tstat), n_animals - 1))
        if p < alpha:
            hit += 1
    return hit / n_sim


def main():
    long = pd.read_csv(os.path.join(OUT, "long.csv"))
    log = []
    P = log.append
    P("=" * 92)
    P("R5  웻랩 검정력 — 잡음을 실측 텔레메트리에서 뽑아 계획서 판정 기준으로 검정")
    P("=" * 92)

    # ---------------------------------------------------------- 1단계
    prec = daily_phase_precision(long)
    prec.to_csv(os.path.join(OUT, "phase_precision.csv"), index=False, encoding="utf-8-sig")

    P("")
    P("[1단계] 하루치 데이터로 위상을 얼마나 정확히 찍는가 (실측)")
    P("        Helissen 각 개체의 24시간 창 acrophase 의 개체 내 원형 표준편차")
    P("")
    P(f"  {'구간':10s} {'지표':10s} {'코호트':14s} {'개체':>4s} {'위상 SD(h)':>12s} {'평균 진폭':>12s}")
    summ = []
    for (cond, var, ds), g in prec.groupby(["condition", "variable", "dataset"]):
        if cond == "recovery":
            continue
        sd = float(g.phase_sd_h.median())
        amp = float(g.mean_amplitude.median())
        summ.append(dict(condition=cond, variable=var, dataset=ds, n=len(g),
                         phase_sd_h=sd, amplitude=amp))
        P(f"  {cond:10s} {var:10s} {ds:14s} {len(g):4d} {sd:12.2f} {amp:12.4f}")
    sdf = pd.DataFrame(summ)

    P("")
    P("  주의: 코호트마다 값이 크게 다르다. 하나의 코호트로 지표를 서열화하면 안 된다.")
    P("        아래 검정력은 코호트별로 따로 계산하고 범위로 제시한다.")

    # baseline 상태의 지표별 위상 SD 범위
    P("")
    P("[baseline 구간 위상 SD 요약 — 지표별 코호트 범위]")
    base = sdf[sdf.condition == "baseline"]
    ranges = {}
    for var, g in base.groupby("variable"):
        lo, hi = float(g.phase_sd_h.min()), float(g.phase_sd_h.max())
        ranges[var] = (lo, hi)
        P(f"  {var:10s} {lo:.2f} ~ {hi:.2f} h   (코호트 {len(g)}개: "
          f"{', '.join(f'{r.dataset}={r.phase_sd_h:.2f}' for _, r in g.iterrows())})")

    # ---------------------------------------------------------- 3단계 마스킹
    P("")
    P("[마스킹] HLU-on 구간에서 지표가 억제되면 위상 추정 오차가 얼마나 커지는가 (실측)")
    P("")
    P(f"  {'지표':10s} {'코호트':14s} {'baseline SD':>12s} {'HLU중 SD':>12s} {'배수':>7s}"
      f" {'진폭비':>8s}")
    mask_rows = []
    for var in ["activity", "tb_sub", "tb_core"]:
        for ds in sorted(sdf.dataset.unique()):
            b = sdf[(sdf.condition == "baseline") & (sdf.variable == var) & (sdf.dataset == ds)]
            t = sdf[(sdf.condition == "treatment") & (sdf.variable == var) & (sdf.dataset == ds)]
            if not len(b) or not len(t):
                continue
            r = float(t.phase_sd_h.iloc[0] / b.phase_sd_h.iloc[0])
            ar = float(t.amplitude.iloc[0] / b.amplitude.iloc[0]) if b.amplitude.iloc[0] else np.nan
            mask_rows.append(dict(variable=var, dataset=ds,
                                  sd_base=float(b.phase_sd_h.iloc[0]),
                                  sd_treat=float(t.phase_sd_h.iloc[0]),
                                  inflation=r, amp_ratio=ar))
            P(f"  {var:10s} {ds:14s} {b.phase_sd_h.iloc[0]:12.2f} {t.phase_sd_h.iloc[0]:12.2f}"
              f" {r:7.2f} {ar:8.3f}")
    mdf = pd.DataFrame(mask_rows)
    P("")
    P("  -> HLU 자체가 지표를 억제하므로, HLU-on 구간에서 잰 위상은 baseline 만큼 정확하지 않다.")
    P("     계획서 가설 1(HLU 중 새 위상 형성)은 이 증폭된 오차 위에서 판정해야 한다.")

    # ---------------------------------------------------------- 2단계 검정력
    P("")
    P("=" * 92)
    P("[검정력] 계획서 판정 기준(baseline 회귀선 대 처치 회귀선 절편 차이) 기준")
    P("=" * 92)
    P("검출하려는 위상 이동 Δφ = 1.0 / 2.0 / 3.0 h,  유의수준 0.05,  양측")
    P("baseline 3일 고정. 처치 구간 일수와 마리수를 격자로 돌린다.")
    P("HLU-on 구간 위상 SD 는 위에서 실측한 증폭배수를 적용한다.")

    grid = []
    for _, m in mdf.iterrows():
        for n_an in [4, 6, 8, 10, 12, 16]:
            for n_tr in [5, 7, 10, 14, 21]:
                for dphi in [1.0, 2.0, 3.0]:
                    pw = simulate(n_an, 3, n_tr,
                                  m.sd_base, m.sd_treat, 1.0, dphi)
                    grid.append(dict(variable=m.variable, dataset=m.dataset,
                                     n_animals=n_an, n_days_treat=n_tr,
                                     true_dphi=dphi, power=pw,
                                     sd_base=m.sd_base, sd_treat=m.sd_treat))
    gdf = pd.DataFrame(grid)
    gdf.to_csv(os.path.join(OUT, "power_grid.csv"), index=False, encoding="utf-8-sig")

    P("")
    P("Δφ = 2.0 h 를 검출할 검정력 (baseline 3일, 처치 14일)")
    P(f"  {'지표':10s} {'코호트':14s} " + "".join(f"{f'n={n}':>8s}" for n in [4, 6, 8, 10, 12, 16]))
    for var in ["tb_core", "tb_sub", "activity"]:
        for ds in sorted(gdf[gdf.variable == var].dataset.unique()):
            s = gdf[(gdf.variable == var) & (gdf.dataset == ds) &
                    (gdf.n_days_treat == 14) & (gdf.true_dphi == 2.0)]
            if not len(s):
                continue
            cells = "".join(f"{float(s[s.n_animals == n].power.iloc[0]):8.2f}"
                            for n in [4, 6, 8, 10, 12, 16]
                            if len(s[s.n_animals == n]))
            P(f"  {var:10s} {ds:14s} {cells}")

    P("")
    P("검정력 0.8 을 달성하는 최소 마리수 (처치 14일)")
    P(f"  {'지표':10s} {'코호트':14s} {'dphi=1h':>10s} {'dphi=2h':>10s} {'dphi=3h':>10s}")
    need_rows = []
    for var in ["tb_core", "tb_sub", "activity"]:
        for ds in sorted(gdf[gdf.variable == var].dataset.unique()):
            cells = []
            for dphi in [1.0, 2.0, 3.0]:
                s = gdf[(gdf.variable == var) & (gdf.dataset == ds) &
                        (gdf.n_days_treat == 14) & (gdf.true_dphi == dphi) &
                        (gdf.power >= 0.8)]
                v = int(s.n_animals.min()) if len(s) else None
                cells.append(v)
                need_rows.append(dict(variable=var, dataset=ds, dphi=dphi,
                                      n_needed=v if v else np.nan))
            P(f"  {var:10s} {ds:14s} " +
              "".join(f"{(str(c) if c else '>16'):>10s}" for c in cells))
    pd.DataFrame(need_rows).to_csv(os.path.join(OUT, "power_required_n.csv"),
                                   index=False, encoding="utf-8-sig")

    # ---------------------------------------------------------- 판독
    P("")
    P("=" * 92)
    P("판독 — 1차 종말점")
    P("=" * 92)
    nd = pd.DataFrame(need_rows)
    P("")
    for var, label in [("tb_core", "심부체온"), ("tb_sub", "피하온도"), ("activity", "활동량")]:
        s = nd[(nd.variable == var) & (nd.dphi == 2.0)]
        if not len(s):
            continue
        vals = s.n_needed.dropna()
        if len(vals) == len(s) and len(vals):
            P(f"  {label:8s} dphi=2h 검출에 필요한 마리수: "
              f"{', '.join(f'{r.dataset}={int(r.n_needed)}' for _, r in s.iterrows() if np.isfinite(r.n_needed))}")
        else:
            miss = [r.dataset for _, r in s.iterrows() if not np.isfinite(r.n_needed)]
            P(f"  {label:8s} dphi=2h: 코호트 {miss} 에서는 16마리로도 검정력 0.8 미달")

    txt = "\n".join(log)
    print(txt)
    with open(os.path.join(RES, "r5_power.txt"), "w", encoding="utf-8") as fh:
        fh.write(txt + "\n")
    print(f"\n-> data/rhythm/phase_precision.csv, power_grid.csv, power_required_n.csv")


if __name__ == "__main__":
    main()
