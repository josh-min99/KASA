"""
Week 6: 웻랩 설계 파라미터 산출 — 표본수 / 시점수 / 시점 간격.

입력 (전부 앞 단계에서 실측·산출된 값이며 임의 가정이 아니다)
  진폭 A   : Zhang 2014 아틀라스에서 코어 clock 유전자의 log2 진폭 실측
  잡음 sd  : OSDR 우주비행 마우스 조직의 군내 SD 실측
             (망막 0.343 / 시신경 0.670 / 가자미근 0.504 log2 — scripts/clock_dose_response.py)
  효과크기 : 진동자 모델 PRC 예측 위상 이동
             (CT16 인가 -4.13 h / CT21 인가 +1.87 h — scripts/oscillator_model.py)

검정
  cosinor 모델  y = M_g + A_g*cos(w(t - phi_g)) + eps
  H0: phi_1 = phi_2   vs   H1: phi_1 != phi_2
  phi 가 주어지면 (M, A) 에 대해 선형이므로 phi 격자 위에서 profile 해 정확히 푼다.
  검정통계량은 F 통계량, 유의수준 0.05.

출력: data/power_curve.csv, data/design_comparison.csv
"""
import os, re, gzip, warnings
import numpy as np
import pandas as pd
from scipy import stats

warnings.filterwarnings("ignore")

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
CACHE = os.path.join(ROOT, "cache")
OUT = os.path.join(ROOT, "data")
os.makedirs(OUT, exist_ok=True)

OMEGA = 2 * np.pi / 24
ALPHA = 0.05
N_SIM = 4000
PHI_GRID = np.arange(0, 24, 0.05)

# OSDR 실측 군내 SD (log2)
SD_MEASURED = {"망막 (OSD-758)": 0.343, "가자미근 (OSD-714)": 0.504, "시신경 (OSD-759)": 0.670}
# 모델 PRC 예측 위상 이동 (h)
DPHI_MODEL = {"CT21 인가 (전진)": 1.87, "CT16 인가 (지연)": 4.13}

CLOCK = ["Arntl", "Npas2", "Per1", "Per2", "Per3", "Cry1", "Cry2",
         "Nr1d1", "Nr1d2", "Dbp", "Tef", "Hlf", "Ciart", "Bhlhe41"]


# ------------------------------------------------ 진폭 추정 (아틀라스 실측)
def estimate_amplitude():
    p = os.path.join(CACHE, "GSE54650_series_matrix.txt.gz")
    ann = os.path.join(CACHE, "GPL6246.annot.gz")
    if not (os.path.exists(p) and os.path.exists(ann)):
        print("  아틀라스 캐시 없음 — scripts/phase_predictor.py 를 먼저 실행할 것")
        return None
    titles = None
    with gzip.open(p, "rt", encoding="utf-8", errors="replace") as f:
        for line in f:
            if line.startswith("!Sample_title"):
                titles = [t.strip('"') for t in line.rstrip("\n").split("\t")[1:]]
            if line.startswith("!series_matrix_table_begin"):
                break
        expr = pd.read_csv(f, sep="\t", index_col=0, comment="!")
    expr.index = expr.index.astype(str).str.strip('"')
    expr.columns = titles[:expr.shape[1]]
    expr = expr.dropna()
    # GSE54650 은 GC-RMA 정규화된 '선형 강도' 다 (중앙값 ~100, 최대 ~45,000).
    # log2 로 바꾸지 않으면 진폭이 강도 단위로 나와 OSDR 의 log2 SD 와 단위가 어긋난다
    # (초기 구현에서 진폭 173.6 'log2' 로 나와 모든 검정력이 1.00 이 되는 형태로 드러났다).
    expr = np.log2(expr.clip(lower=1.0))

    rows = []
    with gzip.open(ann, "rt", encoding="utf-8", errors="replace") as f:
        for line in f:
            if line.startswith("!platform_table_begin"):
                break
        hdr = next(f).rstrip("\n").split("\t")
        gi = hdr.index("Gene symbol")
        for line in f:
            if line.startswith("!platform_table_end"):
                break
            c = line.rstrip("\n").split("\t")
            if len(c) > gi and c[gi] in CLOCK:
                rows.append((c[0], c[gi]))
    probe2sym = dict(rows)

    tissue = pd.Series([c.split("_")[0] for c in expr.columns], index=expr.columns)
    ct = pd.Series([float(re.search(r"CT(\d+)", c).group(1)) % 24 for c in expr.columns],
                   index=expr.columns)

    recs = []
    for probe, symbol in probe2sym.items():
        if probe not in expr.index:
            continue
        for t in tissue.unique():
            cols = tissue[tissue == t].index
            y = expr.loc[probe, cols].to_numpy(float)
            x = ct[cols].to_numpy(float)
            X = np.column_stack([np.ones_like(x), np.cos(OMEGA * x), np.sin(OMEGA * x)])
            b, *_ = np.linalg.lstsq(X, y, rcond=None)
            fit = X @ b
            r2 = 1 - ((y - fit) ** 2).sum() / max(((y - y.mean()) ** 2).sum(), 1e-12)
            recs.append({"symbol": symbol, "tissue": t,
                         "amplitude_log2": float(np.hypot(b[1], b[2])), "r2": r2})
    df = pd.DataFrame(recs)
    strong = df[df.r2 > 0.5]
    print(f"  cosinor R2>0.5 인 (유전자 x 조직) 조합 {len(strong)} / {len(df)}")
    q = strong.amplitude_log2.quantile([0.25, 0.5, 0.75])
    print(f"  clock 유전자 log2 진폭:  25% {q[0.25]:.3f}  중앙값 {q[0.5]:.3f}  75% {q[0.75]:.3f}")
    top = strong.groupby("symbol").amplitude_log2.median().sort_values(ascending=False)
    print("  유전자별 중앙 진폭 (상위):")
    print("   " + ", ".join(f"{k} {v:.2f}" for k, v in top.head(8).items()))
    return float(q[0.5]), float(q[0.25]), float(q[0.75])


# ------------------------------------------------------------ cosinor 검정
def _fit(Y, t):
    """균형 설계(24h 등간격, 시점당 동수)에서 [1, cos, sin] 은 직교한다.
    따라서 lstsq 없이 해석적으로 푼다. Y shape (B, n)."""
    n = t.size
    c, s = np.cos(OMEGA * t), np.sin(OMEGA * t)
    M = Y.mean(axis=1)
    b1 = 2.0 / n * (Y @ c)
    b2 = 2.0 / n * (Y @ s)
    rss = (Y ** 2).sum(axis=1) - n * M ** 2 - (n / 2) * (b1 ** 2 + b2 ** 2)
    return M, b1, b2, np.maximum(rss, 1e-12)


def _check_orthogonal(t):
    c, s = np.cos(OMEGA * t), np.sin(OMEGA * t)
    n = t.size
    assert abs(c.sum()) < 1e-8 and abs(s.sum()) < 1e-8, "설계가 불균형: cos/sin 합이 0이 아님"
    assert abs(c @ s) < 1e-8, "cos·sin 직교 실패"
    assert abs(c @ c - n / 2) < 1e-8 and abs(s @ s - n / 2) < 1e-8, "노름이 n/2 가 아님"


def simulate_power(n_per_tp, n_tp, dphi, sd, amp, n_sim=N_SIM, seed=0):
    """
    H0: phi_1 = phi_2   vs   H1: 자유

    위상 phi 로 제약했을 때의 RSS 는 전체모형 RSS 에
        (n/2) * (b1*sin(w*phi) - b2*cos(w*phi))^2
    만큼 더해진다 (계수벡터 b 를 방향 u(phi) 에 사영하고 남은 직교성분).
    따라서 phi 격자 위에서 산술만으로 profile 할 수 있다.
    """
    rng = np.random.default_rng(seed)
    t = np.repeat(np.arange(n_tp) * (24.0 / n_tp), n_per_tp)
    _check_orthogonal(t)
    n = t.size
    B = n_sim

    mu1 = amp * np.cos(OMEGA * (t - 6.0))
    mu2 = amp * np.cos(OMEGA * (t - 6.0 - dphi))
    Y1 = mu1 + rng.normal(0, sd, (B, n))
    Y2 = mu2 + rng.normal(0, sd, (B, n))

    _, b11, b12, rss1 = _fit(Y1, t)
    _, b21, b22, rss2 = _fit(Y2, t)

    cw = np.cos(OMEGA * PHI_GRID)[None, :]        # (1, G)
    sw = np.sin(OMEGA * PHI_GRID)[None, :]
    pen = (n / 2) * ((b11[:, None] * sw - b12[:, None] * cw) ** 2
                     + (b21[:, None] * sw - b22[:, None] * cw) ** 2)
    rss_h1 = rss1 + rss2
    rss_h0 = rss_h1 + pen.min(axis=1)

    dfd = 2 * n - 6
    F = ((rss_h0 - rss_h1) / 1.0) / (rss_h1 / dfd)
    p = 1 - stats.f.cdf(F, 1, dfd)
    return float((p < ALPHA).mean())


def main():
    print("=" * 78)
    print("[1] 진폭 추정 — Zhang 2014 아틀라스 실측")
    est = estimate_amplitude()
    amp_med, amp_lo, amp_hi = est if est else (0.5, 0.3, 0.9)

    print("\n[2] 검정력 곡선 — 군당 마리 수 x 시점 수")
    print(f"    진폭 {amp_med:.2f} log2, 유의수준 {ALPHA}, 시뮬레이션 {N_SIM}회/조건")
    rows = []
    for sd_label, sd in SD_MEASURED.items():
        for dphi_label, dphi in DPHI_MODEL.items():
            for n_tp in [3, 4, 6]:
                for n_per in [3, 4, 5, 6, 8, 10]:
                    pw = simulate_power(n_per, n_tp, dphi, sd, amp_med)
                    rows.append({"tissue_sd": sd_label, "sd": sd,
                                 "effect": dphi_label, "dphi_h": dphi,
                                 "n_timepoints": n_tp, "n_per_timepoint": n_per,
                                 "total_mice": n_tp * n_per * 2, "power": pw})
    pc = pd.DataFrame(rows)
    pc.to_csv(os.path.join(OUT, "power_curve.csv"), index=False, encoding="utf-8-sig")

    for dphi_label in DPHI_MODEL:
        print(f"\n  --- 효과 {dphi_label} (Δφ={DPHI_MODEL[dphi_label]} h) ---")
        sub = pc[pc.effect == dphi_label]
        piv = sub.pivot_table(index=["tissue_sd", "n_timepoints"],
                              columns="n_per_timepoint", values="power")
        print(piv.round(2).to_string())

    print("\n[3] 80% 검정력 달성에 필요한 최소 설계")
    need = []
    for (sd_label, dphi_label, n_tp), g in pc.groupby(["tissue_sd", "effect", "n_timepoints"]):
        ok = g[g.power >= 0.80].sort_values("total_mice")
        need.append({"tissue": sd_label, "effect": dphi_label, "n_timepoints": n_tp,
                     "min_n_per_tp": int(ok.iloc[0].n_per_timepoint) if len(ok) else None,
                     "total_mice": int(ok.iloc[0].total_mice) if len(ok) else None})
    nd = pd.DataFrame(need)
    print(nd.to_string(index=False))
    nd.to_csv(os.path.join(OUT, "design_comparison.csv"), index=False, encoding="utf-8-sig")

    print("\n[4] 시점 수 축소가 정당한가 — 총 마리 수를 고정하고 배분만 변경")
    print("    검정력이 포화되지 않는 예산 구간에서 비교해야 의미가 있다.")
    rows4 = []
    for budget in [12, 18, 24, 36]:
        print(f"\n    예산 {budget}마리/군, 효과 Δφ=4.13 h:")
        for sd_label, sd in SD_MEASURED.items():
            cells = []
            for n_tp in [3, 4, 6]:
                n_per = budget // n_tp
                if n_per < 2:
                    cells.append(f"{n_tp}시점: n/a")
                    continue
                pw = simulate_power(n_per, n_tp, DPHI_MODEL["CT16 인가 (지연)"], sd, amp_med)
                cells.append(f"{n_tp}시점x{n_per}마리={pw:.2f}")
                rows4.append({"budget_per_group": budget, "tissue_sd": sd_label,
                              "n_timepoints": n_tp, "n_per_timepoint": n_per, "power": pw})
            print(f"      {sd_label:>18}  " + "   ".join(cells))
    pd.DataFrame(rows4).to_csv(os.path.join(OUT, "design_timepoint_tradeoff.csv"),
                               index=False, encoding="utf-8-sig")
    print("\n    ※ 주의: 주기를 24h 로 고정한 cosinor 위상 검정은 (M, A, phi) 식별에")
    print("       최소 3시점이면 충분하다. 시점을 늘려 얻는 것은 위상 검정의 검정력이 아니라")
    print("       주기 오설정과 비정현 파형에 대한 강건성이다. 3시점 설계는 주기가")
    print("       정확히 24h 라는 가정에 전적으로 의존하므로, 자유진동 실험에는 쓸 수 없다.")

    print("\n[5] 진폭 불확실성 민감도 (사분위 범위)")
    for label, a in [("25분위", amp_lo), ("중앙값", amp_med), ("75분위", amp_hi)]:
        pw = simulate_power(6, 6, DPHI_MODEL["CT16 인가 (지연)"], SD_MEASURED["망막 (OSD-758)"], a)
        print(f"    진폭 {label} {a:.2f} log2 → 6시점x6마리 검정력 {pw:.2f}")

    print(f"\n-> {OUT}/power_curve.csv, design_comparison.csv")


if __name__ == "__main__":
    main()
