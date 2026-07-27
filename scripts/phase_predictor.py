"""
Week 5: 단일 샘플 위상 추정기 (single-sample phase inference).

왜 필요한가
  우주비행 전사체는 희생 시각이 기록돼 있지 않다(11개 인공중력 스터디 전수 확인).
  따라서 clock gene 발현의 분산은 처리효과가 아니라 개체별 내부 위상차가 지배한다.
  각 샘플의 내부 위상을 추정할 수 있으면
    (a) 위상을 공변량으로 보정해 검정력을 회복하고
    (b) '인공중력이 개체 간 위상 동기를 유지하는가' 라는 진짜 질문을 직접 검정할 수 있다.

방법 — molecular timetable / cosinor 기반
  학습: Zhang et al. 2014 마우스 circadian atlas (GSE54650)
        12조직 x 24시점(CT18-64, 2h 간격), Affymetrix GPL6246
  각 유전자에 대해 z 표준화 후 cosinor 적합: z_g(CT) = A_g * cos(w(CT - phi_g))
  추정: 관측 z 벡터와 각 후보 CT 의 예측 프로파일 간 상관을 최대화

검증 — leave-one-tissue-out (LOTO)
  OSDR 표적 조직(망막/시신경)은 atlas 에 없다. 따라서 '학습에 없던 조직' 으로
  일반화되는지가 유일하게 의미 있는 검증이다. 같은 조직 내 hold-out 은 과대평가된다.
  유전자 선택도 매 fold 안에서만 수행한다(선택 편향 차단).

사전 규정 기준
  LOTO 중앙값 절대오차 < 3.0 h 를 통과하지 못하면 이 분석은 계획서에서 제외한다.

출력: data/phase_predictor_loto.csv, data/phase_predictor_genes.csv
"""
import os, re, gzip, warnings
import requests
import numpy as np
import pandas as pd

warnings.filterwarnings("ignore")
S = requests.Session()
S.headers.update({"User-Agent": "KASA-circadian-gravity/1.0"})

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
CACHE = os.path.join(ROOT, "cache")
OUT = os.path.join(ROOT, "data")
os.makedirs(CACHE, exist_ok=True)
os.makedirs(OUT, exist_ok=True)

SERIES = ("https://ftp.ncbi.nlm.nih.gov/geo/series/GSE54nnn/GSE54650/matrix/"
          "GSE54650_series_matrix.txt.gz")
PLATFORM = "https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPL6nnn/GPL6246/annot/GPL6246.annot.gz"
N_PANEL = 100
MDE_CRITERION_H = 3.0
OMEGA = 2 * np.pi / 24


def download(url, name):
    p = os.path.join(CACHE, name)
    if not os.path.exists(p):
        r = S.get(url, timeout=3600, stream=True)
        r.raise_for_status()
        with open(p, "wb") as f:
            for ch in r.iter_content(1 << 20):
                f.write(ch)
    return p


def load_atlas():
    p = download(SERIES, "GSE54650_series_matrix.txt.gz")
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
    expr = expr.dropna(how="any")
    # GC-RMA '선형 강도' 이므로 log2 로 변환한다 (중앙값 ~100, 최대 ~45,000).
    expr = np.log2(expr.clip(lower=1.0))
    tissue = pd.Series([c.split("_")[0] for c in expr.columns], index=expr.columns)
    ct = pd.Series([float(re.search(r"CT(\d+)", c).group(1)) % 24 for c in expr.columns],
                   index=expr.columns)
    return expr, tissue, ct


def load_symbols():
    p = download(PLATFORM, "GPL6246.annot.gz")
    rows = []
    with gzip.open(p, "rt", encoding="utf-8", errors="replace") as f:
        for line in f:
            if line.startswith("!platform_table_begin"):
                break
        hdr = next(f).rstrip("\n").split("\t")
        gi = hdr.index("Gene symbol") if "Gene symbol" in hdr else 1
        for line in f:
            if line.startswith("!platform_table_end"):
                break
            c = line.rstrip("\n").split("\t")
            if len(c) > gi:
                rows.append((c[0], c[gi]))
    return pd.Series({k: v for k, v in rows if v})


def zscore_within_tissue(expr, tissue):
    """조직별로 유전자 z 표준화. 조직 간 기저 발현 차이를 제거한다."""
    out = expr.copy().astype(float)
    for t in tissue.unique():
        cols = tissue[tissue == t].index
        blk = out[cols]
        out[cols] = blk.sub(blk.mean(axis=1), axis=0).div(blk.std(axis=1) + 1e-9, axis=0)
    return out


def fit_cosinor(Z, ct):
    """유전자별 cosinor 적합. 반환: 진폭 A, 위상 phi, 설명력 R2."""
    X = np.column_stack([np.ones(len(ct)), np.cos(OMEGA * ct), np.sin(OMEGA * ct)])
    Y = Z.to_numpy(float).T                              # (n_samples, n_genes)
    beta, *_ = np.linalg.lstsq(X, Y, rcond=None)         # (3, n_genes)
    fit = X @ beta
    ss_res = ((Y - fit) ** 2).sum(axis=0)
    ss_tot = ((Y - Y.mean(axis=0)) ** 2).sum(axis=0)
    r2 = 1 - ss_res / np.where(ss_tot > 1e-12, ss_tot, np.nan)
    # A*cos(w t - phi) = A cos(phi) cos(w t) + A sin(phi) sin(w t)
    #   => beta1 = A cos(phi), beta2 = A sin(phi), 따라서 phi = atan2(beta2, beta1).
    # 부호를 뒤집으면 유전자마다 위상이 반사되어 템플릿이 뒤섞이고
    # 예측 정확도가 무작위 수준으로 붕괴한다(실제로 초기 구현에서 발생).
    A = np.hypot(beta[1], beta[2])
    phi = np.arctan2(beta[2], beta[1]) % (2 * np.pi)     # peak 위상 (rad)
    return pd.DataFrame({"A": A, "phi": phi, "r2": r2}, index=Z.index)


def predict_ct(Z_test, panel_fit, grid=None):
    """관측 z 벡터와 후보 CT 예측 프로파일의 상관 최대화."""
    grid = np.arange(0, 24, 0.1) if grid is None else grid
    A = panel_fit["A"].to_numpy()[:, None]
    phi = panel_fit["phi"].to_numpy()[:, None]
    tmpl = A * np.cos(OMEGA * grid[None, :] - phi)       # (n_genes, n_grid)
    tmpl = tmpl - tmpl.mean(axis=0, keepdims=True)
    tmpl = tmpl / (tmpl.std(axis=0, keepdims=True) + 1e-9)
    Obs = Z_test.to_numpy(float)                          # (n_genes, n_samples)
    Obs = Obs - Obs.mean(axis=0, keepdims=True)
    Obs = Obs / (Obs.std(axis=0, keepdims=True) + 1e-9)
    corr = tmpl.T @ Obs / len(A)                          # (n_grid, n_samples)
    return grid[np.argmax(corr, axis=0)], corr.max(axis=0)


def circ_err(pred, true):
    d = (pred - true + 12) % 24 - 12
    return np.abs(d), d


def main():
    print("아틀라스 로드 중...")
    expr, tissue, ct = load_atlas()
    sym = load_symbols()
    print(f"  프로브 {expr.shape[0]:,} x 샘플 {expr.shape[1]}  조직 {tissue.nunique()}종")
    print(f"  조직: {sorted(tissue.unique())}")
    print(f"  CT 분포: {sorted(ct.unique())}")

    Z = zscore_within_tissue(expr, tissue)

    # ---------------- LOTO 검증
    print(f"\n[검증] leave-one-tissue-out  (사전 기준: 중앙값 절대오차 < {MDE_CRITERION_H} h)")
    recs, per_sample = [], []
    for held in sorted(tissue.unique()):
        tr_cols = tissue[tissue != held].index
        te_cols = tissue[tissue == held].index
        fit = fit_cosinor(Z[tr_cols], ct[tr_cols].to_numpy())     # 학습 조직만으로 적합
        panel = fit.nlargest(N_PANEL, "r2")                       # 유전자 선택도 fold 내부
        pred, conf = predict_ct(Z.loc[panel.index, te_cols], panel)
        err, signed = circ_err(pred, ct[te_cols].to_numpy())
        recs.append({"held_out_tissue": held, "n": len(te_cols),
                     "median_abs_err_h": round(float(np.median(err)), 2),
                     "mean_abs_err_h": round(float(err.mean()), 2),
                     "frac_within_3h": round(float((err < 3).mean()), 3),
                     "median_r2_panel": round(float(panel.r2.median()), 3)})
        for s, p_, t_, e_ in zip(te_cols, pred, ct[te_cols], err):
            per_sample.append({"tissue": held, "sample": s, "true_CT": t_,
                               "pred_CT": round(float(p_), 2), "abs_err_h": round(float(e_), 2)})

    res = pd.DataFrame(recs).sort_values("median_abs_err_h")
    print(res.to_string(index=False))
    overall = float(np.median([r["abs_err_h"] for r in per_sample]))
    within = float(np.mean([r["abs_err_h"] < 3 for r in per_sample]))
    print(f"\n  전체 중앙값 절대오차: {overall:.2f} h")
    print(f"  3시간 이내 비율:      {within:.1%}")
    print(f"  무작위 추정 기준선:    6.00 h (균등분포 기대 절대오차)")
    verdict = "통과" if overall < MDE_CRITERION_H else "미달"
    print(f"\n  >>> 사전 기준 {MDE_CRITERION_H} h : {verdict}")
    if verdict == "미달":
        print("      계획서에서 이 분석을 제외하거나, 방법을 교체해야 한다.")

    pd.DataFrame(per_sample).to_csv(os.path.join(OUT, "phase_predictor_loto.csv"),
                                    index=False, encoding="utf-8-sig")

    # ---------------- 현실 조건 검증
    # 위 LOTO 는 낙관적이다. 보류 조직을 z 표준화할 때 24시간에 고루 퍼진 24개 샘플을 썼다.
    # OSDR 실제 상황은 다르다 — 한 스터디의 샘플은 사실상 같은 시각에 채취됐다.
    # 공통 위상이 센터링으로 제거되므로 '절대 CT' 는 원리적으로 복원 불가능하다.
    # 그러나 우리의 estimand 는 절대 CT 가 아니라 군 내 위상 산포다.
    # 따라서 좁은 시간창 안에서 '상대 위상차' 를 복원할 수 있는지가 실제 관건이다.
    print("\n[현실 조건 검증] 좁은 시간창 내 상대 위상 복원")
    print("  (한 스터디의 샘플이 모두 비슷한 시각에 채취된 상황을 모사)")
    rows = []
    rng = np.random.default_rng(0)
    for width_h in [4, 6, 8, 10]:
        errs, corrs = [], []
        for held in sorted(tissue.unique()):
            tr_cols = tissue[tissue != held].index
            fit = fit_cosinor(Z[tr_cols], ct[tr_cols].to_numpy())
            panel = fit.nlargest(N_PANEL, "r2")
            te = tissue[tissue == held].index
            te_ct = ct[te]
            for start in range(0, 24, 2):
                sel = [c for c in te if ((te_ct[c] - start) % 24) < width_h]
                if len(sel) < 4:
                    continue
                # 이 부분집합만으로 센터링 (실제로 가용한 정보만 사용)
                blk = expr[sel].astype(float)
                Zi = blk.sub(blk.mean(axis=1), axis=0).div(blk.std(axis=1) + 1e-9, axis=0)
                pred, _ = predict_ct(Zi.loc[panel.index], panel)
                true = te_ct[sel].to_numpy()
                # 상대 위상: 각자의 평균을 뺀 뒤 비교 (절대 CT 는 복원 불가)
                pr = (pred - pred.mean() + 12) % 24 - 12
                tr_ = (true - true.mean() + 12) % 24 - 12
                errs.append(np.abs(pr - tr_))
                if np.std(tr_) > 1e-9 and np.std(pr) > 1e-9:
                    corrs.append(np.corrcoef(pr, tr_)[0, 1])
        e = np.concatenate(errs)
        rows.append({"window_h": width_h, "n_subsets": len(errs),
                     "median_rel_err_h": round(float(np.median(e)), 2),
                     "frac_within_2h": round(float((e < 2).mean()), 3),
                     "mean_corr_true_vs_pred": round(float(np.mean(corrs)), 3)})
    rel = pd.DataFrame(rows)
    print(rel.to_string(index=False))
    print("  → 시간창이 좁을수록 어렵다. 실제 OSDR 스터디는 창이 사실상 0 에 가깝다.")
    print("  → 절대 CT 는 포기하고, 군 간 '위상 산포 비교' 로만 쓸 것.")
    rel.to_csv(os.path.join(OUT, "phase_predictor_relative.csv"),
               index=False, encoding="utf-8-sig")

    # ---------------- 전체 데이터로 최종 패널
    fit_all = fit_cosinor(Z, ct.to_numpy())
    panel_all = fit_all.nlargest(N_PANEL, "r2").copy()
    panel_all["symbol"] = [sym.get(i, "") for i in panel_all.index]
    panel_all["peak_CT_h"] = (panel_all.phi / OMEGA) % 24
    panel_all.round(4).to_csv(os.path.join(OUT, "phase_predictor_genes.csv"),
                              encoding="utf-8-sig")
    print(f"\n[패널] 전체 학습 상위 {N_PANEL}개 중 clock 유전자:")
    known = ["Bmal1", "Arntl", "Npas2", "Per1", "Per2", "Per3", "Cry1", "Cry2",
             "Nr1d1", "Nr1d2", "Dbp", "Tef", "Hlf", "Ciart", "Bhlhe40", "Bhlhe41"]
    hit = panel_all[panel_all.symbol.isin(known)]
    print(hit[["symbol", "r2", "peak_CT_h"]].round(3).to_string())
    print(f"  (상위 패널 {N_PANEL}개 중 알려진 코어 clock 유전자 {len(hit)}개 — "
          f"패널이 실제로 시계를 잡고 있다는 확인)")
    print(f"\n-> {OUT}/phase_predictor_loto.csv, phase_predictor_genes.csv")


if __name__ == "__main__":
    main()
