"""
예비분석 2 (양성대조): 파이프라인이 실제 신호를 잡을 수 있는지 검증.

논리
  clock 모듈에서 null 이 나왔다(scripts/clock_dose_response.py). 그 null 이
  "신호가 없다" 인지 "우리가 분석을 못했다" 인지 구별하려면, 같은 데이터·같은 파이프라인으로
  원논문(OSD-758/759)이 보고한 효과가 재현되는지 보여야 한다.

  원논문 주장: "Adding artificial gravity on board the ISS can attenuate the
  transcriptomic response to microgravity in a dose-dependent manner."

검정 3종
  A) DEG 개수 사다리 — 지상대조 대비 DEG 수가 uG > 0.33G > 0.66G > 1G 로 감소하는가
  B) 우주비행 반응 점수 — uG-vs-지상 시그니처에 대한 개체별 투영이 G 에 따라 감소하는가
     (uG 군 순환논리를 피하기 위해 leave-one-out 으로 시그니처를 매번 재추정)
  C) 경로 모듈 — 원논문이 지목한 산화스트레스/염증/세포사멸/지질대사 모듈의 G 용량반응
     → 같은 검정을 clock 모듈에도 적용해 직접 비교

출력: data/positive_control_*.csv
"""
import os, re, io, zipfile, warnings
import requests
import numpy as np
import pandas as pd
from scipy import stats

warnings.filterwarnings("ignore")
S = requests.Session()
S.headers.update({"User-Agent": "KASA-circadian-gravity/1.0"})

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
CACHE = os.path.join(ROOT, "cache")
OUT = os.path.join(ROOT, "data")
os.makedirs(CACHE, exist_ok=True)
os.makedirs(OUT, exist_ok=True)

TARGETS = [(758, "GLDS-664", "Retina"), (759, "GLDS-665", "Optic nerve")]

GMAP = {"uG": 0.0, "0.33G by centrifugation": 0.33,
        "0.66G by centrifugation": 0.66, "1G by centrifugation": 1.0}

# 원논문이 명시한 4개 경로 + 비교용 clock 모듈
MODULES = {
    "oxidative_stress": ["Nfe2l2", "Hmox1", "Nqo1", "Gclc", "Gclm", "Txnrd1", "Gpx1", "Gpx4",
                         "Sod1", "Sod2", "Cat", "Prdx1", "Prdx3", "Srxn1", "Gsta3", "Gstp1",
                         "Slc7a11", "Ftl1", "Fth1"],
    "inflammation":     ["Tnf", "Il1b", "Il6", "Nfkb1", "Nfkbia", "Ccl2", "Ccl3", "Cxcl1",
                         "Cxcl10", "Tlr4", "Ptgs2", "Icam1", "Vcam1", "Casp1", "Nlrp3",
                         "Stat3", "Socs3", "Cd68", "Aif1", "Gfap"],
    "apoptosis":        ["Casp3", "Casp8", "Casp9", "Bax", "Bcl2", "Bcl2l1", "Trp53",
                         "Cdkn1a", "Bid", "Apaf1"],
    "lipid_metabolism": ["Srebf1", "Fasn", "Scd1", "Acaca", "Cpt1a", "Ppara", "Pparg",
                         "Lpl", "Cd36", "Hmgcr"],
    "circadian_clock":  ["Bmal1", "Clock", "Npas2", "Per1", "Per2", "Per3", "Cry1", "Cry2",
                         "Nr1d1", "Nr1d2", "Dbp", "Tef", "Hlf", "Ciart", "Bhlhe40",
                         "Bhlhe41", "Rorb"],
}
N_PERM = 2000
FDR = 0.05


def fetch(osd, pattern):
    fs = S.get(f"https://osdr.nasa.gov/osdr/data/osd/files/{osd}", timeout=120).json()
    fs = fs["studies"][f"OSD-{osd}"]["study_files"]
    f = [x for x in fs if re.search(pattern, x["file_name"])][0]
    path = os.path.join(CACHE, f["file_name"])
    if not os.path.exists(path):
        with open(path, "wb") as fh:
            fh.write(S.get("https://osdr.nasa.gov" + f["remote_url"], timeout=1800).content)
    return path


def sample_table(osd):
    z = zipfile.ZipFile(fetch(osd, r"ISA\.zip$"))
    nm = [n for n in z.namelist() if n.split("/")[-1].startswith("s_")][0]
    return pd.read_csv(io.BytesIO(z.read(nm)), sep="\t", dtype=str)


# ---------------------------------------------------------------- A) DEG 사다리
def deg_ladder(de, glds, label):
    rows = []
    for g_label, g_val in [("uG", 0.0), ("0.33G by centrifugation", 0.33),
                           ("0.66G by centrifugation", 0.66), ("1G by centrifugation", 1.0)]:
        col = f"Adj.p.value_(Space Flight & {g_label})v(Ground Control & 1G on Earth)"
        lfc = f"Log2fc_(Space Flight & {g_label})v(Ground Control & 1G on Earth)"
        if col not in de.columns:
            continue
        sig = de[col] < FDR
        sig_lfc = sig & (de[lfc].abs() > 0.5)
        rows.append({"tissue": label, "gravity_g": g_val, "group": g_label,
                     "n_DEG_fdr05": int(sig.sum()),
                     "n_DEG_fdr05_lfc0.5": int(sig_lfc.sum())})
    return pd.DataFrame(rows)


# ------------------------------------------------- B) 우주비행 반응 점수 (LOO)
def signature(X, cols_a, cols_b, n_top=200):
    """cols_a(uG) vs cols_b(ground) Welch t 통계량 상하위 n_top 유전자."""
    a, b = X[cols_a].to_numpy(float), X[cols_b].to_numpy(float)
    t = stats.ttest_ind(a, b, axis=1, equal_var=False).statistic
    t = pd.Series(t, index=X.index).dropna()
    return t.nlargest(n_top).index, t.nsmallest(n_top).index


def response_score(X, sample, up, dn):
    z = X.sub(X.mean(axis=1), axis=0).div(X.std(axis=1) + 1e-9, axis=0)
    return float(z.loc[up, sample].mean() - z.loc[dn, sample].mean())


def loo_scores(X, flight, ground, g):
    """각 비행 샘플에 대해, 그 샘플을 제외하고 시그니처를 재추정한 뒤 점수 산출."""
    ug = [s for s, gv in zip(flight, g) if gv == 0.0]
    out = {}
    for i, s in enumerate(flight):
        ug_i = [u for u in ug if u != s]
        if len(ug_i) < 2:
            continue
        up, dn = signature(X, ug_i, ground)
        out[s] = response_score(X, s, up, dn)
    return out


# ---------------------------------------------------------- C) 모듈 용량반응
def module_trend(Xf, g, genes, rng):
    """Xf 는 비행 샘플 열만 포함해야 한다 (g 와 길이 일치)."""
    assert Xf.shape[1] == len(g), f"열 {Xf.shape[1]} != G {len(g)}"
    X = Xf
    present = [k for k in genes if k in X.index and np.std(X.loc[k]) > 0]
    if len(present) < 5:
        return None
    rhos = np.array([stats.spearmanr(g, X.loc[k].to_numpy(float))[0] for k in present])
    obs = np.abs(rhos).mean()
    null = np.array([np.abs([stats.spearmanr(gp, X.loc[k].to_numpy(float))[0]
                             for k in present]).mean()
                     for gp in (rng.permutation(g) for _ in range(N_PERM))])
    return dict(n_genes=len(present), mean_abs_rho=round(obs, 3),
                null_mean=round(float(null.mean()), 3),
                permutation_p=round(float((null >= obs).mean()), 4))


# -------------------------------------------- D) 경쟁적 귀무분포 (competitive null)
def module_competitive(Xf, g, modules, rng, n_draw=5000):
    """
    [C] 의 permutation 은 self-contained null (라벨 셔플) 이라
    '이 모듈이 반응하는가' 만 본다. 모듈이 전부 null 이면
    '검정력 부족' 과 '특이적 무반응' 을 구별할 수 없다.

    여기서는 발현량을 매칭한 무작위 유전자셋을 뽑아 경쟁적 귀무분포를 만든다.
    → '전사체 평균 대비 이 모듈이 유난히 반응/무반응인가' 를 본다.
    """
    rho_all = np.array([stats.spearmanr(g, Xf.loc[k].to_numpy(float))[0] for k in Xf.index])
    rho_all = pd.Series(rho_all, index=Xf.index).dropna()
    expr = Xf.mean(axis=1).reindex(rho_all.index)
    order = expr.sort_values().index
    rank = pd.Series(np.arange(len(order)), index=order)

    rows = []
    for name, genes in modules.items():
        present = [k for k in genes if k in rho_all.index]
        if len(present) < 5:
            continue
        obs = rho_all[present].abs().mean()
        # 발현량 순위가 비슷한 구간에서 표본추출
        idx = rank[present].to_numpy()
        draws = np.empty(n_draw)
        for i in range(n_draw):
            jitter = idx + rng.integers(-250, 251, size=len(idx))
            jitter = np.clip(jitter, 0, len(order) - 1)
            draws[i] = rho_all[order[jitter]].abs().mean()
        pct = float((draws < obs).mean())
        rows.append(dict(module=name, n_genes=len(present),
                         mean_abs_rho=round(float(obs), 3),
                         random_set_mean=round(float(draws.mean()), 3),
                         percentile_vs_random=round(pct, 3),
                         p_two_sided=round(2 * min(pct, 1 - pct), 4)))
    return pd.DataFrame(rows)


def main():
    rng = np.random.default_rng(0)
    ladders, scores_all, modules_all, comp_all = [], [], [], []

    for osd, glds, label in TARGETS:
        print("\n" + "=" * 78)
        print(f"OSD-{osd}  {label}")
        de = pd.read_csv(fetch(osd, rf"{glds}_rna_seq_differential_expression_GLbulkRNAseq\.csv$"),
                         low_memory=False)
        sym = de.drop_duplicates("ENSEMBL").set_index("ENSEMBL")["SYMBOL"]

        # ---- A
        lad = deg_ladder(de, glds, label)
        print("\n[A] 지상대조 대비 DEG 개수")
        print(lad.to_string(index=False))
        rho_a = stats.spearmanr(lad.gravity_g, lad.n_DEG_fdr05)
        print(f"    Spearman(G, nDEG) rho={rho_a[0]:+.2f} p={rho_a[1]:.3f}  (4점이라 검정력 낮음; 단조성이 핵심)")
        ladders.append(lad)

        # ---- 발현행렬
        st = sample_table(osd)
        st["G"] = st["Factor Value[Altered Gravity]"].map(GMAP)
        cnt = pd.read_csv(fetch(osd, rf"{glds}_rna_seq_Normalized_Counts_GLbulkRNAseq\.csv$"), index_col=0)
        X = np.log2(cnt + 1)
        X.index = X.index.map(lambda e: sym.get(e, e))
        X = X[~X.index.duplicated()]

        flt = st[(st["Factor Value[Spaceflight]"] == "Space Flight") & st["G"].notna()]
        flight = [s for s in flt["Sample Name"] if s in X.columns]
        g = flt.set_index("Sample Name").loc[flight, "G"].to_numpy(float)
        ground = [s for s in st[st["Factor Value[Spaceflight]"] == "Ground Control"]["Sample Name"]
                  if s in X.columns]

        # ---- B
        sc = loo_scores(X, flight, ground, g)
        sdf = pd.DataFrame({"sample": list(sc), "score": list(sc.values())})
        sdf["gravity_g"] = [g[flight.index(s)] for s in sdf["sample"]]
        sdf["tissue"] = label
        print("\n[B] 우주비행 반응 점수 (LOO, 값이 클수록 uG 유사)")
        print(sdf.groupby("gravity_g")["score"].agg(["mean", "std", "count"]).round(3).to_string())
        rb = stats.spearmanr(sdf.gravity_g, sdf.score)
        print(f"    Spearman(G, score) rho={rb[0]:+.3f} p={rb[1]:.4f}")
        ag = sdf[sdf.gravity_g > 0]
        rc = stats.spearmanr(ag.gravity_g, ag.score)
        print(f"    인공중력 3군만(0.33/0.66/1G): rho={rc[0]:+.3f} p={rc[1]:.4f}")
        scores_all.append(sdf)

        # ---- C
        print("\n[C] 모듈별 G 용량반응")
        rows = []
        for name, genes in MODULES.items():
            r = module_trend(X[flight], g, genes, rng)
            if r:
                r.update(module=name, tissue=label)
                rows.append(r)
        mdf = pd.DataFrame(rows)[["tissue", "module", "n_genes", "mean_abs_rho", "null_mean", "permutation_p"]]
        print(mdf.to_string(index=False))
        modules_all.append(mdf)

        # ---- D
        Xf = X[flight]
        Xf = Xf[Xf.mean(axis=1) > 1]          # 저발현 유전자 제거
        cdf = module_competitive(Xf, g, MODULES, rng)
        cdf.insert(0, "tissue", label)
        print("\n[D] 경쟁적 귀무분포 (발현량 매칭 무작위 유전자셋 5,000회)")
        print(cdf.to_string(index=False))
        comp_all.append(cdf)

    pd.concat(ladders).to_csv(os.path.join(OUT, "positive_control_deg_ladder.csv"),
                              index=False, encoding="utf-8-sig")
    pd.concat(scores_all).to_csv(os.path.join(OUT, "positive_control_response_scores.csv"),
                                 index=False, encoding="utf-8-sig")
    pd.concat(modules_all).to_csv(os.path.join(OUT, "positive_control_modules.csv"),
                                  index=False, encoding="utf-8-sig")
    pd.concat(comp_all).to_csv(os.path.join(OUT, "positive_control_competitive.csv"),
                               index=False, encoding="utf-8-sig")
    print(f"\n-> {OUT}/positive_control_*.csv")


if __name__ == "__main__":
    main()
