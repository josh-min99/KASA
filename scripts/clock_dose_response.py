"""
예비분석 1: 인공중력 용량(G)에 대한 clock gene 발현의 단조 반응 검정.

대상: 비행 중 중력 단계가 3개 이상인 스터디 (OSD-758 / 759 / 714)
검정: 유전자별 Spearman(G, log2 expr) + clock 모듈 전체 permutation test
부가: uG vs 1G 최소검출효과크기(MDE) 산출 → 웻랩 power analysis 의 분산 사전값

주의사항 2가지
  1) MGI 공식 심볼 변경: Arntl -> Bmal1. 'Arntl' 로 조회하면 조용히 NaN 이 된다.
  2) 중력 조건은 ISA.zip 의 s_*.txt 에서만 읽을 것 (meta JSON API 는 null 반환).
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

TARGETS = [(758, "GLDS-664", "Retina"), (759, "GLDS-665", "Optic nerve"), (714, "GLDS-638", "Soleus")]
CLOCK = ["Bmal1", "Clock", "Npas2", "Per1", "Per2", "Per3", "Cry1", "Cry2", "Nr1d1",
         "Nr1d2", "Dbp", "Tef", "Hlf", "Ciart", "Bhlhe40", "Bhlhe41", "Rorb"]
GMAP = {"uG": 0.0, "1/6G with centrifugation": 1 / 6, "0.33G by centrifugation": 0.33,
        "0.66G by centrifugation": 0.66, "1G by centrifugation": 1.0, "1G with centrifugation": 1.0}
N_PERM = 2000


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


def ensembl_to_symbol():
    p = fetch(758, r"GLDS-664_rna_seq_differential_expression_GLbulkRNAseq\.csv$")
    m = pd.read_csv(p, usecols=["ENSEMBL", "SYMBOL"])
    return m.drop_duplicates("ENSEMBL").set_index("ENSEMBL")["SYMBOL"]


def load(osd, glds, symbol_map):
    st = sample_table(osd)
    st["G"] = st["Factor Value[Altered Gravity]"].map(GMAP)
    flt = st[(st["Factor Value[Spaceflight]"] == "Space Flight") & st["G"].notna()]
    cnt = pd.read_csv(fetch(osd, rf"{glds}_rna_seq_Normalized_Counts_GLbulkRNAseq\.csv$"), index_col=0)
    ids = [s for s in flt["Sample Name"] if s in cnt.columns]
    X = np.log2(cnt[ids] + 1)
    X.index = X.index.map(lambda g: symbol_map.get(g, g))
    X = X[~X.index.duplicated()]
    return X, flt.set_index("Sample Name").loc[ids, "G"].to_numpy(float)


def mde(X, g, panel, alpha=0.05, power=0.80):
    """uG vs 1G 두 군 t검정 기준 최소검출효과크기 (log2FC)."""
    n0, n1 = int((g == 0).sum()), int((g == 1.0).sum())
    sds = [np.std(X.loc[k].to_numpy(float)[g == lv], ddof=1) for k in panel for lv in (0.0, 1.0)]
    sd = float(np.mean(sds))
    dfree = n0 + n1 - 2
    k = np.sqrt(1 / n0 + 1 / n1)
    unc = (stats.t.ppf(1 - alpha / 2, dfree) + stats.t.ppf(power, dfree)) * k * sd
    bon = (stats.t.ppf(1 - alpha / (2 * len(panel)), dfree) + stats.t.ppf(power, dfree)) * k * sd
    return n0, n1, sd, unc, bon


def main():
    sym = ensembl_to_symbol()
    panel_all = [c for c in CLOCK if c in set(sym.dropna())]
    missing = [c for c in CLOCK if c not in panel_all]
    if missing:
        print(f"[warn] 어노테이션에 없는 심볼: {missing}")
    rng = np.random.default_rng(0)
    records = []

    for osd, glds, label in TARGETS:
        X, g = load(osd, glds, sym)
        panel = [k for k in panel_all if k in X.index and np.std(X.loc[k]) > 0]
        print("\n" + "=" * 78)
        print(f"OSD-{osd}  {label}   n_flight={len(g)}   "
              f"levels={ {round(float(k), 2): int((g == k).sum()) for k in sorted(set(g))} }")

        tab = pd.DataFrame(
            [dict(zip(("gene", "rho", "p"), (k, *stats.spearmanr(g, X.loc[k].to_numpy(float)))))
             for k in panel]).set_index("gene")
        print(tab.round(3).to_string())
        print(f"  p<0.05: {int((tab.p < 0.05).sum())}/{len(tab)}  (우연 기대 {0.05 * len(tab):.1f})")

        obs = tab.rho.abs().mean()
        null = np.array([np.mean([abs(stats.spearmanr(gp, X.loc[k].to_numpy(float))[0]) for k in panel])
                         for gp in (rng.permutation(g) for _ in range(N_PERM))])
        p_mod = float((null >= obs).mean())
        print(f"  모듈 mean|rho|={obs:.3f}  null={null.mean():.3f}  permutation p={p_mod:.3f}")

        n0, n1, sd, unc, bon = mde(X, g, panel)
        print(f"  uG(n={n0}) vs 1G(n={n1})  within-group SD={sd:.3f} log2")
        print(f"  MDE(80% power): 단일검정 {unc:.2f} log2FC ({2 ** unc:.2f}x) / "
              f"Bonferroni {bon:.2f} log2FC ({2 ** bon:.2f}x)")

        records.append(dict(study=f"OSD-{osd}", tissue=label, n_flight=len(g),
                            n_sig_p05=int((tab.p < 0.05).sum()), n_genes=len(tab),
                            module_mean_abs_rho=round(obs, 3), permutation_p=p_mod,
                            within_group_sd_log2=round(sd, 3),
                            mde_log2fc=round(unc, 2), mde_log2fc_bonferroni=round(bon, 2)))
        tab.round(4).to_csv(os.path.join(OUT, f"clock_rho_OSD-{osd}.csv"), encoding="utf-8-sig")

    summary = pd.DataFrame(records)
    summary.to_csv(os.path.join(OUT, "preliminary_summary.csv"), index=False, encoding="utf-8-sig")
    print("\n" + "=" * 78)
    print(summary.to_string(index=False))


if __name__ == "__main__":
    main()
