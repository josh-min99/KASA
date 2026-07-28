"""
Stage 3b — 배치 교란을 제거한 매칭 비교.

s3_hlu_apply.py 의 문제
  FLIGHT 8조직은 전부 RR-1 '단일 미션' 인데 HLU 7건은 '서로 다른 7개 실험' 이었다.
  조직쌍 프로파일 상관 0.857 vs 0.000 의 상당 부분이
  '같은 미션 vs 다른 실험' 이라는 배치 구조에서 올 수 있다. 공정한 비교가 아니다.

해결
  양쪽 모두 '단일 실험의 다조직' 세트로 맞춘다.
    FLIGHT-A : RR-1  (OSD-98,99,101,102,103,104,105,168) 8조직   [Stage 1 에서 확보]
    FLIGHT-B : RR-6  (OSD-243,244,245,246,247,248)       6조직   [독립 미션, 재현용]
    HLU-A    : OSD-334,335,337 (+336 혈장)               동일 실험 다조직
  이러면 '단일 실험 내 조직 간 일관성' 이라는 같은 단위로 비교된다.

출력: results/tables/stage3b_*.csv
"""
import os, re, time, warnings
import requests
import numpy as np
import pandas as pd
from scipy import stats

warnings.filterwarnings("ignore")
S = requests.Session()
S.headers.update({"User-Agent": "KASA-drylab/1.0"})

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RAW = os.path.join(ROOT, "data", "raw")
PROC = os.path.join(ROOT, "data", "processed")
TAB = os.path.join(ROOT, "results", "tables")
for d in (RAW, PROC, TAB):
    os.makedirs(d, exist_ok=True)

CLOCK = ["Bmal1", "Arntl", "Clock", "Npas2", "Per1", "Per2", "Per3", "Cry1", "Cry2",
         "Nr1d1", "Nr1d2", "Dbp", "Tef", "Hlf", "Ciart", "Bhlhe40", "Bhlhe41"]

SERIES = {
    "FLIGHT-B (RR-6)": {243: "dorsal skin", 244: "thymus", 245: "liver",
                        246: "spleen", 247: "colon", 248: "lung"},
    "HLU-A":           {334: "heart", 335: "liver", 336: "plasma", 337: "soleus"},
}
FTREAT = re.compile(r"space ?flight", re.I)
FCTRL = re.compile(r"ground control|vivarium", re.I)
HTREAT = re.compile(r"hindlimb unload", re.I)
HCTRL = re.compile(r"normally loaded|loaded control|cage control|ground control|"
                   r"ambulatory|non[- ]?unloaded|sham", re.I)
DIRTY = re.compile(r"reload|radiation|corticosterone|mKO|knockout|\bKO\b|treated|spike|"
                   r"irradiat|gray|\bGy\b", re.I)


def get(url, tries=3, timeout=1800):
    for k in range(tries):
        try:
            r = S.get(url, timeout=timeout)
            if r.status_code == 200:
                return r
        except Exception:
            pass
        time.sleep(2 + 3 * k)
    return None


def de_path(osd):
    d = os.path.join(RAW, f"OSD-{osd}")
    os.makedirs(d, exist_ok=True)
    loc = [f for f in os.listdir(d) if re.search(r"differential_expression.*\.csv$", f, re.I)]
    if loc:
        loc.sort(key=lambda f: os.path.getsize(os.path.join(d, f)))
        return os.path.join(d, loc[0])
    r = get(f"https://osdr.nasa.gov/osdr/data/osd/files/{osd}", timeout=180)
    if r is None:
        return None
    st = r.json()["studies"].get(f"OSD-{osd}")
    if not st:
        return None
    c = [f for f in st["study_files"]
         if re.search(r"differential_expression.*\.csv$", f["file_name"], re.I)]
    if not c:
        return None
    c.sort(key=lambda f: f["file_size"])
    f = c[0]
    if f["file_size"] / 1e6 > 400:
        return None
    rr = get("https://osdr.nasa.gov" + f["remote_url"])
    if rr is None:
        return None
    p = os.path.join(d, f["file_name"])
    with open(p, "wb") as fh:
        fh.write(rr.content)
    return p


def contrasts(hdr):
    return [(c, m.group(1), m.group(2)) for c in hdr
            if (m := re.match(r"Log2fc_\((.*?)\)v\((.*?)\)$", c))]


def choose(cands, tre, cre):
    sc = []
    for c, a, b in cands:
        for t, k, flip in ((a, b, False), (b, a, True)):
            if tre.search(t) and cre.search(k):
                sc.append((len(DIRTY.findall(t)) * 10 + len(DIRTY.findall(k)) * 10
                           + t.count("&") + k.count("&"), c, flip, t, k))
    if not sc:
        return None
    sc.sort()
    return {"col": sc[0][1], "flip": sc[0][2], "treat": sc[0][3], "ctrl": sc[0][4],
            "dirt": sc[0][0]}


def clock_fc(osd, tissue, tre, cre):
    p = de_path(osd)
    if p is None:
        return None, {"OSD_ID": f"OSD-{osd}", "tissue": tissue, "note": "DE 없음/실패"}
    hdr = pd.read_csv(p, nrows=0).columns.tolist()
    sym = next((c for c in hdr if c.upper() == "SYMBOL"), None)
    ch = choose(contrasts(hdr), tre, cre)
    if not (sym and ch):
        return None, {"OSD_ID": f"OSD-{osd}", "tissue": tissue, "note": "적합 대비 없음",
                      "available": "; ".join(c for c, _, _ in contrasts(hdr))[:200]}
    ca = ch["col"].replace("Log2fc_", "Adj.p.value_")
    use = [sym, ch["col"]] + ([ca] if ca in hdr else [])
    df = pd.read_csv(p, usecols=use, low_memory=False)
    df.columns = ["SYMBOL", "log2fc"] + (["fdr"] if ca in hdr else [])
    if "fdr" not in df.columns:
        df["fdr"] = np.nan
    if ch["flip"]:
        df["log2fc"] = -df["log2fc"]
    n_deg = int((df.fdr < 0.05).sum())
    d = df[df.SYMBOL.isin(CLOCK)].copy()
    d["tissue"] = tissue
    d["OSD_ID"] = f"OSD-{osd}"
    return d, {"OSD_ID": f"OSD-{osd}", "tissue": tissue, "treat": ch["treat"][:60],
               "ctrl": ch["ctrl"][:60], "flipped": ch["flip"], "dirt": ch["dirt"],
               "n_DEG_fdr05": n_deg, "note": ""}


def coherence(piv, min_deg_ok):
    """단일 실험 내 조직 간 clock 반응 일관성."""
    out = {}
    for g in ("Bmal1", "Per2"):
        if g in piv.index:
            v = piv.loc[g].dropna()
            if len(v) >= 2:
                out[f"{g}_mean_log2fc"] = round(float(v.mean()), 3)
                out[f"{g}_sign_concord"] = round(float(max((v > 0).sum(), (v < 0).sum()) / len(v)), 3)
                out[f"{g}_n_tissue"] = int(len(v))
    cols = list(piv.columns)
    rs = []
    for i in range(len(cols)):
        for j in range(i + 1, len(cols)):
            a, b = piv[cols[i]], piv[cols[j]]
            ok = a.notna() & b.notna()
            if ok.sum() >= 6:
                rs.append(stats.pearsonr(a[ok], b[ok])[0])
    if rs:
        out["profile_r_median"] = round(float(np.median(rs)), 3)
        out["profile_r_min"] = round(float(np.min(rs)), 3)
        out["profile_r_max"] = round(float(np.max(rs)), 3)
        out["n_tissue_pairs"] = len(rs)
    out["n_datasets_with_DEG"] = min_deg_ok
    return out


def main():
    allmeta, summary, pivots = [], [], {}

    # FLIGHT-A: Stage 1 결과 재사용
    F = pd.read_csv(os.path.join(PROC, "life2020_clock_fc.csv"))
    pivots["FLIGHT-A (RR-1)"] = F.pivot_table(index="SYMBOL", columns="tissue", values="log2fc")

    for name, spec in SERIES.items():
        tre, cre = (HTREAT, HCTRL) if name.startswith("HLU") else (FTREAT, FCTRL)
        frames = []
        for osd, tis in spec.items():
            d, m = clock_fc(osd, tis, tre, cre)
            m["series"] = name
            allmeta.append(m)
            print(f"[{name}] OSD-{osd} {tis}: "
                  f"{'clock %d개, DEG %d' % (len(d), m['n_DEG_fdr05']) if d is not None else m['note']}",
                  flush=True)
            if d is not None:
                frames.append(d)
        if frames:
            A = pd.concat(frames, ignore_index=True)
            pivots[name] = A.pivot_table(index="SYMBOL", columns="tissue", values="log2fc")
            A.to_csv(os.path.join(PROC, f"clock_fc_{name.split()[0]}.csv"),
                     index=False, encoding="utf-8-sig")

    M = pd.DataFrame(allmeta)
    M.to_csv(os.path.join(TAB, "stage3b_contrasts.csv"), index=False, encoding="utf-8-sig")

    print("\n" + "=" * 96)
    print("단일 실험 내 조직 간 clock 반응 일관성 — 같은 단위로 비교")
    print("=" * 96)
    for name, piv in pivots.items():
        ok = int(M[(M.series == name) & (M.n_DEG_fdr05 > 0)].shape[0]) if "series" in M else len(piv.columns)
        c = coherence(piv, ok)
        c["series"] = name
        c["n_tissues"] = piv.shape[1]
        summary.append(c)
        print(f"\n--- {name} ({piv.shape[1]} 조직) ---")
        print(f"  Bmal1 평균 log2FC {c.get('Bmal1_mean_log2fc')}  부호일치 {c.get('Bmal1_sign_concord')}")
        print(f"  Per2  평균 log2FC {c.get('Per2_mean_log2fc')}   부호일치 {c.get('Per2_sign_concord')}")
        print(f"  조직쌍 프로파일 r 중앙값 {c.get('profile_r_median')} "
              f"(범위 {c.get('profile_r_min')}~{c.get('profile_r_max')}, 쌍 {c.get('n_tissue_pairs')})")
        print(piv.round(2).to_string())

    Sm = pd.DataFrame(summary)
    cols = ["series", "n_tissues", "Bmal1_mean_log2fc", "Bmal1_sign_concord",
            "Per2_mean_log2fc", "Per2_sign_concord", "profile_r_median",
            "profile_r_min", "profile_r_max", "n_tissue_pairs"]
    Sm = Sm[[c for c in cols if c in Sm.columns]]
    Sm.to_csv(os.path.join(TAB, "stage3b_series_coherence.csv"),
              index=False, encoding="utf-8-sig")
    print("\n" + "=" * 96)
    print(Sm.to_string(index=False))


if __name__ == "__main__":
    main()
