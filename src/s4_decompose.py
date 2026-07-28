"""
Stage 4 — 분해와 검정력.

Stage 3c 에서 배치 교란을 제거한 뒤 FLIGHT 와 HLU 의 clock 반응 크기가 크게 달랐다.
그러나 HLU 실험은 전반적으로 DEG 가 적다(7~196 vs RR-1 242~5002).
따라서 두 가지를 구별해야 한다.
  (a) clock 이 '특이적으로' 반응하지 않는다
  (b) 그 실험 전체의 신호가 약해 clock 도 덩달아 작다

검정 1 — 전사체 전역 대비 정규화
  각 데이터셋에서 clock 유전자의 |log2FC| 중앙값을, 같은 데이터셋의
  전체 유전자 |log2FC| 중앙값으로 나눈다(clock enrichment ratio).
  (a) 라면 FLIGHT 는 비율이 크고 HLU 는 작다.
  (b) 라면 두 비율이 비슷하다.

검정 2 — 검정력
  '중력만으로 RR-1 수준의 Bmal1 반응(+1.105 log2FC)이 일어났다면
   HLU 데이터에서 검출됐을 것인가' 를 계산한다.
  GeneLab DE 테이블에는 개체별 값이 없으므로, 같은 데이터셋의
  전체 유전자 log2FC 분포의 표준편차를 잡음 대리값으로 쓴다(보수적).
  검출 불가면 '판정 불가 + 필요한 최소 n' 을 산출한다.

주의: 군당 n=3~5 다. p 값 단독 주장 금지. 효과크기와 구간을 함께 낸다.

출력: results/tables/stage4_*.csv
"""
import os, re, warnings
import numpy as np
import pandas as pd
from scipy import stats

warnings.filterwarnings("ignore")
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RAW = os.path.join(ROOT, "data", "raw")
PROC = os.path.join(ROOT, "data", "processed")
TAB = os.path.join(ROOT, "results", "tables")
os.makedirs(TAB, exist_ok=True)

CLOCK = ["Bmal1", "Arntl", "Clock", "Npas2", "Per1", "Per2", "Per3", "Cry1", "Cry2",
         "Nr1d1", "Nr1d2", "Dbp", "Tef", "Hlf", "Ciart", "Bhlhe40", "Bhlhe41"]
DE_RE = re.compile(r"differential_expression.*\.csv$", re.I)

# (osd, 조직, 시리즈, 처리정규식, 대조정규식)
SETS = [
    (98,  "adrenal gland",     "FLIGHT-A", r"space ?flight", r"ground control"),
    (99,  "EDL",               "FLIGHT-A", r"space ?flight", r"ground control"),
    (101, "gastrocnemius",     "FLIGHT-A", r"space ?flight", r"ground control"),
    (102, "kidney",            "FLIGHT-A", r"space ?flight", r"ground control"),
    (103, "quadriceps",        "FLIGHT-A", r"space ?flight", r"ground control"),
    (104, "soleus",            "FLIGHT-A", r"space ?flight", r"ground control"),
    (105, "tibialis anterior", "FLIGHT-A", r"space ?flight", r"ground control"),
    (168, "liver",             "FLIGHT-A", r"space ?flight", r"ground control"),
    (202, "brain",             "HLU-A",    r"hindlimb unloaded", r"normally loaded"),
    (203, "retina",            "HLU-A",    r"hindlimb unloaded", r"normally loaded"),
    (211, "spleen",            "HLU-A",    r"hindlimb unloaded", r"normally loaded"),
    (237, "dorsal skin",       "HLU-A",    r"hindlimb unloaded", r"normally loaded"),
]
IRRADIATED = re.compile(r"gamma|cobalt|(?<!non-)irradiated", re.I)


def de_path(osd):
    d = os.path.join(RAW, f"OSD-{osd}")
    if not os.path.isdir(d):
        return None
    c = [f for f in os.listdir(d) if DE_RE.search(f)]
    if not c:
        return None
    c.sort(key=lambda f: os.path.getsize(os.path.join(d, f)))
    return os.path.join(d, c[0])


def pick(hdr, tre, cre, avoid_irr):
    best = None
    for c in hdr:
        m = re.match(r"Log2fc_\((.*?)\)v\((.*?)\)$", c)
        if not m:
            continue
        a, b = m.group(1), m.group(2)
        for t, k, flip in ((a, b, False), (b, a, True)):
            if not (re.search(tre, t, re.I) and re.search(cre, k, re.I)):
                continue
            if avoid_irr and (IRRADIATED.search(t) or IRRADIATED.search(k)):
                continue
            sc = t.count("&") + k.count("&")
            if best is None or sc < best[0]:
                best = (sc, c, flip)
    return (best[1], best[2]) if best else (None, False)


def load(osd, tre, cre, avoid_irr):
    p = de_path(osd)
    if p is None:
        return None
    hdr = pd.read_csv(p, nrows=0).columns.tolist()
    sym = next((c for c in hdr if c.upper() == "SYMBOL"), None)
    col, flip = pick(hdr, tre, cre, avoid_irr)
    if not (sym and col):
        return None
    ca = col.replace("Log2fc_", "Adj.p.value_")
    use = [sym, col] + ([ca] if ca in hdr else [])
    df = pd.read_csv(p, usecols=use, low_memory=False)
    df.columns = ["SYMBOL", "log2fc"] + (["fdr"] if ca in hdr else [])
    if "fdr" not in df.columns:
        df["fdr"] = np.nan
    if flip:
        df["log2fc"] = -df["log2fc"]
    return df.dropna(subset=["log2fc"])


def main():
    rows = []
    for osd, tis, ser, tre, cre in SETS:
        df = load(osd, tre, cre, avoid_irr=(ser == "HLU-A"))
        if df is None:
            rows.append({"OSD_ID": f"OSD-{osd}", "tissue": tis, "series": ser,
                         "note": "대비 실패"})
            print(f"[OSD-{osd}] {tis}: 대비 실패", flush=True)
            continue
        allabs = df.log2fc.abs()
        ck = df[df.SYMBOL.isin(CLOCK)]
        ckabs = ck.log2fc.abs()
        sd_all = float(df.log2fc.std())
        rows.append({
            "OSD_ID": f"OSD-{osd}", "tissue": tis, "series": ser,
            "n_genes": len(df),
            "n_DEG_fdr05": int((df.fdr < 0.05).sum()),
            "all_absfc_median": round(float(allabs.median()), 4),
            "clock_absfc_median": round(float(ckabs.median()), 4),
            "clock_enrichment": round(float(ckabs.median() / max(allabs.median(), 1e-9)), 3),
            "clock_pctile": round(float((allabs < ckabs.median()).mean()), 4),
            "sd_log2fc_all": round(sd_all, 4),
            "Bmal1_log2fc": round(float(ck[ck.SYMBOL == "Bmal1"].log2fc.iloc[0]), 3)
                            if (ck.SYMBOL == "Bmal1").any() else np.nan,
        })
        print(f"[OSD-{osd}] {tis} ({ser}): clock/전역 |log2FC| 비 "
              f"{rows[-1]['clock_enrichment']}  DEG {rows[-1]['n_DEG_fdr05']}", flush=True)

    D = pd.DataFrame(rows)
    D.to_csv(os.path.join(TAB, "stage4_normalized.csv"), index=False, encoding="utf-8-sig")

    print("\n" + "=" * 96)
    print("검정 1 — clock 반응을 전사체 전역 반응으로 정규화")
    print("=" * 96)
    print(D.to_string(index=False))
    ok = D.dropna(subset=["clock_enrichment"])
    print("\n시리즈별 요약 (clock |log2FC| 중앙값 / 전역 |log2FC| 중앙값):")
    g = ok.groupby("series").clock_enrichment.agg(["median", "min", "max", "count"])
    print(g.round(3).to_string())
    print("\n시리즈별 clock 백분위 (전역 분포에서 clock 중앙값의 위치):")
    print(ok.groupby("series").clock_pctile.agg(["median", "min", "max"]).round(3).to_string())

    if ok.series.nunique() == 2:
        a = ok[ok.series == "FLIGHT-A"].clock_enrichment
        b = ok[ok.series == "HLU-A"].clock_enrichment
        u = stats.mannwhitneyu(a, b, alternative="two-sided")
        d = float(a.median() - b.median())
        print(f"\n  FLIGHT-A 중앙값 {a.median():.3f} vs HLU-A 중앙값 {b.median():.3f}"
              f"  (차이 {d:+.3f})")
        print(f"  Mann-Whitney U p={u.pvalue:.4f}  [n={len(a)} vs {len(b)} — 참고치]")
        print("  * n 이 작아 p 값 단독으로 주장하지 않는다. 효과크기와 개별값을 함께 본다.")

    # ---------- 검정 2: 검정력
    print("\n" + "=" * 96)
    print("검정 2 — HLU 데이터에서 RR-1 수준의 Bmal1 반응을 검출할 수 있었는가")
    print("=" * 96)
    target = 1.105          # RR-1 Bmal1 평균 log2FC
    pw = []
    for _, r in ok[ok.series == "HLU-A"].iterrows():
        sd = r.sd_log2fc_all
        obs = r.Bmal1_log2fc
        # 전역 log2FC 분포에서 target 이 차지하는 백분위 = 검출 가능성의 대리 지표
        pctl = None
        p = de_path(int(r.OSD_ID.split("-")[1]))
        if p:
            hdr = pd.read_csv(p, nrows=0).columns.tolist()
            sym = next((c for c in hdr if c.upper() == "SYMBOL"), None)
            col, flip = pick(hdr, r"hindlimb unloaded", r"normally loaded", True)
            if sym and col:
                dd = pd.read_csv(p, usecols=[sym, col], low_memory=False).dropna()
                v = dd[col].abs()
                pctl = float((v < target).mean())
        pw.append({"OSD_ID": r.OSD_ID, "tissue": r.tissue,
                   "sd_log2fc_all": sd, "Bmal1_observed": obs,
                   "target_RR1_Bmal1": target,
                   "target_z_vs_sd": round(target / sd, 2) if sd else np.nan,
                   "target_pctile_in_dataset": round(pctl, 4) if pctl is not None else np.nan})
    P = pd.DataFrame(pw)
    P.to_csv(os.path.join(TAB, "stage4_power.csv"), index=False, encoding="utf-8-sig")
    print(P.to_string(index=False))
    print("\n  target_z_vs_sd = RR-1 수준 효과(1.105) / 그 데이터셋 전역 log2FC 표준편차")
    print("  target_pctile  = 그 데이터셋에서 |log2FC| 1.105 이상인 유전자의 비율(1-값)")
    if len(P):
        med_z = P.target_z_vs_sd.median()
        print(f"\n  HLU 데이터셋의 중앙 z = {med_z:.2f}")
        if med_z >= 2:
            print("  >>> RR-1 수준의 효과였다면 HLU 데이터에서 눈에 띄었을 것이다.")
            print("      따라서 '중력만으로는 RR-1 수준 반응이 일어나지 않는다' 로 읽을 수 있다.")
        else:
            print("  >>> RR-1 수준의 효과여도 이 잡음 수준에서는 묻힌다. **판정 불가.**")


if __name__ == "__main__":
    main()
