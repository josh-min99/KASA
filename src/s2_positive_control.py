"""
Stage 2 — 양성대조.

목적: Stage 1·3 과 동일한 파이프라인(GeneLab DE 테이블의 log2FC/FDR)이
      이 데이터에서 '알려진 확실한 효과'를 잡는지 확인한다.
      못 잡으면 이후 모든 '신호 없음' 은 해석 불가다.

대상 효과: 폐용성 근위축의 표준 마커
  Fbxo32 (atrogin-1 / MAFbx), Trim63 (MuRF1)
  → 근육 폐용(HLU, 우주비행) 에서 강하게 상향된다는 것이 확립돼 있다.
    보조: Foxo1/Foxo3(상류), Mstn(myostatin), Ampd3

비교
  (a) FLIGHT  : RR-1 근육 4조직 (OSD-99 EDL, 101 gastroc, 103 quad, 104 soleus, 105 TA)
  (b) HLU     : OSD-876/880 (비복근), 935/949 (가자미근)
  (c) 앵커     : OSD-21 (마이크로어레이, HLU / HLU+Reloaded / Normally Loaded)

판정 (사전 규정)
  Fbxo32 또는 Trim63 이 근육 폐용 조건에서 log2FC > 0 이고 FDR < 0.05 인 데이터셋이
  절반 이상이면 통과. 아니면 실패로 기록하고 FINAL_REPORT 맨 앞에 쓴다.

출력: results/tables/stage2_positive_control.csv
"""
import os, re, glob, warnings
import numpy as np
import pandas as pd

warnings.filterwarnings("ignore")
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RAW = os.path.join(ROOT, "data", "raw")
TAB = os.path.join(ROOT, "results", "tables")
os.makedirs(TAB, exist_ok=True)

ATRO = ["Fbxo32", "Trim63", "Foxo1", "Foxo3", "Mstn", "Ampd3"]

# (osd, 조직, 패러다임)
SETS = [
    (99,  "EDL",               "FLIGHT"),
    (101, "gastrocnemius",     "FLIGHT"),
    (103, "quadriceps",        "FLIGHT"),
    (104, "soleus",            "FLIGHT"),
    (105, "tibialis anterior", "FLIGHT"),
    (876, "gastrocnemius",     "HLU"),
    (880, "gastrocnemius",     "HLU"),
    (935, "soleus",            "HLU"),
    (949, "soleus",            "HLU"),
]

# 각 패러다임에서 '부하 감소' 방향이 되도록 대비를 고르고 부호를 맞춘다
TREAT = re.compile(r"space ?flight|hindlimb unload|unloaded|suspend", re.I)
CTRL = re.compile(r"ground control|normally loaded|loaded control|vivarium|"
                  r"cage control|ambulatory|sham", re.I)


def de_path(osd):
    d = os.path.join(RAW, f"OSD-{osd}")
    if not os.path.isdir(d):
        return None
    c = [f for f in os.listdir(d) if re.search(r"differential_expression.*\.csv$", f, re.I)]
    if not c:
        return None
    c.sort(key=lambda f: os.path.getsize(os.path.join(d, f)))
    return os.path.join(d, c[0])


def pick_contrast(cols):
    """(처리)v(대조) 형태를 찾아 (컬럼, 부호반전여부) 반환."""
    best = None
    for c in cols:
        m = re.match(r"Log2fc_\((.*?)\)v\((.*?)\)$", c)
        if not m:
            continue
        a, b = m.group(1), m.group(2)
        if TREAT.search(a) and CTRL.search(b):
            return c, False
        if TREAT.search(b) and CTRL.search(a):
            best = best or (c, True)
    return best if best else ((cols[0], False) if cols else (None, False))


def run(osd, tissue, para):
    p = de_path(osd)
    if p is None:
        return [{"OSD_ID": f"OSD-{osd}", "tissue": tissue, "paradigm": para,
                 "gene": g, "log2fc": np.nan, "fdr": np.nan,
                 "note": "DE 파일 없음"} for g in ATRO]
    hdr = pd.read_csv(p, nrows=0).columns.tolist()
    sym = next((c for c in hdr if c.upper() == "SYMBOL"), None)
    lfc = [c for c in hdr if c.startswith("Log2fc_")]
    if not (sym and lfc):
        return [{"OSD_ID": f"OSD-{osd}", "tissue": tissue, "paradigm": para,
                 "gene": g, "log2fc": np.nan, "fdr": np.nan,
                 "note": "컬럼 불일치"} for g in ATRO]
    cl, flip = pick_contrast(lfc)
    ca = cl.replace("Log2fc_", "Adj.p.value_")
    use = [sym, cl] + ([ca] if ca in hdr else [])
    df = pd.read_csv(p, usecols=use, low_memory=False)
    df.columns = ["SYMBOL", "log2fc"] + (["fdr"] if ca in hdr else [])
    if "fdr" not in df.columns:
        df["fdr"] = np.nan
    if flip:
        df["log2fc"] = -df["log2fc"]
    out = []
    for g in ATRO:
        r = df[df.SYMBOL == g]
        out.append({"OSD_ID": f"OSD-{osd}", "tissue": tissue, "paradigm": para, "gene": g,
                    "log2fc": float(r.log2fc.iloc[0]) if len(r) else np.nan,
                    "fdr": float(r.fdr.iloc[0]) if len(r) else np.nan,
                    "contrast": cl[:80], "flipped": flip,
                    "note": "" if len(r) else "유전자 없음"})
    return out


def main():
    rows = []
    for osd, tis, para in SETS:
        rows += run(osd, tis, para)
    D = pd.DataFrame(rows)
    D.to_csv(os.path.join(TAB, "stage2_positive_control.csv"),
             index=False, encoding="utf-8-sig")

    print("=== 근위축 마커 log2FC (부하감소 방향 = 양수 기대) ===")
    piv = D.pivot_table(index="gene", columns=["paradigm", "OSD_ID"], values="log2fc")
    print(piv.round(2).to_string())
    print("\n=== FDR ===")
    pf = D.pivot_table(index="gene", columns=["paradigm", "OSD_ID"], values="fdr")
    print(pf.round(4).to_string())

    print("\n대비 컬럼 확인:")
    for _, r in D.drop_duplicates("OSD_ID").iterrows():
        print(f"  {r.OSD_ID:<9} {r.paradigm:<7} {str(r.contrast)[:74]}"
              f"{'  [부호반전]' if r.flipped else ''}")

    # ---------- 판정
    key = D[D.gene.isin(["Fbxo32", "Trim63"])].dropna(subset=["log2fc"])
    hit = key[(key.log2fc > 0) & (key.fdr < 0.05)]
    n_set = key.OSD_ID.nunique()
    n_hit = hit.OSD_ID.nunique()
    print("\n" + "=" * 88)
    print(f"판정: Fbxo32/Trim63 이 상향+FDR<0.05 인 데이터셋 {n_hit}/{n_set}")
    verdict = "통과" if n_set and n_hit / n_set >= 0.5 else "실패"
    print(f"  >>> Stage 2 {verdict}  (기준: 절반 이상)")
    for para, g in key.groupby("paradigm"):
        h = g[(g.log2fc > 0) & (g.fdr < 0.05)]
        print(f"    {para}: {h.OSD_ID.nunique()}/{g.OSD_ID.nunique()} 데이터셋")
    pd.DataFrame([{"n_datasets": n_set, "n_hit": n_hit, "verdict": verdict}]).to_csv(
        os.path.join(TAB, "stage2_verdict.csv"), index=False, encoding="utf-8-sig")


if __name__ == "__main__":
    main()
