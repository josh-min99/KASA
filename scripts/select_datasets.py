"""
전수 카탈로그에서 연구에 사용할 데이터셋을 확정한다.

선별 기준을 코드로 명시한다. 카탈로그(catalog_studies.csv / catalog_columns.csv)는
필터 없이 만들어졌으므로, 여기서의 기준 변경만으로 선별을 재검토할 수 있다.

축 정의
  G  중력 조작이 있는가, 몇 단계인가
  T  시간축 — 샘플 간 '하루 중 시각'이 달라지는가 (기록 여부가 아니라 변이 여부)
  R  circadian 판독 가능한가 (전사체 / 시계열 생리·행동)
  N  표본수

출력: data/dataset_selection.csv
"""
import os, re, json, warnings
import pandas as pd

warnings.filterwarnings("ignore")
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, "data")

# --- 축 T: 값이 '하루 중 시각'인 패턴 (명암주기 표기·희석비·날짜는 제외)
CLOCK = re.compile(r"\b\d{1,2}:\d{2}\b|\b\d{1,2}\s*(?:am|pm)\b|\bZT\s*-?\d{1,2}\b"
                   r"|\bCT\s*-?\d{1,2}\b|zeitgeber", re.I)
NOISE = re.compile(r"12\s*:\s*12|1\s*:\s*\d{2,}|light/?dark|light-dark|\d{2}-[A-Za-z]{3}-\d{4}", re.I)
LIGHTCOL = re.compile(r"light cycle", re.I)

# --- 축 G: 중력 조작 표기 (Factor / Parameter / Characteristics 어디에 있든)
GRAV_COL = re.compile(r"gravity|hypergrav|unloading|centrifug|weight bearing|suspension", re.I)
GRAV_VAL = re.compile(r"\b\d\.?\d*\s*g\b|micrograv|hypergrav|hindlimb|unload|"
                      r"centrifug|weight bearing|lunar|martian", re.I)

# --- 축 R: circadian 판독 가능성
OMICS = re.compile(r"transcription profiling|rna-?seq|microarray|proteom|metabolom", re.I)
BEHAV = re.compile(r"behavior|telemetr|temperature|activity|actigraph", re.I)
RODENT = re.compile(r"Mus musculus|Rattus", re.I)


def main():
    sp = os.path.join(OUT, "catalog_studies.csv")
    cp = os.path.join(OUT, "catalog_columns.csv")
    if not (os.path.exists(sp) and os.path.exists(cp)):
        print("카탈로그가 없다. scripts/build_catalog.py 를 먼저 실행할 것.")
        return
    S = pd.read_csv(sp).fillna("")
    C = pd.read_csv(cp).fillna("")
    print(f"카탈로그: 스터디 {len(S)}건, 컬럼 레코드 {len(C)}건")

    # ---------- 축 T
    tmask = C.values.astype(str)
    t_hits = C[C["values"].apply(lambda v: bool(CLOCK.search(str(v)) and not NOISE.search(str(v))))]
    t_vary = t_hits[(t_hits.n_unique > 1) & (~t_hits.column.str.contains(LIGHTCOL))]
    t_const = t_hits[(t_hits.n_unique <= 1) | (t_hits.column.str.contains(LIGHTCOL))]
    T_vary = set(t_vary.OSD_ID); T_const = set(t_const.OSD_ID)
    print(f"\n[T] 시각이 값에 나타나는 스터디: {t_hits.OSD_ID.nunique()}건")
    print(f"    그중 샘플 간 변이 있음: {len(T_vary)}건 -> {sorted(T_vary)}")
    if len(t_vary):
        print(t_vary[["OSD_ID", "column", "n_unique", "values"]].to_string(index=False, max_colwidth=70))

    # ---------- 축 G
    g_rows = C[C.column.str.contains(GRAV_COL) |
               C["values"].apply(lambda v: bool(GRAV_VAL.search(str(v))))]
    glev = {}
    for osd, g in g_rows.groupby("OSD_ID"):
        lv = set()
        for _, r in g.iterrows():
            if not GRAV_COL.search(str(r.column)):
                continue
            for v in str(r["values"]).split(" | "):
                if GRAV_VAL.search(v):
                    lv.add(v.strip()[:44])
        if lv:
            glev[osd] = sorted(lv)
    print(f"\n[G] 중력 조작 표기가 있는 스터디: {len(glev)}건")

    # ---------- 통합
    rows = []
    for _, s in S.iterrows():
        osd = s.OSD_ID
        lv = glev.get(osd, [])
        rows.append({
            "OSD_ID": osd,
            "organism": s.organism,
            "rodent": bool(RODENT.search(str(s.organism))),
            "material": str(s.material)[:70],
            "n_samples": s.n_samples,
            "has_counts": s.has_geneLab_counts,
            "omics": bool(OMICS.search(str(s.assays) + str(s.factors))),
            "behavior_physio": bool(BEHAV.search(str(s.assays))),
            "G_levels": "; ".join(lv)[:150],
            "n_G_levels": len(lv),
            "T_varies": osd in T_vary,
            "T_recorded_const": osd in T_const,
        })
    D = pd.DataFrame(rows)

    # ---------- 선별
    D["tier"] = ""
    D.loc[D.rodent & (D.n_G_levels >= 3) & D.has_counts, "tier"] = "1_dose_response"
    D.loc[(D.tier == "") & D.rodent & D.T_varies & D.has_counts, "tier"] = "2_phase_controlled"
    D.loc[(D.tier == "") & D.rodent & D.behavior_physio & (D.n_G_levels >= 1), "tier"] = "3_timeseries"
    D.loc[(D.tier == "") & D.rodent & (D.n_G_levels >= 2) & D.has_counts, "tier"] = "4_binary_gravity"

    sel = D[D.tier != ""].sort_values(["tier", "OSD_ID"])
    D.to_csv(os.path.join(OUT, "dataset_selection.csv"), index=False, encoding="utf-8-sig")

    print("\n" + "=" * 100)
    print("선별 결과")
    for t, g in sel.groupby("tier"):
        print(f"\n--- {t}  ({len(g)}건) ---")
        print(g[["OSD_ID", "organism", "material", "n_samples",
                 "n_G_levels", "G_levels", "T_varies"]].to_string(index=False, max_colwidth=52))
    print(f"\n-> {OUT}/dataset_selection.csv")


if __name__ == "__main__":
    main()
