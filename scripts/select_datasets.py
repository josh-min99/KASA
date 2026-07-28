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
import requests

# main() 안에서 S 를 DataFrame 으로 재사용하므로 세션은 다른 이름을 쓴다
SESSION = requests.Session()
SESSION.headers.update({"User-Agent": "KASA-circadian-gravity/1.0"})

# 발현 행렬 파일 패턴. RNA-seq counts 뿐 아니라 마이크로어레이도 포함해야 한다.
# counts 만 찾으면 구형 array 스터디(예: OSD-63 의 쥐 2G 과중력)를
# '표현형 전용' 으로 잘못 분류한다.
EXPR_FILE = re.compile(
    r"_rna_seq_Normalized_Counts_GLbulkRNAseq\.csv$"
    r"|_array_normalized_expression_probeset_GLmicroarray\.csv$"
    r"|_array_normalized_intensities_probe_GLmicroarray\.csv$", re.I)


def has_expression(osd_id):
    """files API 로 발현 행렬 보유 여부를 확인 (카탈로그가 RNA-seq 만 기록하므로 보강)."""
    n = osd_id.split("-")[1]
    try:
        j = SESSION.get(f"https://osdr.nasa.gov/osdr/data/osd/files/{n}", timeout=90).json()
        st = j.get("studies", {}).get(osd_id)
        if not st:
            return False
        return any(EXPR_FILE.search(f["file_name"]) for f in st["study_files"])
    except Exception:
        return False

warnings.filterwarnings("ignore")
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, "data")

# --- 축 T: 값이 '하루 중 시각'인 패턴 (명암주기 표기·희석비·날짜는 제외)
CLOCK = re.compile(r"\b\d{1,2}:\d{2}\b|\b\d{1,2}\s*(?:am|pm)\b|\bZT\s*-?\d{1,2}\b"
                   r"|\bCT\s*-?\d{1,2}\b|zeitgeber", re.I)
# 주의: 오탐 패턴에 반드시 \b 를 붙일 것.
#   경계 없는 `1\s*:\s*\d{2,}` 은 희석비 1:100 뿐 아니라 시각 11:18 의 '1:18' 에도 매칭돼
#   실제 시각 기록(OSD-612)을 통째로 걸러낸다. 초기 구현에서 실제로 발생했다.
NOISE = re.compile(r"\b12\s*:\s*12\b|\b1\s*:\s*\d{2,}\b|light/?dark|light-dark"
                   r"|\d{2}-[A-Za-z]{3}-\d{4}", re.I)
LIGHTCOL = re.compile(r"light cycle", re.I)


def has_clock_value(values_field):
    """값을 이어붙인 문자열이 아니라 개별 값 단위로 검사한다.
    이어붙여 검사하면 한 값의 오탐이 나머지 전체를 오염시킨다."""
    for v in str(values_field).split(" | "):
        if CLOCK.search(v) and not NOISE.search(v):
            return True
    return False


# --- 축 G: 중력 조작 표기 (Factor / Parameter / Characteristics 어디에 있든)
#   Spaceflight 도 중력 조건이다 (궤도 μg vs 지상 1g). 초기 구현이 이를 누락했다.
GRAV_COL = re.compile(r"gravity|hypergrav|unloading|centrifug|weight bearing|suspension"
                      r"|spaceflight|space flight", re.I)
GRAV_VAL = re.compile(r"\b\d\.?\d*\s*g\b|micrograv|hypergrav|hindlimb|unload|"
                      r"centrifug|weight bearing|lunar|martian|space ?flight|"
                      r"ground control|vivarium", re.I)

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
    t_hits = C[C["values"].apply(has_clock_value)]
    t_vary = t_hits[(t_hits.n_unique > 1) & (~t_hits.column.str.contains(LIGHTCOL))]
    t_const = t_hits[(t_hits.n_unique <= 1) | (t_hits.column.str.contains(LIGHTCOL))]
    T_vary = set(t_vary.OSD_ID); T_const = set(t_const.OSD_ID)
    print(f"\n[T] 시각이 값에 나타나는 스터디: {t_hits.OSD_ID.nunique()}건")
    print(f"    그중 샘플 간 변이 있음: {len(T_vary)}건 -> {sorted(T_vary)}")
    if len(t_vary):
        print(t_vary[["OSD_ID", "column", "n_unique", "values"]].to_string(index=False, max_colwidth=70))

    # ---------- 축 G
    # 문자열이 아니라 '실제 중력 크기(g)' 로 환산해 센다.
    # Ground Control 과 Vivarium Control 은 표기가 달라도 둘 다 1 g 이므로
    # 서로 다른 중력 단계가 아니다. 문자열을 그대로 세면 3단계로 잘못 집계된다.
    def to_g(v):
        s = str(v).strip().lower()
        if re.search(r"not applicable|^nan$|^$", s):
            return None
        if re.search(r"\bug\b|microgravity|space ?flight|hindlimb unload|"
                     r"0g |^0g|0% partial|antiorthostatic|suspension", s):
            return 0.0
        m = re.search(r"(\d+(?:\.\d+)?)\s*/\s*(\d+)\s*g", s)      # 1/6G
        if m:
            return round(float(m.group(1)) / float(m.group(2)), 3)
        m = re.search(r"(\d+(?:\.\d+)?)\s*g\b", s)                 # 0.33G, 1G
        if m:
            return round(float(m.group(1)), 3)
        # 주의: 'mars' 를 화성으로 읽으면 안 된다. JAXA 의 인공중력 장비 이름이
        # MARS(Multiple Artificial-gravity Research System) 라서, 장비명이 적힌
        # 자유 텍스트가 0.38 g 로 오인된다(OSD-714 에서 실제로 발생).
        if re.search(r"\blunar\b", s):
            return 0.167
        if re.search(r"\bmartian\b", s):
            return 0.38
        if re.search(r"ground control|vivarium|loaded control|sham|reload|"
                     r"basal|1 ?g on earth|habitat ground", s):
            return 1.0
        return None

    # 중력 '단계' 는 Factor Value 컬럼에서만 센다.
    # Parameter Value 의 자유 텍스트(장비명·반경 등)까지 넣으면 오탐이 섞인다.
    g_rows = C[C.column.str.startswith("Factor Value") & C.column.str.contains(GRAV_COL, na=False)]
    glev, glab = {}, {}
    for osd, g in g_rows.groupby("OSD_ID"):
        vals, labs = set(), set()
        for _, r in g.iterrows():
            for v in str(r["values"]).split(" | "):
                gv = to_g(v)
                if gv is not None:
                    vals.add(gv); labs.add(v.strip()[:40])
        if vals:
            glev[osd] = sorted(vals); glab[osd] = sorted(labs)
    print(f"\n[G] 중력 크기로 환산된 스터디: {len(glev)}건")
    multi = {k: v for k, v in glev.items() if len(v) >= 3}
    print(f"    3단계 이상: {len(multi)}건")

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
            "G_values": "; ".join(map(str, lv)),
            "G_labels": "; ".join(glab.get(osd, []))[:150],
            "n_G_levels": len(lv),
            "T_varies": osd in T_vary,
            "T_recorded_const": osd in T_const,
        })
    D = pd.DataFrame(rows)

    # ---------- 발현 행렬 보유 여부 보강 (마이크로어레이 포함)
    cand = D[D.rodent & ((D.n_G_levels >= 2) | D.T_varies)].OSD_ID.tolist()
    print(f"\n[보강] 발현 행렬 확인 대상 {len(cand)}건 (마이크로어레이 포함)")
    expr = {}
    for i, osd in enumerate(cand, 1):
        expr[osd] = has_expression(osd)
        if i % 40 == 0:
            print(f"    {i}/{len(cand)}")
    D["has_expr"] = D.OSD_ID.map(expr).fillna(D.has_counts).astype(bool)
    added = [k for k, v in expr.items()
             if v and not bool(D.set_index("OSD_ID").loc[k, "has_counts"])]
    if added:
        print(f"    counts 로는 놓쳤다가 array 로 확인된 스터디: {added}")

    # ---------- 선별
    D["tier"] = ""
    # 1: 비행 중 중력이 3단계 이상 + 전사체 -> 용량반응 검정 가능
    D.loc[D.rodent & (D.n_G_levels >= 3) & D.has_expr, "tier"] = "1_dose_response"
    # 2: 샘플별 시각 변이 + 전사체 -> 위상을 공변량으로 통제 가능
    D.loc[(D.tier == "") & D.rodent & D.T_varies & D.has_expr, "tier"] = "2_phase_controlled"
    # 3: 시계열 생리·행동 판독
    D.loc[(D.tier == "") & D.rodent & D.behavior_physio, "tier"] = "3_timeseries"
    # 4: 중력 2단계 + 전사체
    D.loc[(D.tier == "") & D.rodent & (D.n_G_levels >= 2) & D.has_expr, "tier"] = "4_binary_gravity"
    # 5: 중력 3단계 이상이나 전사체 없음 (표현형 전용)
    D.loc[(D.tier == "") & D.rodent & (D.n_G_levels >= 3), "tier"] = "5_phenotype_only"

    sel = D[D.tier != ""].sort_values(["tier", "OSD_ID"])
    D.to_csv(os.path.join(OUT, "dataset_selection.csv"), index=False, encoding="utf-8-sig")

    print("\n" + "=" * 100)
    print("선별 결과")
    for t, g in sel.groupby("tier"):
        print(f"\n--- {t}  ({len(g)}건) ---")
        print(g[["OSD_ID", "material", "n_samples", "has_expr",
                 "n_G_levels", "G_values", "T_varies"]].to_string(index=False, max_colwidth=46))
    print(f"\n-> {OUT}/dataset_selection.csv")


if __name__ == "__main__":
    main()
