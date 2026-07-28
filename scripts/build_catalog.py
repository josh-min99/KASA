"""
OSDR 전수 카탈로그 구축 — 검색 인덱스에 의존하지 않는다.

배경
  초기 분석은 탐색 범위를 세 번 잘못 잡았다.
    1) organism=Mus musculus 로만 검색 -> 쥐(Rattus) 스터디 누락
    2) 검색 API 인덱스 필드만 정규식 검사 -> ISA 전용 컬럼 누락
    3) Factor Value[Altered Gravity] 가 있는 스터디만 대상 -> 다른 표기 누락
  세 오류 모두 '무엇을 검색할지 미리 정한' 데서 나왔다.

  따라서 여기서는 OSD ID 를 1 부터 MAX_ID 까지 전수 순회하고,
  각 스터디의 ISA 아카이브를 열어 모든 컬럼과 값을 기록한다.
  필터링은 카탈로그를 다 만든 뒤에 한다.

출력: data/catalog_studies.csv   스터디 1행
      data/catalog_columns.csv   스터디 x 컬럼 1행 (값 샘플 포함)
"""
import os, re, io, json, zipfile, warnings, time
import requests
import pandas as pd

warnings.filterwarnings("ignore")
S = requests.Session()
S.headers.update({"User-Agent": "KASA-circadian-gravity/1.0"})

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
CACHE = os.path.join(ROOT, "cache")
ISADIR = os.path.join(CACHE, "isa")
OUT = os.path.join(ROOT, "data")
STATE = os.path.join(CACHE, "catalog_state.json")
for d in (CACHE, ISADIR, OUT):
    os.makedirs(d, exist_ok=True)

MAX_ID = 1100
SKIP_COL = re.compile(r"^(Term Source REF|Term Accession Number|Unit)", re.I)


def get_json(url, timeout=90, tries=3):
    for k in range(tries):
        try:
            r = S.get(url, timeout=timeout)
            if r.status_code == 200:
                return r.json()
            if r.status_code in (404, 400):
                return None
        except Exception:
            time.sleep(1 + k)
    return None


def study_files(osd):
    j = get_json(f"https://osdr.nasa.gov/osdr/data/osd/files/{osd}")
    if not j:
        return None
    st = j.get("studies", {}).get(f"OSD-{osd}")
    return st["study_files"] if st else None


def isa_tables(osd, files):
    isa = [f for f in files if f["file_name"].endswith("ISA.zip")]
    if not isa:
        return {}
    p = os.path.join(ISADIR, isa[0]["file_name"])
    if not os.path.exists(p):
        try:
            b = S.get("https://osdr.nasa.gov" + isa[0]["remote_url"], timeout=900).content
        except Exception:
            return {}
        with open(p, "wb") as fh:
            fh.write(b)
    try:
        z = zipfile.ZipFile(p)
    except Exception:
        return {}
    out = {}
    for n in z.namelist():
        base = n.split("/")[-1]
        if base.startswith(("s_", "a_")):
            try:
                out[base] = pd.read_csv(io.BytesIO(z.read(n)), sep="\t", dtype=str)
            except Exception:
                pass
    return out


def uniq(df, col, k=8):
    v = df[col].dropna().unique()
    return len(v), [str(x)[:90] for x in v[:k]]


def main():
    done = {}
    if os.path.exists(STATE):
        done = json.load(open(STATE, encoding="utf-8"))
    studies, columns = [], []

    for osd in range(1, MAX_ID + 1):
        key = str(osd)
        if osd % 25 == 0:
            print(f"  {osd}/{MAX_ID}  (수집 {len(studies)})", flush=True)
            json.dump(done, open(STATE, "w", encoding="utf-8"))

        files = study_files(osd)
        if not files:
            done[key] = "none"
            continue
        tabs = isa_tables(osd, files)
        if not tabs:
            done[key] = "no-isa"
            continue

        has_counts = any(re.search(r"Normalized_Counts_GLbulkRNAseq\.csv$", f["file_name"])
                         for f in files)
        rec = {"OSD_ID": f"OSD-{osd}", "n_isa_tables": len(tabs),
               "has_geneLab_counts": has_counts, "n_files": len(files)}
        org, mat, assay, facs, nsamp = set(), set(), set(), {}, 0

        for base, df in tabs.items():
            if base.startswith("s_"):
                nsamp = max(nsamp, len(df))
            for c in df.columns:
                if SKIP_COL.match(c):
                    continue
                n_u, vals = uniq(df, c)
                columns.append({"OSD_ID": f"OSD-{osd}", "table": base, "column": c,
                                "n_unique": n_u, "values": " | ".join(vals)})
                lc = c.lower()
                if "organism]" in lc:
                    org |= set(df[c].dropna().unique())
                if "material type" in lc or "organism part" in lc:
                    mat |= set(df[c].dropna().unique())
                if c.startswith("Factor Value"):
                    facs[c] = vals
            if base.startswith("a_"):
                assay.add(base)

        rec.update({"organism": "; ".join(sorted(map(str, org)))[:120],
                    "material": "; ".join(sorted(map(str, mat)))[:160],
                    "n_samples": nsamp,
                    "n_assays": len(assay),
                    "assays": "; ".join(sorted(assay))[:300],
                    "factors": "; ".join(sorted(facs))[:400],
                    "factor_values": json.dumps(facs, ensure_ascii=False)[:900]})
        studies.append(rec)
        done[key] = "ok"

    json.dump(done, open(STATE, "w", encoding="utf-8"))
    sdf = pd.DataFrame(studies)
    cdf = pd.DataFrame(columns)
    sdf.to_csv(os.path.join(OUT, "catalog_studies.csv"), index=False, encoding="utf-8-sig")
    cdf.to_csv(os.path.join(OUT, "catalog_columns.csv"), index=False, encoding="utf-8-sig")
    print(f"\n스터디 {len(sdf)}건, 컬럼 레코드 {len(cdf)}건")
    print(f"-> {OUT}/catalog_studies.csv, catalog_columns.csv")


if __name__ == "__main__":
    main()
