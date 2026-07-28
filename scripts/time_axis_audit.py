"""
시간축 감사 — OSDR 설치류 스터디에 '희생 시각'이 기록돼 있는가.

이 스크립트는 초기 분석의 오류를 바로잡기 위해 작성됐다.
초기에는 검색 API 인덱스에 노출된 필드(Parameter Value, Characteristics,
Study Protocol Description 등)만 정규식으로 훑고 "시각 기록 0건"이라고 결론지었다.
검색 인덱스는 ISA 파일의 모든 컬럼을 노출하지 않으므로 그 방법으로는 불완전하다.

여기서는 설치류 스터디 전수의 ISA 아카이브를 실제로 열어
**모든 컬럼의 모든 값**을 시각 패턴으로 검사한다.

추가로 두 가지 방법론 오류를 함께 수정한다.
  (1) 초기 인벤토리는 organism=Mus musculus 로만 검색해 쥐(Rattus) 스터디를 놓쳤다.
      OSD-616/617/652/653 은 Altered Gravity 가 factor 인 쥐 부분중력 연구다.
  (2) 검색 인덱스 기반 지상 스터디 추출은 payload API 보다 불완전하다.

출력: data/time_axis_audit.csv, data/gravity_studies_all_rodent.csv
"""
import os, re, io, zipfile, json, warnings
import requests
import pandas as pd

warnings.filterwarnings("ignore")
S = requests.Session()
S.headers.update({"User-Agent": "KASA-circadian-gravity/1.0"})

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
CACHE = os.path.join(ROOT, "cache")
OUT = os.path.join(ROOT, "data")
os.makedirs(CACHE, exist_ok=True)
os.makedirs(OUT, exist_ok=True)

ORGANISMS = ["Mus musculus", "Rattus norvegicus"]

# 값이 '하루 중 시각'으로 보이는 패턴
CLOCK = re.compile(r"\b\d{1,2}:\d{2}\b|\b\d{1,2}\s*(?:am|pm)\b|\bZT\s*-?\d{1,2}\b"
                   r"|\bCT\s*-?\d{1,2}\b|zeitgeber", re.I)
# 오탐: 명암주기 표기(12:12), 희석비(1:100), 날짜
NOISE = re.compile(r"12\s*:\s*12|1\s*:\s*\d{2,}|light/?dark|light-dark|\d{2}-[A-Za-z]{3}-\d{4}", re.I)
SKIP_COL = re.compile(r"^(Term Source REF|Term Accession Number|Unit)", re.I)
GRAV = re.compile(r"gravity|hypergrav|microgravity simulation|centrifug|unloading", re.I)


def search_all(params, page=100, cap=9000):
    hits, frm = [], 0
    while frm < cap:
        p = dict(params); p.update(size=page, **{"from": frm})
        r = S.get("https://osdr.nasa.gov/osdr/data/search", params=p, timeout=180).json()
        h = r["hits"]["hits"]
        if not h:
            break
        hits += h; frm += page
        if frm >= r["hits"]["total"]:
            break
    return hits


def rodent_accessions():
    accs = set()
    for org in ORGANISMS:
        for h in search_all({"term": "study", "ffield": "organism", "fvalue": org}):
            a = h["_source"].get("Accession")
            if isinstance(a, str) and a.startswith("OSD-"):
                accs.add(int(a.split("-")[1]))
    return sorted(accs)


def isa_tables(osd):
    r = S.get(f"https://osdr.nasa.gov/osdr/data/osd/files/{osd}", timeout=120).json()
    st = r["studies"].get(f"OSD-{osd}")
    if not st:
        return {}
    isa = [f for f in st["study_files"] if f["file_name"].endswith("ISA.zip")]
    if not isa:
        return {}
    p = os.path.join(CACHE, isa[0]["file_name"])
    if not os.path.exists(p):
        with open(p, "wb") as fh:
            fh.write(S.get("https://osdr.nasa.gov" + isa[0]["remote_url"], timeout=900).content)
    z = zipfile.ZipFile(p)
    return {n: pd.read_csv(io.BytesIO(z.read(n)), sep="\t", dtype=str)
            for n in z.namelist() if n.split("/")[-1].startswith(("s_", "a_"))}


def main():
    accs = rodent_accessions()
    print(f"설치류 OSD 스터디 {len(accs)}건 (Mus + Rattus)")

    hits, grav = [], []
    for i, osd in enumerate(accs, 1):
        if i % 50 == 0:
            print(f"  {i}/{len(accs)}")
        try:
            tabs = isa_tables(osd)
        except Exception:
            continue
        gcols = set()
        for nm, df in tabs.items():
            for c in df.columns:
                if SKIP_COL.match(c):
                    continue
                if c.startswith("Factor Value") and GRAV.search(c):
                    gcols.add(c)
                vals = df[c].dropna().unique()
                for v in vals:
                    sv = str(v)
                    if CLOCK.search(sv) and not NOISE.search(sv):
                        hits.append({"OSD_ID": f"OSD-{osd}", "column": c,
                                     "n_unique": len(vals),
                                     "varies_across_samples": len(vals) > 1,
                                     "sample_values": " | ".join(map(str, vals[:6]))})
                        break
        if gcols:
            grav.append({"OSD_ID": f"OSD-{osd}", "gravity_factors": "; ".join(sorted(gcols))})

    hdf = pd.DataFrame(hits)
    hdf.to_csv(os.path.join(OUT, "time_axis_audit.csv"), index=False, encoding="utf-8-sig")
    gdf = pd.DataFrame(grav)
    gdf.to_csv(os.path.join(OUT, "gravity_studies_all_rodent.csv"), index=False, encoding="utf-8-sig")

    print("\n" + "=" * 88)
    print(f"시각 패턴이 값에 나타나는 항목: {len(hdf)}건")
    if len(hdf):
        print(hdf.to_string(index=False, max_colwidth=60))

    varying = hdf[hdf.varies_across_samples] if len(hdf) else hdf
    # 명암주기 기준점 표기는 시각이지만 샘플 간 변이가 없다
    real = varying[~varying.column.str.contains("light cycle", case=False)] if len(varying) else varying
    print("\n" + "=" * 88)
    print(f"샘플 간 시각이 '달라지는' 스터디: {real.OSD_ID.nunique() if len(real) else 0}건")
    if len(real):
        print(real.to_string(index=False, max_colwidth=60))
    print("\n→ 리듬 분석에 필요한 것은 '기록 여부'가 아니라 '샘플 간 변이'다.")
    print(f"\n중력이 factor 인 설치류 스터디: {len(gdf)}건")
    print(f"-> {OUT}/time_axis_audit.csv, gravity_studies_all_rodent.csv")


if __name__ == "__main__":
    main()
