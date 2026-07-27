"""
Week 1 잔여 + Week 2: 착륙 교란이 없는 지상 데이터 + 리듬 시계열 데이터 탐색.

동기
  비행 데이터는 전부 착륙 후 채취라 중력 효과와 재적응 효과가 분리되지 않는다
  (OSD-758: splashdown 6시간 이내, OSD-714: 착륙 2일 후).
  지상 모사(hindlimb unloading, 지상 원심분리)는 그 교란이 없다.

  또한 전사체는 위상 정보를 주지 못하므로, 실제 시계열이 존재하는
  생리·행동 텔레메트리(체온/활동량)가 있는지 확인한다. 있으면 리듬 분석의 본진이 된다.

출력: data/ground_studies.csv, data/telemetry_candidates.csv
"""
import os, re, io, zipfile, json, warnings, collections
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


def search_all(params, page=100, cap=8000):
    hits, frm = [], 0
    while frm < cap:
        p = dict(params); p.update(size=page, **{"from": frm})
        r = S.get("https://osdr.nasa.gov/osdr/data/search", params=p, timeout=180).json()
        h = r["hits"]["hits"]
        if not h:
            break
        hits += h
        frm += page
        if frm >= r["hits"]["total"]:
            break
    return hits


def uniq_by_accession(hits):
    out = {}
    for h in hits:
        out.setdefault(h["_source"].get("Accession"), h["_source"])
    return out


# ------------------------------------------------------- 1) 지상 중력 모사 스터디
print("=" * 78)
print("[1] 지상 중력 모사 / 지상 원심분리 마우스 스터디")

mouse = uniq_by_accession(search_all({"term": "study", "ffield": "organism", "fvalue": "Mus musculus"}))
print(f"    마우스 레코드 {len(mouse)}건 스캔")

GROUND_PAT = re.compile(r"hindlimb unloading|hindlimb suspension|microgravity simulation|"
                        r"tail suspension|clinostat|random positioning|rotating wall|"
                        r"hypergravity|centrifugation", re.I)
FLIGHT_PAT = re.compile(r"spaceflight|space flight", re.I)

rows = []
for acc, s in mouse.items():
    if not str(acc).startswith("OSD-"):
        continue
    fac = str(s.get("Study Factor Name", ""))
    title = str(s.get("Study Title", ""))
    desc = str(s.get("Study Description", ""))
    blob = f"{fac} {title}"
    if not GROUND_PAT.search(blob):
        continue
    rows.append({
        "OSD_ID": acc,
        "title": title[:120],
        "factors": re.sub(r"\s{2,}", " | ", fac),
        "material": str(s.get("Material Type", ""))[:60],
        "assay": str(s.get("Study Assay Technology Type", ""))[:60],
        "has_flight_factor": bool(FLIGHT_PAT.search(fac)),
        "ground_cue": GROUND_PAT.search(blob).group(0).lower(),
        "brain_tissue": bool(re.search(r"brain|hypothalam|suprachiasm|SCN|cortex|hippocamp",
                                       str(s.get("Material Type", "")) + title, re.I)),
    })

gdf = pd.DataFrame(rows).drop_duplicates("OSD_ID").sort_values("OSD_ID")
gdf.to_csv(os.path.join(OUT, "ground_studies.csv"), index=False, encoding="utf-8-sig")
print(f"    지상 모사 후보 {len(gdf)}건")
print(gdf["ground_cue"].value_counts().to_string())
print("\n    조직 분포 상위:")
print(gdf["material"].value_counts().head(12).to_string())
print(f"\n    뇌 조직 포함: {int(gdf.brain_tissue.sum())}건")
if gdf.brain_tissue.any():
    print(gdf[gdf.brain_tissue][["OSD_ID", "material", "ground_cue", "title"]].to_string(index=False))


# ------------------------------------------------- 2) 리듬 시계열 / 텔레메트리 탐색
print("\n" + "=" * 78)
print("[2] 생리·행동 텔레메트리 (체온 / 활동량 / 심박) 탐색")

# 2-1. assay measurement type 전수
mt = collections.Counter()
for s in mouse.values():
    v = s.get("Study Assay Measurement Type", "")
    for tok in re.split(r"\s{2,}", str(v)):
        if tok.strip():
            mt[tok.strip()] += 1
TELE_PAT = re.compile(r"telemetr|temperature|activity|behavio|locomot|heart|physiolog|"
                      r"actigraph|circadian|rhythm|environmental", re.I)
print("\n    [2-1] 측정 유형 중 시계열 가능성 있는 것:")
found_mt = {k: v for k, v in mt.items() if TELE_PAT.search(k)}
print("        " + (json.dumps(found_mt, ensure_ascii=False) if found_mt else "없음"))
print(f"    (전체 측정 유형 {len(mt)}종. 상위 10: {[k for k, _ in mt.most_common(10)]})")

# 2-2. 비-omics 데이터 소스 유형
dst = collections.Counter(str(s.get("Data Source Type", "")) for s in mouse.values())
print("\n    [2-2] Data Source Type 분포:")
for k, v in dst.most_common(10):
    print(f"        [{v}] {k[:70]}")

# 2-3. 자유 텍스트 검색
print("\n    [2-3] 자유 텍스트 검색:")
tele_rows = []
for term in ["telemetry", "core body temperature", "locomotor activity",
             "circadian rhythm", "actigraphy", "heart rate"]:
    hits = uniq_by_accession(search_all({"term": term}, cap=300))
    osd = [a for a in hits if str(a).startswith("OSD-")]
    print(f"        '{term}': {len(hits)}건 (OSD {len(osd)}건)")
    for a in osd:
        s = hits[a]
        tele_rows.append({"term": term, "OSD_ID": a,
                          "title": str(s.get("Study Title", ""))[:100],
                          "organism": str(s.get("organism", ""))[:40],
                          "measurement": str(s.get("Study Assay Measurement Type", ""))[:60],
                          "source_type": str(s.get("Data Source Type", ""))[:40]})

tdf = pd.DataFrame(tele_rows).drop_duplicates(["OSD_ID", "term"])
tdf.to_csv(os.path.join(OUT, "telemetry_candidates.csv"), index=False, encoding="utf-8-sig")
print(f"\n    후보 {tdf.OSD_ID.nunique()}건 -> data/telemetry_candidates.csv")
if len(tdf):
    print(tdf.drop_duplicates("OSD_ID").head(25).to_string(index=False))

# 2-4. ALSDA / 비-omics 저장소 확인
print("\n    [2-4] ALSDA(생명과학 데이터 아카이브) 엔드포인트 확인:")
for u in ["https://osdr.nasa.gov/osdr/data/search?term=ALSDA",
          "https://osdr.nasa.gov/geode-py/ws/api/payloads",
          "https://osdr.nasa.gov/geode-py/ws/api/biospecimens"]:
    try:
        r = S.get(u, timeout=90)
        n = ""
        if r.status_code == 200:
            try:
                j = r.json()
                n = f" hits={j.get('hits', {}).get('total')}" if "hits" in j else f" keys={list(j)[:3]}"
            except Exception:
                n = ""
        print(f"        {r.status_code}{n}  {u}")
    except Exception as e:
        print(f"        ERR {u}: {e}")

print(f"\n-> {OUT}/ground_studies.csv, {OUT}/telemetry_candidates.csv")
