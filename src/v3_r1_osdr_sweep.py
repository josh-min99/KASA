"""
R1-a: OSDR/ALSDA 전수 스윕 — 일주기 해상도를 가질 가능성이 있는 파일 탐색.

기존 `scripts/ground_and_telemetry.py` 와의 차이
  그 스크립트는 (1) 마우스로 organism 을 한정했고 (2) 검색 API 의 인덱스 필드만 봤다.
  데이터_지도.md §6.1 에 그 두 가지가 범위 오류로 기록돼 있다.
  여기서는 종을 한정하지 않고, **각 스터디의 실제 파일 목록까지 열어** 파일명 수준에서
  시계열 신호를 찾는다. 메타데이터에 'transcription profiling' 만 적혀 있어도
  부속 파일에 텔레메트리가 들어 있을 수 있기 때문이다(OSD-952 가 실제로 그런 사례였다).

산출
  data/rhythm/osdr_files_index.csv   전 스터디 파일 목록 (캐시)
  data/rhythm/osdr_candidates.csv    리듬 신호가 걸린 스터디/파일

실행: python src/v3_r1_osdr_sweep.py
"""
import os
import re
import sys
import json
import time
import warnings

import requests
import pandas as pd

warnings.filterwarnings("ignore")
if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
CACHE = os.path.join(ROOT, "cache", "rhythm")
OUT = os.path.join(ROOT, "data", "rhythm")
os.makedirs(CACHE, exist_ok=True)
os.makedirs(OUT, exist_ok=True)

S = requests.Session()
S.headers.update({"User-Agent": "KASA-circadian-gravity/3.0"})

MAX_ID = 1100  # scripts/build_catalog.py 와 동일

# 파일명·어세이명에서 시계열 가능성을 시사하는 토큰.
# 넓게 잡고 나중에 사람이 걸러낸다. 좁게 잡아 놓치는 쪽이 훨씬 비싸다.
RHYTHM_PAT = re.compile(
    r"telemetr|biotelemetr|actigraph|actogram|activity|locomot|behavio|ethogram|"
    r"circadian|rhythm|acrophase|cosinor|"
    r"temperatur|\bTb\b|thermo|"
    r"heart[\s_-]?rate|\bECG\b|\bEKG\b|\bEEG\b|polysomn|sleep|"
    r"accelerom|wheel|infrared|video|segmentation|pose|"
    r"time[\s_-]?series|timeseries|continuous|hourly|per[\s_-]?minute",
    re.I)

# 중력 조건. Altered Gravity factor 가 없는 스터디도 잡히도록 자유 텍스트까지 본다.
GRAV_PAT = re.compile(
    r"space[\s_-]?flight|microgravity|micro-g|\buG\b|"
    r"hindlimb (unloading|suspension)|tail suspension|antiorthostatic|"
    r"bed[\s_-]?rest|head[\s_-]?down|dry immersion|"
    r"clinostat|random positioning|rotating wall|"
    r"hypergravity|centrifug|parabolic|"
    r"partial gravity|artificial gravity|altered gravity|lunar gravity|mars gravity",
    re.I)

# 명백히 시계열이 아닌데 토큰만 걸리는 것들 (omics 파일명에 'activity' 등이 흔하다).
# 배제 패턴은 과거에 정탐을 지운 전력이 있으므로, 지운 목록을 반드시 남긴다.
OMICS_PAT = re.compile(
    r"\.fastq|\.fq\.gz|\.bam\b|\.bai\b|\.cel\b|\.idat|\.vcf|"
    r"counts?[_\-\.]|normalized[_\-]counts|differential[_\-]expression|"
    r"multiqc|rRNA|trimmed|aligned\.out|genebody|rsem|salmon|kallisto",
    re.I)


def get_json(url, params=None, tries=3, timeout=180):
    for i in range(tries):
        try:
            r = S.get(url, params=params, timeout=timeout)
            if r.status_code == 200:
                return r.json()
            if r.status_code == 404:
                return None
        except Exception:
            pass
        time.sleep(2 * (i + 1))
    return None


def all_study_ids():
    """전 스터디 메타데이터.

    주의 — 검색 API 로만 열거하면 안 된다.
    초판에서 `term=study` 페이지네이션으로 열거했더니 411건만 잡혔는데,
    같은 저장소의 `data/catalog_studies.csv` 에는 626건이 있었다. 222건 누락이다.
    `scripts/build_catalog.py` 가 ID 를 1..MAX_ID 로 전수 순회하는 이유가 이것이다.
    따라서 여기서도 전수 순회를 기준으로 삼고, 검색 API 결과는 메타데이터 보강에만 쓴다.
    """
    cache = os.path.join(CACHE, "osdr_all_studies.json")
    if os.path.exists(cache):
        with open(cache, encoding="utf-8") as fh:
            return json.load(fh)

    # (1) 검색 API — 메타데이터(설명·factor)가 풍부하므로 보강용으로 먼저 받는다
    meta, frm = {}, 0
    while frm < 20000:
        j = get_json("https://osdr.nasa.gov/osdr/data/search",
                     {"term": "study", "size": 100, "from": frm})
        if not j:
            break
        hits = j.get("hits", {}).get("hits", [])
        if not hits:
            break
        for h in hits:
            s = h["_source"]
            acc = s.get("Accession")
            if acc and str(acc).startswith("OSD-"):
                meta[acc] = {
                    "title": str(s.get("Study Title", ""))[:200],
                    "organism": str(s.get("organism", ""))[:80],
                    "factors": str(s.get("Study Factor Name", ""))[:300],
                    "measurement": str(s.get("Study Assay Measurement Type", ""))[:200],
                    "source_type": str(s.get("Data Source Type", ""))[:60],
                    "description": str(s.get("Study Description", ""))[:1500],
                }
        frm += 100
        if frm >= j.get("hits", {}).get("total", 0):
            break
    print(f"    검색 API: {len(meta)}건")

    # (2) 기존 카탈로그의 ID 도 합친다 (같은 저장소가 이미 확인한 전수 목록)
    known = {}
    cat = os.path.join(ROOT, "data", "catalog_studies.csv")
    if os.path.exists(cat):
        c = pd.read_csv(cat)
        for _, r in c.iterrows():
            known[str(r["OSD_ID"])] = {
                "organism": str(r.get("organism", ""))[:80],
                "factors": str(r.get("factors", ""))[:300],
                "measurement": str(r.get("assays", ""))[:200],
                "factor_values": str(r.get("factor_values", ""))[:1500],
            }
        print(f"    기존 카탈로그: {len(known)}건")

    # (3) 전수 ID 순회로 합집합을 만든다
    ids = set(meta) | set(known)
    print(f"    합집합 {len(ids)}건 (누락 방지용 전수 순회는 MAX_ID={MAX_ID} 까지)")

    out = []
    for n in range(1, MAX_ID + 1):
        oid = f"OSD-{n}"
        m = meta.get(oid, {})
        k = known.get(oid, {})
        if not m and not k and oid not in ids:
            continue
        out.append({
            "OSD_ID": oid,
            "title": m.get("title", ""),
            "organism": m.get("organism") or k.get("organism", ""),
            "factors": m.get("factors") or k.get("factors", ""),
            "measurement": m.get("measurement") or k.get("measurement", ""),
            "source_type": m.get("source_type", ""),
            "description": (m.get("description", "") + " " + k.get("factor_values", ""))[:2000],
        })
    with open(cache, "w", encoding="utf-8") as fh:
        json.dump(out, fh, ensure_ascii=False)
    return out


def files_for(osd_id):
    """스터디 파일 목록. 캐시."""
    num = osd_id.split("-")[1]
    cache = os.path.join(CACHE, f"files_{num}.json")
    if os.path.exists(cache):
        try:
            with open(cache, encoding="utf-8") as fh:
                return json.load(fh)
        except Exception:
            pass
    j = get_json(f"https://osdr.nasa.gov/osdr/data/osd/files/{num}", timeout=120)
    files = []
    if j:
        try:
            st = j["studies"][osd_id]["study_files"]
            files = [{"file_name": f.get("file_name", ""),
                      "remote_url": f.get("remote_url", ""),
                      "file_size": f.get("file_size", 0),
                      "category": str(f.get("category", "")),
                      "subcategory": str(f.get("subcategory", ""))} for f in st]
        except Exception:
            files = []
    with open(cache, "w", encoding="utf-8") as fh:
        json.dump(files, fh, ensure_ascii=False)
    return files


def main():
    print("=" * 78)
    print("R1-a  OSDR 전수 스윕")
    print("=" * 78)

    studies = all_study_ids()
    print(f"\n전 스터디 {len(studies)}건 (종 무관)")
    org = pd.Series([s["organism"] for s in studies]).value_counts()
    print("\n상위 종:")
    print(org.head(12).to_string())

    # 1) 중력 조건이 걸리는 스터디로 1차 축소 (제목+factor+description+measurement)
    grav = []
    for s in studies:
        blob = " ".join([s["title"], s["factors"], s["description"], s["measurement"]])
        m = GRAV_PAT.search(blob)
        if m:
            s = dict(s)
            s["gravity_cue"] = m.group(0).lower()
            grav.append(s)
    print(f"\n중력 조건이 걸린 스터디: {len(grav)}건")

    # 2) 파일 목록은 **전 스터디** 를 연다.
    #    중력 단서로 미리 거르지 않는다. 초판에서 그렇게 했다가, 메타데이터에
    #    'transcription profiling' 만 적혀 있고 부속 파일에 텔레메트리가 들어 있는
    #    구조(OSD-952 가 그런 사례)를 놓칠 위험이 있음을 확인했다.
    #    중력 조건 판정은 파일을 다 본 뒤에 붙인다.
    gcue = {s["OSD_ID"]: s.get("gravity_cue", "") for s in grav}
    targets = {}
    for s in studies:
        s = dict(s)
        s["gravity_cue"] = gcue.get(s["OSD_ID"], "")
        targets[s["OSD_ID"]] = s
    print(f"파일 목록을 열 대상: {len(targets)}건 (전 스터디)\n")

    idx_rows, cand_rows, dropped = [], [], []
    for i, (oid, s) in enumerate(sorted(targets.items(), key=lambda x: int(x[0].split("-")[1])), 1):
        fs = files_for(oid)
        if i % 25 == 0:
            print(f"    [{i}/{len(targets)}] {oid}  누적 후보 {len(cand_rows)}")
        for f in fs:
            fn = f["file_name"]
            idx_rows.append({"OSD_ID": oid, "file_name": fn, "file_size": f["file_size"],
                             "category": f["category"], "subcategory": f["subcategory"]})
            hit = RHYTHM_PAT.search(fn + " " + f["category"] + " " + f["subcategory"])
            if not hit:
                continue
            if OMICS_PAT.search(fn):
                dropped.append({"OSD_ID": oid, "file_name": fn, "token": hit.group(0)})
                continue
            cand_rows.append({
                "OSD_ID": oid, "organism": s["organism"], "title": s["title"],
                "gravity_cue": s.get("gravity_cue", ""), "source_type": s["source_type"],
                "file_name": fn, "file_size": f["file_size"],
                "category": f["category"], "subcategory": f["subcategory"],
                "rhythm_token": hit.group(0).lower(),
                "remote_url": f["remote_url"],
            })

    idx = pd.DataFrame(idx_rows)
    cand = pd.DataFrame(cand_rows)
    drop = pd.DataFrame(dropped)

    idx.to_csv(os.path.join(OUT, "osdr_files_index.csv"), index=False, encoding="utf-8-sig")
    cand.to_csv(os.path.join(OUT, "osdr_candidates.csv"), index=False, encoding="utf-8-sig")
    drop.to_csv(os.path.join(OUT, "osdr_dropped_by_omics_filter.csv"), index=False, encoding="utf-8-sig")

    print("\n" + "=" * 78)
    print(f"파일 총 {len(idx):,}개 색인")
    print(f"리듬 후보 파일 {len(cand):,}개 / 스터디 {cand.OSD_ID.nunique() if len(cand) else 0}건")
    print(f"omics 필터로 제외 {len(drop):,}개 (전량 osdr_dropped_by_omics_filter.csv 에 보존)")

    if len(cand):
        print("\n스터디별 후보 파일 수 상위 25:")
        top = (cand.groupby(["OSD_ID", "organism", "gravity_cue"])
               .agg(n_files=("file_name", "size"),
                    tokens=("rhythm_token", lambda x: ",".join(sorted(set(x))[:5])))
               .reset_index().sort_values("n_files", ascending=False))
        print(top.head(25).to_string(index=False))

        print("\n중력 단서가 있는 스터디만:")
        g = top[top.gravity_cue != ""]
        print(f"  {len(g)}건")
        print(g.head(30).to_string(index=False))

    print(f"\n-> {OUT}/osdr_candidates.csv")


if __name__ == "__main__":
    main()
