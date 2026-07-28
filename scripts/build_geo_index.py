"""
외부(NCBI GEO) 참조 데이터셋 색인.

OSDR 은 우주생명과학 데이터를 모으지만, 지상 circadian 기준 데이터는 GEO 에 있다.
주제가 바뀌어도 재사용할 수 있도록 검색 조건을 코드로 남긴다.

수집 대상
  - SCN(시교차상핵) circadian 시계열
  - 조직별 circadian atlas
  - 중력/과중력 관련 설치류 데이터

출력: data/geo_index.csv
"""
import os, time, json, warnings
import requests
import pandas as pd

warnings.filterwarnings("ignore")
S = requests.Session()
S.headers.update({"User-Agent": "KASA-circadian-gravity/1.0"})

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, "data")
os.makedirs(OUT, exist_ok=True)

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"

QUERIES = {
    "SCN_circadian": "suprachiasmatic AND circadian AND Mus musculus[ORGN] AND gse[ETYP]",
    "circadian_atlas": "circadian AND (atlas OR time course OR timepoints) AND "
                       "Mus musculus[ORGN] AND gse[ETYP]",
    "hypergravity_rodent": "(hypergravity OR centrifugation) AND "
                           "(Mus musculus[ORGN] OR Rattus norvegicus[ORGN]) AND gse[ETYP]",
    "spaceflight_rodent": "(spaceflight OR microgravity) AND "
                          "(Mus musculus[ORGN] OR Rattus norvegicus[ORGN]) AND gse[ETYP]",
    "vestibular_circadian": "(vestibular OR otolith) AND circadian AND gse[ETYP]",
}
RETMAX = 60


def esearch(term, retmax=RETMAX):
    r = S.get(f"{EUTILS}/esearch.fcgi", timeout=120,
              params={"db": "gds", "term": term, "retmax": retmax, "retmode": "json"})
    j = r.json()["esearchresult"]
    return j.get("idlist", []), int(j.get("count", 0))


def esummary(ids):
    if not ids:
        return {}
    out = {}
    for i in range(0, len(ids), 40):
        chunk = ids[i:i + 40]
        r = S.get(f"{EUTILS}/esummary.fcgi", timeout=180,
                  params={"db": "gds", "id": ",".join(chunk), "retmode": "json"})
        out.update({k: v for k, v in r.json().get("result", {}).items() if k != "uids"})
        time.sleep(0.4)
    return out


def main():
    rows = []
    for tag, term in QUERIES.items():
        ids, total = esearch(term)
        print(f"[{tag}] 전체 {total}건, 상위 {len(ids)}건 조회")
        for uid, d in esummary(ids).items():
            rows.append({
                "query": tag,
                "accession": d.get("accession"),
                "title": str(d.get("title", ""))[:160],
                "organism": d.get("taxon"),
                "type": str(d.get("gdstype", ""))[:60],
                "n_samples": int(d.get("n_samples", 0) or 0),
                "platform": d.get("gpl"),
                "pubmed": ";".join(map(str, d.get("pubmedids", []) or [])),
                "summary": str(d.get("summary", ""))[:300],
            })
        time.sleep(0.5)

    df = pd.DataFrame(rows).drop_duplicates("accession").sort_values(
        ["query", "n_samples"], ascending=[True, False])
    df.to_csv(os.path.join(OUT, "geo_index.csv"), index=False, encoding="utf-8-sig")

    print(f"\n총 {len(df)}건 색인")
    for q, g in df.groupby("query"):
        print(f"\n--- {q} (상위 8건) ---")
        print(g.head(8)[["accession", "n_samples", "organism", "title"]]
              .to_string(index=False, max_colwidth=76))
    print(f"\n-> {OUT}/geo_index.csv")


if __name__ == "__main__":
    main()
