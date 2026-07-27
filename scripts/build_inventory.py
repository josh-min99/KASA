"""
OSDR 인공중력(Altered Gravity) 마우스 스터디 전수 인벤토리 생성.

핵심 주의: 중력 조건은 /osdr/data/osd/meta/{id} (JSON) 에서 읽으면 안 된다.
그 API는 Factor Value[Altered Gravity] 를 null 로 반환한다(2026-07 확인).
반드시 ISA.zip 안의 s_*.txt 를 파싱할 것.

출력: data/studies_gravity.csv
"""
import os, re, io, zipfile, warnings
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

# Study Factor Name 에 "Altered Gravity" 가 들어간 마우스 스터디 (search API 로 도출)
AG_STUDIES = [758, 759, 714, 686, 288, 289, 532, 426, 238, 239, 29]

GMAP = {
    "uG": 0.0,
    "1/6G with centrifugation": 1 / 6,
    "0.33G by centrifugation": 0.33,
    "0.66G by centrifugation": 0.66,   # 논문 본문은 0.67G 로 표기. 메타데이터는 0.66.
    "1G by centrifugation": 1.0,
    "1G with centrifugation": 1.0,
    "1G on Earth": 1.0,
}


def study_files(osd):
    r = S.get(f"https://osdr.nasa.gov/osdr/data/osd/files/{osd}", timeout=120).json()
    st = r["studies"].get(f"OSD-{osd}")
    return st["study_files"] if st else []


def fetch(osd, pattern):
    cand = [f for f in study_files(osd) if re.search(pattern, f["file_name"])]
    if not cand:
        return None
    f = cand[0]
    path = os.path.join(CACHE, f["file_name"])
    if not os.path.exists(path):
        with open(path, "wb") as fh:
            fh.write(S.get("https://osdr.nasa.gov" + f["remote_url"], timeout=1800).content)
    return path


def sample_table(osd):
    p = fetch(osd, r"ISA\.zip$")
    if not p:
        return None
    z = zipfile.ZipFile(p)
    nm = [n for n in z.namelist() if n.split("/")[-1].startswith("s_")]
    return pd.read_csv(io.BytesIO(z.read(nm[0])), sep="\t", dtype=str) if nm else None


def col(df, *needles):
    for c in df.columns:
        if all(n.lower() in c.lower() for n in needles):
            return c
    return None


rows = []
for osd in AG_STUDIES:
    df = sample_table(osd)
    if df is None:
        continue
    gc = col(df, "Altered Gravity")
    tc = col(df, "Material Type") or col(df, "Organism Part")
    lc = col(df, "light cycle")
    ec = col(df, "Euthanasia Date")
    fc = col(df, "Factor Value[Spaceflight]") or col(df, "Spaceflight")

    counts = df[gc].value_counts().to_dict() if gc else {}
    g_levels = sorted({GMAP[k] for k in counts if k in GMAP})
    flight_levels = sorted({GMAP[k] for k in counts if k in GMAP and k != "1G on Earth"})

    rows.append({
        "OSD_ID": f"OSD-{osd}",
        "tissue": df[tc].dropna().unique()[0] if tc is not None and df[tc].notna().any() else "",
        "n_total": len(df),
        "n_by_gravity": "; ".join(f"{k}={v}" for k, v in sorted(counts.items())),
        "gravity_levels_g": ", ".join(f"{g:.2f}" for g in g_levels),
        "n_inflight_levels": len(flight_levels),
        "dose_response_usable": len(flight_levels) >= 3,
        "light_cycle": df[lc].dropna().unique()[0] if lc is not None and df[lc].notna().any() else "NOT RECORDED",
        "euthanasia_date": ", ".join(sorted(df[ec].dropna().unique())) if ec is not None else "NOT RECORDED",
        "time_of_day_recorded": False,   # 11개 스터디 전수 확인 결과 예외 없음
        "controls": ", ".join(sorted(df[fc].dropna().unique())) if fc is not None else "",
    })

out = pd.DataFrame(rows).sort_values(["n_inflight_levels", "n_total"], ascending=False)
path = os.path.join(OUT, "studies_gravity.csv")
out.to_csv(path, index=False, encoding="utf-8-sig")
print(out.to_string(index=False))
print(f"\n-> {path}")
