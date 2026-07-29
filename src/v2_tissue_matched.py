"""
v2 — 조직을 맞춘 FLIGHT vs HLU 비교.

문제의식
  어제 Stage 3c/4 는 FLIGHT-A(EDL·부신·비복근·신장·간·사두근·가자미근·전경골근)와
  HLU-A(뇌·망막·비장·등쪽피부)를 비교했다. 두 세트의 조직이 하나도 겹치지 않는다.
  따라서 '비행 3.11 vs HLU 0.78' 이 처리 효과인지 조직 효과인지 구별할 수 없다.
  I4 는 데이터셋 내부에서 자기 정규화되므로 배치 효과에 원래 강하다.
  배치를 피하려고 조직 매칭을 포기할 이유가 없었다. 조직을 맞춰 다시 계산한다.

대비 선택 규칙 (어제 정규식 오탐 사고를 구조적으로 막기 위해 방식을 바꿨다)
  블랙리스트로 '나쁜 대비를 지우는' 대신,
  **중력/하중 이외의 모든 조건이 양쪽에서 동일할 것**을 요구한다.
    1) 대비 문자열 "(A)v(B)" 를 ' & ' 로 분해한다
    2) 각 성분을 중력 성분 / 배경 성분으로 나눈다
    3) 배경 성분의 다중집합이 A 와 B 에서 완전히 같아야 채택
    4) 중력 성분은 A 가 처리(비행·현수), B 가 대조여야 한다
  이러면 유전형(Zfp697 mKO, AMPKalpha mKO), 병용처치(x-ray, corticosterone),
  주사(tetanus toxoid, CpG), 시점, 계통이 자동으로 매칭된다.
  매칭되는 대비가 없으면 '대비 불가' 로 기록하고 제외한다. 억지로 만들지 않는다.

품질 게이트
  DEG(FDR<0.05) 개수가 임계 미만인 데이터셋은 '신호 없음' 으로 제외한다.
  어제 RR-6 처럼 DEG 한 자릿수인 것을 '재현 실패' 로 오독하지 않기 위함이다.

출력: results/v2/tables/*.csv
"""
import os, re, time, warnings
import requests
import numpy as np
import pandas as pd

warnings.filterwarnings("ignore")
S = requests.Session()
S.headers.update({"User-Agent": "KASA-drylab-v2/1.0"})

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RAW = os.path.join(ROOT, "data", "raw")
V2 = os.path.join(ROOT, "results", "v2")
TAB = os.path.join(V2, "tables")
for d in (TAB, os.path.join(V2, "figures")):
    os.makedirs(d, exist_ok=True)

CLOCK = ["Bmal1", "Arntl", "Clock", "Npas2", "Per1", "Per2", "Per3", "Cry1", "Cry2",
         "Nr1d1", "Nr1d2", "Dbp", "Tef", "Hlf", "Ciart", "Bhlhe40", "Bhlhe41"]
DEG_MIN = 30                    # 품질 게이트
MAX_DE_MB = 500                 # DE 파일 다운로드 상한

# 조직 | OSD | 처리유형
TARGETS = [
    ("비복근",   876, "HLU"), ("비복근",   880, "HLU"),
    ("비복근",   101, "FLIGHT"), ("비복근",   419, "FLIGHT"),
    ("가자미근", 935, "HLU"), ("가자미근", 949, "HLU"),
    ("가자미근", 104, "FLIGHT"), ("가자미근", 714, "FLIGHT"), ("가자미근", 770, "FLIGHT"),
    ("비장",     201, "HLU"), ("비장",     211, "HLU"),
    ("비장",     246, "FLIGHT"), ("비장",     288, "FLIGHT"), ("비장",     506, "FLIGHT"),
    ("망막",     203, "HLU"),
    ("망막",      87, "FLIGHT"), ("망막",    194, "FLIGHT"),
    ("망막",     255, "FLIGHT"), ("망막",    758, "FLIGHT"),
    ("등쪽피부", 237, "HLU"),
    ("등쪽피부", 238, "FLIGHT"), ("등쪽피부", 240, "FLIGHT"),
    ("등쪽피부", 243, "FLIGHT"), ("등쪽피부", 254, "FLIGHT"),
]

# 중력/하중 성분으로 인식할 토큰
GRAV_TOKEN = re.compile(
    r"space ?flight|ground control|vivarium|basal control"
    r"|hindlimb unloaded|hindlimb loaded control|normally loaded"
    r"|^\s*u\s*g\s*$|micrograv|\d\s*/?\s*\d*\s*g\b|centrifugation|on earth", re.I)
# 처리(=중력 감소) 판정
TREAT_TOKEN = re.compile(r"space ?flight|hindlimb unloaded|^\s*u\s*g\s*$|micrograv", re.I)
CTRL_TOKEN = re.compile(r"ground control|vivarium|hindlimb loaded control"
                        r"|normally loaded|1\s*g on earth", re.I)
# 처리군에 붙어서는 안 되는 것 (중력 외 추가 처치)
DISQUALIFY = re.compile(r"reloaded", re.I)
# 궤도에서 원심분리로 1G 를 '복원' 한 군. Space Flight 이지만 중력 감소가 아니다.
# (OSD-758·288·714·238 에 존재. 규칙 검증에서 잘못 채택되는 것을 발견해 추가)
RESTORED_G = re.compile(r"1\s*g\s*(by|with)\s*centrifugation", re.I)
# 배경이 양쪽에서 일치해도 '깨끗한' 배경이 낫다. 동률일 때 이것으로 우선순위를 정한다.
# 주의: 문자열 전체에 정규식을 걸면 'not CpG injected' 가 'CpG injected' 와 같은 점수를 받는다.
# 어제 non-irradiated 사고와 같은 유형이므로 성분 단위로 판정한다.
_DIRTY_WORD = re.compile(r"gamma|cobalt|x-ray|irradiat|corticosterone|injected|"
                         r"treated|mKO|knockout", re.I)
_NEGATED = re.compile(r"^\s*(not|non[- ]|sham|un)\b", re.I)


def dirtiness(background_terms):
    """배경 성분 중 실제 처치에 해당하는 것의 개수.
    'not CpG injected', 'non-irradiated', 'sham irradiation' 은 처치가 아니다."""
    n = 0
    for t in background_terms:
        if _NEGATED.search(t.strip()):
            continue
        if _DIRTY_WORD.search(t):
            n += 1
    return n
DE_RE = re.compile(r"differential_expression.*\.csv$", re.I)


def get(url, tries=3, timeout=2400):
    for k in range(tries):
        try:
            r = S.get(url, timeout=timeout)
            if r.status_code == 200:
                return r
        except Exception:
            pass
        time.sleep(2 + 3 * k)
    return None


def de_path(osd):
    """DE 테이블 확보. 없으면 내려받는다(상한 MAX_DE_MB)."""
    d = os.path.join(RAW, f"OSD-{osd}")
    os.makedirs(d, exist_ok=True)
    loc = [f for f in os.listdir(d) if DE_RE.search(f)]
    if loc:
        loc.sort(key=lambda f: os.path.getsize(os.path.join(d, f)))
        return os.path.join(d, loc[0]), ""
    r = get(f"https://osdr.nasa.gov/osdr/data/osd/files/{osd}", timeout=180)
    if r is None:
        return None, "files API 실패"
    st = r.json()["studies"].get(f"OSD-{osd}")
    if not st:
        return None, "스터디 없음"
    cand = sorted([f for f in st["study_files"] if DE_RE.search(f["file_name"])],
                  key=lambda f: f["file_size"])
    if not cand:
        return None, "DE 테이블 없음(플랫폼상 미제공)"
    f = cand[0]
    mb = f["file_size"] / 1e6
    if mb > MAX_DE_MB:
        return None, f"DE 파일 용량초과 {mb:.0f}MB"
    print(f"    DE 다운로드 중 {f['file_name'][:52]} ({mb:.0f}MB)", flush=True)
    rr = get("https://osdr.nasa.gov" + f["remote_url"])
    if rr is None:
        return None, "DE 다운로드 실패"
    p = os.path.join(d, f["file_name"])
    with open(p, "wb") as fh:
        fh.write(rr.content)
    return p, ""


def split_terms(side):
    return [t.strip() for t in side.split("&") if t.strip()]


def classify(terms):
    grav = [t for t in terms if GRAV_TOKEN.search(t)]
    back = [t for t in terms if not GRAV_TOKEN.search(t)]
    return grav, back


def find_contrast(hdr):
    """중력 이외 배경이 완전히 일치하는 (처리)v(대조) 대비를 찾는다.
    반환: dict 또는 None, 그리고 검토한 후보 전체 로그."""
    cands, log = [], []
    for c in hdr:
        m = re.match(r"Log2fc_\((.*?)\)v\((.*?)\)$", c)
        if not m:
            continue
        for A, B, flip in ((m.group(1), m.group(2), False),
                           (m.group(2), m.group(1), True)):
            ga, ba = classify(split_terms(A))
            gb, bb = classify(split_terms(B))
            reason = None
            if not ga or not gb:
                reason = "중력 성분 없음"
            elif not any(TREAT_TOKEN.search(t) for t in ga):
                reason = "처리측이 중력 감소 조건 아님"
            elif not any(CTRL_TOKEN.search(t) for t in gb):
                reason = "대조측이 정상 중력 조건 아님"
            elif any(TREAT_TOKEN.search(t) for t in gb):
                reason = "대조측도 중력 감소 조건"
            elif DISQUALIFY.search(A) or DISQUALIFY.search(B):
                reason = "재하중 등 추가 처치 포함"
            elif RESTORED_G.search(A):
                reason = "처리측이 궤도 원심분리로 1G 복원 (중력 감소 아님)"
            elif sorted(x.lower() for x in ba) != sorted(x.lower() for x in bb):
                reason = f"배경 불일치 {sorted(ba)} vs {sorted(bb)}"
            if reason is None:
                # 배경이 양쪽에서 일치하더라도 '무엇이 일치하는가' 가 중요하다.
                # 방사선 조사 개체 안에서의 HLU 효과는 중력 단독 효과가 아니다.
                # 또 부분중력(0.33G)보다 완전 미세중력(uG)이 대비로서 깨끗하다.
                dirty = dirtiness(ba)
                partial = 0 if any(re.search(r"^\s*u\s*g\s*$|micrograv", t, re.I)
                                   for t in ga) else 1
                # 지상 HLU 는 uG 표기가 없으므로 partial 패널티를 주지 않는다
                if not any(re.search(r"space ?flight", t, re.I) for t in ga):
                    partial = 0
                cands.append({"col": c, "flip": flip, "treat": A, "ctrl": B,
                              "n_background": len(ba), "dirty_bg": dirty,
                              "partial_g": partial})
            else:
                log.append({"col": c, "dir": "rev" if flip else "fwd", "reason": reason})
    if not cands:
        return None, log
    # 우선순위: 배경 오염 적음 > 완전 미세중력 > 배경 성분 수 적음
    cands.sort(key=lambda d: (d["dirty_bg"], d["partial_g"], d["n_background"]))
    return cands[0], log


def analyse(osd, tissue, para):
    p, err = de_path(osd)
    base = {"조직": tissue, "OSD_ID": f"OSD-{osd}", "처리": para}
    if p is None:
        return {**base, "채택": False, "제외사유": err}
    hdr = pd.read_csv(p, nrows=0).columns.tolist()
    sym = next((c for c in hdr if c.upper() == "SYMBOL"), None)
    ch, rejlog = find_contrast(hdr)
    if sym is None:
        return {**base, "채택": False, "제외사유": "SYMBOL 컬럼 없음"}
    if ch is None:
        top = "; ".join(f"{r['reason']}" for r in rejlog[:3])
        return {**base, "채택": False,
                "제외사유": f"대비 불가 (후보 {len(rejlog)}개 전부 탈락: {top})"}
    ca = ch["col"].replace("Log2fc_", "Adj.p.value_")
    use = [sym, ch["col"]] + ([ca] if ca in hdr else [])
    df = pd.read_csv(p, usecols=use, low_memory=False)
    df.columns = ["SYMBOL", "log2fc"] + (["fdr"] if ca in hdr else [])
    if "fdr" not in df.columns:
        df["fdr"] = np.nan
    if ch["flip"]:
        df["log2fc"] = -df["log2fc"]
    df = df.dropna(subset=["log2fc"])
    n_deg = int((df.fdr < 0.05).sum())
    ck = df[df.SYMBOL.isin(CLOCK)]
    all_med = float(df.log2fc.abs().median())
    ck_med = float(ck.log2fc.abs().median()) if len(ck) else np.nan
    i4 = ck_med / all_med if all_med > 0 else np.nan
    rec = {**base,
           "대비_처리": ch["treat"][:70], "대비_대조": ch["ctrl"][:70],
           "부호반전": ch["flip"], "n_유전자": len(df), "DEG": n_deg,
           "clock_n": len(ck),
           "전역_absfc_중앙값": round(all_med, 4),
           "clock_absfc_중앙값": round(ck_med, 4) if ck_med == ck_med else np.nan,
           "I4": round(i4, 3) if i4 == i4 else np.nan,
           "clock_백분위": round(float((df.log2fc.abs() < ck_med).mean()), 4)
                           if ck_med == ck_med else np.nan,
           "Bmal1_log2fc": round(float(ck[ck.SYMBOL == "Bmal1"].log2fc.iloc[0]), 3)
                           if (ck.SYMBOL == "Bmal1").any() else np.nan,
           "채택": n_deg >= DEG_MIN,
           "제외사유": "" if n_deg >= DEG_MIN else f"신호 없음 (DEG {n_deg} < {DEG_MIN})"}
    return rec


def main():
    rows, rejlogs = [], []
    for tissue, osd, para in TARGETS:
        print(f"[{tissue}] OSD-{osd} ({para})", flush=True)
        r = analyse(osd, tissue, para)
        rows.append(r)
        msg = (f"    I4={r.get('I4')} DEG={r.get('DEG')} 채택={r.get('채택')}"
               if r.get("채택") is not None else "")
        print(f"    {r.get('제외사유') or msg}", flush=True)
        if r.get("대비_처리"):
            print(f"    ({r['대비_처리'][:56]}) v ({r['대비_대조'][:52]})", flush=True)

    D = pd.DataFrame(rows)
    D.to_csv(os.path.join(TAB, "v2_per_dataset.csv"), index=False, encoding="utf-8-sig")

    print("\n" + "=" * 100)
    print("데이터셋별 결과")
    print("=" * 100)
    cols = ["조직", "OSD_ID", "처리", "DEG", "I4", "clock_백분위", "Bmal1_log2fc", "채택", "제외사유"]
    print(D[cols].to_string(index=False, max_colwidth=46))

    ok = D[D.채택 == True]
    print(f"\n채택 {len(ok)} / 전체 {len(D)}")

    print("\n" + "=" * 100)
    print("조직 내 비교 (FLIGHT vs HLU) — 평균 내지 않고 개별 보고")
    print("=" * 100)
    summ = []
    for tis, g in ok.groupby("조직", sort=False):
        f = g[g.처리 == "FLIGHT"]; h = g[g.처리 == "HLU"]
        print(f"\n--- {tis} ---")
        for _, r in g.sort_values("처리").iterrows():
            print(f"    {r.처리:<7} {r.OSD_ID:<9} I4={r.I4:<6} DEG={r.DEG:<6} "
                  f"clock백분위={r.clock_백분위}")
        if len(f) and len(h):
            verdict = ("FLIGHT>HLU" if f.I4.min() > h.I4.max()
                       else "HLU>FLIGHT" if h.I4.min() > f.I4.max() else "겹침")
            print(f"    → FLIGHT {f.I4.min():.2f}~{f.I4.max():.2f} vs "
                  f"HLU {h.I4.min():.2f}~{h.I4.max():.2f}  [{verdict}]")
            summ.append({"조직": tis, "n_FLIGHT": len(f), "n_HLU": len(h),
                         "FLIGHT_I4_min": round(f.I4.min(), 3),
                         "FLIGHT_I4_max": round(f.I4.max(), 3),
                         "FLIGHT_I4_median": round(f.I4.median(), 3),
                         "HLU_I4_min": round(h.I4.min(), 3),
                         "HLU_I4_max": round(h.I4.max(), 3),
                         "HLU_I4_median": round(h.I4.median(), 3),
                         "판정": verdict})
        else:
            miss = "HLU 없음" if not len(h) else "FLIGHT 없음"
            print(f"    → 비교 불가 ({miss})")
            summ.append({"조직": tis, "n_FLIGHT": len(f), "n_HLU": len(h), "판정": f"비교 불가({miss})"})

    Sm = pd.DataFrame(summ)
    Sm.to_csv(os.path.join(TAB, "v2_tissue_summary.csv"), index=False, encoding="utf-8-sig")
    print("\n" + "=" * 100)
    print(Sm.to_string(index=False))
    print(f"\n-> {TAB}")


if __name__ == "__main__":
    main()
