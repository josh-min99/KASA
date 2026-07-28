"""
Stage 1 — Life 2020 (doi:10.3390/life10090196, PMC7555136) 재현.

원논문 방법 (전문 §2.1, §3.3 에서 확인. data/raw/_refs/life2020_methods.txt 참조)
  - 데이터: GLDS-98, -99, -101, -102, -103, -104, -105, -168 (RR-1 동일 미션 8조직)
  - 조직: adrenal gland, EDL, gastrocnemius, kidney, quadriceps, soleus,
          tibialis anterior, liver
  - n: 각 5~6 FLT / 5~6 GC
  - 지표: GeneLab 이 제공하는 `rna_seq_differential_expression.csv` 의
          조직별 log2FC 와 FDR 을 그대로 사용 (FDR < 0.05)
  - 비동기화의 정의: Arntl(Bmal1) 과 Per2 는 역위상이므로 함께 움직여야 하는데,
          Arntl 은 전 조직에서 일관되게 상향인 반면
          Per2 는 근육에서만 하향이고 부신·간에서는 유의하지 않다.

재현 판정 기준 (사전 규정)
  A. Arntl 이 8조직 중 몇 곳에서 상향인가 (원논문: 전부)
  B. Per2 가 근육 조직에서 하향인가 (원논문: 근육 전부)
  C. Per2 가 부신·간에서 비유의인가 (원논문: 그렇다)
  A·B·C 모두 일치 -> 재현 성공 / 방향만 일치 -> 부분 / 그 외 -> 실패

출력: data/processed/life2020_clock_fc.csv, results/tables/stage1_replication.csv
"""
import os, re, time, warnings
import requests
import numpy as np
import pandas as pd

warnings.filterwarnings("ignore")
S = requests.Session()
S.headers.update({"User-Agent": "KASA-drylab/1.0"})

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RAW = os.path.join(ROOT, "data", "raw")
PROC = os.path.join(ROOT, "data", "processed")
TAB = os.path.join(ROOT, "results", "tables")
for d in (RAW, PROC, TAB):
    os.makedirs(d, exist_ok=True)

# 원논문이 명시한 8개 데이터셋. GLDS-N == OSD-N (RR-1 계열)
LIFE2020 = {98: "adrenal gland", 99: "EDL", 101: "gastrocnemius", 102: "kidney",
            103: "quadriceps", 104: "soleus", 105: "tibialis anterior", 168: "liver"}
MUSCLE = {"EDL", "gastrocnemius", "quadriceps", "soleus", "tibialis anterior"}

CLOCK = ["Arntl", "Bmal1", "Clock", "Npas2", "Per1", "Per2", "Per3",
         "Cry1", "Cry2", "Nr1d1", "Nr1d2", "Dbp", "Tef", "Hlf", "Ciart",
         "Bhlhe40", "Bhlhe41"]
FDR_TH = 0.05


def get(url, tries=3, timeout=1800):
    for k in range(tries):
        try:
            r = S.get(url, timeout=timeout)
            if r.status_code == 200:
                return r
        except Exception:
            pass
        time.sleep(2 + 3 * k)
    return None


def fetch_de(osd):
    """DE 테이블 경로 반환. 없으면 None."""
    d = os.path.join(RAW, f"OSD-{osd}")
    os.makedirs(d, exist_ok=True)
    local = [f for f in os.listdir(d) if re.search(r"differential_expression", f, re.I)]
    if local:
        return os.path.join(d, local[0])
    r = get(f"https://osdr.nasa.gov/osdr/data/osd/files/{osd}", timeout=180)
    if r is None:
        return None
    st = r.json()["studies"].get(f"OSD-{osd}")
    if not st:
        return None
    cand = [f for f in st["study_files"]
            if re.search(r"differential_expression.*\.csv$", f["file_name"], re.I)]
    if not cand:
        return None
    cand.sort(key=lambda f: f["file_size"])           # 작은 것부터
    f = cand[0]
    p = os.path.join(d, f["file_name"])
    rr = get("https://osdr.nasa.gov" + f["remote_url"])
    if rr is None:
        return None
    with open(p, "wb") as fh:
        fh.write(rr.content)
    return p


def extract(path, tissue, osd):
    hdr = pd.read_csv(path, nrows=0).columns.tolist()
    sym = next((c for c in hdr if c.upper() == "SYMBOL"), None)
    lfc = [c for c in hdr if c.startswith("Log2fc_")]
    adj = [c for c in hdr if c.startswith("Adj.p.value_")]
    if not (sym and lfc and adj):
        return None, {"osd": osd, "tissue": tissue, "error": f"컬럼 불일치: sym={sym} lfc={len(lfc)} adj={len(adj)}"}
    # Space Flight vs Ground Control 대비 선택
    def pick(cols):
        for c in cols:
            if re.search(r"space ?flight", c, re.I) and re.search(r"ground control", c, re.I):
                return c
        return cols[0]
    cl, ca = pick(lfc), pick(adj)
    df = pd.read_csv(path, usecols=[sym, cl, ca], low_memory=False)
    df.columns = ["SYMBOL", "log2fc", "fdr"]
    df = df[df.SYMBOL.isin(CLOCK)].dropna(subset=["log2fc"])
    # 대비 방향 확인: (Space Flight)v(Ground Control) 이면 FLT/GC
    flt_first = bool(re.match(r".*\(Space ?Flight[^)]*\)v\(Ground", cl, re.I))
    if not flt_first:
        df["log2fc"] = -df["log2fc"]
    df["tissue"] = tissue
    df["OSD_ID"] = f"OSD-{osd}"
    df["contrast_col"] = cl
    df["flipped"] = not flt_first
    return df, {"osd": osd, "tissue": tissue, "contrast": cl, "flipped": not flt_first,
                "n_clock_found": len(df)}


def main():
    frames, log = [], []
    for osd, tissue in LIFE2020.items():
        p = fetch_de(osd)
        if p is None:
            log.append({"osd": osd, "tissue": tissue, "error": "DE 파일 없음/다운로드 실패"})
            print(f"[OSD-{osd}] {tissue}: 실패", flush=True)
            continue
        df, meta = extract(p, tissue, osd)
        log.append(meta)
        if df is None:
            print(f"[OSD-{osd}] {tissue}: {meta['error']}", flush=True)
            continue
        frames.append(df)
        print(f"[OSD-{osd}] {tissue}: clock {len(df)}개 / 대비 {meta['contrast'][:70]}"
              f"{' [부호반전]' if meta['flipped'] else ''}", flush=True)

    if not frames:
        print("\n### 재현 불가 — DE 테이블을 하나도 확보하지 못했다.")
        pd.DataFrame(log).to_csv(os.path.join(TAB, "stage1_download_log.csv"),
                                 index=False, encoding="utf-8-sig")
        return

    A = pd.concat(frames, ignore_index=True)
    A.to_csv(os.path.join(PROC, "life2020_clock_fc.csv"), index=False, encoding="utf-8-sig")
    pd.DataFrame(log).to_csv(os.path.join(TAB, "stage1_download_log.csv"),
                             index=False, encoding="utf-8-sig")

    print(f"\n확보 조직 {A.tissue.nunique()}/8: {sorted(A.tissue.unique())}")

    piv = A.pivot_table(index="SYMBOL", columns="tissue", values="log2fc")
    pfdr = A.pivot_table(index="SYMBOL", columns="tissue", values="fdr")
    piv.round(3).to_csv(os.path.join(TAB, "stage1_clock_log2fc.csv"), encoding="utf-8-sig")
    pfdr.round(4).to_csv(os.path.join(TAB, "stage1_clock_fdr.csv"), encoding="utf-8-sig")

    print("\n=== clock gene log2FC (Space Flight vs Ground Control) ===")
    print(piv.round(2).to_string())
    print("\n=== FDR ===")
    print(pfdr.round(3).to_string())

    # ---------- 사전 규정한 재현 판정
    print("\n" + "=" * 90)
    print("재현 판정")
    print("=" * 90)
    res = {}

    def row(g):
        return (piv.loc[g] if g in piv.index else None,
                pfdr.loc[g] if g in pfdr.index else None)

    bm = "Arntl" if "Arntl" in piv.index else ("Bmal1" if "Bmal1" in piv.index else None)
    if bm:
        v, f = row(bm)
        up = (v > 0)
        res["A_Arntl_up_all"] = f"{int(up.sum())}/{int(v.notna().sum())} 조직에서 상향"
        print(f"A. {bm} 상향: {res['A_Arntl_up_all']}   (원논문: 전 조직 상향)")
        print(f"   조직별 log2FC: {v.round(2).to_dict()}")
    else:
        res["A_Arntl_up_all"] = "Arntl/Bmal1 없음"
        print("A. Arntl/Bmal1 을 DE 테이블에서 찾지 못함")

    if "Per2" in piv.index:
        v, f = row("Per2")
        mus = [t for t in v.index if t in MUSCLE]
        non = [t for t in v.index if t not in MUSCLE]
        dn_m = (v[mus] < 0)
        res["B_Per2_down_muscle"] = f"{int(dn_m.sum())}/{len(mus)} 근육 조직에서 하향"
        print(f"\nB. Per2 근육 하향: {res['B_Per2_down_muscle']}   (원논문: 근육 전부 하향)")
        print(f"   근육 log2FC: {v[mus].round(2).to_dict()}")
        sig_non = (f[non] < FDR_TH)
        res["C_Per2_ns_nonmuscle"] = (f"비근육 {len(non)}곳 중 유의 {int(sig_non.sum())}곳 "
                                      f"(FDR<{FDR_TH})")
        print(f"\nC. Per2 비근육 비유의: {res['C_Per2_ns_nonmuscle']}   (원논문: 부신·간 비유의)")
        print(f"   비근육 log2FC: {v[non].round(2).to_dict()}")
        print(f"   비근육 FDR   : {f[non].round(3).to_dict()}")
    else:
        res["B_Per2_down_muscle"] = res["C_Per2_ns_nonmuscle"] = "Per2 없음"
        print("\nPer2 를 DE 테이블에서 찾지 못함")

    pd.DataFrame([res]).to_csv(os.path.join(TAB, "stage1_replication.csv"),
                               index=False, encoding="utf-8-sig")
    print(f"\n-> {TAB}/stage1_replication.csv")


if __name__ == "__main__":
    main()
