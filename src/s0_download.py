"""
Stage 0 — 데이터 확보.

docs/proposal_draft.md §3.2(앵커 OSD-21) 및 §3.3(조직 매칭 6종)의 스터디를 내려받는다.
각 스터디마다 발현행렬 + 샘플 메타데이터(ISA) + 프로토콜(i_Investigation)을 확보한다.

규칙: 다운로드 실패는 실패로 기록. 재시도 3회. 합성 데이터 금지.
출력: data/raw/<OSD-nnn>/ , data_inventory.md
"""
import os, re, io, json, time, zipfile, warnings
import requests
import pandas as pd

warnings.filterwarnings("ignore")
S = requests.Session()
S.headers.update({"User-Agent": "KASA-drylab/1.0"})

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RAW = os.path.join(ROOT, "data", "raw")
os.makedirs(RAW, exist_ok=True)

# §3.2 앵커 + §3.3 조직 매칭
TARGETS = {
    21:  ("gastrocnemius", "ANCHOR"),
    # 등쪽 피부
    237: ("dorsal skin", "HLU"), 238: ("dorsal skin", "FLIGHT"),
    240: ("dorsal skin", "FLIGHT"), 243: ("dorsal skin", "FLIGHT"),
    254: ("dorsal skin", "FLIGHT"),
    # 망막
    203: ("retina", "HLU"), 87: ("retina", "FLIGHT"), 194: ("retina", "FLIGHT"),
    255: ("retina", "FLIGHT"), 758: ("retina", "FLIGHT"),
    # 비복근
    876: ("gastrocnemius", "HLU"), 880: ("gastrocnemius", "HLU"),
    101: ("gastrocnemius", "FLIGHT"), 419: ("gastrocnemius", "FLIGHT"),
    # 가자미근
    935: ("soleus", "HLU"), 949: ("soleus", "HLU"),
    104: ("soleus", "FLIGHT"), 714: ("soleus", "FLIGHT"), 770: ("soleus", "FLIGHT"),
    # 비장
    201: ("spleen", "HLU"), 211: ("spleen", "HLU"),
    246: ("spleen", "FLIGHT"), 288: ("spleen", "FLIGHT"), 506: ("spleen", "FLIGHT"),
    # 골수
    214: ("bone marrow", "HLU"), 690: ("bone marrow", "FLIGHT"),
}

EXPR_PAT = re.compile(
    r"_rna_seq_Normalized_Counts_GLbulkRNAseq\.csv$"
    r"|_array_normalized_expression_probeset_GLmicroarray\.csv$"
    r"|_rna_seq_differential_expression_GLbulkRNAseq\.csv$"
    r"|_array_differential_expression_GLmicroarray\.csv$", re.I)
ISA_PAT = re.compile(r"ISA\.zip$", re.I)
MAX_EXPR_MB = 120          # 초대형 DE 파일은 건너뛴다


def get(url, tries=3, **kw):
    for k in range(tries):
        try:
            r = S.get(url, timeout=kw.pop("timeout", 180), **kw)
            if r.status_code == 200:
                return r
        except Exception:
            pass
        time.sleep(2 + 3 * k)
    return None


def download(osd):
    d = os.path.join(RAW, f"OSD-{osd}")
    os.makedirs(d, exist_ok=True)
    rec = {"OSD_ID": f"OSD-{osd}", "isa": False, "expr_files": [], "errors": []}

    r = get(f"https://osdr.nasa.gov/osdr/data/osd/files/{osd}")
    if r is None:
        rec["errors"].append("files API 실패"); return rec
    try:
        st = r.json()["studies"].get(f"OSD-{osd}")
    except Exception as e:
        rec["errors"].append(f"files JSON 파싱 실패: {e}"); return rec
    if not st:
        rec["errors"].append("스터디 없음"); return rec
    files = st["study_files"]

    for f in files:
        want_isa = ISA_PAT.search(f["file_name"])
        want_expr = EXPR_PAT.search(f["file_name"])
        if not (want_isa or want_expr):
            continue
        if want_expr and f["file_size"] / 1e6 > MAX_EXPR_MB:
            rec["errors"].append(f"용량초과 건너뜀: {f['file_name']} "
                                 f"({f['file_size']/1e6:.0f}MB)")
            continue
        p = os.path.join(d, f["file_name"])
        if not os.path.exists(p) or os.path.getsize(p) == 0:
            rr = get("https://osdr.nasa.gov" + f["remote_url"], timeout=1800)
            if rr is None:
                rec["errors"].append(f"다운로드 실패: {f['file_name']}"); continue
            with open(p, "wb") as fh:
                fh.write(rr.content)
        if want_isa:
            rec["isa"] = True
        else:
            rec["expr_files"].append(f["file_name"])
    return rec


def parse_isa(osd):
    """샘플 테이블에서 군·n·플랫폼·희생시각 기록여부를 추출."""
    d = os.path.join(RAW, f"OSD-{osd}")
    zs = [x for x in os.listdir(d) if x.endswith("ISA.zip")] if os.path.isdir(d) else []
    if not zs:
        return None
    z = zipfile.ZipFile(os.path.join(d, zs[0]))
    out = {"groups": {}, "n": 0, "time_recorded": False, "time_cols": [],
           "assays": [], "tissue_col": ""}
    CLOCK = re.compile(r"\b\d{1,2}:\d{2}\b|\b\d{1,2}\s*(?:am|pm)\b|\bZT\s*\d|\bCT\s*\d", re.I)
    NOISE = re.compile(r"\b12\s*:\s*12\b|\b1\s*:\s*\d{2,}\b|light/?dark|light-dark", re.I)
    for nm in z.namelist():
        base = nm.split("/")[-1]
        if base.startswith("a_"):
            out["assays"].append(base)
        if not base.startswith("s_"):
            continue
        df = pd.read_csv(io.BytesIO(z.read(nm)), sep="\t", dtype=str)
        out["n"] = len(df)
        fv = [c for c in df.columns if c.startswith("Factor Value")]
        if fv:
            key = df[fv].fillna("NA").agg(" | ".join, axis=1)
            out["groups"] = key.value_counts().to_dict()
        for c in df.columns:
            if c.startswith(("Term ", "Unit")):
                continue
            for v in df[c].dropna().unique():
                sv = str(v)
                if CLOCK.search(sv) and not NOISE.search(sv):
                    out["time_recorded"] = True
                    out["time_cols"].append(c)
                    break
        # 프로토콜 저장
        for inm in z.namelist():
            if inm.split("/")[-1].startswith("i_"):
                with open(os.path.join(d, "protocol_investigation.txt"), "wb") as fh:
                    fh.write(z.read(inm))
        df.to_csv(os.path.join(d, "sample_table.tsv"), sep="\t", index=False)
    out["time_cols"] = sorted(set(out["time_cols"]))
    return out


def main():
    rows, log = [], []
    order = [21] + [k for k in TARGETS if k != 21]
    for osd in order:
        tissue, para = TARGETS[osd]
        print(f"[OSD-{osd}] {tissue} / {para} ...", flush=True)
        rec = download(osd)
        isa = parse_isa(osd) if rec["isa"] else None
        plat = "RNA-seq" if any("rna_seq" in f for f in rec["expr_files"]) else \
               ("microarray" if any("array" in f for f in rec["expr_files"]) else "?")
        rows.append({
            "OSD_ID": rec["OSD_ID"], "tissue": tissue, "paradigm": para,
            "n_samples": isa["n"] if isa else 0,
            "groups": "; ".join(f"{k}={v}" for k, v in (isa["groups"] if isa else {}).items())[:220],
            "platform": plat,
            "expr_downloaded": len(rec["expr_files"]) > 0,
            "expr_files": "; ".join(rec["expr_files"])[:160],
            "isa_downloaded": rec["isa"],
            "time_of_day_recorded": bool(isa and isa["time_recorded"]),
            "time_cols": "; ".join(isa["time_cols"]) if isa else "",
            "errors": "; ".join(rec["errors"])[:200],
        })
        log.append(f"OSD-{osd}: isa={rec['isa']} expr={len(rec['expr_files'])} err={rec['errors']}")

    inv = pd.DataFrame(rows)
    inv.to_csv(os.path.join(ROOT, "data", "processed", "data_inventory.csv"),
               index=False, encoding="utf-8-sig")

    md = ["# data_inventory", "",
          f"생성: Stage 0 (`src/s0_download.py`). 총 {len(inv)}건 시도.", "",
          "| 스터디 | 조직 | 패러다임 | n | 플랫폼 | 희생시각 기록 | 발현행렬 | ISA | 비고 |",
          "|---|---|---|---|---|---|---|---|---|"]
    for _, r in inv.iterrows():
        md.append(f"| {r.OSD_ID} | {r.tissue} | {r.paradigm} | {r.n_samples} | {r.platform} "
                  f"| {'O' if r.time_of_day_recorded else 'X'} "
                  f"| {'O' if r.expr_downloaded else '**X**'} "
                  f"| {'O' if r.isa_downloaded else '**X**'} | {r.errors or ''} |")
    md += ["", "## 군 구성", ""]
    for _, r in inv.iterrows():
        if r.groups:
            md.append(f"- **{r.OSD_ID}** ({r.tissue}): {r.groups}")
    open(os.path.join(ROOT, "data_inventory.md"), "w", encoding="utf-8").write("\n".join(md))

    ok = inv[inv.expr_downloaded & inv.isa_downloaded]
    print(f"\n성공 {len(ok)}/{len(inv)}")
    print(f"OSD-21 확보: {bool(len(inv[(inv.OSD_ID=='OSD-21') & inv.expr_downloaded]))}")
    print("\n".join(log))


if __name__ == "__main__":
    main()
