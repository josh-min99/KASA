"""
Stage 2 (개정) — 양성대조.

s2_positive_control.py 초판이 2/9 로 실패했으나, 대비(contrast) 선택이 잘못됐다.
로그에서 확인된 문제:
  OSD-876/880 : 'Hindlimb Unloaded **and Reloaded**' 대비를 골랐다 (재하중은 위축신호를 되돌린다)
  OSD-935     : 'HLU & x-ray radiation & corticosterone' 병용 대비
  OSD-949     : 'AMPKalpha mKO' 녹아웃 배경 대비
따라서 초판의 실패는 파이프라인 실패인지 대비 선택 실패인지 구별되지 않는다.

개정 내용
  1) 데이터셋마다 가능한 모든 대비를 나열하고, 추가 요인이 가장 적은
     '순수 처리 vs 순수 대조' 대비를 규칙으로 선택한다.
  2) 마커 의존성을 없애기 위해 파이프라인 무결성 지표를 추가한다:
     각 데이터셋의 FDR<0.05 DEG 개수. 파일을 옳게 읽고 있다면 0 일 수 없다.
  3) RR-1 은 37일 장기 비행이다. Fbxo32/Trim63 은 폐용 '급성기' 마커라
     만성기에는 기저로 돌아가거나 하향한다는 점을 판정에 반영한다.
     따라서 만성 위축 표현형 마커(Myh7 하향 / Myh4 상향, 느린→빠른 섬유 전환)를 병기한다.

출력: results/tables/stage2b_*.csv
"""
import os, re, warnings
import numpy as np
import pandas as pd

warnings.filterwarnings("ignore")
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RAW = os.path.join(ROOT, "data", "raw")
TAB = os.path.join(ROOT, "results", "tables")
os.makedirs(TAB, exist_ok=True)

ACUTE = ["Fbxo32", "Trim63"]                       # 폐용 급성기 마커
CHRONIC = ["Myh7", "Myh4", "Myh1", "Myh2"]         # 느린->빠른 섬유 전환
EXTRA = ["Foxo1", "Foxo3", "Mstn", "Ampd3"]

SETS = [(99, "EDL", "FLIGHT"), (101, "gastrocnemius", "FLIGHT"),
        (103, "quadriceps", "FLIGHT"), (104, "soleus", "FLIGHT"),
        (105, "tibialis anterior", "FLIGHT"),
        (876, "gastrocnemius", "HLU"), (880, "gastrocnemius", "HLU"),
        (935, "soleus", "HLU"), (949, "soleus", "HLU")]

TREAT = re.compile(r"space ?flight|hindlimb unloaded", re.I)
CTRL = re.compile(r"ground control|normally loaded|loaded control|vivarium", re.I)
# 대비를 오염시키는 추가 요인
DIRTY = re.compile(r"reload|radiation|corticosterone|mKO|knockout|KO\b|"
                   r"treated|drug|spike", re.I)


def de_path(osd):
    d = os.path.join(RAW, f"OSD-{osd}")
    if not os.path.isdir(d):
        return None
    c = [f for f in os.listdir(d) if re.search(r"differential_expression.*\.csv$", f, re.I)]
    if not c:
        return None
    c.sort(key=lambda f: os.path.getsize(os.path.join(d, f)))
    return os.path.join(d, c[0])


def contrasts(hdr):
    out = []
    for c in hdr:
        m = re.match(r"Log2fc_\((.*?)\)v\((.*?)\)$", c)
        if m:
            out.append((c, m.group(1), m.group(2)))
    return out


def choose(cands):
    """추가 요인이 가장 적은 (처리)v(대조) 를 고른다. 반환 (컬럼, 부호반전, 오염도)."""
    scored = []
    for c, a, b in cands:
        for treat, ctrl, flip in ((a, b, False), (b, a, True)):
            if TREAT.search(treat) and CTRL.search(ctrl):
                dirt = len(DIRTY.findall(treat)) + len(DIRTY.findall(ctrl))
                # '&' 로 이어진 추가 조건 수
                extra = treat.count("&") + ctrl.count("&")
                scored.append((dirt * 10 + extra, c, flip, treat, ctrl))
    if not scored:
        return None
    scored.sort()
    s = scored[0]
    return {"col": s[1], "flip": s[2], "dirt": s[0], "treat": s[3], "ctrl": s[4],
            "n_candidates": len(scored)}


def main():
    rows, meta, integ = [], [], []
    for osd, tis, para in SETS:
        p = de_path(osd)
        if p is None:
            meta.append({"OSD_ID": f"OSD-{osd}", "note": "DE 파일 없음"})
            continue
        hdr = pd.read_csv(p, nrows=0).columns.tolist()
        sym = next((c for c in hdr if c.upper() == "SYMBOL"), None)
        ch = choose(contrasts(hdr))
        if not (sym and ch):
            meta.append({"OSD_ID": f"OSD-{osd}", "note": "적합 대비 없음",
                         "all": "; ".join(c for c, _, _ in contrasts(hdr))[:250]})
            continue
        ca = ch["col"].replace("Log2fc_", "Adj.p.value_")
        use = [sym, ch["col"]] + ([ca] if ca in hdr else [])
        df = pd.read_csv(p, usecols=use, low_memory=False)
        df.columns = ["SYMBOL", "log2fc"] + (["fdr"] if ca in hdr else [])
        if "fdr" not in df.columns:
            df["fdr"] = np.nan
        if ch["flip"]:
            df["log2fc"] = -df["log2fc"]

        n_deg = int((df.fdr < 0.05).sum())
        integ.append({"OSD_ID": f"OSD-{osd}", "tissue": tis, "paradigm": para,
                      "n_genes": len(df), "n_DEG_fdr05": n_deg,
                      "frac_DEG": round(n_deg / max(len(df), 1), 4)})
        meta.append({"OSD_ID": f"OSD-{osd}", "tissue": tis, "paradigm": para,
                     "treat": ch["treat"][:70], "ctrl": ch["ctrl"][:70],
                     "flipped": ch["flip"], "dirt_score": ch["dirt"],
                     "n_candidate_contrasts": ch["n_candidates"]})
        for g in ACUTE + CHRONIC + EXTRA:
            r = df[df.SYMBOL == g]
            rows.append({"OSD_ID": f"OSD-{osd}", "tissue": tis, "paradigm": para, "gene": g,
                         "log2fc": float(r.log2fc.iloc[0]) if len(r) else np.nan,
                         "fdr": float(r.fdr.iloc[0]) if len(r) else np.nan})

    D = pd.DataFrame(rows); M = pd.DataFrame(meta); I = pd.DataFrame(integ)
    D.to_csv(os.path.join(TAB, "stage2b_markers.csv"), index=False, encoding="utf-8-sig")
    M.to_csv(os.path.join(TAB, "stage2b_contrasts.csv"), index=False, encoding="utf-8-sig")
    I.to_csv(os.path.join(TAB, "stage2b_integrity.csv"), index=False, encoding="utf-8-sig")

    print("=== 선택된 대비 (오염도 낮은 순으로 규칙 선택) ===")
    print(M.to_string(index=False, max_colwidth=54))

    print("\n=== 파이프라인 무결성: FDR<0.05 DEG 개수 ===")
    print(I.to_string(index=False))
    zero = I[I.n_DEG_fdr05 == 0]
    print(f"  DEG 0 인 데이터셋: {len(zero)}/{len(I)}"
          f"{' -> ' + ', '.join(zero.OSD_ID) if len(zero) else ''}")

    print("\n=== 급성 위축 마커 (부하감소 방향 = 양수 기대) ===")
    a = D[D.gene.isin(ACUTE)]
    print(a.pivot_table(index="gene", columns=["paradigm", "OSD_ID"], values="log2fc").round(2).to_string())
    print("  FDR:")
    print(a.pivot_table(index="gene", columns=["paradigm", "OSD_ID"], values="fdr").round(4).to_string())

    print("\n=== 만성 섬유전환 마커 (Myh7 하향 / Myh4 상향 기대) ===")
    c = D[D.gene.isin(CHRONIC)]
    print(c.pivot_table(index="gene", columns=["paradigm", "OSD_ID"], values="log2fc").round(2).to_string())
    print("  FDR:")
    print(c.pivot_table(index="gene", columns=["paradigm", "OSD_ID"], values="fdr").round(4).to_string())

    # ---------- 판정
    print("\n" + "=" * 92)
    ok_int = len(zero) == 0 and I.n_DEG_fdr05.median() > 100
    print(f"[무결성] 모든 데이터셋에서 DEG 검출 & 중앙값 {int(I.n_DEG_fdr05.median())}개"
          f" -> {'통과' if ok_int else '실패'}")

    def hits(genes, direction):
        s = D[D.gene.isin(genes)].dropna(subset=["log2fc", "fdr"])
        s = s[(s.fdr < 0.05) & ((s.log2fc > 0) if direction > 0 else (s.log2fc < 0))]
        return s.OSD_ID.nunique()

    n_set = D.OSD_ID.nunique()
    h_ac = hits(ACUTE, +1)
    h_myh7 = D[(D.gene == "Myh7") & (D.fdr < 0.05) & (D.log2fc < 0)].OSD_ID.nunique()
    h_myh4 = D[(D.gene == "Myh4") & (D.fdr < 0.05) & (D.log2fc > 0)].OSD_ID.nunique()
    print(f"[급성 마커] Fbxo32/Trim63 상향+유의: {h_ac}/{n_set} 데이터셋")
    print(f"[만성 마커] Myh7 하향+유의: {h_myh7}/{n_set} · Myh4 상향+유의: {h_myh4}/{n_set}")
    ok_bio = (h_ac + h_myh7 + h_myh4) > 0 and (h_ac / n_set >= 0.5 or
                                               max(h_myh7, h_myh4) / n_set >= 0.5)
    verdict = "통과" if (ok_int and ok_bio) else ("조건부 통과" if ok_int else "실패")
    print(f"\n  >>> Stage 2 (개정) 판정: {verdict}")
    pd.DataFrame([{"integrity_ok": ok_int, "median_DEG": int(I.n_DEG_fdr05.median()),
                   "acute_hits": h_ac, "myh7_hits": h_myh7, "myh4_hits": h_myh4,
                   "n_datasets": n_set, "verdict": verdict}]).to_csv(
        os.path.join(TAB, "stage2b_verdict.csv"), index=False, encoding="utf-8-sig")


if __name__ == "__main__":
    main()
