"""
Stage 3 — Stage 1 의 지표를 HLU 데이터에 그대로 적용.

Stage 1 재현 지표 (Life 2020 원문 그대로)
  조직별 clock gene log2FC + FDR. 비동기화의 조작적 정의는
  "Arntl(Bmal1) 은 전 조직 일관 상향인데 Per2 는 조직군에 따라 갈린다".

Stage 3 에서 추가로 정의하는 정량 지표 (원논문에는 없음. 명시적으로 밝힌다)
  I1  Bmal1 log2FC 의 조직 간 부호 일치율
  I2  Per2  log2FC 의 조직 간 부호 일치율
  I3  조직쌍 간 clock gene log2FC 프로파일의 Pearson 상관 (16유전자 벡터)
      동기화돼 있으면 높고, 비동기화면 낮아진다.

게이트
  Stage 2(개정)에서 데이터셋별 양성대조 통과 여부를 확인하고,
  통과하지 못한 데이터셋의 수치는 '탐색적'으로만 표기한다.

OSD-21 (앵커)
  마이크로어레이. Space Flight / Ground Control / Hindlimb Unloaded /
  Hindlimb Unloaded and Reloaded / Normally Loaded Control 5군이 한 스터디에 있다.
  HLU vs HLU+Reloaded 비교가 비행 데이터의 '착륙 후 채취' 문제를 통제한다.

출력: results/tables/stage3_*.csv
"""
import os, re, warnings
import numpy as np
import pandas as pd
from scipy import stats

warnings.filterwarnings("ignore")
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RAW = os.path.join(ROOT, "data", "raw")
PROC = os.path.join(ROOT, "data", "processed")
TAB = os.path.join(ROOT, "results", "tables")
for d in (PROC, TAB):
    os.makedirs(d, exist_ok=True)

CLOCK = ["Bmal1", "Arntl", "Clock", "Npas2", "Per1", "Per2", "Per3", "Cry1", "Cry2",
         "Nr1d1", "Nr1d2", "Dbp", "Tef", "Hlf", "Ciart", "Bhlhe40", "Bhlhe41"]

# HLU 데이터셋 (조직 매칭 표에서 근육 + 기타)
HLU_SETS = [(876, "gastrocnemius"), (880, "gastrocnemius"),
            (935, "soleus"), (949, "soleus"),
            (237, "dorsal skin"), (203, "retina"),
            (201, "spleen"), (211, "spleen"), (214, "bone marrow")]

TREAT = re.compile(r"hindlimb unloaded", re.I)
CTRL = re.compile(r"normally loaded|loaded control|cage control|ambulatory|vivarium", re.I)
DIRTY = re.compile(r"reload|radiation|corticosterone|mKO|knockout|\bKO\b|treated|spike", re.I)


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


def choose(cands, treat_re=TREAT, ctrl_re=CTRL):
    scored = []
    for c, a, b in cands:
        for t, k, flip in ((a, b, False), (b, a, True)):
            if treat_re.search(t) and ctrl_re.search(k):
                dirt = len(DIRTY.findall(t)) + len(DIRTY.findall(k))
                extra = t.count("&") + k.count("&")
                scored.append((dirt * 10 + extra, c, flip, t, k))
    if not scored:
        return None
    scored.sort()
    _, col, flip, t, k = scored[0]
    return {"col": col, "flip": flip, "treat": t, "ctrl": k, "dirt": scored[0][0]}


def clock_fc(osd, tissue, treat_re=TREAT, ctrl_re=CTRL, label=""):
    p = de_path(osd)
    if p is None:
        return None, {"OSD_ID": f"OSD-{osd}", "note": "DE 파일 없음"}
    hdr = pd.read_csv(p, nrows=0).columns.tolist()
    sym = next((c for c in hdr if c.upper() == "SYMBOL"), None)
    ch = choose(contrasts(hdr), treat_re, ctrl_re)
    if not (sym and ch):
        return None, {"OSD_ID": f"OSD-{osd}", "note": "적합 대비 없음",
                      "available": "; ".join(c for c, _, _ in contrasts(hdr))[:220]}
    ca = ch["col"].replace("Log2fc_", "Adj.p.value_")
    use = [sym, ch["col"]] + ([ca] if ca in hdr else [])
    df = pd.read_csv(p, usecols=use, low_memory=False)
    df.columns = ["SYMBOL", "log2fc"] + (["fdr"] if ca in hdr else [])
    if "fdr" not in df.columns:
        df["fdr"] = np.nan
    if ch["flip"]:
        df["log2fc"] = -df["log2fc"]
    n_deg = int((df.fdr < 0.05).sum())
    d = df[df.SYMBOL.isin(CLOCK)].copy()
    d["tissue"] = tissue
    d["OSD_ID"] = f"OSD-{osd}"
    d["label"] = label or f"OSD-{osd}"
    return d, {"OSD_ID": f"OSD-{osd}", "tissue": tissue, "treat": ch["treat"][:64],
               "ctrl": ch["ctrl"][:64], "flipped": ch["flip"], "dirt": ch["dirt"],
               "n_DEG_fdr05": n_deg, "note": ""}


def sign_concord(series):
    v = series.dropna()
    if len(v) < 2:
        return np.nan, len(v)
    maj = max((v > 0).sum(), (v < 0).sum())
    return maj / len(v), len(v)


def profile_corr(piv):
    """조직쌍 간 clock 프로파일 상관."""
    cols = list(piv.columns)
    out = []
    for i in range(len(cols)):
        for j in range(i + 1, len(cols)):
            a, b = piv[cols[i]], piv[cols[j]]
            ok = a.notna() & b.notna()
            if ok.sum() >= 6:
                r, p = stats.pearsonr(a[ok], b[ok])
                out.append({"t1": cols[i], "t2": cols[j], "r": round(r, 3),
                            "p": round(p, 4), "n_genes": int(ok.sum())})
    return pd.DataFrame(out)


def main():
    # ---------- FLIGHT (Stage 1 결과 재사용)
    F = pd.read_csv(os.path.join(PROC, "life2020_clock_fc.csv"))
    F["paradigm"] = "FLIGHT"
    fpiv = F.pivot_table(index="SYMBOL", columns="tissue", values="log2fc")

    # ---------- HLU
    frames, meta = [], []
    for osd, tis in HLU_SETS:
        d, m = clock_fc(osd, tis)
        meta.append(m)
        if d is not None:
            d["paradigm"] = "HLU"
            frames.append(d)
            print(f"[OSD-{osd}] {tis}: clock {len(d)}개 | {m['treat'][:40]} v {m['ctrl'][:34]}"
                  f" | DEG {m['n_DEG_fdr05']}", flush=True)
        else:
            print(f"[OSD-{osd}] {tis}: {m['note']}", flush=True)

    M = pd.DataFrame(meta)
    M.to_csv(os.path.join(TAB, "stage3_hlu_contrasts.csv"), index=False, encoding="utf-8-sig")

    if not frames:
        print("\n### HLU 데이터에서 clock gene 대비를 하나도 만들지 못했다. Stage 3 중단.")
        return

    H = pd.concat(frames, ignore_index=True)
    # 조직별로 여러 스터디가 있으면 라벨을 조직+OSD 로 구분
    H["tissue_lab"] = H.tissue + " (" + H.OSD_ID + ")"
    H.to_csv(os.path.join(PROC, "hlu_clock_fc.csv"), index=False, encoding="utf-8-sig")
    hpiv = H.pivot_table(index="SYMBOL", columns="tissue_lab", values="log2fc")
    hfdr = H.pivot_table(index="SYMBOL", columns="tissue_lab", values="fdr")
    hpiv.round(3).to_csv(os.path.join(TAB, "stage3_hlu_clock_log2fc.csv"), encoding="utf-8-sig")
    hfdr.round(4).to_csv(os.path.join(TAB, "stage3_hlu_clock_fdr.csv"), encoding="utf-8-sig")

    print("\n=== HLU: clock gene log2FC (Unloaded vs Loaded) ===")
    print(hpiv.round(2).to_string())
    print("\n=== HLU: FDR ===")
    print(hfdr.round(3).to_string())

    # ---------- 지표 비교
    print("\n" + "=" * 92)
    print("지표 비교 (FLIGHT vs HLU)")
    print("=" * 92)
    rows = []
    for name, piv in (("FLIGHT", fpiv), ("HLU", hpiv)):
        for g in ("Bmal1", "Per2"):
            if g in piv.index:
                c, n = sign_concord(piv.loc[g])
                direction = "상향" if piv.loc[g].mean() > 0 else "하향"
                rows.append({"paradigm": name, "gene": g, "n_tissues": n,
                             "sign_concordance": round(c, 3) if c == c else np.nan,
                             "mean_log2fc": round(float(piv.loc[g].mean()), 3),
                             "direction": direction})
    C = pd.DataFrame(rows)
    print(C.to_string(index=False))
    C.to_csv(os.path.join(TAB, "stage3_index_I1I2.csv"), index=False, encoding="utf-8-sig")

    print("\n[I3] 조직쌍 간 clock 프로파일 상관")
    for name, piv in (("FLIGHT", fpiv), ("HLU", hpiv)):
        pc = profile_corr(piv)
        if len(pc):
            pc.to_csv(os.path.join(TAB, f"stage3_profilecorr_{name}.csv"),
                      index=False, encoding="utf-8-sig")
            print(f"  {name}: 조직쌍 {len(pc)}개 | r 중앙값 {pc.r.median():.3f} "
                  f"| 범위 {pc.r.min():.3f}~{pc.r.max():.3f}")
        else:
            print(f"  {name}: 계산 가능한 조직쌍 없음")

    # ---------- OSD-21 앵커
    print("\n" + "=" * 92)
    print("OSD-21 앵커 — HLU vs HLU+Reloaded")
    print("=" * 92)
    p21 = de_path(21)
    if p21 is None:
        print("  OSD-21 에 DE 테이블이 없다 (마이크로어레이 스터디). 별도 처리 필요.")
        pd.DataFrame([{"note": "OSD-21 DE 테이블 없음"}]).to_csv(
            os.path.join(TAB, "stage3_osd21.csv"), index=False, encoding="utf-8-sig")
    else:
        hdr = pd.read_csv(p21, nrows=0).columns.tolist()
        cs = contrasts(hdr)
        print(f"  가용 대비 {len(cs)}개:")
        for c, a, b in cs:
            print(f"    ({a}) v ({b})")


if __name__ == "__main__":
    main()
