"""
Stage 3d — OSD-21 앵커. HLU vs HLU+Reloaded 비교 (명세 §Stage3 필수 항목).

OSD-21 (STS-108, 마이크로어레이) 은 한 스터디에 5군이 있다.
  Space Flight / Ground Control / Hindlimb Unloaded /
  Hindlimb Unloaded and Reloaded / Normally Loaded Control
Reloaded 군은 원저자가 '착륙부터 희생까지의 3.5시간' 을 모사하려고 넣은 것으로,
비행 데이터의 '착륙 후 채취' 문제를 통제한다.

문제: OSD-21 에는 GeneLab DE 테이블이 없다(구형 마이크로어레이).
따라서 정규화 발현행렬 + 샘플 테이블에서 직접 군 평균차를 계산한다.
n=3~5 이므로 p 값 단독 주장은 하지 않는다. 효과크기와 개체값을 함께 낸다.

출력: results/tables/stage3d_osd21_*.csv
"""
import os, re, warnings
import numpy as np
import pandas as pd
from scipy import stats

warnings.filterwarnings("ignore")
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RAW = os.path.join(ROOT, "data", "raw", "OSD-21")
TAB = os.path.join(ROOT, "results", "tables")
PROC = os.path.join(ROOT, "data", "processed")
for d in (TAB, PROC):
    os.makedirs(d, exist_ok=True)

CLOCK = ["Bmal1", "Arntl", "Clock", "Npas2", "Per1", "Per2", "Per3", "Cry1", "Cry2",
         "Nr1d1", "Nr1d2", "Dbp", "Tef", "Hlf", "Ciart", "Bhlhe40", "Bhlhe41"]


def main():
    expr_f = [f for f in os.listdir(RAW) if "normalized_expression" in f]
    if not expr_f:
        print("### OSD-21 정규화 발현행렬 없음. 앵커 분석 불가.")
        pd.DataFrame([{"note": "발현행렬 없음"}]).to_csv(
            os.path.join(TAB, "stage3d_osd21_status.csv"), index=False, encoding="utf-8-sig")
        return
    E = pd.read_csv(os.path.join(RAW, expr_f[0]), low_memory=False)
    st = pd.read_csv(os.path.join(RAW, "sample_table.tsv"), sep="\t", dtype=str)

    print(f"발현행렬 {E.shape} | 샘플표 {st.shape}")
    sym_col = next((c for c in E.columns[:14] if re.search(r"symbol|gene", c, re.I)), None)
    print(f"심볼 컬럼: {sym_col}")
    print(f"발현행렬 앞 컬럼: {list(E.columns[:10])}")

    # 샘플명 매칭
    samp = [c for c in E.columns if c in set(st["Sample Name"])]
    if not samp:
        # 부분일치 시도
        names = list(st["Sample Name"])
        samp = [c for c in E.columns if any(n in c or c in n for n in names)]
    print(f"매칭된 샘플 열: {len(samp)}/{len(st)}")
    if not samp or sym_col is None:
        print("### 발현행렬과 샘플표를 연결하지 못했다. 앵커 분석 불가.")
        pd.DataFrame([{"note": "샘플명 매칭 실패",
                       "expr_cols": "; ".join(map(str, E.columns[:12])),
                       "sample_names": "; ".join(names[:8])}]).to_csv(
            os.path.join(TAB, "stage3d_osd21_status.csv"), index=False, encoding="utf-8-sig")
        return

    # 군 라벨
    fv = [c for c in st.columns if c.startswith("Factor Value")]
    st["grp"] = st[fv].fillna("NA").agg(" | ".join, axis=1)
    lab = st.set_index("Sample Name")["grp"].to_dict()

    def which(s):
        g = lab.get(s, "")
        if re.search(r"Hindlimb Unloaded and Reloaded", g, re.I):
            return "HLU+Reloaded"
        if re.search(r"Hindlimb Unloaded", g, re.I):
            return "HLU"
        if re.search(r"Normally Loaded", g, re.I):
            return "NormalLoaded"
        if re.search(r"Space Flight", g, re.I):
            return "Flight"
        if re.search(r"Ground Control", g, re.I):
            return "GroundControl"
        return "other"

    grp = pd.Series({s: which(s) for s in samp})
    print("\n군 구성:"); print(grp.value_counts().to_string())

    X = E[[sym_col] + samp].copy()
    X[sym_col] = X[sym_col].astype(str)
    X = X[X[sym_col].isin(CLOCK)]
    num = X[samp].apply(pd.to_numeric, errors="coerce")
    # 마이크로어레이 정규화값이 선형이면 log2
    if np.nanmedian(num.values) > 30:
        num = np.log2(num.clip(lower=1))
        print("  (선형 강도로 판단해 log2 변환)")
    X = pd.concat([X[[sym_col]].reset_index(drop=True), num.reset_index(drop=True)], axis=1)
    X = X.groupby(sym_col).mean()
    print(f"\nclock 유전자 {len(X)}개 확보: {sorted(X.index)}")

    COMPS = [("HLU", "NormalLoaded"), ("HLU+Reloaded", "NormalLoaded"),
             ("HLU+Reloaded", "HLU"), ("Flight", "GroundControl")]
    rows = []
    for a, b in COMPS:
        ca = [s for s in samp if grp[s] == a]
        cb = [s for s in samp if grp[s] == b]
        if len(ca) < 2 or len(cb) < 2:
            rows.append({"comparison": f"{a} vs {b}", "n_a": len(ca), "n_b": len(cb),
                         "note": "표본 부족"})
            continue
        for g in X.index:
            va, vb = X.loc[g, ca].astype(float), X.loc[g, cb].astype(float)
            d = float(va.mean() - vb.mean())
            sp = np.sqrt(((len(va) - 1) * va.var(ddof=1) + (len(vb) - 1) * vb.var(ddof=1))
                         / max(len(va) + len(vb) - 2, 1))
            se = sp * np.sqrt(1 / len(va) + 1 / len(vb)) if sp > 0 else np.nan
            tcrit = stats.t.ppf(0.975, len(va) + len(vb) - 2)
            rows.append({"comparison": f"{a} vs {b}", "gene": g,
                         "n_a": len(ca), "n_b": len(cb),
                         "diff_log2": round(d, 3),
                         "ci_low": round(d - tcrit * se, 3) if se == se else np.nan,
                         "ci_high": round(d + tcrit * se, 3) if se == se else np.nan,
                         "cohens_d": round(d / sp, 2) if sp > 0 else np.nan,
                         "p_uncorrected": round(float(stats.ttest_ind(va, vb,
                                                equal_var=True).pvalue), 4)})
    R = pd.DataFrame(rows)
    R.to_csv(os.path.join(TAB, "stage3d_osd21_clock.csv"), index=False, encoding="utf-8-sig")

    print("\n" + "=" * 96)
    print("OSD-21 clock gene 군간 차이 (log2, 95% CI)  — n=3~5, p는 참고치")
    print("=" * 96)
    for comp, g in R[R.gene.notna()].groupby("comparison"):
        print(f"\n--- {comp}  (n={g.n_a.iloc[0]} vs {g.n_b.iloc[0]}) ---")
        print(g[["gene", "diff_log2", "ci_low", "ci_high", "cohens_d",
                 "p_uncorrected"]].to_string(index=False))
        big = g[g.diff_log2.abs() > 0.5]
        print(f"    |차이|>0.5 인 유전자: {len(big)}/{len(g)}"
              f"{' -> ' + ', '.join(big.gene) if len(big) else ''}")

    # 개체별 산점도용 원자료
    long = []
    for g in X.index:
        for s in samp:
            long.append({"gene": g, "sample": s, "group": grp[s],
                         "value_log2": float(X.loc[g, s])})
    pd.DataFrame(long).to_csv(os.path.join(PROC, "osd21_clock_long.csv"),
                              index=False, encoding="utf-8-sig")
    print(f"\n-> {TAB}/stage3d_osd21_clock.csv, {PROC}/osd21_clock_long.csv")


if __name__ == "__main__":
    main()
