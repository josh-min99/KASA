"""
Stage 3c — 배치 교란을 제거한 최종 매칭 비교.

s3b 에서 시도한 HLU-A(OSD-334~337)는 발현 데이터가 아예 없었다(miRNA/대사체 계열).
대신 진짜 매칭 세트를 찾았다.

  HLU-A : OSD-202(뇌) / 203(망막) / 211(비장) / 237(등쪽 피부)
          네 스터디가 동일한 2x2 설계를 공유한다.
            Factor Value[Hindlimb Unloading] : Hindlimb Unloaded / Normally Loaded Control
            Factor Value[Ionizing Radiation] : non-irradiated / cobalt-57 gamma radiation
          '비조사 HLU vs 비조사 정상하중' 대비를 쓰면 하중 변화만 분리된다.

비교 대상
  FLIGHT-A (RR-1, 8조직)  — Stage 1 에서 확보. 단일 미션
  HLU-A    (4조직)        — 단일 실험, 비조사 조건만

지표는 Stage 1 과 동일: 조직별 clock gene log2FC, 그리고
  I1/I2 Bmal1·Per2 의 조직 간 부호 일치율과 평균 크기
  I3    조직쌍 프로파일 상관

주의: DEG 개수가 0 에 가까운 데이터셋은 '신호 없음' 이지 '비동기화' 가 아니다.
      반드시 함께 보고한다.

출력: results/tables/stage3c_*.csv
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
HLU_A = {202: "brain", 203: "retina", 211: "spleen", 237: "dorsal skin"}

# 'non-irradiated' 는 조사군이 아니다. 뒤돌아보기로 'non-' 을 배제한다.
IRRADIATED = re.compile(r"gamma|cobalt|(?<!non-)irradiated", re.I)
DE_RE = re.compile(r"differential_expression.*\.csv$", re.I)


def de_path(osd):
    d = os.path.join(RAW, f"OSD-{osd}")
    os.makedirs(d, exist_ok=True)
    c = [f for f in os.listdir(d) if DE_RE.search(f)]
    if not c:
        import requests
        sess = requests.Session()
        sess.headers.update({"User-Agent": "KASA-drylab/1.0"})
        try:
            j = sess.get(f"https://osdr.nasa.gov/osdr/data/osd/files/{osd}",
                         timeout=180).json()
            fs = j["studies"][f"OSD-{osd}"]["study_files"]
            cand = sorted([f for f in fs if DE_RE.search(f["file_name"])],
                          key=lambda f: f["file_size"])
            if cand and cand[0]["file_size"] / 1e6 < 400:
                rr = sess.get("https://osdr.nasa.gov" + cand[0]["remote_url"], timeout=2400)
                if rr.status_code == 200:
                    with open(os.path.join(d, cand[0]["file_name"]), "wb") as fh:
                        fh.write(rr.content)
        except Exception:
            pass
        c = [f for f in os.listdir(d) if DE_RE.search(f)]
    if not c:
        return None
    c.sort(key=lambda f: os.path.getsize(os.path.join(d, f)))
    return os.path.join(d, c[0])


def contrasts(hdr):
    return [(c, m.group(1), m.group(2)) for c in hdr
            if (m := re.match(r"Log2fc_\((.*?)\)v\((.*?)\)$", c))]


def pick_pure_hlu(cands):
    """비조사 HLU vs 비조사 정상하중. 시간 요인이 있으면 가장 짧은 것 하나."""
    best = None
    for c, a, b in cands:
        for t, k, flip in ((a, b, False), (b, a, True)):
            if not (re.search(r"hindlimb unloaded", t, re.I)
                    and re.search(r"normally loaded", k, re.I)):
                continue
            # 조사군 제외.
            # 주의: 'irradiat' 로 거르면 'non-irradiated' 까지 함께 지워진다.
            # 초기 구현이 그 때문에 필요한 대비를 전부 잃었다.
            # 'non-' 이 앞에 붙지 않은 조사 표기만 배제한다.
            if IRRADIATED.search(t) or IRRADIATED.search(k):
                continue
            score = t.count("&") + k.count("&")
            if best is None or score < best[0]:
                best = (score, c, flip, t, k)
    if best is None:
        return None
    return {"col": best[1], "flip": best[2], "treat": best[3], "ctrl": best[4]}


def clock_fc(osd, tissue):
    p = de_path(osd)
    if p is None:
        return None, {"OSD_ID": f"OSD-{osd}", "tissue": tissue, "note": "DE 파일 없음"}
    hdr = pd.read_csv(p, nrows=0).columns.tolist()
    sym = next((c for c in hdr if c.upper() == "SYMBOL"), None)
    ch = pick_pure_hlu(contrasts(hdr))
    if not (sym and ch):
        return None, {"OSD_ID": f"OSD-{osd}", "tissue": tissue,
                      "note": "순수 HLU 대비 없음",
                      "available": "; ".join(c for c, _, _ in contrasts(hdr))[:300]}
    ca = ch["col"].replace("Log2fc_", "Adj.p.value_")
    use = [sym, ch["col"]] + ([ca] if ca in hdr else [])
    df = pd.read_csv(p, usecols=use, low_memory=False)
    df.columns = ["SYMBOL", "log2fc"] + (["fdr"] if ca in hdr else [])
    if "fdr" not in df.columns:
        df["fdr"] = np.nan
    if ch["flip"]:
        df["log2fc"] = -df["log2fc"]
    d = df[df.SYMBOL.isin(CLOCK)].copy()
    d["tissue"] = tissue
    d["OSD_ID"] = f"OSD-{osd}"
    return d, {"OSD_ID": f"OSD-{osd}", "tissue": tissue, "treat": ch["treat"][:70],
               "ctrl": ch["ctrl"][:70], "flipped": ch["flip"],
               "n_DEG_fdr05": int((df.fdr < 0.05).sum()), "note": ""}


def coherence(piv, label):
    out = {"series": label, "n_tissues": piv.shape[1]}
    for g in ("Bmal1", "Per2"):
        if g in piv.index:
            v = piv.loc[g].dropna()
            if len(v) >= 2:
                out[f"{g}_mean"] = round(float(v.mean()), 3)
                out[f"{g}_absmean"] = round(float(v.abs().mean()), 3)
                out[f"{g}_concord"] = round(float(max((v > 0).sum(), (v < 0).sum()) / len(v)), 3)
    cols = list(piv.columns); rs = []
    for i in range(len(cols)):
        for j in range(i + 1, len(cols)):
            a, b = piv[cols[i]], piv[cols[j]]
            ok = a.notna() & b.notna()
            if ok.sum() >= 6:
                rs.append(stats.pearsonr(a[ok], b[ok])[0])
    if rs:
        out.update({"profile_r_median": round(float(np.median(rs)), 3),
                    "profile_r_min": round(float(np.min(rs)), 3),
                    "profile_r_max": round(float(np.max(rs)), 3),
                    "n_pairs": len(rs)})
    # 전체 clock 반응 크기
    out["clock_absmean_log2fc"] = round(float(piv.abs().mean().mean()), 3)
    return out


def main():
    frames, meta = [], []
    for osd, tis in HLU_A.items():
        d, m = clock_fc(osd, tis)
        meta.append(m)
        if d is not None:
            frames.append(d)
            print(f"[OSD-{osd}] {tis}: clock {len(d)}개 | DEG {m['n_DEG_fdr05']}")
            print(f"          ({m['treat']}) v ({m['ctrl']})")
        else:
            print(f"[OSD-{osd}] {tis}: {m['note']}")
            if "available" in m:
                print(f"          가용: {m['available'][:280]}")

    M = pd.DataFrame(meta)
    M.to_csv(os.path.join(TAB, "stage3c_contrasts.csv"), index=False, encoding="utf-8-sig")
    if not frames:
        print("\n### 순수 HLU 대비를 하나도 만들지 못했다.")
        return

    H = pd.concat(frames, ignore_index=True)
    H.to_csv(os.path.join(PROC, "hlu_matched_clock_fc.csv"), index=False, encoding="utf-8-sig")
    hpiv = H.pivot_table(index="SYMBOL", columns="tissue", values="log2fc")
    hfdr = H.pivot_table(index="SYMBOL", columns="tissue", values="fdr")
    hpiv.round(3).to_csv(os.path.join(TAB, "stage3c_hlu_log2fc.csv"), encoding="utf-8-sig")
    hfdr.round(4).to_csv(os.path.join(TAB, "stage3c_hlu_fdr.csv"), encoding="utf-8-sig")

    F = pd.read_csv(os.path.join(PROC, "life2020_clock_fc.csv"))
    fpiv = F.pivot_table(index="SYMBOL", columns="tissue", values="log2fc")

    print("\n=== HLU-A: clock gene log2FC (비조사 HLU vs 비조사 정상하중) ===")
    print(hpiv.round(2).to_string())
    print("\n=== HLU-A: FDR ===")
    print(hfdr.round(3).to_string())

    S = pd.DataFrame([coherence(fpiv, "FLIGHT-A (RR-1, 8조직, 단일미션)"),
                      coherence(hpiv, "HLU-A (4조직, 단일실험)")])
    S.to_csv(os.path.join(TAB, "stage3c_coherence.csv"), index=False, encoding="utf-8-sig")
    print("\n" + "=" * 96)
    print("단일 실험 내 조직 간 clock 반응 — 같은 단위 비교")
    print("=" * 96)
    print(S.to_string(index=False))

    print("\n데이터셋별 DEG 수 (신호 유무 확인):")
    print(M[["OSD_ID", "tissue", "n_DEG_fdr05"]].to_string(index=False))


if __name__ == "__main__":
    main()
