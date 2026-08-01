"""
R6: 역문제 폐쇄 — 웻랩이 Δφ 를 재 오면 기존 비행 스냅샷이 해석 가능해지는가.

문제 구조
  v1/v2 에서 확정한 사실: 비행 전사체는 단일 시점이라 **진폭 변화와 위상 이동을 구별하지
  못한다**. 그래서 OSDR 자료만으로는 중력의 clock 효과를 판정할 수 없다.
  그런데 이것은 '데이터가 나쁘다' 가 아니라 '미지수가 하나 더 있다' 는 뜻이다.
  Δφ 를 밖에서 넣어 주면 남는 미지수가 줄어 스냅샷이 해석 가능해진다.

여기서 하는 것
  1. 위상만 Δφ 만큼 이동했을 때 단일 시점 clock 유전자 벡터가 어떻게 보이는지 순방향 계산
  2. 그 예측 크기를 v2 에서 실측한 잡음 바닥(results/v2/noise_floor.txt)과 대조해
     **검출 가능한 Δφ 하한**을 구한다
  3. 반대로 RR-1 / OSD-21 에서 관측된 벡터 크기가 어떤 Δφ 구간과 정합적인지 역추정한다
  4. 사전등록형 예측 진술을 생성한다

모델
  clock 유전자 g 의 단일시점 log2 발현을  x_g(t) = M_g + A_g*cos(w(t - phi_g))  로 둔다.
  위상만 Δφ 이동하면 관측되는 log2FC 는
      d_g(Δφ) = A_g * [cos(w(t - phi_g - Δφ)) - cos(w(t - phi_g))]
  진폭 A_g 는 지상 아틀라스(GSE54650) 실측값을 쓴다 — `scripts/power_analysis.py` 와 동일 출처.
  희생 시각 t 는 OSDR 에 기록이 없으므로(v1 확정 사실) t 를 모르는 채로 다뤄야 한다.
  따라서 t 를 0-24 전 구간에 균등하게 두고 |d_g| 의 분포를 구한다.

산출: data/rhythm/inverse_dphi.csv, results/v3/r6_inverse.txt
"""
import os
import re
import sys
import gzip
import warnings

import numpy as np
import pandas as pd

warnings.filterwarnings("ignore")
if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
CACHE = os.path.join(ROOT, "cache")
OUT = os.path.join(ROOT, "data", "rhythm")
RES = os.path.join(ROOT, "results", "v3")
os.makedirs(RES, exist_ok=True)

OMEGA = 2 * np.pi / 24
# 아틀라스는 12개 조직 x 24시점이다. 조직마다 위상이 달라 섞으면 진폭이 뭉개진다.
# 간(Liv)으로 한정한다 — 비행 전사체에서 clock 분석이 가장 많이 이뤄진 조직이고
# v1 의 Life 2020 재현에도 포함돼 있다.
TISSUE = "Liv"
CLOCK = ["Arntl", "Bmal1", "Npas2", "Clock", "Per1", "Per2", "Per3",
         "Cry1", "Cry2", "Nr1d1", "Nr1d2", "Dbp", "Tef", "Hlf", "Ciart", "Bhlhe41"]

# v2 에서 실측한 전사체 잡음 바닥 (results/v2/noise_floor.txt)
NOISE_FLOOR = {
    "HLU+재하중 vs 정상하중": 0.063,
    "우주비행 vs 지상대조": 0.088,
    "HLU vs 정상하중": 0.120,
    "HLU+재하중 vs HLU": 0.139,
}
OBSERVED = {  # 같은 표의 관측 |차이| 중앙값
    "HLU+재하중 vs 정상하중": 0.224,
    "우주비행 vs 지상대조": 0.215,
    "HLU vs 정상하중": 0.182,
    "HLU+재하중 vs HLU": 0.138,
}


def atlas_amplitudes():
    """지상 아틀라스(GSE54650)에서 clock 유전자 log2 진폭·위상 실측.
    `scripts/phase_predictor.py` 가 받아 둔 캐시를 그대로 쓴다."""
    p = os.path.join(CACHE, "GSE54650_series_matrix.txt.gz")
    ann = os.path.join(CACHE, "GPL6246.annot.gz")
    if not (os.path.exists(p) and os.path.exists(ann)):
        return None

    # 주의: GPL6246.annot 의 열 순서는 ID / Gene title / Gene symbol / ... 이다.
    # 초판에서 f[1](Gene title)을 심볼로 읽어 매칭이 전멸했고 대체값으로 넘어갔었다.
    # 헤더에서 'Gene symbol' 열 위치를 찾아 쓴다.
    sym = {}
    with gzip.open(ann, "rt", errors="replace") as fh:
        started, gcol = False, None
        for line in fh:
            if line.startswith("!platform_table_begin"):
                hdr = next(fh).rstrip("\n").split("\t")
                gcol = hdr.index("Gene symbol") if "Gene symbol" in hdr else 2
                started = True
                continue
            if line.startswith("!platform_table_end"):
                break
            if started:
                f = line.rstrip("\n").split("\t")
                if len(f) > gcol:
                    sym[f[0]] = f[gcol]

    rows, samples = [], None
    with gzip.open(p, "rt", errors="replace") as fh:
        for line in fh:
            if line.startswith("!Sample_title"):
                samples = [x.strip('"') for x in line.rstrip("\n").split("\t")[1:]]
            if line.startswith("!series_matrix_table_begin"):
                hdr = next(fh)
                for l2 in fh:
                    if l2.startswith("!series_matrix_table_end"):
                        break
                    rows.append(l2.rstrip("\n").split("\t"))
                break
    if not rows or samples is None:
        return None

    # 샘플명은 'Adr_CT18' 형태이고 12개 조직 x 24시점 = 288개다.
    # 조직마다 위상이 다르므로 섞으면 안 된다. 간(Liv)으로 한정한다.
    ct = []
    for s in samples:
        m = re.match(r"^(\w+?)_CT(\d+)", s)
        ct.append(float(m.group(2)) if (m and m.group(1) == TISSUE) else np.nan)
    ct = np.array(ct, float)
    ok = np.isfinite(ct)
    if ok.sum() < 8:
        return None
    t = ct[ok] % 24

    X = np.column_stack([np.ones(len(t)), np.cos(OMEGA * t), np.sin(OMEGA * t)])
    out = []
    for r in rows:
        pid = r[0].strip('"')
        g = sym.get(pid, "")
        if g not in CLOCK:
            continue
        try:
            y = np.array([float(v) for v in r[1:]], float)[ok]
        except ValueError:
            continue
        if not np.isfinite(y).all():
            continue
        # GSE54650 series matrix 는 선형 강도값이다(중앙값 수백~수천).
        # 비행 대비는 log2FC 단위이므로 log2 로 맞춰야 한다.
        # 초판에서 이 변환을 빠뜨려 진폭이 344 로 나왔었다.
        if np.nanmedian(y) > 32:
            y = np.log2(np.clip(y, 1, None))
        beta, *_ = np.linalg.lstsq(X, y, rcond=None)
        amp = float(np.hypot(beta[1], beta[2]))
        phi = float((np.arctan2(beta[2], beta[1]) / OMEGA) % 24)
        out.append(dict(gene=g, amplitude=amp, phase=phi))
    if not out:
        return None
    df = pd.DataFrame(out).groupby("gene").agg(amplitude=("amplitude", "median"),
                                              phase=("phase", "median")).reset_index()
    return df


# 아틀라스를 못 읽을 때 쓰는 대체값 (마우스 간 clock 유전자 log2 진폭의 통상 범위)
FALLBACK = pd.DataFrame({
    "gene": CLOCK,
    "amplitude": [0.90, 0.90, 0.45, 0.35, 1.10, 1.05, 0.85, 0.45, 0.55,
                  1.30, 1.00, 1.40, 0.80, 0.75, 1.10, 0.70],
    "phase": [23.0, 23.0, 22.0, 22.5, 12.0, 16.0, 15.0, 18.0, 20.0,
              6.0, 8.0, 10.0, 10.5, 11.0, 9.5, 10.0],
})


def predicted_shift_magnitude(amp, phi, dphi, n_t=240):
    """희생 시각을 모르므로 t 를 0-24 균등으로 두고 |d_g| 의 유전자 중앙값 분포를 구한다."""
    ts = np.linspace(0, 24, n_t, endpoint=False)
    med = []
    for t in ts:
        d = amp * (np.cos(OMEGA * (t - phi - dphi)) - np.cos(OMEGA * (t - phi)))
        med.append(np.median(np.abs(d)))
    return np.array(med)


def main():
    log = []
    P = log.append
    P("=" * 92)
    P("R6  역문제 — 웻랩의 Δφ 가 기존 비행 스냅샷을 해석 가능하게 만드는가")
    P("=" * 92)

    at = atlas_amplitudes()
    if at is None or len(at) < 6:
        at = FALLBACK
        src = "대체값 (아틀라스 캐시 없음)"
    else:
        src = f"GSE54650 아틀라스 실측 ({len(at)} 유전자)"
    P(f"\nclock 유전자 진폭 출처: {src}")
    P(f"  log2 진폭 중앙값 {at.amplitude.median():.3f} (범위 {at.amplitude.min():.3f}~{at.amplitude.max():.3f})")

    amp, phi = at.amplitude.values, at.phase.values

    # ---------------------------------------------- 1) 순방향
    P("")
    P("[1] 위상만 Δφ 이동했을 때 단일 시점에서 보이는 |log2FC| (유전자 중앙값)")
    P("    희생 시각이 기록돼 있지 않으므로 t 를 0-24h 균등으로 두고 분포로 제시한다.")
    P("")
    P(f"  {'Δφ(h)':>7s} {'중앙값':>9s} {'25%':>9s} {'75%':>9s} {'최소':>9s} {'최대':>9s}")
    rows = []
    for dphi in [0.5, 1, 1.5, 2, 3, 4, 6, 8, 12]:
        m = predicted_shift_magnitude(amp, phi, dphi)
        rows.append(dict(dphi=dphi, med=np.median(m), q25=np.percentile(m, 25),
                         q75=np.percentile(m, 75), lo=m.min(), hi=m.max()))
        P(f"  {dphi:7.1f} {np.median(m):9.3f} {np.percentile(m, 25):9.3f} "
          f"{np.percentile(m, 75):9.3f} {m.min():9.3f} {m.max():9.3f}")
    fwd = pd.DataFrame(rows)

    # ---------------------------------------------- 2) 검출 가능한 Δφ 하한
    P("")
    P("[2] 검출 가능한 Δφ 하한 — v2 에서 실측한 전사체 잡음 바닥과 대조")
    P("    잡음 바닥은 results/v2/noise_floor.txt (CI 폭에서 SE 를 역산한 반정규 중앙값)")
    P("")
    P(f"  {'대비':26s} {'잡음바닥':>9s} {'관측':>9s} {'검출가능 Δφ 하한':>18s}")
    fine = np.arange(0.1, 12.01, 0.05)
    med_by_dphi = np.array([np.median(predicted_shift_magnitude(amp, phi, d, n_t=120))
                            for d in fine])
    lim_rows = []
    for k, floor in NOISE_FLOOR.items():
        idx = np.where(med_by_dphi >= floor)[0]
        dmin = float(fine[idx[0]]) if len(idx) else np.nan
        idx2 = np.where(med_by_dphi >= OBSERVED[k])[0]
        dobs = float(fine[idx2[0]]) if len(idx2) else np.nan
        lim_rows.append(dict(contrast=k, noise_floor=floor, observed=OBSERVED[k],
                             dphi_min_detectable=dmin, dphi_matching_observed=dobs))
        P(f"  {k:26s} {floor:9.3f} {OBSERVED[k]:9.3f} {dmin:18.2f}")
    lim = pd.DataFrame(lim_rows)

    # ---------------------------------------------- 3) 역추정
    P("")
    P("[3] 관측된 벡터 크기와 정합적인 Δφ (위상 이동만으로 설명한다면)")
    P("")
    P(f"  {'대비':26s} {'관측 |차이|':>10s} {'정합 Δφ(h)':>12s} {'잡음대비':>9s}")
    for _, r in lim.iterrows():
        mult = r.observed / r.noise_floor
        P(f"  {r.contrast:26s} {r.observed:10.3f} {r.dphi_matching_observed:12.2f} {mult:8.2f}x")
    P("")
    P("  주의 — 이 표는 '위상 이동만으로 설명한다면' 이라는 가정 아래의 값이다.")
    P("  진폭 변화가 섞이면 같은 관측이 다른 Δφ 와도 정합적이다. 그 축퇴가 바로")
    P("  단일 시점 자료의 한계이고, 웻랩이 Δφ 를 독립적으로 재 와야 하는 이유다.")

    # ---------------------------------------------- 4) 사전등록형 예측
    P("")
    P("=" * 92)
    P("[4] 사전등록형 예측 진술 (웻랩 결과가 나오기 전에 확정해 둘 것)")
    P("=" * 92)
    d_lo = float(lim.dphi_min_detectable.min())
    P("")
    P(f"  P1. 웻랩에서 측정된 Δφ 가 {d_lo:.1f} h 미만이면,")
    P(f"      그 크기의 위상 이동은 기존 비행 전사체의 잡음 바닥 아래에 있다.")
    P(f"      -> 비행 스냅샷에서 clock 신호가 안 잡히는 것은 '효과 없음' 의 증거가 아니라")
    P(f"         분해능의 한계임이 **정량적으로** 확정된다.")
    P("")
    P(f"  P2. 웻랩 Δφ 가 {d_lo:.1f} h 이상이면, 비행 대비의 clock |log2FC| 중앙값은")
    P(f"      위 [1] 표의 해당 행 이상이어야 한다. 실측(0.215)이 그보다 작으면")
    P(f"      비행에서는 지상 아날로그만큼의 위상 이동이 일어나지 않았다는 뜻이고,")
    P(f"      그 차이가 곧 '중력 이외 요인의 기여분' 의 하한이 된다.")
    P("")
    P("  P3. 어느 쪽이 나오든 결론이 나온다. 이것이 웻랩 없이는 닫히지 않는 고리다.")

    fwd.to_csv(os.path.join(OUT, "inverse_forward.csv"), index=False, encoding="utf-8-sig")
    lim.to_csv(os.path.join(OUT, "inverse_dphi.csv"), index=False, encoding="utf-8-sig")

    txt = "\n".join(log)
    print(txt)
    with open(os.path.join(RES, "r6_inverse.txt"), "w", encoding="utf-8") as fh:
        fh.write(txt + "\n")
    print(f"\n-> data/rhythm/inverse_dphi.csv, inverse_forward.csv, results/v3/r6_inverse.txt")


if __name__ == "__main__":
    main()
