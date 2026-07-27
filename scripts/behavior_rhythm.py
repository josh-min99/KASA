"""
Week 2: 우주비행 중 행동 리듬 시계열 분석 (OSD-952 / RR-1).

왜 중요한가
  전사체는 단일시점이라 위상 정보를 주지 못한다(scripts/clock_dose_response.py 참조).
  OSDR 의 생체원격 데이터(OSD-681 등)는 기간 평균으로 요약돼 일주기 해상도가 없다.
  OSD-952 는 예외다 — ISS Rodent Habitat 영상을 DeepEthogram 으로 프레임 단위 분류한
  188만 행 데이터에 Day(0-32) 와 Light Cycle Phase(Light/Dark) 가 붙어 있다.
  즉 **미세중력에서 행동 리듬을 직접 측정할 수 있는 유일한 공개 데이터**다.

한계 (계획서에 반드시 명시)
  - RR-1 은 인공중력 군이 없다. 전부 uG. 따라서 중력 용량반응은 볼 수 없다.
  - 영상 표본이 시간대별로 균일하지 않다 → 노출시간 보정 필수.
  - 개체 식별이 안 되므로 케이지 단위 집계다.

출력: data/behavior_rhythm_daily.csv, data/behavior_rhythm_summary.csv
"""
import os, re, warnings
import requests
import numpy as np
import pandas as pd
from scipy import stats

warnings.filterwarnings("ignore")
S = requests.Session()
S.headers.update({"User-Agent": "KASA-circadian-gravity/1.0"})

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
CACHE = os.path.join(ROOT, "cache")
OUT = os.path.join(ROOT, "data")
os.makedirs(CACHE, exist_ok=True)
os.makedirs(OUT, exist_ok=True)

FNAME = "LSDS-168_Behavioral Segmentation_inferences_aligned.csv"
# 활동성 행동 / 비활동성 행동
ACTIVE = ["crawling", "eating", "grooming", "circling", "backflips", "drifting"]
INACTIVE = ["huddling"]


def get_file():
    path = os.path.join(CACHE, FNAME)
    if os.path.exists(path):
        return path
    fs = S.get("https://osdr.nasa.gov/osdr/data/osd/files/952", timeout=120).json()
    fs = fs["studies"]["OSD-952"]["study_files"]
    f = [x for x in fs if x["file_name"] == FNAME][0]
    with open(path, "wb") as fh:
        fh.write(S.get("https://osdr.nasa.gov" + f["remote_url"], timeout=3600).content)
    return path


def parse_clock(fn):
    """'269_18-02-15_1389_Feeder_1.1' -> 18:02:15 (녹화 시작 시각)."""
    m = re.match(r"^\d+_(\d{2})-(\d{2})-(\d{2})_", str(fn))
    if not m:
        return np.nan
    h, mi, s = map(int, m.groups())
    if h > 23 or mi > 59 or s > 59:
        return np.nan
    return h + mi / 60 + s / 3600


def main():
    df = pd.read_csv(get_file())
    print(f"원본 {df.shape[0]:,} 행 / 영상 {df.original_file.nunique()} 개 / Day {df.Day.min()}-{df.Day.max()}")

    # 프레임은 다중 라벨(한 프레임에 eating+crawling 동시). 프레임 단위로 재구성.
    df["frame_key"] = df.original_file.astype(str) + "#" + df.original_frame.astype(str)
    df["clock_h"] = df.original_file.map(parse_clock)
    print(f"파일명에서 시각 파싱 성공: {df.clock_h.notna().mean():.1%}")

    frames = df.drop_duplicates("frame_key")[["frame_key", "Day", "Light Cycle Phase",
                                              "original_file", "clock_h"]]
    print(f"고유 프레임 {len(frames):,}")

    act = df[df.predicted_behavior.isin(ACTIVE)].frame_key.unique()
    frames["active"] = frames.frame_key.isin(act)

    # ---------------- 1) 명/암 활동 대비
    print("\n[1] 명암주기별 활동 프레임 비율")
    t = frames.groupby("Light Cycle Phase")["active"].agg(["mean", "count"])
    print(t.round(4).to_string())
    dark = frames[frames["Light Cycle Phase"] == "Dark"].active
    light = frames[frames["Light Cycle Phase"] == "Light"].active
    chi = stats.chi2_contingency([[dark.sum(), len(dark) - dark.sum()],
                                  [light.sum(), len(light) - light.sum()]])
    print(f"    Dark {dark.mean():.4f} vs Light {light.mean():.4f}   chi2 p={chi.pvalue:.3e}")
    print("    (설치류는 야행성 → Dark 활동이 높아야 정상. 이것이 리듬 존재의 1차 증거)")

    # ---------------- 2) 행동별 명암 편향
    print("\n[2] 행동별 명암 편향 (Dark 비율 - Light 비율)")
    rows = []
    for b in sorted(df.predicted_behavior.unique()):
        k = set(df[df.predicted_behavior == b].frame_key)
        f = frames.assign(has=frames.frame_key.isin(k))
        d = f[f["Light Cycle Phase"] == "Dark"].has.mean()
        l = f[f["Light Cycle Phase"] == "Light"].has.mean()
        rows.append({"behavior": b, "dark_frac": round(d, 4), "light_frac": round(l, 4),
                     "dark_minus_light": round(d - l, 4)})
    bdf = pd.DataFrame(rows).sort_values("dark_minus_light", ascending=False)
    print(bdf.to_string(index=False))

    # ---------------- 3) 일자별 리듬 진폭 추이
    print("\n[3] 비행일수에 따른 리듬 진폭 (Dark-Light 활동 차이) 추이")
    daily = (frames.groupby(["Day", "Light Cycle Phase"])["active"]
             .agg(["mean", "count"]).reset_index())
    piv = daily.pivot(index="Day", columns="Light Cycle Phase", values="mean")
    cnt = daily.pivot(index="Day", columns="Light Cycle Phase", values="count")
    piv = piv.join(cnt, rsuffix="_n")
    # 양쪽 위상 모두 최소 500 프레임 있는 날만
    ok = piv.dropna()
    ok = ok[(ok.get("Dark_n", 0) >= 500) & (ok.get("Light_n", 0) >= 500)]
    ok["amplitude"] = ok["Dark"] - ok["Light"]
    print(ok.round(4).to_string())
    print(f"\n    유효 일수 {len(ok)} / 전체 {piv.shape[0]}")
    if len(ok) >= 4:
        r = stats.spearmanr(ok.index, ok.amplitude)
        lr = stats.linregress(ok.index.astype(float), ok.amplitude.astype(float))
        print(f"    Spearman(Day, 진폭) rho={r[0]:+.3f} p={r[1]:.4f}")
        print(f"    선형 기울기 {lr.slope:+.5f}/day  (p={lr.pvalue:.4f}, R2={lr.rvalue**2:.3f})")
        print("    → 음의 기울기면 미세중력 체류에 따른 행동 리듬 진폭 감소를 시사")

    # ---------------- 4) 시각별 프로파일 (파일명 타임스탬프)
    print("\n[4] 녹화 시각별 활동 프로파일")
    fh = frames.dropna(subset=["clock_h"]).copy()
    if len(fh):
        fh["hour_bin"] = fh.clock_h.astype(int)
        prof = fh.groupby("hour_bin")["active"].agg(["mean", "count"])
        prof = prof[prof["count"] >= 300]
        print(prof.round(4).to_string())
        if len(prof) >= 6:
            # 24시간 코사인 적합
            h = prof.index.to_numpy(float)
            y = prof["mean"].to_numpy(float)
            Xd = np.column_stack([np.ones_like(h), np.cos(2 * np.pi * h / 24), np.sin(2 * np.pi * h / 24)])
            beta, *_ = np.linalg.lstsq(Xd, y, rcond=None)
            fit = Xd @ beta
            ss = 1 - ((y - fit) ** 2).sum() / ((y - y.mean()) ** 2).sum()
            amp = float(np.hypot(beta[1], beta[2]))
            acro = float((-np.arctan2(beta[2], beta[1])) % (2 * np.pi) / (2 * np.pi) * 24)
            print(f"    코사인 적합: MESOR={beta[0]:.4f}  진폭={amp:.4f}  acrophase={acro:.2f} h  R2={ss:.3f}")
            print("    (주의: 녹화 시각이 균일 표집이 아니므로 해석에 제약. 계획서에 명시할 것)")

    # ---------------- 5) 민감도 검정: 'background' 교란
    # background 가 Light 에서 54% 다. 이것이 '수면' 인지 '탐지 실패' 인지 구별되지 않는다.
    # 만약 탐지 실패라면 [1] 의 명암 차이는 리듬이 아니라 영상 품질의 인공물일 수 있다.
    print("\n[5] 민감도 검정 — background 교란")
    bg = set(df[df.predicted_behavior == "background"].frame_key)
    frames["is_bg"] = frames.frame_key.isin(bg)
    print("    background 프레임 비율:")
    print(frames.groupby("Light Cycle Phase")["is_bg"].mean().round(4).to_string())

    # (a) background 프레임을 아예 제외하고 재계산
    sub = frames[~frames.is_bg]
    ta = sub.groupby("Light Cycle Phase")["active"].agg(["mean", "count"])
    print("\n    (a) background 제외 후 활동 비율:")
    print(ta.round(4).to_string())
    print("        ※ huddling 이 사실상 0 이라 이 지표는 정의상 1.0 이 된다 — 비정보적.")
    print("        따라서 (b) 가 실질적인 민감도 검정이다.")

    # (b) 동물이 탐지된 프레임 안에서, 활동 vs 비활동 행동의 상대 비율
    #     탐지 여부와 무관한 지표 → 이게 살아남아야 리듬 주장이 성립한다
    eat = set(df[df.predicted_behavior == "eating"].frame_key)
    crawl = set(df[df.predicted_behavior == "crawling"].frame_key)
    frames["eating"] = frames.frame_key.isin(eat)
    frames["crawling"] = frames.frame_key.isin(crawl)
    det = frames[frames.crawling]          # 이동이 탐지된 프레임 = 확실히 관측된 개체
    tb = det.groupby("Light Cycle Phase")["eating"].agg(["mean", "count"])
    print("\n    (b) 이동(crawling) 관측 프레임 중 섭식 동반 비율 — 탐지 편향에 강건:")
    print(tb.round(4).to_string())
    if len(tb) == 2:
        d2 = det[det["Light Cycle Phase"] == "Dark"].eating
        l2 = det[det["Light Cycle Phase"] == "Light"].eating
        c2 = stats.chi2_contingency([[d2.sum(), len(d2) - d2.sum()],
                                     [l2.sum(), len(l2) - l2.sum()]])
        print(f"        Dark {d2.mean():.4f} vs Light {l2.mean():.4f}  chi2 p={c2.pvalue:.3e}")
        print("        → 이 대비가 유지되면 명암 리듬은 탐지 인공물이 아니다")

    ok.to_csv(os.path.join(OUT, "behavior_rhythm_daily.csv"), encoding="utf-8-sig")
    bdf.to_csv(os.path.join(OUT, "behavior_rhythm_summary.csv"), index=False, encoding="utf-8-sig")
    print(f"\n-> {OUT}/behavior_rhythm_daily.csv, behavior_rhythm_summary.csv")


if __name__ == "__main__":
    main()
