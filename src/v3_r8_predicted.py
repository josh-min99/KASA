"""
R8: 웻랩 예상 결과 시뮬레이션 그림.

양식 요구사항
  KASA 연구계획서 「3 예상되는 결과」는 "연구 수행과정에서 산출될 것으로 예상되는
  중간 결과물(설계도, 모형, 시뮬레이션 자료 등)을 포함" 하라고 명시한다.
  따라서 예상 결과를 말로만 쓰지 않고 시뮬레이션 그림으로 제시한다.

이 그림이 '그럴듯한 만화'가 아닌 이유
  오차막대를 임의로 그리지 않는다. R5 에서 Helissen 원자료로 실측한
  구간별 위상 추정 오차를 그대로 쓴다.
      HLU 전 0.31 h / HLU 중 2.77 h / 해제 후 1.67 h  (심부체온, 개체 내)
  동조 여부를 판정할 수 있는 구간이 어디인지가 이 오차에서 그림으로 드러난다.

모델
  암전(DD)에서 마우스 고유 주기 tau = 23.7 h 로 둔다(문헌 통상값).
  acrophase 는 하루에 (period - 24) 시간만큼 시계시각 기준으로 이동한다.
    - 동조 시나리오 : HLU 구간에서 period = T (10/12/14시간 on-off -> T = 20/24/28 h)
    - 마스킹 시나리오: 내부 시계는 계속 tau. 관측만 HLU 주기에 끌려간다
  두 시나리오는 HLU 구간에서 관측이 거의 같다. HLU 를 제거해야 갈라진다.
  이것이 free-run 구간을 주 판정 경로로 두는 이유이며, 그림의 핵심 메시지다.

산출: results/v3/figures/그림1_웻랩예상결과.png
      data/rhythm/predicted_wetlab.csv
"""
import os
import sys
import textwrap

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager

if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA = os.path.join(ROOT, "data", "rhythm")
FIG = os.path.join(ROOT, "results", "v3", "figures")
os.makedirs(FIG, exist_ok=True)

for f in ["Malgun Gothic", "NanumGothic", "Gulim"]:
    if any(x.name == f for x in font_manager.fontManager.ttflist):
        plt.rcParams["font.family"] = f
        break
plt.rcParams.update({
    "figure.facecolor": "white", "axes.facecolor": "white",
    "font.size": 11, "axes.labelsize": 11, "xtick.labelsize": 10,
    "ytick.labelsize": 10, "legend.fontsize": 9.5,
    "axes.unicode_minus": False, "savefig.dpi": 300,
    "savefig.bbox": "tight", "savefig.facecolor": "white",
})

RNG = np.random.default_rng(20260804)

TAU = 23.7           # 암전에서 마우스 고유 주기 (문헌 통상값)
N_BASE, N_HLU, N_FREE = 3, 14, 14
N_ANIMALS = 8        # R5b 에서 심부체온 절편 검정으로 dphi=2h 를 검출하는 데 필요한 수

C_ENTRAIN = "#c8452e"
C_MASK = "#1f6fb4"
C_GREY = "#8a8a8a"
BAND = "#f0f0f0"


def phase_sd_by_phase():
    """R5 산출물에서 심부체온의 구간별 위상 추정 오차를 읽는다(실측값)."""
    pr = pd.read_csv(os.path.join(DATA, "phase_precision.csv"))
    s = (pr[pr.variable == "tb_core"]
         .groupby("condition").phase_sd_h.median())
    return float(s["baseline"]), float(s["treatment"]), float(s["recovery"])


def trajectory(T):
    """두 시나리오의 일별 acrophase(시계시각 기준, 누적 표기).

    반환: days, 동조, 마스킹
    """
    d = np.arange(N_BASE + N_HLU + N_FREE)
    ent = np.zeros(len(d), float)
    msk = np.zeros(len(d), float)

    drift_free = TAU - 24.0        # -0.3 h/day
    drift_ent = T - 24.0           # 동조 시 HLU 구간의 하루 이동량

    for i, day in enumerate(d):
        if day < N_BASE:                       # baseline (DD, 자유진행)
            ent[i] = msk[i] = drift_free * day
        elif day < N_BASE + N_HLU:             # HLU 구간
            k = day - (N_BASE - 1)
            base = drift_free * (N_BASE - 1)
            # 동조: 내부 주기가 T 로 끌려간다
            ent[i] = base + drift_ent * k
            # 마스킹: 내부는 tau 그대로지만 관측이 HLU 주기에 끌려가 거의 같아 보인다
            msk[i] = base + drift_ent * k
        else:                                  # free-run (HLU 제거)
            k = day - (N_BASE + N_HLU - 1)
            ent_end = drift_free * (N_BASE - 1) + drift_ent * N_HLU
            ent[i] = ent_end + drift_free * k
            # 마스킹이었다면 내부 시계는 내내 tau 였으므로 baseline 연장선으로 복귀
            msk[i] = drift_free * day
    return d, ent, msk


def wrap12(x):
    """관측 가능한 위상차는 24시간 순환이므로 -12~+12 로 접는다."""
    return (x + 12) % 24 - 12


def main():
    sd_b, sd_t, sd_r = phase_sd_by_phase()
    sem = lambda sd: sd / np.sqrt(N_ANIMALS)
    print(f"실측 위상 추정 오차(심부체온) : HLU전 {sd_b:.2f} / HLU중 {sd_t:.2f} / 해제후 {sd_r:.2f} h")
    print(f"n={N_ANIMALS} 일 때 군평균 표준오차 : {sem(sd_b):.2f} / {sem(sd_t):.2f} / {sem(sd_r):.2f} h")

    at = pd.read_csv(os.path.join(DATA.replace("rhythm", ""), "model_arnold_tongue.csv"))

    groups = [(24, "12시간 on/off  (T = 24 h)"),
              (20, "10시간 on/off  (T = 20 h)"),
              (28, "14시간 on/off  (T = 28 h)")]

    fig, axes = plt.subplots(1, 3, figsize=(10.4, 3.9))
    rows = []

    # 실제로 관측되는 값은 '하루 중 몇 시에 정점인가'(0-24h)이다. 누적 이동량이 아니다.
    # 초판에서 누적값으로 그렸더니 패널마다 y 범위가 2h~60h 로 달라져,
    # 판정 간격 화살표의 길이가 실제 크기와 반대로 읽혔다. 관측값으로 바꿔 축을 통일한다.
    PHI0 = 14.0  # baseline 정점 시각 (임의 기준점)

    for ax, (T, title) in zip(axes[:3], groups):
        d, ent, msk = trajectory(T)
        b_end, h_end = N_BASE - 1, N_BASE + N_HLU - 1
        ent_o = (PHI0 + ent) % 24
        msk_o = (PHI0 + msk) % 24

        ax.axvspan(N_BASE - 0.5, N_BASE + N_HLU - 0.5, color=BAND, lw=0)
        ax.text(N_BASE + N_HLU / 2 - 0.5, 0.972, "HLU 주기 인가",
                transform=ax.get_xaxis_transform(), ha="center", va="top",
                fontsize=8, color="#555555")
        ax.text(N_BASE + N_HLU + N_FREE / 2 - 0.5, 0.972, "free-run (판정 구간)",
                transform=ax.get_xaxis_transform(), ha="center", va="top",
                fontsize=8, color="#2c6b2c", fontweight="bold")

        # 예측선 — 24h 를 넘어가며 접히므로 끊어 그린다
        def seg_plot(x, y, **kw):
            b = np.where(np.abs(np.diff(y)) > 12)[0]
            for xs, ys in zip(np.split(x, b + 1), np.split(y, b + 1)):
                if len(xs) > 1:
                    ax.plot(xs, ys, **kw)
                    kw.pop("label", None)

        seg_plot(d[:h_end + 1], ent_o[:h_end + 1], ls="-", color="#444444", lw=1.4)
        seg_plot(d[h_end:], ent_o[h_end:], ls="-", color=C_ENTRAIN, lw=2.0,
                 label="동조 예측")
        seg_plot(d[h_end:], msk_o[h_end:], ls="--", color=C_MASK, lw=2.0,
                 label="마스킹 예측")

        # 관측 모사 (실측 오차로 잡음 부여)
        for seg, sd, col in [(slice(0, b_end + 1), sd_b, "#444444"),
                             (slice(b_end, h_end + 1), sd_t, "#444444"),
                             (slice(h_end, len(d)), sd_r, C_ENTRAIN)]:
            dd = d[seg]
            yy = (ent_o[seg] + RNG.normal(0, sem(sd), len(dd))) % 24
            ax.errorbar(dd, yy, yerr=sem(sd), fmt="o", ms=3.0, lw=0.9,
                        color=col, alpha=0.85, capsize=1.6, zorder=3, ls="none")

        gap = wrap12(ent[-1] - msk[-1])
        rows.append(dict(T_hours=T, group=title, gap_h=gap,
                         sd_base=sd_b, sd_hlu=sd_t, sd_free=sd_r,
                         sem_free=sem(sd_r), detectable=abs(gap) > 3 * sem(sd_r)))

        y1, y2 = ent_o[-1], msk_o[-1]
        ax.annotate("", xy=(d[-1], y1), xytext=(d[-1], y2),
                    arrowprops=dict(arrowstyle="<->", color="#2c6b2c", lw=1.3))
        ax.text(d[-1] - 0.8, (y1 + y2) / 2, f"판정 간격\n{abs(gap):.1f} h",
                ha="right", va="center", fontsize=8, color="#2c6b2c",
                fontweight="bold")

        ax.set_title(title, fontsize=10, loc="left")
        ax.set_xlabel("실험 일수")
        ax.set_xlim(-1, len(d) + 1.5)
        ax.set_ylim(0, 24)
        ax.set_yticks([0, 6, 12, 18, 24])
        for s in ("top", "right"):
            ax.spines[s].set_visible(False)
    axes[0].set_ylabel("관측 정점 시각 acrophase (h)")
    axes[0].legend(loc="lower left", frameon=False)

    # (진동자 모델 패널은 5단계를 계획서에서 빼면서 함께 제거)

    cap = ("그림 1. 웻랩 시뮬레이션 결과. 음영은 후지현수 인가 14일, 오른쪽은 제거한 자유진행 14일이다. "
           "오차막대는 실측 심부체온 위상 오차(인가 전 0.31 / 중 2.77 / 해제 후 1.67 h)의 8마리 군평균이다. "
           "인가 중에는 동조와 마스킹이 같은 궤적을 그려 구별되지 않고, 제거 후 갈라지는 간격(초록 화살표)이 판정 신호다.")

    # 캡션을 이미지 안에 넣는다(사용자 요청). 문서 본문에서 캡션 문단이 차지하던
    # 줄을 없애 분량을 줄이기 위해서다. 대신 글씨가 그림과 함께 축소되므로
    # 그림 폭을 너무 줄이면 읽히지 않는다.
    fig.subplots_adjust(bottom=0.30, left=0.06, right=0.99, top=0.93)
    fig.text(0.012, 0.055, "\n".join(textwrap.wrap(cap, 132)),
             fontsize=9.0, va="top", ha="left", color="#222222")
    with open(os.path.join(FIG, "caption_fig1.txt"), "w", encoding="utf-8") as fh:
        fh.write(cap + "\n")
    fig.savefig(os.path.join(FIG, "그림1_웻랩예상결과.png"))
    plt.close(fig)

    df = pd.DataFrame(rows)
    # Arnold tongue 요약도 함께 저장
    need = []
    for T in [20, 24, 28]:
        sub = at[(at.drive_period_h.between(T - 0.4, T + 0.4)) & (at.entrained)]
        need.append(dict(T_hours=T,
                         min_strength=float(sub.strength.min()) if len(sub) else np.nan))
    df = df.merge(pd.DataFrame(need), on="T_hours")
    df.to_csv(os.path.join(DATA, "predicted_wetlab.csv"), index=False, encoding="utf-8-sig")

    print("\n=== 군별 예상 판정 간격 ===")
    print(df[["T_hours", "gap_h", "sem_free", "detectable", "min_strength"]]
          .round(3).to_string(index=False))
    print(f"\n-> {FIG}/그림1_웻랩예상결과.png, {DATA}/predicted_wetlab.csv")


if __name__ == "__main__":
    main()
