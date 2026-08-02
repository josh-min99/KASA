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

산출: results/v3/figures/그림5_웻랩예상결과.png
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
    "font.size": 9, "axes.labelsize": 9, "xtick.labelsize": 8.5,
    "ytick.labelsize": 8.5, "legend.fontsize": 8,
    "axes.unicode_minus": False, "savefig.dpi": 300,
    "savefig.bbox": "tight", "savefig.facecolor": "white",
})

RNG = np.random.default_rng(20260804)

TAU = 23.7           # 암전에서 마우스 고유 주기 (문헌 통상값)
N_BASE, N_HLU, N_FREE = 3, 14, 14
N_ANIMALS = 8        # R5 에서 dphi=1h 검출에 필요한 심부체온 마리 수

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

    fig, axes = plt.subplots(1, 4, figsize=(15.4, 4.1),
                             gridspec_kw={"width_ratios": [1, 1, 1, 0.92]})
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

    # ---- 4번째 패널 : Arnold tongue (동조 가능 영역)
    ax = axes[3]
    ent_pts = at[at.entrained]
    ax.scatter(ent_pts.drive_period_h, ent_pts.strength, s=9, color="#8fb14a",
               alpha=0.85, lw=0, label="동조 가능")
    for T, c in [(20, C_MASK), (24, "#444444"), (28, C_ENTRAIN)]:
        sub = at[(at.drive_period_h.between(T - 0.4, T + 0.4)) & (at.entrained)]
        smin = float(sub.strength.min()) if len(sub) else np.nan
        ax.axvline(T, color=c, lw=1.0, ls=":")
        if np.isfinite(smin):
            ax.plot([T], [smin], "v", color=c, ms=7, zorder=4)
            ax.text(T, smin - 0.012, f"{smin:.2f}", ha="center", va="top",
                    fontsize=8, color=c, fontweight="bold")
    ax.set_xlabel("HLU 주기 T (h)")
    ax.set_ylabel("필요한 자극 세기 (모델 단위)")
    ax.set_xlim(18.5, 29.5)
    ax.set_ylim(-0.03, 0.21)
    ax.set_title("동조에 필요한 자극 세기", fontsize=10, loc="left")
    ax.legend(loc="upper left", frameon=False)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)

    cap = ("그림 E. 웻랩 예상 결과 시뮬레이션. 세로 음영은 HLU 주기를 인가하는 14일, 그 오른쪽은 HLU 를 "
           "완전히 제거한 free-run 14일이다. 점은 예상 관측값이며, 오차막대는 임의값이 아니라 Helissen "
           f"원자료에서 실측한 심부체온 위상 추정 오차(HLU 전 {sd_b:.2f} h / HLU 중 {sd_t:.2f} h / 해제 후 "
           f"{sd_r:.2f} h)를 {N_ANIMALS}마리 군평균으로 환산한 값이다. HLU 인가 구간에서는 동조와 마스킹이 "
           "같은 궤적을 그리고 오차도 가장 크므로 두 가설을 구별할 수 없다. 두 예측은 HLU 를 제거한 뒤에야 "
           "갈라지며, 그 간격(초록 화살표)이 판정 신호다. 세 군 모두 간격이 군평균 표준오차의 3배를 넘어 "
           "검출 가능하다. 오른쪽 패널은 본 연구의 진동자 모델이 산출한 동조 가능 영역이다. 세 군 모두 동조 "
           "가능 범위 안에 있으나 필요한 자극 세기가 다르며, T = 28 h(14시간 on/off)가 T = 20 h(10시간 "
           "on/off)보다 뚜렷하게 강한 자극을 요구한다(모델 단위 0.13 대 0.05). 이 모델은 파라미터가 "
           "실측값이 아니므로 상대적 난이도라는 정성적 구조에만 쓰고 배수를 결론에 쓰지 않는다.")
    # 캡션은 이미지에 굽지 않는다. 한글에서 그림 크기를 조절하면 이미지에 박힌 글씨도
    # 같이 확대·축소되어 본문 글꼴과 어긋나기 때문이다. 원문만 파일로 내보낸다.
    with open(os.path.join(FIG, "caption_fig5.txt"), "w", encoding="utf-8") as fh:
        fh.write(cap + "\n")

    fig.tight_layout()
    fig.savefig(os.path.join(FIG, "그림5_웻랩예상결과.png"))
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
    print(f"\n-> {FIG}/그림5_웻랩예상결과.png, {DATA}/predicted_wetlab.csv")


if __name__ == "__main__":
    main()
