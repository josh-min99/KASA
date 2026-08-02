"""
R9: 그림 1(웻랩 시뮬레이션)과 그림 2(역문제)를 한 장에 담는다.

왜 합치는가
  계획서 분량을 2쪽으로 맞춰야 한다. 그림을 두 블록으로 나누면 블록 사이 여백과
  캡션이 두 벌 들어가 세로로 10cm 가까이 먹는다. 한 줄 4패널로 합치면 4.2cm 다.
  그래프 자체를 합치는 것이 아니라(축이 전혀 달라 합칠 수 없다) 배치만 합친다.
  본문에서는 '그림 1', '그림 2' 로 그대로 나눠 부른다. 패널 위 제목이 그 구분이다.

글씨 크기
  이 그림은 폭 16cm(본문 전폭)로 삽입할 것을 전제한다. 그림 폭이 13.6in 이므로
  문서에서 0.46배로 줄어든다. 따라서 그림 안 글씨는 목표 크기의 약 2.2배로 키워
  둔다(예: 문서에서 7pt 로 보이려면 15pt). 이전 판이 읽히지 않았던 원인이 이것이다.

산출: results/v3/figures/그림_통합.png
"""
import os
import sys

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from v3_r8_predicted import (BAND, C_ENTRAIN, C_MASK, N_ANIMALS, N_BASE,
                             N_FREE, N_HLU, RNG, TAU, phase_sd_by_phase,
                             trajectory, wrap12)

if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA = os.path.join(ROOT, "data", "rhythm")
FIG = os.path.join(ROOT, "results", "v3", "figures")

for f in ["Malgun Gothic", "NanumGothic", "Gulim"]:
    if any(x.name == f for x in font_manager.fontManager.ttflist):
        plt.rcParams["font.family"] = f
        break

# 문서 삽입 폭 16cm 기준 환산 배율 (13.6in -> 6.3in)
S = 2.16
TICK, LABEL, TITLE, ANNO = 15, 16, 17, 13.5

plt.rcParams.update({
    "figure.facecolor": "white", "axes.facecolor": "white",
    "axes.labelsize": LABEL, "xtick.labelsize": TICK, "ytick.labelsize": TICK,
    "legend.fontsize": ANNO, "axes.unicode_minus": False,
    "savefig.dpi": 300, "savefig.facecolor": "white",
})

C_GROUND = "#2f6f9f"
C_FLIGHT = "#c8452e"
PHI0 = 14.0          # baseline 정점 시각 (임의 기준점)


def panel_wetlab(ax, T, title, sd_b, sd_t, sd_r, show_y):
    sem = lambda sd: sd / np.sqrt(N_ANIMALS)
    d, ent, msk = trajectory(T)
    b_end, h_end = N_BASE - 1, N_BASE + N_HLU - 1
    ent_o = (PHI0 + ent) % 24
    msk_o = (PHI0 + msk) % 24

    ax.axvspan(N_BASE - 0.5, N_BASE + N_HLU - 0.5, color=BAND, lw=0)

    def seg_plot(x, y, **kw):
        b = np.where(np.abs(np.diff(y)) > 12)[0]
        for xs, ys in zip(np.split(x, b + 1), np.split(y, b + 1)):
            if len(xs) > 1:
                ax.plot(xs, ys, **kw)
                kw.pop("label", None)

    seg_plot(d[:h_end + 1], ent_o[:h_end + 1], ls="-", color="#444444", lw=2.2)
    seg_plot(d[h_end:], ent_o[h_end:], ls="-", color=C_ENTRAIN, lw=3.0,
             label="동조 예측")

    # 마스킹이면 내부 시계는 내내 고유 주기 tau 로 돈다. 그 궤적을 처음부터 끝까지
    # 한 직선으로 그린다.
    #   초판은 해제 시점부터만 그렸는데, 그 시작점이 '마스킹된 관측값' 이라
    #   해제 직후 하루 만에 급히 꺾이는 가짜 구간이 생겼다(10/14시간 군에서 눈에 띈다).
    #   인가 중 구간은 관측되지 않는 값이므로 흐리게 둔다.
    hidden = (PHI0 + (TAU - 24.0) * d) % 24
    seg_plot(d[:h_end + 1], hidden[:h_end + 1], ls="--", color=C_MASK, lw=1.8,
             alpha=0.45, label="  (인가 중 · 미관측)")
    seg_plot(d[h_end:], hidden[h_end:], ls="--", color=C_MASK, lw=3.0,
             label="마스킹 예측")

    # 경계 날짜를 두 구간에 모두 넣으면 점이 두 번 찍힌다(초판 문제).
    for seg, sd, col in [(slice(0, b_end + 1), sd_b, "#444444"),
                         (slice(b_end + 1, h_end + 1), sd_t, "#444444"),
                         (slice(h_end + 1, len(d)), sd_r, C_ENTRAIN)]:
        dd = d[seg]
        yy = (ent_o[seg] + RNG.normal(0, sem(sd), len(dd))) % 24
        ax.errorbar(dd, yy, yerr=sem(sd), fmt="o", ms=4.5, lw=1.4,
                    color=col, alpha=0.85, capsize=2.4, zorder=3, ls="none")

    gap = wrap12(ent[-1] - msk[-1])
    y1, y2 = ent_o[-1], msk_o[-1]
    ax.annotate("", xy=(d[-1], y1), xytext=(d[-1], y2),
                arrowprops=dict(arrowstyle="<->", color="#2c6b2c", lw=2.0))

    # 간격 수치를 화살표 옆에 두면 패널마다 위치가 달라 데이터와 겹친다(초판 문제).
    # 제목 둘째 줄로 올려 고정한다.
    # 표준오차 대비 배수(7.1 / 6.4 / 20.0)는 넣지 않는다. 패널 폭에 안 들어가고,
    # 본문에 '표준오차의 6.4~20배' 로 범위가 이미 적혀 있다. 이걸로 표 2 를 없앴다.
    ax.set_title(f"{title}\n판정 간격 {abs(gap):.1f} h", fontsize=TITLE,
                 loc="left", pad=3, linespacing=1.25)
    ax.set_xlabel("실험 일수")
    ax.set_xlim(-1, len(d) + 1.5)
    ax.set_xticks([0, 10, 20, 30])
    ax.set_ylim(0, 24)
    ax.set_yticks([0, 6, 12, 18, 24])
    if show_y:
        ax.set_ylabel("정점 시각 (h)")
        ax.legend(loc="lower left", frameon=False, handlelength=1.3,
                  borderpad=0.1, labelspacing=0.2)
    else:
        ax.set_yticklabels([])
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    return gap


def panel_inverse(ax):
    fw = pd.read_csv(os.path.join(DATA, "inverse_forward.csv"))
    ax.plot(fw.dphi, fw["med"], "-", color=C_GROUND, lw=3.0)
    ax.fill_between(fw.dphi, fw.q25, fw.q75, color=C_GROUND, alpha=0.15, lw=0)
    ax.axhline(0.088, color=C_FLIGHT, lw=2.0, ls="--")
    ax.axhline(0.215, color="#333333", lw=1.6)
    # 두 선 사이 간격(0.127)이 글자 높이보다 좁다. 선 옆에 붙이면 반드시 겹치므로
    # 곡선 아래 빈 영역에 띄워 둔다.
    ax.text(11.8, 0.50, "비행 실측 0.215", ha="right", va="bottom",
            fontsize=ANNO, color="#333333")
    ax.text(11.8, 0.30, "잡음 바닥 0.088", ha="right", va="bottom",
            fontsize=ANNO, color=C_FLIGHT)
    ax.axvspan(1.0, 12, color="#8fb14a", alpha=0.10, lw=0)
    ax.text(1.5, 1.42, "웻랩 검출 가능 (1 h 이상)", fontsize=ANNO, color="#4a7020")
    ax.plot([0.55, 0.55], [0, 0.088], ls=":", color=C_FLIGHT, lw=1.8)
    ax.annotate("0.55 h", xy=(0.55, 0.088), xytext=(0.9, 0.62),
                fontsize=ANNO, color=C_FLIGHT,
                arrowprops=dict(arrowstyle="-", color=C_FLIGHT, lw=1.2))
    # 본문에 '0.55시간과 1시간이 겹친다' 가 그대로 있어 제목에서는 뺀다.
    ax.set_title("위상 이동 Δφ 대\n단일 시점 발현 변화", fontsize=TITLE,
                 loc="left", pad=3, linespacing=1.25)
    ax.set_xlabel("위상 이동 Δφ (h)")
    ax.set_ylabel("clock |log2FC|")
    ax.set_xlim(0, 12)
    ax.set_ylim(0, 1.6)
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)


def main():
    sd_b, sd_t, sd_r = phase_sd_by_phase()
    print(f"실측 위상 오차(심부체온): 전 {sd_b:.2f} / 중 {sd_t:.2f} / 후 {sd_r:.2f} h")

    fig = plt.figure(figsize=(13.6, 3.6))
    # 4열 사이에 빈 열을 하나 넣어 그림 1 묶음과 그림 2 를 눈으로 갈라 놓는다.
    gs = fig.add_gridspec(1, 5, width_ratios=[1, 1, 1, 0.22, 1.30], wspace=0.20,
                          left=0.048, right=0.995, top=0.715, bottom=0.150)
    axes = [fig.add_subplot(gs[0, i]) for i in (0, 1, 2, 4)]

    groups = [(24, "12시간 on/off (T = 24 h)"),
              (20, "10시간 on/off (T = 20 h)"),
              (28, "14시간 on/off (T = 28 h)")]
    gaps = [panel_wetlab(ax, T, t, sd_b, sd_t, sd_r, show_y=(i == 0))
            for i, (ax, (T, t)) in enumerate(zip(axes[:3], groups))]
    panel_inverse(axes[3])

    # 그림 번호는 본문에서 '그림 1', '그림 2' 로 부르므로 묶음 위에 크게 붙인다.
    p0, p2, p3 = (a.get_position() for a in (axes[0], axes[2], axes[3]))
    fig.text((p0.x0 + p2.x1) / 2, 0.988, "그림 1.  웻랩 시뮬레이션 결과",
             ha="center", va="top", fontsize=18, fontweight="bold")
    fig.text((p3.x0 + p3.x1) / 2, 0.988, "그림 2.  역문제",
             ha="center", va="top", fontsize=18, fontweight="bold")

    out = os.path.join(FIG, "그림_통합.png")
    fig.savefig(out)
    plt.close(fig)

    for (T, _), g in zip(groups, gaps):
        print(f"  T={T}h  판정 간격 {abs(g):.1f} h")
    print(f"-> {out}")


if __name__ == "__main__":
    main()
