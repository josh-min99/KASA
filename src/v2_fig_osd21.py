"""
v2 — OSD-21 앵커 그림 (전면 재설계).

재설계 사유
  이전 판은 clock 유전자 15개 중 4개(Arntl, Per1, Nr1d1, Cry2)만 그렸다.
  그 4개는 전부 '비행에서 p<0.05 이고 상향인' 유전자였다.
  유의한 5개 중 유일하게 하향인 Hlf 는 빠져 있었고, 나머지 10개도 빠져 있었다.
  즉 "비행은 움직이고 HLU 는 안 움직인다" 는 그림의 메시지가
  비행이 움직인 유전자를 골라 그린 결과였다. 순환 논증이다.

  특히 Per2 는 1단계 재현에서 '사전에' 판정 기준으로 정한 유전자인데
  (Life 2020: 근육 5/5 하향), OSD-21 에서는 +0.04 (d=0.33) 로 아무 일도
  일어나지 않는다. 사전 지정 유전자는 재현되지 않고 사후 선택 유전자만
  그려져 있었다.

  또한 이전 판은 절대 발현량을 그렸는데, 두 arm 의 기준선 차이(Arntl 0.49)가
  주장하는 효과(0.46)보다 커서, 캡션이 '보지 말라'고 한 세로 간격이
  그림에서 가장 눈에 띄는 요소였다.

수정 내용
  (A) 15개 전부를 대비(군간 차이)로 그린다.
      - 절대값이 아닌 차이를 그리므로 arm 기준선 차이가 개입하지 않는다.
      - 유전자 순서는 결과가 아니라 시계 회로의 기능 구획으로 고정한다.
        (효과크기 순으로 정렬하면 시각적으로 선택 효과가 재발한다)
      - 착륙 교란 통제군(HLU+재하중)도 같이 그린다.
  (B) 사전 지정 유전자 2개(Arntl, Per2)만 개체값으로 보여준다.
      어느 것이 사전 지정이고 어느 것이 탐색인지 그림 안에서 구별된다.

출력: results/v2/figures/fig3_v2_osd21_anchor.png
"""
import os, warnings
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.lines import Line2D

warnings.filterwarnings("ignore")
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PROC = os.path.join(ROOT, "data", "processed")
TAB = os.path.join(ROOT, "results", "tables")
FIG = os.path.join(ROOT, "results", "v2", "figures")
os.makedirs(FIG, exist_ok=True)

for f in ["Malgun Gothic", "NanumGothic", "Gulim"]:
    if any(x.name == f for x in font_manager.fontManager.ttflist):
        plt.rcParams["font.family"] = f
        break
plt.rcParams.update({
    "figure.facecolor": "white", "axes.facecolor": "white",
    "font.size": 9, "axes.labelsize": 9, "xtick.labelsize": 8,
    "ytick.labelsize": 8.5, "legend.fontsize": 8,
    "axes.unicode_minus": False, "savefig.dpi": 300,
    "savefig.bbox": "tight", "savefig.facecolor": "white",
})
INK = "#2b2b2b"

# 유전자 순서는 시계 회로의 기능 구획으로 고정한다.
# 효과크기 순으로 정렬하면 큰 값이 위로 모여 선택 효과가 시각적으로 재발한다.
GENE_BLOCKS = [
    ("활성화 축",   ["Arntl", "Clock", "Npas2"]),
    ("억제 축",     ["Per1", "Per2", "Per3", "Cry1", "Cry2"]),
    ("보조 루프",   ["Nr1d1", "Nr1d2", "Bhlhe40", "Bhlhe41"]),
    ("출력 유전자", ["Dbp", "Tef", "Hlf"]),
]
# Life 2020 재현 단계에서 분석 전에 지정한 유전자
PRESPEC = {"Arntl", "Per2"}

CONTRASTS = [
    ("Flight vs GroundControl",      "우주비행 vs 지상대조",   "o", INK,     1.1),
    ("HLU vs NormalLoaded",          "HLU vs 정상하중",        "s", "white", 1.1),
    ("HLU+Reloaded vs NormalLoaded", "HLU+재하중 vs 정상하중", "^", "#b0b0b0", 0.9),
]
# 패널 C: 처리가 같고 희생 시각만 3.5시간 다른 대비를 포함해 넷을 나란히 본다.
# 라벨은 짧게 유지한다. 길면 왼쪽 패널 영역까지 침범한다.
SUMMARY = CONTRASTS + [
    ("HLU+Reloaded vs HLU", "HLU+재하중 vs HLU †", "D", "white", 0.9)]

GROUPS = [("GroundControl", "지상대조"), ("Flight", "우주비행"),
          ("NormalLoaded", "정상하중"), ("HLU", "HLU"),
          ("HLU+Reloaded", "HLU+재하중")]


def beeswarm_offsets(values, y_range, thresh_frac=0.025, step=0.075):
    """겹치는 점만 중앙 기준 대칭으로 벌린다. 가로 위치에 정보는 없다."""
    v = np.asarray(values, dtype=float)
    order = np.argsort(v)
    thresh = y_range * thresh_frac
    off = np.zeros(len(v))
    cluster, last, clusters = [], None, []
    for idx in order:
        if last is not None and abs(v[idx] - v[last]) > thresh:
            clusters.append(cluster); cluster = []
        cluster.append(idx); last = idx
    if cluster:
        clusters.append(cluster)
    for cl in clusters:
        if len(cl) == 1:
            continue
        seq = sorted((k + 1) // 2 * (1 if k % 2 else -1) for k in range(len(cl)))
        for idx, s in zip(cl, seq):
            off[idx] = s * step
    return off


def panel_forest(ax, R):
    """A: 15개 유전자 전부, 세 대비의 효과크기와 95% 신뢰구간."""
    ypos, ylab, blocks = {}, [], {}
    y = 0.0
    for bname, genes in GENE_BLOCKS:
        top = y
        for g in genes:
            ypos[g] = y
            ylab.append((y, g))
            y -= 1.0
        blocks[bname] = (top, y + 1.0)
        y -= 0.8

    ax.axvline(0, color="#444444", lw=1.0, ls=(0, (4, 3)), zorder=1)
    off = {c[0]: d for c, d in zip(CONTRASTS, (0.24, 0.0, -0.24))}

    for cname, lab, mk, fc, lw in CONTRASTS:
        sub = R[R.comparison == cname].set_index("gene")
        for g in ypos:
            if g not in sub.index:
                continue
            r = sub.loc[g]
            yy = ypos[g] + off[cname]
            ax.plot([r.ci_low, r.ci_high], [yy, yy], color=INK, lw=lw, zorder=2,
                    alpha=0.85 if fc != "#b0b0b0" else 0.45)
            ax.scatter(r.diff_log2, yy, marker=mk, s=34, facecolor=fc,
                       edgecolor=INK, lw=0.9, zorder=3,
                       alpha=1.0 if fc != "#b0b0b0" else 0.7)

    ax.set_yticks([p for p, _ in ylab])
    ax.set_yticklabels(
        [f"{g}*" if g in PRESPEC else g for _, g in ylab])
    for t, (_, g) in zip(ax.get_yticklabels(), ylab):
        if g in PRESPEC:
            t.set_fontweight("bold")
    ax.tick_params(axis="y", length=0, pad=3)
    for sp in ("top", "right", "left"):
        ax.spines[sp].set_visible(False)
    ax.grid(True, axis="x", color="#eeeeee", lw=0.8)
    ax.set_axisbelow(True)
    ax.set_ylim(y + 0.4, 0.9)
    ax.set_xlabel("군간 차이 (log2)")

    # 유전자명이 길어(Bhlhe40/41) 축 라벨 영역이 넓다. 구획 라벨은 그보다 더 바깥에 둔다.
    for bname, (top, bot) in blocks.items():
        ax.text(-0.345, (top + bot) / 2, bname, transform=ax.get_yaxis_transform(),
                ha="center", va="center", fontsize=8, color="#555555",
                rotation=90, clip_on=False)
        ax.plot([-0.255, -0.255], [bot - 0.3, top + 0.3],
                transform=ax.get_yaxis_transform(), color="#bbbbbb", lw=0.9,
                clip_on=False, zorder=2)

    handles = [Line2D([], [], marker=mk, ls="none", ms=6.2, mfc=fc, mec=INK,
                      mew=0.9, label=lab) for _, lab, mk, fc, _ in CONTRASTS]
    ax.legend(handles=handles, loc="lower right", frameon=False,
              handletextpad=0.5, borderpad=0.2)
    ax.set_title("(A) clock 유전자 15개 전부 — 효과크기와 95% 신뢰구간",
                 loc="left", fontsize=9.5, pad=8)


def panel_points(ax, L, gene, show_xlab):
    """B: 사전 지정 유전자의 개체값."""
    sub = L[L.gene == gene]
    yr = sub.value_log2.max() - sub.value_log2.min()
    for i, (key, lab) in enumerate(GROUPS):
        v = sub[sub.group == key].value_log2.values
        if not len(v):
            continue
        ax.scatter(np.full(len(v), i) + beeswarm_offsets(v, yr), v,
                   s=20, facecolor=INK, edgecolor="white", lw=0.6,
                   alpha=0.9, zorder=3)
        ax.plot([i - 0.26, i + 0.26], [v.mean()] * 2, color="black", lw=1.7, zorder=4)
    ax.axvline(1.5, color="#cccccc", lw=0.9, ls=(0, (3, 3)), zorder=1)
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)
    ax.grid(True, axis="y", color="#eeeeee", lw=0.7)
    ax.set_axisbelow(True)
    ax.set_xticks(range(len(GROUPS)))
    ax.set_xlim(-0.6, len(GROUPS) - 0.4)
    if show_xlab:
        ax.set_xticklabels([x[1] for x in GROUPS], fontsize=7, rotation=52, ha="right")
    else:
        ax.set_xticklabels([])
    ax.set_ylabel("정규화 발현 (log2)", fontsize=8)
    ax.tick_params(labelsize=7.5)
    ax.set_title(gene, loc="left", fontsize=9, pad=4)


def panel_precision(ax, R):
    """C: 효과크기는 비슷한데 '유의한 개수'는 신뢰구간 폭이 정한다."""
    rows = []
    for cname, lab, mk, fc, _ in SUMMARY:
        s = R[R.comparison == cname]
        rows.append(dict(lab=lab, mk=mk, fc=fc,
                         med=s.diff_log2.abs().median(),
                         ciw=(s.ci_high - s.ci_low).median(),
                         nsig=int(((s.ci_low > 0) | (s.ci_high < 0)).sum())))
    S = pd.DataFrame(rows)

    y = np.arange(len(S))[::-1].astype(float)
    ax.barh(y, S.ciw, height=0.34, color="#e2e2e2", edgecolor="#9a9a9a",
            lw=0.7, zorder=2, label="신뢰구간 폭 중앙값")
    ax.scatter(S.med, y, marker="o", s=40, facecolor=INK, edgecolor="white",
               lw=0.8, zorder=4, label="|군간 차이| 중앙값")

    # 개수는 축 바깥 오른쪽에 별도 열로 둔다 (축 좌표계 사용).
    yt = ax.get_yaxis_transform()
    for yy, n in zip(y, S.nsig):
        ax.text(1.10, yy, f"{n}/15", transform=yt, fontsize=8.2, va="center",
                ha="right", clip_on=False,
                fontweight="bold" if n else "normal",
                color=INK if n else "#888888")
    ax.text(1.10, y.max() + 0.78, "CI가 0을\n배제한 수", transform=yt, fontsize=7.2,
            va="center", ha="right", color="#555555", linespacing=1.25, clip_on=False)

    ax.set_yticks(y)
    ax.set_yticklabels(S.lab, fontsize=7.4)
    ax.tick_params(axis="y", length=0, pad=3)
    for sp in ("top", "right", "left"):
        ax.spines[sp].set_visible(False)
    ax.grid(True, axis="x", color="#eeeeee", lw=0.8)
    ax.set_axisbelow(True)
    ax.set_xlim(0, 1.0)
    ax.set_ylim(-0.55, y.max() + 1.1)
    ax.set_xlabel("log2", fontsize=8)
    ax.tick_params(axis="x", labelsize=7.5)
    # 범례는 막대와 겹치지 않도록 축 아래로 뺀다.
    ax.legend(loc="upper right", bbox_to_anchor=(1.16, -0.30), ncol=2,
              frameon=False, fontsize=7.2, handletextpad=0.5,
              borderpad=0.15, columnspacing=1.4)
    ax.set_title("(C) 효과크기는 비슷하고, 유의한 개수는 정밀도가 정한다",
                 loc="left", fontsize=9.5, pad=8)


def main():
    R = pd.read_csv(os.path.join(TAB, "stage3d_osd21_clock.csv"))
    R = R[R.gene.notna()]
    L = pd.read_csv(os.path.join(PROC, "osd21_clock_long.csv"))

    fig = plt.figure(figsize=(10.2, 7.6))
    # B2 의 회전된 x라벨과 C 의 제목이 겹치지 않도록 행 간격을 넉넉히 준다.
    gs = fig.add_gridspec(3, 2, width_ratios=[1.46, 1.0],
                          height_ratios=[1.0, 1.0, 1.15],
                          hspace=0.95, wspace=0.40)
    axA = fig.add_subplot(gs[:, 0])
    axB1 = fig.add_subplot(gs[0, 1])
    axB2 = fig.add_subplot(gs[1, 1])
    axC = fig.add_subplot(gs[2, 1])

    panel_forest(axA, R)
    panel_points(axB1, L, "Arntl", False)
    panel_points(axB2, L, "Per2", True)
    axB1.set_title("(B) 사전 지정 유전자의 개체값 — Arntl", loc="left",
                   fontsize=9.5, pad=8)
    axB2.set_title("Per2", loc="left", fontsize=9, pad=4)
    panel_precision(axC, R)

    fig.suptitle("그림 2. 공개 데이터 최선의 설계(OSD-21, STS-108)에서도 판정이 불가능한 이유",
                 x=0.052, ha="left", fontsize=11, y=0.975)

    fig.text(0.052, -0.055,
        "OSD-21 은 우주비행·지상대조·HLU·HLU+재하중·정상하중 5군이 한 실험 안에 모두 있는 유일한 "
        "스터디다. 조직(비복근)이 하나로 고정되어 조직 교락이 없고,\n"
        "실험 배치가 같아 배치 교락이 없으며, 미션이 하나여서 미션 교락도 없다. 즉 본 분석이 다른 곳에서 "
        "지적한 교란이 여기에는 모두 부재하다. 군당 n=4~5.\n\n"
        "(A) 유전자 선택 없이 clock 유전자 15개 전부를 표시하였다. 절대 발현량이 아니라 군간 차이를 "
        "그렸으므로 두 실험의 기준선 차이가 개입하지 않는다. 가로선은\n"
        "95% 신뢰구간이다. 유전자 순서는 결과가 아니라 시계 회로의 기능 구획으로 고정하였다. 15개 유전자 "
        "× 3개 대비의 다중검정이므로 개별 p값은 표시하지 않았다.\n"
        "(B) 별표를 붙인 Arntl 과 Per2 는 1단계 재현에서 분석 전에 지정한 유전자다"
        "(Life 2020, doi:10.3390/life10090196). 나머지 13개는 탐색 대상이다. Arntl 은 재현되나\n"
        "Per2 는 +0.04(d=0.33)로 재현되지 않는다.\n"
        "(C) † 원저자는 HLU 군을 12일째에, HLU+재하중 군을 3.5시간 뒤에 희생시켰다"
        "(STS-108 의 착륙~희생 간격을 맞추기 위한 설계). 따라서 맨 아래 대비는\n"
        "처리가 동일하고 희생 시각만 3.5시간 다르다.\n\n"
        "읽는 법. 네 대비의 |군간 차이| 중앙값은 0.138~0.224 로 서로 비슷하다. 그런데 신뢰구간이 0 을 "
        "배제하는 유전자 수는 8, 5, 0, 0 으로 갈린다.\n"
        "그 순서는 효과크기가 아니라 신뢰구간 폭(0.43, 0.64, 0.82, 0.95)의 역순과 정확히 일치한다. "
        "즉 군당 n=4~5 에서 '유의한 clock 유전자가 몇 개인가'는\n"
        "처리가 아니라 군내 분산이 정한다. 비행에서 5개가 유의하고 HLU 에서 0개인 것도 효과크기 차이가 "
        "아니다(0.215 대 0.182). 여기에 더해 모든 군이\n"
        "단일 시점 채취이므로 관측된 차이가 진폭 변화인지 위상 이동인지 구별할 수 없다. 교란을 모두 "
        "제거해도 남는 이 두 한계는 데이터를 더 고르거나\n"
        "더 잘 골라서 해결되지 않으며, 충분한 표본과 연속 기록을 요구하는 직접적 근거다.",
        fontsize=7.3, va="top", ha="left")

    out = os.path.join(FIG, "fig3_v2_osd21_anchor.png")
    fig.savefig(out)
    plt.close(fig)
    print(f"저장: {out}")


if __name__ == "__main__":
    main()
