"""
Stage 5 — 계획서 4페이지 첨부용 그림.

요구사항: 흑백 인쇄 가능, 폰트 8pt 이상, 캡션 포함.
색이 아니라 명도·마커·해칭으로 구분한다.

Fig 1  Stage 1 재현 — 비행 시 조직별 clock gene log2FC (Bmal1 vs Per2)
Fig 2  Stage 3c/4 핵심 — 비행 vs HLU 의 clock 특이 반응 (정규화)
Fig 3  Stage 3d 앵커 — OSD-21 단일 스터디 내 군간 비교 (개체값 산점도 포함)
"""
import os, warnings
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager

warnings.filterwarnings("ignore")
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
TAB = os.path.join(ROOT, "results", "tables")
PROC = os.path.join(ROOT, "data", "processed")
FIG = os.path.join(ROOT, "results", "figures")
os.makedirs(FIG, exist_ok=True)

for f in ["Malgun Gothic", "NanumGothic", "Gulim"]:
    if any(x.name == f for x in font_manager.fontManager.ttflist):
        plt.rcParams["font.family"] = f
        break
plt.rcParams.update({
    "figure.facecolor": "white", "axes.facecolor": "white",
    "font.size": 9, "axes.labelsize": 9, "axes.titlesize": 10,
    "xtick.labelsize": 8.5, "ytick.labelsize": 8.5, "legend.fontsize": 8.5,
    "axes.unicode_minus": False, "savefig.dpi": 300,
    "savefig.bbox": "tight", "savefig.facecolor": "white",
})
# 흑백 인쇄 대비: 명도 계열
K0, K1, K2, K3 = "#000000", "#4d4d4d", "#8c8c8c", "#c8c8c8"


def tidy(ax):
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    ax.grid(True, axis="y", color="#dddddd", lw=0.7)
    ax.set_axisbelow(True)


# ============================================================ Fig 1
def fig1():
    fc = pd.read_csv(os.path.join(TAB, "stage1_clock_log2fc.csv"), index_col=0)
    tissues = list(fc.columns)
    x = np.arange(len(tissues))
    fig, ax = plt.subplots(figsize=(6.6, 3.2))
    w = 0.38
    b = fc.loc["Bmal1"].values
    p = fc.loc["Per2"].values
    ax.bar(x - w / 2, b, w, color=K1, edgecolor="black", lw=0.6, label="Bmal1")
    ax.bar(x + w / 2, p, w, color="white", edgecolor="black", lw=0.6,
           hatch="////", label="Per2")
    ax.axhline(0, color="black", lw=0.8)
    tidy(ax)
    ax.set_xticks(x)
    ax.set_xticklabels([t.replace(" ", "\n") for t in tissues], fontsize=8)
    ax.set_ylabel("log2 fold change\n(우주비행 / 지상대조)")
    ax.set_title("Fig 1. 우주비행 시 조직 간 clock gene 반응 (RR-1, 8조직)", loc="left")
    ax.legend(loc="lower left", frameon=False)
    ax.text(0.99, 0.96, "Bmal1 8/8 조직 상향 · Per2 근육 5/5 하향\n"
                        "Per2 는 부신·간에서 비유의 (FDR 0.998 / 0.971)",
            transform=ax.transAxes, ha="right", va="top", fontsize=8, color="#333333")
    fig.text(0.0, -0.10,
             "Life 2020 (doi:10.3390/life10090196) 의 조직 간 clock gene 비동기화를 재현한 결과.\n"
             "GeneLab 처리 DE 테이블(OSD-98·99·101·102·103·104·105·168)의 log2FC 를 그대로 사용했다.\n"
             "Bmal1 은 모든 조직에서 일관되게 상향인 반면 Per2 는 근육에서만 하향이며 부신·간에서는\n"
             "변화가 없다. 원논문이 보고한 방향·크기와 일치한다.",
             fontsize=8, va="top")
    fig.savefig(os.path.join(FIG, "fig1_replication.png"))
    plt.close(fig)
    print("  fig1_replication.png")


# ============================================================ Fig 2
def fig2():
    D = pd.read_csv(os.path.join(TAB, "stage4_normalized.csv"))
    D = D.dropna(subset=["clock_enrichment"])
    fig, axes = plt.subplots(1, 2, figsize=(7.0, 3.3),
                             gridspec_kw={"width_ratios": [1.25, 1]})

    ax = axes[0]
    order = {"FLIGHT-A": 0, "HLU-A": 1}
    for i, (ser, g) in enumerate(D.groupby("series")):
        xs = np.full(len(g), order[ser]) + np.linspace(-0.16, 0.16, len(g))
        ax.scatter(xs, g.clock_enrichment, s=42,
                   facecolor=K1 if ser == "FLIGHT-A" else "white",
                   edgecolor="black", lw=0.8, zorder=3)
        ax.plot([order[ser] - 0.28, order[ser] + 0.28],
                [g.clock_enrichment.median()] * 2, color="black", lw=2, zorder=4)
    ax.axhline(1.0, color="#999999", lw=1.0, ls=(0, (4, 3)))
    ax.text(1.52, 1.03, "1.0 = 전역과 동일", fontsize=8, color="#555555", ha="right")
    tidy(ax)
    ax.set_xticks([0, 1])
    ax.set_xticklabels(["우주비행\n(RR-1, 8조직)", "HLU\n(지상, 4조직)"])
    ax.set_ylabel("clock 반응 / 전사체 전역 반응\n(|log2FC| 중앙값 비)")
    ax.set_ylim(0, 6)
    ax.set_title("A. clock 특이 반응", loc="left")

    ax = axes[1]
    lab = ["Bmal1\n평균 log2FC", "Bmal1\n부호일치", "조직쌍\n프로파일 r"]
    fl = [1.105, 1.000, 0.856]
    hl = [0.098, 0.500, -0.199]
    x = np.arange(3); w = 0.36
    ax.bar(x - w / 2, fl, w, color=K1, edgecolor="black", lw=0.6, label="우주비행")
    ax.bar(x + w / 2, hl, w, color="white", edgecolor="black", lw=0.6,
           hatch="////", label="HLU")
    ax.axhline(0, color="black", lw=0.8)
    tidy(ax)
    ax.set_xticks(x); ax.set_xticklabels(lab, fontsize=8)
    ax.set_ylim(-0.45, 1.3)
    ax.set_title("B. 조직 간 일관성", loc="left")
    ax.legend(loc="upper right", frameon=False)

    fig.suptitle("Fig 2. 하중 변화만으로는 우주비행의 clock gene 반응이 재현되지 않는다",
                 x=0.0, ha="left", fontsize=10, y=1.04)
    fig.text(0.0, -0.14,
             "A. 각 데이터셋에서 clock 유전자 |log2FC| 중앙값을 같은 데이터셋의 전체 유전자 |log2FC| 중앙값으로 나눈 값.\n"
             "   점은 조직(데이터셋), 가로선은 중앙값. 우주비행 3.11 (1.90–5.28) vs HLU 0.78 (0.63–1.38).\n"
             "   이 정규화는 'HLU 실험의 신호가 전반적으로 약해서' 라는 대안 설명을 배제한다.\n"
             "B. 비교는 양쪽 모두 단일 실험의 다조직 세트로 맞춰 배치 효과를 제거했다.\n"
             "   HLU 세트는 동일한 HLU×방사선 2×2 설계에서 비조사 조건만 사용했다(OSD-202·203·211·237).\n"
             "주의: 조직당 군 n=3–6. 단일 유전자 수준의 검정력은 없다(본문 참조).",
             fontsize=8, va="top")
    fig.savefig(os.path.join(FIG, "fig2_flight_vs_hlu.png"))
    plt.close(fig)
    print("  fig2_flight_vs_hlu.png")


# ============================================================ Fig 3
def fig3():
    L = pd.read_csv(os.path.join(PROC, "osd21_clock_long.csv"))
    R = pd.read_csv(os.path.join(TAB, "stage3d_osd21_clock.csv"))
    genes = ["Arntl", "Per1", "Nr1d1", "Cry2"]
    groups = ["GroundControl", "Flight", "NormalLoaded", "HLU", "HLU+Reloaded"]
    short = {"GroundControl": "지상\n대조", "Flight": "우주\n비행",
             "NormalLoaded": "정상\n하중", "HLU": "HLU", "HLU+Reloaded": "HLU+\n재하중"}
    fig, axes = plt.subplots(1, len(genes), figsize=(8.2, 3.0), sharex=True)
    fig.subplots_adjust(wspace=0.45)
    for ax, g in zip(axes, genes):
        sub = L[L.gene == g]
        for i, grp in enumerate(groups):
            v = sub[sub.group == grp].value_log2.values
            if not len(v):
                continue
            ax.scatter(np.full(len(v), i) + np.linspace(-0.13, 0.13, len(v)), v,
                       s=26, facecolor="white" if i >= 2 else K1,
                       edgecolor="black", lw=0.7, zorder=3)
            ax.plot([i - 0.26, i + 0.26], [v.mean()] * 2, color="black", lw=1.8, zorder=4)
        tidy(ax)
        ax.set_xticks(range(len(groups)))
        ax.set_xticklabels([short[x] for x in groups], fontsize=7.0, rotation=45, ha="right")
        ax.set_title(g, fontsize=9.5)
        if ax is axes[0]:
            ax.set_ylabel("정규화 발현 (log2)")
    fig.suptitle("Fig 3. 단일 스터디(OSD-21, STS-108) 내 5군 비교 — 개체값",
                 x=0.0, ha="left", fontsize=10, y=1.06)
    fig.text(0.0, -0.20,
             "동일 스터디에 우주비행·지상대조·HLU·HLU+재하중·정상하중 5군이 모두 있어 배치 효과가 없다.\n"
             "점은 개체, 가로선은 군 평균. 군당 n=4–5.\n"
             "우주비행에서 Arntl +0.46 (Cohen d=2.49), Per1 +1.59 (d=3.21), Nr1d1 +0.37 (d=3.13) 이 관측되나\n"
             "HLU 에서는 같은 방향의 변화가 나타나지 않는다 (Arntl -0.30 으로 오히려 반대).\n"
             "HLU+재하중 군은 원저자가 착륙–희생 간 3.5시간을 모사하려 설계한 군으로,\n"
             "비행 데이터의 '착륙 후 채취' 교란을 통제한다.",
             fontsize=8, va="top")
    fig.savefig(os.path.join(FIG, "fig3_osd21_anchor.png"))
    plt.close(fig)
    print("  fig3_osd21_anchor.png")


if __name__ == "__main__":
    print("그림 생성:")
    fig1(); fig2(); fig3()
    print(f"-> {FIG}")
