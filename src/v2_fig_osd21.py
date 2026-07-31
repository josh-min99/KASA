"""
v2 — OSD-21 앵커 그림 수정판.

수정 사유
  기존 results/figures/fig3_osd21_anchor.png 은 마커의 채움 여부로
  '비행 arm'(지상대조·우주비행)과 'HLU arm'(정상하중·HLU·HLU+재하중)을 구분했으나,
  범례도 캡션 설명도 없어 독자가 그 부호화를 알 수 없었다.
  캡션이 'arm 내부 대비만 사용했다'고 서술하는데 정작 어느 점이 어느 arm 인지
  표시되지 않아, 그림을 올바르게 읽을 수 없는 상태였다.

수정 내용
  1) 범례 추가 — 두 arm 을 명시
  2) arm 경계에 구분선과 라벨 추가
  3) 캡션에 부호화 설명 명시
  4) 군 이름을 가로로 배치해 회전 라벨 겹침 해소

기존 파일은 덮어쓰지 않고 results/v2/figures/ 에 새로 저장한다.
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
    "ytick.labelsize": 8.5, "legend.fontsize": 8.5,
    "axes.unicode_minus": False, "savefig.dpi": 300,
    "savefig.bbox": "tight", "savefig.facecolor": "white",
})
INK = "#2b2b2b"

GENES = ["Arntl", "Per1", "Nr1d1", "Cry2"]
# (군, 표시명, arm)
GROUPS = [("GroundControl", "지상대조", "flight"),
          ("Flight",        "우주비행", "flight"),
          ("NormalLoaded",  "정상하중", "hlu"),
          ("HLU",           "HLU",      "hlu"),
          ("HLU+Reloaded",  "HLU+재하중", "hlu")]


def main():
    L = pd.read_csv(os.path.join(PROC, "osd21_clock_long.csv"))
    fig, axes = plt.subplots(1, len(GENES), figsize=(8.6, 3.4), sharex=True)

    for ax, g in zip(axes, GENES):
        sub = L[L.gene == g]
        for i, (key, lab, arm) in enumerate(GROUPS):
            v = sub[sub.group == key].value_log2.values
            if not len(v):
                continue
            face = INK if arm == "flight" else "white"
            ax.scatter(np.full(len(v), i) + np.linspace(-0.14, 0.14, len(v)), v,
                       s=30, facecolor=face, edgecolor="black", lw=0.8, zorder=3)
            ax.plot([i - 0.28, i + 0.28], [v.mean()] * 2, color="black", lw=1.9, zorder=4)
        # arm 경계선
        ax.axvline(1.5, color="#aaaaaa", lw=1.0, ls=(0, (3, 3)), zorder=1)
        for sp in ("top", "right"):
            ax.spines[sp].set_visible(False)
        ax.grid(True, axis="y", color="#e8e8e8", lw=0.7)
        ax.set_axisbelow(True)
        ax.set_xticks(range(len(GROUPS)))
        ax.set_xticklabels([x[1] for x in GROUPS], fontsize=7.4, rotation=52, ha="right")
        ax.set_xlim(-0.6, len(GROUPS) - 0.4)
        ax.set_title(g, fontsize=10)
        if ax is axes[0]:
            ax.set_ylabel("정규화 발현 (log2)")

    # arm 구분은 범례와 세로 점선으로 전달한다.
    # 패널 안에 arm 이름을 넣으면 유전자 제목과 겹친다.

    handles = [
        Line2D([], [], marker="o", linestyle="none", markersize=7,
               markerfacecolor=INK, markeredgecolor="black",
               label="비행 arm — 지상대조 · 우주비행"),
        Line2D([], [], marker="o", linestyle="none", markersize=7,
               markerfacecolor="white", markeredgecolor="black",
               label="HLU arm — 정상하중 · HLU · HLU+재하중"),
        Line2D([], [], color="black", lw=1.9, label="군 평균"),
    ]
    fig.legend(handles=handles, loc="lower center", ncol=3, frameon=False,
               bbox_to_anchor=(0.5, -0.16), fontsize=8.5)

    fig.suptitle("그림 2. 단일 스터디(OSD-21, STS-108) 내 5군 비교 — 개체값",
                 x=0.0, ha="left", fontsize=10.5, y=1.08)
    fig.text(0.0, -0.30,
             "동일 스터디에 우주비행·지상대조·HLU·HLU+재하중·정상하중 5군이 모두 있어 배치 효과가 없다.\n"
             "점은 개체, 가로선은 군 평균이다. 군당 n=4-5.\n"
             "마커의 채움은 실험 arm 을 나타낸다. 채워진 점은 비행 arm(STS-108 비행 실험),\n"
             "빈 점은 지상 HLU arm 이다. 두 arm 은 사육 조건과 기준선이 다르므로\n"
             "arm 을 가로질러 비교하지 않고 각 arm 내부의 대비만 사용하였다(세로 점선이 경계).\n"
             "우주비행에서 Arntl +0.46 (Cohen d=2.49), Per1 +1.59 (d=3.21), Nr1d1 +0.37 (d=3.13) 이\n"
             "관측되나 HLU 에서는 같은 방향의 변화가 나타나지 않는다 (Arntl -0.30 으로 오히려 반대).\n"
             "HLU+재하중 군은 원저자가 착륙-희생 간 3.5시간을 모사하려 설계한 군으로,\n"
             "비행 데이터의 '착륙 후 채취' 교란을 통제한다.",
             fontsize=7.6, va="top")

    out = os.path.join(FIG, "fig3_v2_osd21_anchor.png")
    fig.savefig(out)
    plt.close(fig)
    print(f"저장: {out}")


if __name__ == "__main__":
    main()
