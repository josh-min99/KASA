"""
v2 — OSD-21 앵커 그림 수정판.

수정 사유
  기존 results/figures/fig3_osd21_anchor.png 은 마커의 채움 여부로
  '비행 arm'(지상대조·우주비행)과 'HLU arm'(정상하중·HLU·HLU+재하중)을 구분했으나,
  범례도 캡션 설명도 없어 독자가 그 부호화를 알 수 없었다.
  캡션이 'arm 내부 대비만 사용했다'고 서술하는데 정작 어느 점이 어느 arm 인지
  표시되지 않아, 그림을 올바르게 읽을 수 없는 상태였다.

수정 내용
  x축에 군 이름이 이미 적혀 있으므로 마커 채움은 같은 정보를 중복 부호화한 것이었다.
  중복을 없애고 arm 경계는 세로 점선 하나로만 표시한다.
  마커를 통일하면 범례 자체가 불필요해진다.
  (첫 수정에서는 범례를 추가했으나, 그것은 불필요한 부호화를 설명하려고
   요소를 더 늘린 것이었다. 부호화를 제거하는 쪽이 맞다.)

기존 파일은 덮어쓰지 않고 results/v2/figures/ 에 새로 저장한다.
"""
import os, warnings
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager
# (범례를 없앴으므로 Line2D 불필요)

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

def beeswarm_offsets(values, y_range, thresh_frac=0.025, step=0.075):
    """겹치는 점만 중앙 기준 대칭으로 벌린다.

    x 축이 범주형이므로 가로 위치에는 정보가 없다. 따라서
      - 겹치지 않는 점은 정확히 중앙에 둔다
      - 겹치는 점만, 그것도 좌우 대칭으로 배치한다 (0, +1, -1, +2, -2 ...)
    이전 판은 np.linspace 로 모든 점을 샘플 순서대로 좌->우 배치해
    가로 위치가 무언가를 뜻하는 것처럼 보였다.
    """
    v = np.asarray(values, dtype=float)
    order = np.argsort(v)
    thresh = y_range * thresh_frac
    off = np.zeros(len(v))
    cluster, last = [], None
    clusters = []
    for idx in order:
        if last is not None and abs(v[idx] - v[last]) > thresh:
            clusters.append(cluster); cluster = []
        cluster.append(idx); last = idx
    if cluster:
        clusters.append(cluster)
    for cl in clusters:
        if len(cl) == 1:
            continue
        # 0, +1, -1, +2, -2 ... 순으로 대칭 배치
        seq = []
        for k in range(len(cl)):
            m = (k + 1) // 2
            seq.append(m if k % 2 == 1 else -m)
        seq = sorted(seq)
        for idx, s in zip(cl, seq):
            off[idx] = s * step
    return off


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
        yr = sub.value_log2.max() - sub.value_log2.min()
        for i, (key, lab, arm) in enumerate(GROUPS):
            v = sub[sub.group == key].value_log2.values
            if not len(v):
                continue
            ax.scatter(np.full(len(v), i) + beeswarm_offsets(v, yr), v,
                       s=26, facecolor=INK, edgecolor="white", lw=0.7,
                       alpha=0.9, zorder=3)
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

    fig.suptitle("그림 2. 단일 스터디(OSD-21, STS-108) 내 5군 비교 — 개체값",
                 x=0.0, ha="left", fontsize=10.5, y=1.04)
    fig.text(0.0, -0.34,
             "동일 스터디에 우주비행·지상대조·HLU·HLU+재하중·정상하중 5군이 모두 있어 배치 효과가 없다.\n"
             "점은 개체, 가로선은 군 평균이다. 군당 n=4-5.\n"
             "가로축은 범주형이며 점의 가로 위치에는 의미가 없다.\n"
             "값이 서로 가까워 겹치는 개체만 중앙을 기준으로 좌우 대칭으로 벌려 표시하였다.\n"
             "세로 점선 왼쪽 두 군은 STS-108 비행 실험, 오른쪽 세 군은 지상 후지현수 실험이다.\n"
             "두 실험은 사육 조건과 기준선이 다르므로 점선을 가로질러 비교하지 않고\n"
             "각 실험 내부의 대비만 사용하였다.\n"
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
