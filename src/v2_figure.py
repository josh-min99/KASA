"""
v2 — 조직별 I4 forest plot.

요구: x축 I4, 1.0 기준선, FLIGHT/HLU 다른 마커, 흑백 인쇄 가능, 캡션 포함.
채택된 데이터셋만 표시하고 제외된 것은 캡션에 사유와 함께 적는다.

출력: results/v2/figures/fig2_v2_tissue_matched.png
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
TAB = os.path.join(ROOT, "results", "v2", "tables")
FIG = os.path.join(ROOT, "results", "v2", "figures")
os.makedirs(FIG, exist_ok=True)

for f in ["Malgun Gothic", "NanumGothic", "Gulim"]:
    if any(x.name == f for x in font_manager.fontManager.ttflist):
        plt.rcParams["font.family"] = f
        break
plt.rcParams.update({
    "figure.facecolor": "white", "axes.facecolor": "white",
    "font.size": 9, "axes.labelsize": 9, "xtick.labelsize": 8.5,
    "ytick.labelsize": 8, "legend.fontsize": 8.5,
    "axes.unicode_minus": False, "savefig.dpi": 300,
    "savefig.bbox": "tight", "savefig.facecolor": "white",
})

TISSUE_ORDER = ["비복근", "가자미근", "비장", "망막", "등쪽피부"]


def main():
    D = pd.read_csv(os.path.join(TAB, "v2_per_dataset.csv"))
    ok = D[D["채택"] == True].copy()
    ex = D[D["채택"] != True].copy()

    rows, y = [], 0.0
    bounds = {}
    for tis in TISSUE_ORDER:
        g = ok[ok["조직"] == tis]
        if not len(g):
            continue
        top = y
        for para in ("FLIGHT", "HLU"):
            for _, r in g[g["처리"] == para].sort_values("OSD_ID").iterrows():
                rows.append({"y": y, "I4": r.I4, "para": para,
                             "lab": f"{r.OSD_ID}  (DEG {int(r.DEG):,})",
                             "tissue": tis})
                y -= 1.0
        bounds[tis] = (top, y + 1.0)
        y -= 0.7
    R = pd.DataFrame(rows)

    fig, ax = plt.subplots(figsize=(7.0, 0.36 * len(R) + 1.6))
    ax.axvline(1.0, color="#444444", lw=1.1, ls=(0, (4, 3)), zorder=1)

    for para, mk, fc, lb in (("FLIGHT", "o", "#2b2b2b", "우주비행"),
                             ("HLU", "s", "white", "지상 HLU")):
        s = R[R.para == para]
        if len(s):
            ax.scatter(s.I4, s.y, marker=mk, s=66, facecolor=fc, edgecolor="black",
                       lw=1.1, zorder=3, label=lb)

    ax.set_yticks(R.y.tolist())
    ax.set_yticklabels(R.lab.tolist())
    ax.tick_params(axis="y", length=0, pad=3)
    for sp in ("top", "right", "left"):
        ax.spines[sp].set_visible(False)
    ax.grid(True, axis="x", color="#eeeeee", lw=0.8)
    ax.set_axisbelow(True)
    xmax = max(3.0, R.I4.max() * 1.10)
    ax.set_xlim(0, xmax)
    ax.set_ylim(R.y.min() - 0.9, R.y.max() + 0.9)

    # 조직 이름은 축 바깥 왼쪽에 별도 배치 (라벨과 겹치지 않게)
    for tis, (top, bot) in bounds.items():
        ax.text(-0.30, (top + bot) / 2, tis, transform=ax.get_yaxis_transform(),
                ha="center", va="center", fontsize=9.5, fontweight="bold", rotation=90)
        ax.plot([-0.24, -0.24], [bot - 0.28, top + 0.28],
                transform=ax.get_yaxis_transform(), color="#999999", lw=1.0,
                clip_on=False, zorder=2)

    # 축 라벨 용어는 본문과 통일한다 ('I4' 는 스크립트 내부 이름이다).
    ax.set_xlabel("clock 특이 반응 지수  =  clock 유전자 |log2FC| 중앙값 ÷ 전체 유전자 |log2FC| 중앙값")
    ax.legend(loc="lower right", frameon=False)
    # 제목은 계획서의 그림 번호를 따른다.
    # (v2 작업 당시 파일명 기준으로 'Fig 2(v2)' 라 붙였으나 계획서에서는 그림 1이다.)
    ax.set_title("그림 1. 조직을 맞춘 clock 특이 반응 비교", loc="left",
                 fontsize=10.5, pad=10)
    ax.text(1.04, R.y.max() + 0.6, "1.0 = 전체 유전자와 동일", fontsize=8, color="#555555")

    # 제외 사유는 잘라내지 않고 줄바꿈한다.
    # ([:74] 로 자르면 OSD-935 의 '배경 불일치 [...] vs [Non-' 에서 끊겨 뜻이 사라진다)
    import textwrap
    exlines = []
    for _, r in ex.iterrows():
        head = f"{r['OSD_ID']} ({r['조직']}, {r['처리']}): "
        wrapped = textwrap.wrap(str(r["제외사유"]), width=104 - len(head)) or [""]
        exlines.append(head + wrapped[0])
        exlines.extend(" " * (len(head) + 2) + w for w in wrapped[1:])
    cap = ("각 점은 데이터셋 하나이며, 같은 조직 안에서 우주비행과 지상 HLU 를 나란히 배치했다.\n"
           "이 지수는 데이터셋 내부에서 자기 정규화되므로 실험 간 신호 규모 차이에 영향을 받지 않는다.\n"
           "대비는 중력 이외의 모든 배경 조건이 양쪽에서 동일할 것을 요구해 선택했다.\n"
           "품질 게이트로 FDR<0.05 인 DEG 가 30개 미만인 데이터셋은 신호 없음으로 제외했다.\n"
           "괄호 안은 해당 대비의 DEG 수다.\n\n"
           "제외된 데이터셋 9건\n  " + "\n  ".join(exlines))
    fig.text(0.02, -0.015, cap, fontsize=7.4, va="top", ha="left",
             transform=fig.transFigure)

    fig.savefig(os.path.join(FIG, "fig2_v2_tissue_matched.png"))
    plt.close(fig)
    print(f"fig2_v2_tissue_matched.png 생성 (채택 {len(ok)} / 제외 {len(ex)})")


if __name__ == "__main__":
    main()
