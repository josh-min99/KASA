"""
Week 7-8: 계획서 첨부용 핵심 그림 5장.

  Fig 1  인공중력 데이터 지도 — 무엇이 있고 무엇이 없는가
  Fig 2  파이프라인은 신호를 잡는다 / 그러나 모듈 규모에서는 아무것도 안 보인다
  Fig 3  위상 추정 — 되는 조건과 안 되는 조건
  Fig 4  진동자 모델 — PRC 와 동조 영역
  Fig 5  검정력 — 웻랩 표본수 근거

색상은 검증된 기본 팔레트를 그대로 쓴다 (validate_palette.js 통과).
대비가 3:1 미만인 슬롯(aqua)은 직접 라벨로 보완한다.
"""
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA = os.path.join(ROOT, "data")
FIG = os.path.join(ROOT, "figures")
os.makedirs(FIG, exist_ok=True)

# ---- 팔레트 (references/palette.md 값 그대로)
SURFACE = "#fcfcfb"
INK, INK2, MUTED = "#0b0b0b", "#52514e", "#898781"
GRID, AXIS = "#e1e0d9", "#c3c2b7"
S1, S2, S3 = "#2a78d6", "#eb6834", "#1baf7a"      # categorical 1,2,3
SEQ = ["#86b6ef", "#5598e7", "#2a78d6", "#1c5cab"]  # ordinal blue 250/350/450/550
POS, NEG, MIDG = "#2a78d6", "#e34948", "#f0efec"    # diverging blue<->red

for f in ["Malgun Gothic", "NanumGothic", "Gulim"]:
    if any(x.name == f for x in font_manager.fontManager.ttflist):
        plt.rcParams["font.family"] = f
        break
plt.rcParams.update({
    "figure.facecolor": SURFACE, "axes.facecolor": SURFACE,
    "axes.edgecolor": AXIS, "axes.linewidth": 0.8,
    "axes.labelcolor": INK2, "axes.titlecolor": INK,
    "xtick.color": MUTED, "ytick.color": MUTED,
    "xtick.labelcolor": INK2, "ytick.labelcolor": INK2,
    "grid.color": GRID, "grid.linewidth": 0.8,
    "font.size": 9, "axes.titlesize": 10.5, "axes.titleweight": "bold",
    "legend.frameon": False, "axes.unicode_minus": False,
    "savefig.facecolor": SURFACE, "savefig.bbox": "tight", "savefig.dpi": 200,
})


def tidy(ax, grid_axis="y"):
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    ax.grid(True, axis=grid_axis, alpha=0.9, zorder=0)
    ax.set_axisbelow(True)


# ============================================================== Fig 1
def fig1():
    df = pd.read_csv(os.path.join(DATA, "studies_gravity.csv"))
    rows = []
    for _, r in df.iterrows():
        for tok in str(r.n_by_gravity).split("; "):
            lab, n = tok.rsplit("=", 1)
            if lab == "1G on Earth":
                continue                      # 지상대조는 비행 중 중력 단계가 아님
            g = {"uG": 0.0, "1/6G with centrifugation": 1 / 6,
                 "0.33G by centrifugation": .33, "0.66G by centrifugation": .66,
                 "1G by centrifugation": 1.0, "1G with centrifugation": 1.0}.get(lab)
            if g is None:
                continue
            tis = r.tissue if isinstance(r.tissue, str) and r.tissue.strip() else "(조직 미기재)"
            rows.append({"study": r.OSD_ID, "tissue": tis, "g": g, "n": int(n),
                         "usable": bool(r.dose_response_usable)})
    d = pd.DataFrame(rows)
    order = (d.groupby("study").g.nunique().sort_values(ascending=True).index.tolist())
    ypos = {s: i for i, s in enumerate(order)}

    fig, ax = plt.subplots(figsize=(8.4, 5.0))
    for _, r in d.iterrows():
        ci = min(int(r.g * 3.99), 3) if r.g > 0 else 0
        ax.scatter(r.g, ypos[r.study], s=150 + r.n * 22, color=SEQ[ci],
                   edgecolor=SURFACE, linewidth=1.6, zorder=3)
        ax.text(r.g, ypos[r.study], str(r.n), ha="center", va="center",
                fontsize=8.5, color="white", zorder=4, fontweight="bold")
    for s, y in ypos.items():
        ax.plot([-0.04, 1.04], [y, y], color=GRID, lw=0.8, zorder=1)
        t = d[d.study == s].tissue.iloc[0]
        lab = f"{s}  ·  {t}"
        usable = d[d.study == s].usable.iloc[0]
        ax.text(-0.09, y, lab, ha="right", va="center", fontsize=8.5,
                color=INK if usable else MUTED,
                fontweight="bold" if usable else "normal")
    ax.set_yticks([]); ax.set_ylim(-0.8, len(order) - 0.2)
    ax.set_xlim(-0.55, 1.12)
    ax.set_xticks([0, 1 / 6, .33, .66, 1.0])
    ax.set_xticklabels(["0 g\n(무중력)", "1/6 g\n(달)", "0.33 g", "0.66 g", "1 g\n(지구)"])
    ax.set_xlabel("비행 중 인공중력 수준  (원 안 숫자 = 개체 수)")
    for sp in ("top", "right", "left"):
        ax.spines[sp].set_visible(False)
    ax.set_title("인공중력 마우스 스터디 11건 — 중력 단계가 3개 이상인 것은 3건뿐",
                 loc="left", pad=14)
    ax.text(-0.55, len(order) - 0.45,
            "굵은 글씨 = 용량반응 분석 가능 (OSD-758/759/714)",
            fontsize=8, color=INK2, ha="left")
    fig.savefig(os.path.join(FIG, "fig1_data_map.png"))
    plt.close(fig)
    print("  fig1_data_map.png")


# ============================================================== Fig 2
def fig2():
    sc = pd.read_csv(os.path.join(DATA, "positive_control_response_scores.csv"))
    cp = pd.read_csv(os.path.join(DATA, "positive_control_competitive.csv"))
    fig, axes = plt.subplots(1, 2, figsize=(10.6, 4.3))

    ax = axes[0]
    for tis, col in [("Retina", S1), ("Optic nerve", S2)]:
        s = sc[sc.tissue == tis]
        ax.scatter(s.gravity_g + (0.012 if tis == "Retina" else -0.012), s.score,
                   s=34, color=col, alpha=0.55, edgecolor="none", zorder=3)
        m = s.groupby("gravity_g").score.mean()
        ax.plot(m.index, m.values, color=col, lw=2.0, marker="o", ms=8,
                markeredgecolor=SURFACE, markeredgewidth=1.5, zorder=4)
        ax.text(1.09, m.iloc[-1] + (0.10 if tis == "Retina" else -0.10),
                {"Retina": "망막", "Optic nerve": "시신경"}[tis],
                color=col, fontsize=9, fontweight="bold", va="center")
    ax.axhline(0, color=AXIS, lw=0.8, zorder=1)
    tidy(ax)
    ax.set_xticks([0, .33, .66, 1.0]); ax.set_xlim(-0.09, 1.42)
    ax.set_xlabel("인공중력 수준 (g)")
    ax.set_ylabel("우주비행 반응 점수", labelpad=6)
    ax.set_title("A. 전사체 전역 — 신호가 뚜렷하다", loc="left")
    ax.text(-0.09, ax.get_ylim()[1], "↑ 무중력에 가까움", fontsize=8, color=MUTED,
            va="top", ha="left")
    # 주의: U+2212(−) 는 Malgun Gothic 에 없다. ASCII 하이픈을 쓸 것.
    ax.text(0.02, 0.04, "인공중력 3군 Spearman\n망막 rho = -0.80 (p=0.0001)\n시신경 rho = -0.77 (p=0.0002)",
            transform=ax.transAxes, fontsize=8.2, color=INK2, va="bottom")

    ax = axes[1]
    names = {"oxidative_stress": "산화스트레스", "inflammation": "염증",
             "apoptosis": "세포사멸", "lipid_metabolism": "지질대사",
             "circadian_clock": "생체시계"}
    cp["lab"] = cp.module.map(names)
    yl = list(names.values())
    ypos = {k: i for i, k in enumerate(yl)}
    ax.axvspan(0.05, 0.95, color=MIDG, zorder=0)
    for tis, col, mk in [("Retina", S1, "o"), ("Optic nerve", S2, "s")]:
        s = cp[cp.tissue == tis]
        ax.scatter(s.percentile_vs_random, [ypos[l] for l in s.lab],
                   s=64, color=col, marker=mk, edgecolor=SURFACE, linewidth=1.4,
                   zorder=3, label={"Retina": "망막", "Optic nerve": "시신경"}[tis])
    for x in (0.05, 0.95):
        ax.axvline(x, color=AXIS, lw=1.0, ls=(0, (4, 3)), zorder=2)
    ax.set_yticks(range(len(yl))); ax.set_yticklabels(yl)
    ax.set_xlim(-0.03, 1.03); ax.set_ylim(-0.9, len(yl) - 0.4)
    ax.set_xlabel("무작위 유전자셋 대비 백분위")
    tidy(ax, grid_axis="x")
    ax.set_title("B. 모듈 규모 — 아무것도 안 보인다", loc="left")
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.16), ncol=2,
              fontsize=8.5, labelcolor=INK2)
    ax.text(0.5, -0.82, "회색 구간 = 무작위 유전자셋과 구별 불가",
            ha="center", fontsize=8, color=INK2)
    fig.suptitle("공개 데이터의 생체시계 무신호는 '효과 없음'이 아니라 '분해능 없음'이다",
                 x=0.008, ha="left", fontsize=11.5, fontweight="bold", color=INK)
    fig.tight_layout(rect=[0, 0, 1, 0.94])
    fig.savefig(os.path.join(FIG, "fig2_positive_control.png"))
    plt.close(fig)
    print("  fig2_positive_control.png")


# ============================================================== Fig 3
def fig3():
    lo = pd.read_csv(os.path.join(DATA, "phase_predictor_loto.csv"))
    rel = pd.read_csv(os.path.join(DATA, "phase_predictor_relative.csv"))
    fig, axes = plt.subplots(1, 2, figsize=(10.2, 4.3))

    ax = axes[0]
    ax.plot([0, 24], [0, 24], color=AXIS, lw=1.2, ls=(0, (4, 3)), zorder=1)
    ax.scatter(lo.true_CT, lo.pred_CT, s=30, color=S1, alpha=0.5,
               edgecolor="none", zorder=3)
    tidy(ax, grid_axis="both")
    ax.set_xlim(-1, 25); ax.set_ylim(-1, 25)
    ax.set_xticks(range(0, 25, 6)); ax.set_yticks(range(0, 25, 6))
    ax.set_xlabel("실제 시각 (CT, 시)"); ax.set_ylabel("예측 시각 (CT, 시)")
    ax.set_title("A. 24시간이 고루 표집된 경우 — 잘 맞는다", loc="left")
    ax.text(0.97, 0.06, "학습에 없던 조직으로 검증\n중앙값 오차 0.70시간\n96.9%가 3시간 이내\n\n"
                        "※ 24시 = 0시 (원형 시간축)\n   좌상·우하 점들은 실제로는 가까움",
            transform=ax.transAxes, fontsize=8.2, color=INK2, va="bottom", ha="right")

    ax = axes[1]
    x = np.arange(len(rel))
    ax.bar(x, rel.mean_corr_true_vs_pred, width=0.55, color=S2, zorder=3)
    ax.axhline(0.99, color=S1, lw=2.0, zorder=4)
    ax.text(-0.45, 1.03, "A 패널의 성능 수준", color=S1, fontsize=8.5,
            va="bottom", ha="left", fontweight="bold")
    for xi, v in zip(x, rel.mean_corr_true_vs_pred):
        ax.text(xi, v + 0.03, f"{v:.2f}", ha="center", fontsize=8.5, color=INK2)
    ax.set_xticks(x); ax.set_xticklabels([f"{w}시간" for w in rel.window_h])
    ax.set_xlim(-0.6, len(rel) - 0.4)
    ax.set_ylim(0, 1.16); ax.set_yticks([0, 0.25, 0.5, 0.75, 1.0])
    ax.set_xlabel("샘플이 채취된 시간 범위")
    ax.set_ylabel("위상 예측 정확도 (상관계수)", labelpad=6)
    tidy(ax)
    ax.set_title("B. 같은 시각에 몰려 있으면 — 무너진다", loc="left")
    ax.text(0.5, 0.62, "OSDR 우주비행 스터디는\n이 그래프의 왼쪽 끝보다도 좁다",
            transform=ax.transAxes, fontsize=8.4, color=INK2, va="top", ha="center")
    fig.suptitle("위상 추정은 원리적으로 OSDR 후향 데이터에 적용할 수 없다",
                 x=0.008, ha="left", fontsize=11.5, fontweight="bold", color=INK)
    fig.tight_layout(rect=[0, 0, 1, 0.94])
    fig.savefig(os.path.join(FIG, "fig3_phase_inference.png"))
    plt.close(fig)
    print("  fig3_phase_inference.png")


# ============================================================== Fig 4
def fig4():
    prc = pd.read_csv(os.path.join(DATA, "model_prc.csv"))
    at = pd.read_csv(os.path.join(DATA, "model_arnold_tongue.csv"))
    base = prc[prc.condition.str.contains("약한")].sort_values("CT_h")
    fig, axes = plt.subplots(1, 2, figsize=(10.6, 4.4))

    ax = axes[0]
    ax.axhline(0, color=AXIS, lw=1.0, zorder=2)
    ax.fill_between(base.CT_h, 0, base.phase_shift_h,
                    where=base.phase_shift_h >= 0, color=POS, alpha=0.22, zorder=1)
    ax.fill_between(base.CT_h, 0, base.phase_shift_h,
                    where=base.phase_shift_h < 0, color=NEG, alpha=0.22, zorder=1)
    ax.plot(base.CT_h, base.phase_shift_h, color=INK, lw=2.0, zorder=4)
    for ct, col, lab in [(16, NEG, "CT16\n최대 지연"), (21, POS, "CT21\n최대 전진"),
                         (8, MUTED, "CT8\n무반응(음성대조)")]:
        row = base.iloc[(base.CT_h - ct).abs().argmin()]
        ax.scatter([row.CT_h], [row.phase_shift_h], s=90, color=col,
                   edgecolor=SURFACE, linewidth=2, zorder=5)
        ax.annotate(lab, (row.CT_h, row.phase_shift_h),
                    textcoords="offset points",
                    xytext=(0, 16 if row.phase_shift_h >= 0 else -34),
                    ha="center", fontsize=8.3, color=col, fontweight="bold")
    tidy(ax, grid_axis="both")
    ax.set_xlim(-0.5, 24); ax.set_xticks(range(0, 25, 6))
    ax.set_ylim(-6.2, 4.2)
    ax.set_xlabel("중력 펄스를 준 내부 시각 (CT, 시)")
    ax.set_ylabel("위상 이동 (시간)", labelpad=6)
    ax.text(0.02, 0.97, "위쪽 = 시계가 앞당겨짐\n아래쪽 = 늦춰짐",
            transform=ax.transAxes, fontsize=8, color=MUTED, va="top")
    ax.set_title("A. 언제 자극하느냐가 효과의 방향을 바꾼다", loc="left")

    ax = axes[1]
    piv = at.pivot_table(index="strength", columns="period_ratio", values="entrained")
    ax.pcolormesh(piv.columns.values, piv.index.values,
                  piv.values.astype(float), cmap=matplotlib.colors.ListedColormap([SURFACE, S1]),
                  vmin=0, vmax=1, shading="nearest", zorder=2)
    ax.axvline(1.0, color=AXIS, lw=1.0, ls=(0, (4, 3)), zorder=3)
    tidy(ax, grid_axis="both")
    ax.set_xlabel("자극 주기 / 생체시계 고유주기")
    ax.set_ylabel("자극 세기 (모델 단위)", labelpad=6)
    ax.set_title("B. 동조되는 영역 (파란 부분)", loc="left")
    ax.text(0.98, 0.96, "세기를 올리면 넓어지다가\n0.14에서 최대, 0.16에서 소실",
            transform=ax.transAxes, ha="right", va="top", fontsize=8.3, color=INK2)
    fig.suptitle("진동자 모델이 웻랩 원심분리 스케줄을 직접 지정한다",
                 x=0.008, ha="left", fontsize=11.5, fontweight="bold", color=INK)
    fig.tight_layout(rect=[0, 0, 1, 0.94])
    fig.savefig(os.path.join(FIG, "fig4_model.png"))
    plt.close(fig)
    print("  fig4_model.png")


# ============================================================== Fig 5
def fig5():
    pc = pd.read_csv(os.path.join(DATA, "power_curve.csv"))
    pc = pc[pc.n_timepoints == 6]
    fig, axes = plt.subplots(1, 2, figsize=(10.4, 4.3), sharey=True)
    cols = {"망막 (OSD-758)": S1, "가자미근 (OSD-714)": S2, "시신경 (OSD-759)": S3}
    for i, (ax, eff) in enumerate(zip(axes, ["CT16 인가 (지연)", "CT21 인가 (전진)"])):
        sub = pc[pc.effect == eff]
        for tis, col in cols.items():
            s = sub[sub.tissue_sd == tis].sort_values("total_mice")
            ax.plot(s.total_mice, s.power, color=col, lw=2.0, marker="o", ms=8,
                    markeredgecolor=SURFACE, markeredgewidth=1.5, zorder=4,
                    label=tis.split(" (")[0])
            # A 패널은 곡선이 상단에서 겹치므로 직접 라벨을 달지 않는다 (범례로 식별)
            if i == 1:
                ax.text(s.total_mice.iloc[-1] + 2.5, s.power.iloc[-1],
                        tis.split(" (")[0], color=col, fontsize=8.5,
                        fontweight="bold", va="center")
        ax.axhline(0.8, color=AXIS, lw=1.4, ls=(0, (4, 3)), zorder=2)
        tidy(ax, grid_axis="both")
        ax.set_ylim(0, 1.06); ax.set_xlim(30, 148 if i == 0 else 168)
        ax.set_xlabel("총 마리 수 (두 군 합계)")
        dphi = sub.dphi_h.iloc[0]
        ax.set_title(f"{'A' if i == 0 else 'B'}. {eff} — 예상 위상 이동 {dphi} 시간",
                     loc="left")
    axes[0].set_ylabel("효과를 검출할 확률 (검정력)", labelpad=6)
    axes[0].text(146, 0.815, "검정력 80% 기준선", fontsize=8, color=INK2, ha="right")
    axes[0].legend(loc="lower right", fontsize=8.5, labelcolor=INK2,
                   title="조직별 잡음 수준", title_fontsize=8.5)
    axes[0].text(0.03, 0.06, "6시점 채취 기준", transform=axes[0].transAxes,
                 fontsize=8.4, color=INK2)
    axes[1].text(0.03, 0.94, "가자미근·시신경은 어떤 설계로도\n80%에 못 미친다",
                 transform=axes[1].transAxes, fontsize=8.4, color=INK2, va="top")
    fig.suptitle("1차 종말점은 CT16 지연이어야 한다 — CT21 전진은 검출 불가",
                 x=0.008, ha="left", fontsize=11.5, fontweight="bold", color=INK)
    fig.tight_layout(rect=[0, 0, 1, 0.94])
    fig.savefig(os.path.join(FIG, "fig5_power.png"))
    plt.close(fig)
    print("  fig5_power.png")


if __name__ == "__main__":
    print("그림 생성:")
    fig1(); fig2(); fig3(); fig4(); fig5()
    print(f"-> {FIG}")
