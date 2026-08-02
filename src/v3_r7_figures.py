"""
R7: 계획서용 그림 4장.

그림 규칙 (PROGRESS.md 에 재발 기록이 있는 것들)
  - 유니코드 마이너스는 폰트에서 두부로 깨진다 -> axes.unicode_minus=False, 문자열은 ASCII 하이픈
  - 캡션은 textwrap 으로 줄바꿈 (잘림 사고 있었음)
  - 같은 정보를 두 번 부호화하지 않는다 (마커 채움 + x축 라벨 중복 사고 있었음)
  - 범주형 축에 가로로 흩뿌리지 않는다

산출: results/v3/figures/figA~figD.png
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
CAPTIONS = []
os.makedirs(FIG, exist_ok=True)

for f in ["Malgun Gothic", "NanumGothic", "Gulim"]:
    if any(x.name == f for x in font_manager.fontManager.ttflist):
        plt.rcParams["font.family"] = f
        break
plt.rcParams.update({
    "figure.facecolor": "white", "axes.facecolor": "white",
    "font.size": 9, "axes.labelsize": 9, "xtick.labelsize": 8.5,
    "ytick.labelsize": 8.5, "legend.fontsize": 8.5,
    "axes.unicode_minus": False, "savefig.dpi": 300,
    "savefig.bbox": "tight", "savefig.facecolor": "white",
})

C_GROUND = "#1f6fb4"
C_FLIGHT = "#c8452e"
C_GREY = "#7a7a7a"
C_LIGHT = "#c9d6e3"


def caption(fig, text, width=118):
    """캡션을 이미지에 그리지 않고 모아 두기만 한다. dump_captions() 참조."""
    CAPTIONS.append(text)


# ------------------------------------------------------------------ 그림 A
def fig_a():
    au = pd.read_csv(os.path.join(DATA, "audit.csv"))
    order = ["PASS", "CONDITIONAL", "LEVEL_ONLY", "DIRECTION_ONLY", "LITERATURE", "REJECT"]
    cnt = au.verdict.value_counts().reindex(order).fillna(0).astype(int)
    label = {"PASS": "PASS\n전 축 정량", "CONDITIONAL": "CONDITIONAL\n파라미터 1점",
             "LEVEL_ONLY": "LEVEL_ONLY\n수준 축만", "DIRECTION_ONLY": "DIRECTION_ONLY\n방향만",
             "LITERATURE": "LITERATURE\n논문 수치만", "REJECT": "REJECT\n사용 불가"}
    colors = [C_FLIGHT, "#e08a35", "#d8c13a", "#8fb14a", C_GREY, C_LIGHT]

    fig, (ax, ax2) = plt.subplots(1, 2, figsize=(11.2, 4.0),
                                  gridspec_kw={"width_ratios": [1.0, 1.25]})
    y = np.arange(len(order))[::-1]
    ax.barh(y, cnt.values, color=colors, height=0.62)
    for yy, v in zip(y, cnt.values):
        ax.text(v + 0.15, yy, str(v), va="center", fontsize=9)
    ax.set_yticks(y)
    ax.set_yticklabels([label[o] for o in order], fontsize=8)
    ax.set_xlabel("후보 데이터셋 수")
    ax.set_xlim(0, max(cnt.values) + 1.6)
    ax.set_title("(A) 채택 판정", fontsize=10.5, loc="left")
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)

    # 탐색 규모
    steps = ["OSDR 스터디\n633건", "OSDR 파일\n157,801개", "리듬 후보 파일\n52개",
             "리듬 후보 스터디\n5건", "24h 커버리지\n0건"]
    vals = [633, 157801, 52, 5, 0]
    xs = np.arange(len(steps))
    ax2.plot(xs, np.log10(np.array(vals) + 1), "o-", color=C_GROUND, lw=1.6, ms=6)
    for x, v in zip(xs, vals):
        ax2.annotate(f"{v:,}", (x, np.log10(v + 1)), textcoords="offset points",
                     xytext=(0, 9), ha="center", fontsize=8.5)
    ax2.set_xticks(xs)
    ax2.set_xticklabels(steps, fontsize=7.6)
    ax2.set_ylabel("log10 (건수 + 1)")
    ax2.set_ylim(-0.35, 5.9)
    ax2.set_title("(B) OSDR 전수 스윕에서 남은 것", fontsize=10.5, loc="left")
    for s in ("top", "right"):
        ax2.spines[s].set_visible(False)

    caption(fig,
            "그림 A. 리듬 데이터셋 감사. (A) 채택 기준 C1-C7 을 적용한 판정. 기준과 판정 순서는 "
            "본 분석 전에 고정하고, 판정을 이미 아는 6건으로 역방향 검증해 6/6 일치를 확인한 뒤 적용했다. "
            "(B) OSDR 전 스터디 633건의 파일 목록 157,801개를 파일명 수준까지 열어 리듬 신호를 찾은 결과. "
            "후보로 남은 5개 스터디 중 하루 4점 이상을 24시간에 걸쳐 확보한 것은 하나도 없다. "
            "일주기를 목적으로 설계된 유일한 비행 실험(Neurolab STS-90)은 NASA 아카이브에 "
            "'No data submitted by PI' 로 기록돼 원자료가 존재하지 않는다.")
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, "그림1_데이터감사.png"))
    plt.close(fig)


# ------------------------------------------------------------------ 그림 B
def fig_b():
    co = pd.read_csv(os.path.join(DATA, "concordance.csv"))
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11.2, 4.2))

    # 좌: 활동량 수준
    lv = co[(co.axis == "활동량 수준") & co.value.notna()].copy()
    lv["short"] = lv.arm.str.replace("지상 언로딩 설치류 ", "지상 HLU ", regex=False)
    lv["short"] = lv["short"].str.replace("실제 비행 (중력만 분리, 초파리)",
                                          "궤도 uG vs 궤도 1G\n(초파리)", regex=False)
    lv = lv.iloc[::-1]
    y = np.arange(len(lv))
    cols = [C_FLIGHT if "궤도" in s else C_GROUND for s in lv["short"]]
    ax1.barh(y, lv.value.values * 100, color=cols, height=0.55)
    for yy, r in zip(y, lv.itertuples()):
        if np.isfinite(r.ci_lo):
            ax1.plot([r.ci_lo * 100, r.ci_hi * 100], [yy, yy], color="#333333", lw=1.2)
    ax1.axvline(0, color="#555555", lw=0.9)
    ax1.set_yticks(y)
    ax1.set_yticklabels(lv["short"], fontsize=7.8)
    ax1.set_xlabel("대조 대비 활동량 수준 변화 (%)")
    ax1.set_title("(A) 활동량 수준 -- 지상의 큰 감소가 궤도에서는 재현되지 않는다",
                  fontsize=10.5, loc="left")
    for s in ("top", "right"):
        ax1.spines[s].set_visible(False)

    # 우: 체온 acrophase
    tb = co[(co.axis == "체온 acrophase") & co.value.notna()].copy()
    tb["short"] = (tb.arm.str.replace("지상 언로딩 설치류 ", "지상 HLU 마우스\n", regex=False)
                   .str.replace("지상 언로딩 인간 ", "지상 HDBR 인간\n", regex=False))
    tb = tb.iloc[::-1]
    y = np.arange(len(tb))
    ax2.barh(y, tb.value.values, color=C_GROUND, height=0.55)
    for yy, r in zip(y, tb.itertuples()):
        if np.isfinite(r.ci_lo):
            ax2.plot([r.ci_lo, r.ci_hi], [yy, yy], color="#333333", lw=1.2)
        if np.isfinite(r.noise_floor):
            ax2.plot([r.noise_floor, r.noise_floor], [yy - 0.28, yy + 0.28],
                     color=C_FLIGHT, lw=1.4)
    ax2.axvline(0, color="#555555", lw=0.9)
    ax2.set_yticks(y)
    ax2.set_yticklabels(tb["short"], fontsize=7.8)
    ax2.set_xlabel("acrophase 이동 (h, 양수 = 지연)")
    ax2.set_title("(B) 체온 위상 -- 점추정은 5/5 지연, 구간이 0을 배제하는 것은 3/5",
                  fontsize=10.5, loc="left")
    ax2.text(0.98, 0.04, "세로 빨간 선 = 실측 잡음 바닥", transform=ax2.transAxes,
             ha="right", fontsize=7.6, color=C_FLIGHT)
    for s in ("top", "right"):
        ax2.spines[s].set_visible(False)

    caption(fig,
            "그림 B. 축별 정합성. 가로 막대는 개체 자기대조 평균, 가는 선은 부트스트랩 95% 구간이다. "
            "(A) 중력을 낮췄을 때 활동량 '수준' 은 지상 언로딩에서 크게 내려간다(-49 ~ -80%, 세 코호트 "
            "모두 구간이 0을 배제). 그러나 궤도에서 중력만 뺀 비교(같은 habitat 안 1G 원심 대 uG, "
            "방사선·격리·발사·온도가 양쪽 동일)에서는 +2.7% [-0.2, +6.0] 으로 그만한 감소가 없다. "
            "따라서 지상의 큰 감소를 중력 효과로 보기 어렵다. 단 종이 다르므로(마우스 대 초파리) "
            "이 대비는 방향 논증으로만 쓴다. (B) 체온 리듬 위상은 점추정이 5개 비교 모두 지연 방향이고, "
            "실제 비행 영장류(Cosmos 2044/2229)에서 보고된 방향과 일치한다. 다만 95% 구간이 0을 "
            "배제하는 것은 5개 중 3개(심부 코호트, BBR2, Cocktail)다. 세로 빨간 선은 R3.5 에서 실측한 "
            "잡음 바닥이다. 종·측정방식·기간이 교락돼 있으므로 방향만 주장하고 크기는 주장하지 않는다.")
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, "그림2_축별정합성.png"))
    plt.close(fig)


# ------------------------------------------------------------------ 그림 C
def fig_c():
    pr = pd.read_csv(os.path.join(DATA, "phase_precision.csv"))
    s = (pr[pr.condition.isin(["baseline", "treatment"])]
         .groupby(["variable", "dataset", "condition"]).phase_sd_h.median().reset_index())
    varlab = {"tb_core": "심부체온", "tb_sub": "피하온도", "activity": "활동량"}
    order = ["tb_core", "tb_sub", "activity"]

    fig, ax = plt.subplots(figsize=(8.6, 4.2))
    xs, labels, i = [], [], 0
    for v in order:
        for ds in sorted(s[s.variable == v].dataset.unique()):
            b = s[(s.variable == v) & (s.dataset == ds) & (s.condition == "baseline")]
            t = s[(s.variable == v) & (s.dataset == ds) & (s.condition == "treatment")]
            if not len(b) or not len(t):
                continue
            bv, tv = float(b.phase_sd_h.iloc[0]), float(t.phase_sd_h.iloc[0])
            ax.plot([i, i], [bv, tv], color="#999999", lw=1.0, zorder=1)
            ax.scatter([i], [bv], color=C_GROUND, s=42, zorder=3)
            ax.scatter([i], [tv], color=C_FLIGHT, s=42, marker="s", zorder=3)
            ax.annotate(f"x{tv/bv:.0f}", (i, tv), textcoords="offset points",
                        xytext=(9, 0), fontsize=7.6, color=C_FLIGHT, va="center")
            xs.append(i)
            labels.append(f"{varlab[v]}\n{ds.replace('helissen', '')}")
            i += 1
    ax.set_xticks(xs)
    ax.set_xticklabels(labels, fontsize=7.8)
    ax.set_ylabel("하루치 데이터의 위상 추정 표준편차 (h)")
    ax.set_ylim(0, 6.2)
    ax.scatter([], [], color=C_GROUND, s=42, label="HLU 전 (baseline)")
    ax.scatter([], [], color=C_FLIGHT, s=42, marker="s", label="HLU 중")
    ax.legend(loc="upper left", frameon=False)
    ax.set_title("위상 추정 정밀도는 HLU 중에 모든 지표에서 나빠진다", fontsize=10.5, loc="left")
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)

    caption(fig,
            "그림 C. 마스킹의 정량화. Helissen 2023 의 개체별 원자료에서 24시간 창을 하루씩 밀며 "
            "acrophase 를 추정하고, 그 개체 내 원형 표준편차를 구했다. HLU 를 걸면 지표 자체가 억제되어 "
            "(활동 진폭 -91%, 체온 진폭 -37 ~ -79%) 위상 추정 오차가 3.6배에서 38배까지 커진다. "
            "이것은 활동량만의 문제가 아니라 모든 지표에 해당한다. 따라서 'HLU 중에 새 위상이 형성되는가' "
            "(가설 1) 는 원리적으로 판정이 어렵고, HLU 를 제거한 free-run 구간(가설 2) 이 주 판정 "
            "경로가 되어야 한다. 심부체온이 저하 폭이 가장 작다(9.0배).")
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, "그림3_마스킹.png"))
    plt.close(fig)


# ------------------------------------------------------------------ 그림 D
def fig_d():
    g = pd.read_csv(os.path.join(DATA, "power_regression.csv"))
    fw = pd.read_csv(os.path.join(DATA, "inverse_forward.csv"))
    lim = pd.read_csv(os.path.join(DATA, "inverse_dphi.csv"))

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11.2, 4.2))

    varlab = {"tb_core": "심부체온", "tb_sub": "피하온도", "activity": "활동량"}
    style = {"tb_core": ("-", C_FLIGHT), "tb_sub": ("--", C_GROUND), "activity": (":", C_GREY)}
    sub = g[(g.dphi == 2.0) & (g.T_hours == 24)]
    for v in ["tb_core", "tb_sub", "activity"]:
        d = (sub[sub.variable == v].groupby("n_animals")
             .power_intercept.agg(["min", "max", "mean"]))
        ls, c = style[v]
        ax1.plot(d.index, d["mean"], ls, color=c, lw=1.8, label=varlab[v])
        ax1.fill_between(d.index, d["min"], d["max"], color=c, alpha=0.13, lw=0)
    ax1.axhline(0.8, color="#555555", lw=0.9, ls="-.")
    ax1.text(24, 0.815, "검정력 0.8", ha="right", fontsize=7.8, color="#555555")
    ax1.set_xlabel("군당 마리 수")
    ax1.set_ylabel("검정력")
    ax1.set_ylim(0.1, 1.03)
    ax1.legend(loc="lower right", frameon=False)
    ax1.set_title("(A) 절편 검정으로 위상 이동 2 h 검출 (처치 14일)",
                  fontsize=10.5, loc="left")
    for s in ("top", "right"):
        ax1.spines[s].set_visible(False)

    ax2.plot(fw.dphi, fw["med"], "-", color=C_GROUND, lw=1.8, label="예상 clock |log2FC| 중앙값")
    ax2.fill_between(fw.dphi, fw.q25, fw.q75, color=C_GROUND, alpha=0.15, lw=0,
                     label="희생 시각 미상에 따른 범위")
    for _, r in lim.iterrows():
        ax2.axhline(r.noise_floor, color=C_GREY, lw=0.7, ls=":")
    ax2.axhline(0.088, color=C_FLIGHT, lw=1.2, ls="--")
    ax2.text(6.2, 0.10, "비행 대비 잡음 바닥 0.088", ha="left", va="bottom",
             fontsize=7.6, color=C_FLIGHT)
    ax2.axhline(0.215, color="#333333", lw=1.0)
    ax2.text(6.2, 0.235, "비행 대비 실측 0.215", ha="left", va="bottom",
             fontsize=7.6, color="#333333")
    ax2.axvspan(1.0, 12, color="#8fb14a", alpha=0.10, lw=0)
    ax2.text(1.15, 1.50, "웻랩이 검출 가능한 영역 (Δφ >= 1 h)", fontsize=7.8, color="#4a7020")
    ax2.set_xlabel("위상 이동 Δφ (h)")
    ax2.set_ylabel("단일 시점에서 보이는 clock |log2FC|")
    ax2.set_xlim(0, 12)
    ax2.set_ylim(0, 1.62)
    ax2.legend(loc="upper left", bbox_to_anchor=(0.02, 0.93), frameon=False, fontsize=7.8)
    ax2.set_title("(B) 역문제 -- 웻랩 Δφ 가 비행 스냅샷을 해석 가능하게 만든다",
                  fontsize=10.5, loc="left")
    for s in ("top", "right"):
        ax2.spines[s].set_visible(False)

    caption(fig,
            "그림 D. 설계 사양. (A) 잡음을 Helissen 원자료의 실측 위상 추정 오차에서 뽑고, 계획서의 판정 "
            "기준인 회귀선 절편 차이 검정을 그대로 적용한 검정력. 음영은 코호트 간 범위다. "
            "2 h 이동은 심부체온 8마리, 피하온도 10-16마리, 활동량 12-16마리를 요구한다. "
            "주기가 24 h 와 다른 군에서는 기울기 검정이 훨씬 강력해 필요 마리수가 4-12마리로 줄어든다. "
            "(B) 위상만 Δφ 이동했을 때 단일 시점 전사체에서 보이는 clock |log2FC|. 진폭은 지상 아틀라스 "
            "(GSE54650, 간) 실측값이고, 희생 시각이 기록돼 있지 않으므로 0-24 h 전 구간을 두어 범위로 "
            "제시했다. 비행 대비의 잡음 바닥(0.088)을 넘으려면 Δφ >= 0.55 h 여야 하고, 웻랩이 검출할 수 "
            "있는 하한(약 1 h)이 그 위에 있다. 즉 웻랩이 Δφ 를 재 오면 기존 비행 자료의 해석이 닫힌다.")
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, "그림4_설계사양.png"))
    plt.close(fig)


def dump_captions():
    """캡션 원문을 파일로 내보낸다.

    캡션을 이미지에 굽지 않는 이유: 한글 문서에서 그림 크기를 조절하면 이미지에 박힌
    글씨도 같이 확대·축소되어 본문 글꼴과 어긋난다. 캡션은 문서 본문 텍스트로 두고,
    여기서는 원문만 내보내 문서와 그림이 어긋나지 않게 한다.
    """
    p = os.path.join(FIG, "captions.txt")
    with open(p, "w", encoding="utf-8") as fh:
        for i, c in enumerate(CAPTIONS, 1):
            fh.write("[그림 %d]\n%s\n\n" % (i, c))
    print("캡션 원문 -> %s" % p)


# ------------------------------------------------- 그림 2 (마스킹 + 검정력 합본)
def fig_masking_power():
    """계획서 축약본용. 그림 C(마스킹)와 그림 D-A(검정력)를 한 장으로 합친다.

    5단계(진동자 모델)와 예상결과 (2)(5)를 계획서에서 뺐으므로, 남은 그림은
    데이터 감사 / 마스킹+검정력 / 웻랩 예상결과 세 장이다.
    """
    pr = pd.read_csv(os.path.join(DATA, "phase_precision.csv"))
    g = pd.read_csv(os.path.join(DATA, "power_regression.csv"))
    varlab = {"tb_core": "심부체온", "tb_sub": "피하온도", "activity": "활동량"}

    fig, (ax, ax2) = plt.subplots(1, 2, figsize=(11.4, 4.2),
                                  gridspec_kw={"width_ratios": [1.15, 1.0]})

    # (A) 마스킹
    s = (pr[pr.condition.isin(["baseline", "treatment"])]
         .groupby(["variable", "dataset", "condition"]).phase_sd_h.median().reset_index())
    xs, labels, i = [], [], 0
    for v in ["tb_core", "tb_sub", "activity"]:
        for ds in sorted(s[s.variable == v].dataset.unique()):
            b = s[(s.variable == v) & (s.dataset == ds) & (s.condition == "baseline")]
            t = s[(s.variable == v) & (s.dataset == ds) & (s.condition == "treatment")]
            if not len(b) or not len(t):
                continue
            bv, tv = float(b.phase_sd_h.iloc[0]), float(t.phase_sd_h.iloc[0])
            ax.plot([i, i], [bv, tv], color="#999999", lw=1.0, zorder=1)
            ax.scatter([i], [bv], color=C_GROUND, s=40, zorder=3)
            ax.scatter([i], [tv], color=C_FLIGHT, s=40, marker="s", zorder=3)
            ax.annotate(f"x{tv/bv:.0f}", (i, tv), textcoords="offset points",
                        xytext=(8, 0), fontsize=7.4, color=C_FLIGHT, va="center")
            xs.append(i)
            labels.append(varlab[v] + "\n" + ds.replace("helissen", ""))
            i += 1
    ax.set_xticks(xs); ax.set_xticklabels(labels, fontsize=7.6)
    ax.set_ylabel("하루치 데이터의 위상 추정 표준편차 (h)")
    ax.set_ylim(0, 6.2)
    ax.scatter([], [], color=C_GROUND, s=40, label="후지현수 전")
    ax.scatter([], [], color=C_FLIGHT, s=40, marker="s", label="후지현수 중")
    ax.legend(loc="upper left", frameon=False)
    ax.set_title("(A) 후지현수 중 위상 추정 정밀도가 나빠진다", fontsize=10.5, loc="left")
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)

    # (B) 검정력
    style = {"tb_core": ("-", C_FLIGHT), "tb_sub": ("--", C_GROUND), "activity": (":", C_GREY)}
    sub = g[(g.dphi == 2.0) & (g.T_hours == 24) & (g.n_days_base == 3)]
    for v in ["tb_core", "tb_sub", "activity"]:
        d = (sub[sub.variable == v].groupby("n_animals")
             .power_intercept.agg(["min", "max", "mean"]))
        ls, c = style[v]
        ax2.plot(d.index, d["mean"], ls, color=c, lw=1.8, label=varlab[v])
        ax2.fill_between(d.index, d["min"], d["max"], color=c, alpha=0.13, lw=0)
    ax2.axhline(0.8, color="#555555", lw=0.9, ls="-.")
    ax2.text(24, 0.815, "검정력 0.8", ha="right", fontsize=7.8, color="#555555")
    ax2.set_xlabel("군당 마리 수")
    ax2.set_ylabel("검정력")
    ax2.set_ylim(0.1, 1.03)
    ax2.legend(loc="lower right", frameon=False)
    ax2.set_title("(B) 절편 검정으로 위상 이동 2 h 검출", fontsize=10.5, loc="left")
    for sp in ("top", "right"):
        ax2.spines[sp].set_visible(False)

    caption(fig,
            "그림 2. 마스킹의 정량화와 그에 따른 표본 수. (A) 개체별 원자료에서 24시간 창을 하루씩 "
            "밀며 최고시각을 추정하고 개체 내 표준편차를 구했다. 후지현수를 걸면 지표가 억제되어 "
            "위상 추정 오차가 3.6배에서 38배까지 커진다. 활동량만의 문제가 아니라 모든 지표에 "
            "해당하며, 저하 폭은 심부체온이 가장 작다. (B) 이 오차를 잡음으로 두고 계획서의 판정 "
            "기준인 회귀선 절편 차이 검정을 적용한 검정력. 음영은 코호트 간 범위다.")
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, "그림2_마스킹과검정력.png"))
    plt.close(fig)


if __name__ == "__main__":
    # 계획서에 들어가는 것은 두 장. fig_b(축별 정합성)와 fig_d(설계 사양 2패널)는
    # 예상결과 (2)(5)를 뺐으므로 본문에서 제외했으나, 발표자료용으로 코드는 남겨 둔다.
    fig_a(); print("그림1_데이터감사.png")
    fig_masking_power(); print("그림2_마스킹과검정력.png")
    dump_captions()
    print(f"-> {FIG}")
