"""
R4: 축별 정합성 — 지상 아날로그와 실제 비행이 어느 축에서 일치하는가.

논지
  "HLU 가 우주와 같다" 를 통째로 증명하는 것은 불가능하고, v1/v2 에서 전사체 축으로는
  판정조차 불가능함을 이미 보였다. 대신 **전이되는 축과 전이되지 않는 축을 분리**한다.
  1차 종말점을 전이되는 축 위에 올리는 것이 목적이다.

R1 감사 결과에 따른 제약 (게이트 G1: 원자료 정량 사용 가능 2건)
  원자료 기반 칸과 논문 보고 수치 기반 칸을 반드시 구분해 표기한다.
  교락(종/측정방식/조명/노출기간/구속)이 남은 비교는 **방향만** 주장하고 크기는 주장하지 않는다.
  이 규칙은 PROGRESS.md 의 v1 결론 3건 철회 사고에서 나왔다.

산출: data/rhythm/concordance.csv, results/v3/r4_concordance.txt
"""
import os
import sys

import numpy as np
import pandas as pd

if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, "data", "rhythm")
RES = os.path.join(ROOT, "results", "v3")
os.makedirs(RES, exist_ok=True)


def boot_ci(vals, n=4000, seed=7):
    v = np.asarray([x for x in vals if np.isfinite(x)], float)
    if len(v) < 2:
        return (np.nan, np.nan)
    rng = np.random.default_rng(seed)
    m = rng.choice(v, size=(n, len(v)), replace=True).mean(axis=1)
    return (float(np.percentile(m, 2.5)), float(np.percentile(m, 97.5)))


def main():
    par = pd.read_csv(os.path.join(OUT, "params.csv"))
    nf = pd.read_csv(os.path.join(OUT, "noise_floor.csv"))
    long = pd.read_csv(os.path.join(OUT, "long.csv"))
    nfi = {(r.dataset, r.variable): r for _, r in nf.iterrows()}

    cells = []

    def add(axis, arm, tier, value, ci, noise, unit, source, confound, note):
        ratio = (abs(value) / noise) if (np.isfinite(value) and noise and np.isfinite(noise)
                                         and noise > 0) else np.nan
        cells.append(dict(axis=axis, arm=arm, tier=tier, value=value,
                          ci_lo=ci[0] if ci else np.nan, ci_hi=ci[1] if ci else np.nan,
                          noise_floor=noise, noise_multiple=ratio, unit=unit,
                          source=source, confound=confound, note=note))

    # ============================================================ 원자료 기반
    # ---- 지상 언로딩 (설치류) : Helissen. 개체별 baseline 대비 treatment
    for ds, cohort in [("helissen2020", "피하 코호트 A"), ("helissen2021", "피하 코호트 B"),
                       ("helissen2022", "심부 코호트")]:
        for var, axis_amp, axis_acr, axis_lvl in [
                ("activity", "활동량 리듬 진폭", "활동량 acrophase", "활동량 수준"),
                ("tb_sub", "체온 리듬 진폭", "체온 acrophase", "체온 수준"),
                ("tb_core", "체온 리듬 진폭", "체온 acrophase", "체온 수준")]:
            sub = par[(par.dataset == ds) & (par.variable == var)]
            if not len(sub):
                continue
            b = sub[sub.condition == "baseline"].set_index("subject_id")
            t = sub[sub.condition == "treatment"].set_index("subject_id")
            common = sorted(set(b.index) & set(t.index))
            if len(common) < 3:
                continue
            n = nfi.get((ds, var))
            # 진폭비 - 1  (자기대조 정규화)
            ampr = [(t.loc[s, "amplitude"] / b.loc[s, "amplitude"]) - 1 for s in common
                    if b.loc[s, "amplitude"] > 0]
            add(axis_amp, f"지상 언로딩 설치류 ({cohort})", "PASS(원자료)",
                float(np.mean(ampr)), boot_ci(ampr),
                float(n.amp_ratio_floor) if n is not None else np.nan,
                "baseline 대비 비 - 1", f"Helissen {ds}",
                "없음 (동일 개체 자기대조)",
                f"n={len(common)}")
            # acrophase 이동
            dphi = []
            for s in common:
                d = (t.loc[s, "acrophase"] - b.loc[s, "acrophase"] + 12) % 24 - 12
                dphi.append(d)
            add(axis_acr, f"지상 언로딩 설치류 ({cohort})", "PASS(원자료)",
                float(np.mean(dphi)), boot_ci(dphi),
                float(n.acrophase_floor_h) if n is not None else np.nan,
                "시간(+ 지연)", f"Helissen {ds}",
                "HLU 중 지표 억제로 인한 마스킹",
                f"n={len(common)}")
            # 수준 (MESOR 비)
            lvl = [(t.loc[s, "mesor"] / b.loc[s, "mesor"]) - 1 for s in common
                   if b.loc[s, "mesor"] != 0]
            add(axis_lvl, f"지상 언로딩 설치류 ({cohort})", "PASS(원자료)",
                float(np.mean(lvl)), boot_ci(lvl),
                float(n.mesor_floor / abs(b.mesor.mean())) if (n is not None and b.mesor.mean()) else np.nan,
                "baseline 대비 비 - 1", f"Helissen {ds}",
                "없음 (동일 개체 자기대조)", f"n={len(common)}")

    # ---- 지상 언로딩 (인간) : HDBR 3주차 -> 8주차
    for ds in ["hdbr_BBR2", "hdbr_Cocktail"]:
        sub = par[(par.dataset == ds) & (par.variable == "tb_core")]
        b = sub[sub.condition == "baseline"].set_index("subject_id")
        t = sub[sub.condition == "treatment"].set_index("subject_id")
        common = sorted(set(b.index) & set(t.index))
        if len(common) < 3:
            continue
        dphi = [((t.loc[s, "acrophase"] - b.loc[s, "acrophase"] + 12) % 24 - 12) for s in common]
        add("체온 acrophase", f"지상 언로딩 인간 ({ds.split('_')[1]})", "CONDITIONAL(원자료)",
            float(np.mean(dphi)), boot_ci(dphi), np.nan, "시간(+ 지연)",
            f"HDBR figshare 13633790", "3주차 자체가 이미 침상안정 중 — 절대 기준 아님",
            f"n={len(common)}, 3주차 대비 8주차")
        ampr = [(t.loc[s, "amplitude"] / b.loc[s, "amplitude"]) - 1 for s in common
                if b.loc[s, "amplitude"] > 0]
        add("체온 리듬 진폭", f"지상 언로딩 인간 ({ds.split('_')[1]})", "CONDITIONAL(원자료)",
            float(np.mean(ampr)), boot_ci(ampr), np.nan, "비 - 1",
            "HDBR figshare 13633790", "동일", f"n={len(common)}")

    # ---- 실제 비행, 중력만 분리 : OSD-595 (초파리, 활동량 '수준' 만)
    o = long[long.dataset == "osd595"]
    g1 = o[o.gravity == "1G"].groupby("subject_id").value.mean()
    ug = o[o.gravity == "uG"].groupby("subject_id").value.mean()
    rel = float(ug.mean() / g1.mean() - 1)
    rng = np.random.default_rng(11)
    bs = [float(rng.choice(ug.values, len(ug)).mean() / rng.choice(g1.values, len(g1)).mean() - 1)
          for _ in range(4000)]
    add("활동량 수준", "실제 비행 (중력만 분리, 초파리)", "LEVEL_ONLY(원자료)",
        rel, (float(np.percentile(bs, 2.5)), float(np.percentile(bs, 97.5))), np.nan,
        "궤도 1G 대비 비 - 1", "OSD-595 / Cell Rep 40:111279",
        "종(초파리), 명기 전용 측정",
        f"궤도상 1G 원심 {len(g1)}모듈 vs uG {len(ug)}모듈, 16일. "
        f"방사선·격리·발사가 양쪽 동일")

    # ============================================================ 문헌 기반
    LIT = [
        ("활동량 리듬 진폭", "실제 비행 설치류 (RR-1)", "DIRECTION_ONLY",
         np.nan, "유지 (Dark 0.931 vs Light 0.455, chi2 p<1e-300)",
         "OSD-952", "케이지 단위, baseline 없음, 측정방식(영상분류)",
         "미세중력에서 야행성 행동 리듬이 탐지편향 통제 후에도 유지"),
        ("활동량 리듬 진폭", "실제 비행 설치류 (MHU-1)", "DIRECTION_ONLY",
         np.nan, "지상대조에서 주간 32.5% / 야간 63.2%",
         "Shimbo 2021 Sci Rep 11:2827", "영상 2개씩뿐, 24h 연속 없음",
         "궤도상 1g 대조가 있으나 시계열이 끊겨 리듬 비교 불가"),
        ("체온 acrophase", "실제 비행 영장류 (Cosmos 2044/2229)", "LITERATURE",
         np.nan, "위상 지연 (방향만 보고)",
         "Fuller 1996 J Appl Physiol 81:188", "n=4, 1990년대, 원자료 없음",
         "비행 중 체온 리듬 위상 지연 + 평균 심박 감소"),
        ("체온 acrophase", "지상 과중력 설치류 (2G 펄스)", "LITERATURE",
         np.nan, "재동조 가속 16.0+-3.1일 vs 20.4+-0.8일 (p=0.002)",
         "Martin 2020 Sci Rep 10:8646", "과중력이라 미세중력과 방향 반대, 원자료 없음",
         "전정손상 쥐에서 효과 소실 -> 중력 입력이 위상 조절에 관여"),
        ("체온 리듬 진폭", "지상 과중력 설치류 (만성 2G)", "LITERATURE",
         np.nan, "초기 7-10일 진폭 감소 후 회복",
         "Fuller 1994 계열", "원자료 없음, 1990년대",
         "활동 리듬은 같은 기간 일시 소실"),
        ("리듬 존재 여부", "실제 비행 인간 (ISS 6개월+)", "LITERATURE",
         np.nan, "24/12/8시간 성분 유지",
         "PMC5192238 / Sci Rep 8:10381", "HRV 지표, 24h Holter 5시점",
         "우주인 7명에서 일주기 성분이 비행 중에도 유지됨"),
        ("활동량 리듬 진폭", "지상 언로딩 설치류 (스트레스 분해)", "LITERATURE",
         np.nan, "초기 저활동은 사회격리·구속 탓, 후기 활성화가 언로딩 고유 효과",
         "PMC3388052", "쥐, 원자료 없음",
         "5군(쌍사육/개별/언로딩케이지/수평구속/HDT) 분해. sham 설계의 근거"),
    ]
    for axis, arm, tier, val, desc, src, conf, note in LIT:
        cells.append(dict(axis=axis, arm=arm, tier=tier, value=val,
                          ci_lo=np.nan, ci_hi=np.nan, noise_floor=np.nan,
                          noise_multiple=np.nan, unit=desc, source=src,
                          confound=conf, note=note))

    df = pd.DataFrame(cells)
    df.to_csv(os.path.join(OUT, "concordance.csv"), index=False, encoding="utf-8-sig")

    # ------------------------------------------------------------- 보고
    log = []
    P = log.append
    P("=" * 100)
    P("R4  축별 정합성")
    P("=" * 100)
    P("")
    P("표기 규칙")
    P("  tier 가 (원자료)인 칸만 효과 크기를 주장한다.")
    P("  LITERATURE / DIRECTION_ONLY 칸은 방향만 쓰고 크기는 주장하지 않는다.")
    P("  '잡음배수' 는 R3.5 에서 실측한 잡음 바닥 대비 몇 배인가. 1 이하면 잡음과 구별 안 됨.")
    P("")

    for axis in ["활동량 리듬 진폭", "체온 리듬 진폭", "활동량 acrophase",
                 "체온 acrophase", "활동량 수준", "체온 수준", "리듬 존재 여부"]:
        sub = df[df.axis == axis]
        if not len(sub):
            continue
        P("-" * 100)
        P(f"[{axis}]")
        for _, r in sub.iterrows():
            if np.isfinite(r.value):
                ci = (f" [{r.ci_lo:+.3f}, {r.ci_hi:+.3f}]"
                      if np.isfinite(r.ci_lo) else "")
                nm = (f"  잡음대비 {r.noise_multiple:.1f}x"
                      if np.isfinite(r.noise_multiple) else "  잡음바닥 미산출")
                P(f"  {r.arm:42s} {r.value:+8.3f}{ci}{nm}")
                P(f"      {r.tier} / {r.unit} / {r.source}")
            else:
                P(f"  {r.arm:42s} {r.unit}")
                P(f"      {r.tier} / {r.source}")
            P(f"      교락: {r.confound}")
            if r.note:
                P(f"      비고: {r.note}")
        P("")

    P("=" * 100)
    P("판독")
    P("=" * 100)

    # 활동량 축: 지상과 비행의 방향
    act_amp = df[(df.axis == "활동량 리듬 진폭") & df.value.notna()]
    P("")
    P("1) 활동량 축 — 지상 아날로그와 실제 비행이 어긋난다")
    for _, r in act_amp.iterrows():
        P(f"   지상: {r.arm} 진폭비 {r.value:+.3f} (잡음대비 {r.noise_multiple:.1f}x)")
    P("   비행: RR-1 마우스는 야행성 리듬을 유지 (OSD-952, 방향만)")
    lvl = df[(df.axis == "활동량 수준") & (df.tier == "LEVEL_ONLY(원자료)")]
    if len(lvl):
        r = lvl.iloc[0]
        P(f"   비행(중력만 분리): 초파리 uG 활동 수준이 궤도 1G 대비 {r.value:+.1%} "
          f"[{r.ci_lo:+.1%}, {r.ci_hi:+.1%}]")
    hl = df[(df.axis == "활동량 수준") & df.arm.str.contains("지상 언로딩")]
    for _, r in hl.iterrows():
        P(f"   지상 언로딩 활동 수준: {r.value:+.1%}  ({r.arm})")
    P("   -> 중력을 낮췄을 때 활동 '수준' 이 지상에서는 크게 내려가고 궤도에서는 오히려 올라간다.")
    P("      부호가 반대다. 활동량은 아날로그에서 비행으로 전이되는 축이 아니다.")

    P("")
    P("2) 체온 위상 축 — 방향이 일치한다")
    tb = df[(df.axis == "체온 acrophase")]
    for _, r in tb.iterrows():
        if np.isfinite(r.value):
            P(f"   {r.arm}: {r.value:+.2f} h" +
              (f" (잡음바닥 {r.noise_floor:.2f} h, {r.noise_multiple:.1f}x)"
               if np.isfinite(r.noise_multiple) else ""))
        else:
            P(f"   {r.arm}: {r.unit}")
    P("   -> 지상 언로딩(인간·설치류)과 실제 비행(영장류)이 모두 '지연' 방향이다.")
    P("      단 종·측정방식·기간이 교락돼 있으므로 크기는 주장하지 않는다.")

    txt = "\n".join(log)
    print(txt)
    with open(os.path.join(RES, "r4_concordance.txt"), "w", encoding="utf-8") as fh:
        fh.write(txt + "\n")
    print(f"\n-> data/rhythm/concordance.csv, results/v3/r4_concordance.txt")


if __name__ == "__main__":
    main()
