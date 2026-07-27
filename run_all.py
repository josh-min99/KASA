"""
전체 파이프라인 재현 실행기.

  python run_all.py            전체 실행
  python run_all.py --check    산출물 존재 여부와 핵심 수치만 검증 (재실행 없음)

각 스크립트는 독립 실행 가능하며 cache/ 를 공유한다.
최초 실행 시 약 400MB 를 내려받는다(OSDR + GEO). 이후에는 캐시를 쓴다.
"""
import os
import subprocess
import sys
import time

ROOT = os.path.dirname(os.path.abspath(__file__))
SCRIPTS = os.path.join(ROOT, "scripts")
DATA = os.path.join(ROOT, "data")
FIGS = os.path.join(ROOT, "figures")

# (스크립트, 설명, 이 스크립트가 만드는 대표 산출물)
STEPS = [
    ("build_inventory.py",     "Week1 인공중력 스터디 인벤토리",   ["studies_gravity.csv"]),
    ("clock_dose_response.py", "Week1 clock 용량반응 + MDE",        ["preliminary_summary.csv"]),
    ("positive_control.py",    "Week1 양성대조",                    ["positive_control_deg_ladder.csv",
                                                                     "positive_control_competitive.csv"]),
    ("ground_and_telemetry.py","Week1-2 지상 스터디 / 텔레메트리",   ["ground_studies.csv",
                                                                     "telemetry_candidates.csv"]),
    ("behavior_rhythm.py",     "Week2 ISS 행동 리듬",               ["behavior_rhythm_summary.csv"]),
    ("oscillator_model.py",    "Week3-4 진동자 모델 / PRC / 동조",   ["model_prc.csv",
                                                                     "model_arnold_tongue.csv"]),
    ("phase_predictor.py",     "Week5 위상 추정기",                 ["phase_predictor_loto.csv",
                                                                     "phase_predictor_relative.csv"]),
    ("power_analysis.py",      "Week6 검정력 분석",                 ["power_curve.csv"]),
    ("make_figures.py",        "Week7-8 그림 5장",                  []),
]

FIG_FILES = ["fig1_data_map.png", "fig2_positive_control.png", "fig3_phase_inference.png",
             "fig4_model.png", "fig5_power.png"]


def run_all():
    fails = []
    for script, desc, _ in STEPS:
        print("=" * 74)
        print(f"[실행] {script}  —  {desc}")
        t0 = time.time()
        r = subprocess.run([sys.executable, os.path.join(SCRIPTS, script)],
                           capture_output=True, text=True, encoding="utf-8", errors="replace")
        dt = time.time() - t0
        if r.returncode == 0:
            print(f"  완료 ({dt:.1f}s)")
        else:
            print(f"  실패 (exit {r.returncode})")
            print((r.stderr or "")[-1500:])
            fails.append(script)
    print("=" * 74)
    print("실패 없음" if not fails else f"실패: {fails}")
    return 1 if fails else 0


def check():
    ok = True
    print("[산출물 확인]")
    for _, _, outs in STEPS:
        for o in outs:
            p = os.path.join(DATA, o)
            exists = os.path.exists(p)
            ok &= exists
            print(f"  {'OK ' if exists else 'MISSING'}  data/{o}")
    for f in FIG_FILES:
        p = os.path.join(FIGS, f)
        exists = os.path.exists(p)
        ok &= exists
        print(f"  {'OK ' if exists else 'MISSING'}  figures/{f}")

    print("\n[핵심 수치 검증]")
    import pandas as pd

    def expect(label, actual, target, tol):
        good = abs(actual - target) <= tol
        print(f"  {'OK ' if good else 'FAIL'}  {label}: {actual} (기대 {target}±{tol})")
        return good

    inv = pd.read_csv(os.path.join(DATA, "studies_gravity.csv"))
    ok &= expect("인공중력 스터디 수", len(inv), 11, 0)
    ok &= expect("용량반응 가능 스터디 수", int(inv.dose_response_usable.sum()), 3, 0)
    ok &= expect("희생 시각 기록 스터디 수", int(inv.time_of_day_recorded.sum()), 0, 0)

    pre = pd.read_csv(os.path.join(DATA, "preliminary_summary.csv")).set_index("study")
    ok &= expect("OSD-758 clock permutation p", float(pre.loc["OSD-758", "permutation_p"]), 0.61, 0.08)
    ok &= expect("OSD-758 군내 SD (log2)", float(pre.loc["OSD-758", "within_group_sd_log2"]), 0.343, 0.005)

    lad = pd.read_csv(os.path.join(DATA, "positive_control_deg_ladder.csv"))
    ret = lad[lad.tissue == "Retina"].set_index("gravity_g").n_DEG_fdr05
    ok &= expect("망막 DEG uG", int(ret.loc[0.0]), 693, 0)
    ok &= expect("망막 DEG 1G", int(ret.loc[1.0]), 56, 0)

    loto = pd.read_csv(os.path.join(DATA, "phase_predictor_loto.csv"))
    ok &= expect("위상추정 LOTO 중앙값 오차(h)", round(float(loto.abs_err_h.median()), 2), 0.70, 0.15)
    rel = pd.read_csv(os.path.join(DATA, "phase_predictor_relative.csv"))
    ok &= expect("현실조건 최고 상관", round(float(rel.mean_corr_true_vs_pred.max()), 2), 0.16, 0.06)

    prc = pd.read_csv(os.path.join(DATA, "model_prc.csv"))
    b = prc[prc.condition.str.contains("약한")]
    ok &= expect("PRC 최대 지연(h)", round(float(b.phase_shift_h.min()), 2), -4.13, 0.30)
    ok &= expect("PRC 최대 전진(h)", round(float(b.phase_shift_h.max()), 2), 1.87, 0.30)

    print("\n" + ("전체 검증 통과" if ok else "검증 실패 항목 있음"))
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(check() if "--check" in sys.argv else run_all())
