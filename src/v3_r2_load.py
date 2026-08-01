"""
R2: 확보한 원자료를 공통 스키마로 변환 + SHA256 고정.

왜 해시를 고정하는가
  PROGRESS.md [2026-07-31 18:05] 에 표 3 이 오염된 사고가 기록돼 있다.
  원인은 `de_path()` 가 '가장 작은 로컬 파일'을 고르는 방식이라, 다운로드가 진행 중일 때
  실행하면 다른 파일이 선택된 것이었다. 여기서는 파일을 해시로 못 박고,
  이후 모든 스크립트가 manifest 를 검증한 뒤에만 진행한다.

공통 스키마 (long)
  dataset      데이터셋 키
  subject_id   개체(또는 실험단위) 식별자
  species      종
  condition    baseline | treatment | recovery  (자기대조 정규화의 기준)
  gravity      1G | uG | HLU | HDT | 0.001G ...
  t_hours      실험 시작 기준 경과 시간
  clock_hours  하루 중 시각 (0-24). 없으면 NaN
  variable     activity | tb_core | tb_sub | hr | bp_mean
  value        측정값

대상 (R1 감사에서 원자료를 실제로 받은 것만)
  helissen   PASS         HLU 마우스 활동+체온, 2h bin, 개체별
  hdbr       CONDITIONAL  인간 심부체온, 6분 간격 32h, w3/w8 자기대조
  osd595     LEVEL_ONLY   궤도 초파리 활동, 명기 전용 -> 리듬 축 제외

산출: data/rhythm/manifest.csv, data/rhythm/long.csv
"""
import os
import re
import sys
import hashlib
import zipfile
import io

import numpy as np
import pandas as pd

if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
FILES = os.path.join(ROOT, "cache", "rhythm", "files")
OUT = os.path.join(ROOT, "data", "rhythm")
os.makedirs(OUT, exist_ok=True)


def sha256(p):
    h = hashlib.sha256()
    with open(p, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


# ---------------------------------------------------------------- Helissen
# ReadMe: 2h bin, 실험기간 10-11일. 논문(Life 13:844): 대조 3일 / HU 5일 / 회복 2일.
# 파일에 군 라벨이 없으므로 날짜에서 유도해야 한다.
# 유도가 맞는지는 데이터 자체로 검증한다 (HU 시작에서 활동량이 급락해야 한다).
HELISSEN_FILES = {
    2020: "physiological_param_ground_model_microgravity_2020.tab",
    2021: "physiological_param_ground_model_microgravity_2021.tab",
    2022: "physiological_param_ground_model_microgravity_2022.tab",
}


def load_helissen():
    frames = []
    for year, fn in HELISSEN_FILES.items():
        p = os.path.join(FILES, fn)
        # Dataverse 의 ?format=original 은 쉼표 구분 CSV 를 준다 (.tab 확장자지만)
        d = pd.read_csv(p, sep=",")
        d.columns = [c.strip() for c in d.columns]
        d["dt"] = pd.to_datetime(d["Date"], format="%d/%m/%Y %H:%M:%S", errors="coerce")
        d = d.dropna(subset=["dt"])
        d["year"] = year
        frames.append(d)
    raw = pd.concat(frames, ignore_index=True)

    rows = []
    for (year, implant), g in raw.groupby(["year", "Implant"]):
        g = g.sort_values("dt")
        t0 = g.dt.min()
        exp_day = ((g.dt - t0).dt.total_seconds() / 86400).values  # 0 부터
        # 실험 1일차 = day 0. 논문 설계: C-3..C-1 (3일) / HU1..HU5 (5일) / R+1,R+2 (2일)
        day_idx = np.floor(exp_day).astype(int)
        cond = np.where(day_idx < 3, "baseline",
                        np.where(day_idx < 8, "treatment", "recovery"))
        grav = np.where(cond == "treatment", "HLU", "1G")
        var_map = {"A_NPMN(Counts):Activity": "activity",
                   "T_NPMN(Celsius):Temp": "tb_core" if year == 2022 else "tb_sub",
                   "HR(bpm):ECG": "hr",
                   "Mean(mmHg):Pressure": "bp_mean"}
        for col, var in var_map.items():
            if col not in g.columns:
                continue
            v = pd.to_numeric(g[col], errors="coerce").values
            ok = ~np.isnan(v)
            rows.append(pd.DataFrame({
                "dataset": f"helissen{year}", "subject_id": f"{year}_{implant}",
                "species": "mouse", "condition": cond[ok], "gravity": grav[ok],
                "t_hours": exp_day[ok] * 24,
                "clock_hours": g.dt.dt.hour.values[ok] + g.dt.dt.minute.values[ok] / 60,
                "variable": var, "value": v[ok]}))
    return pd.concat(rows, ignore_index=True)


# ---------------------------------------------------------------- HDBR
def load_hdbr():
    z = zipfile.ZipFile(os.path.join(FILES, "HDBR_BBR2_Cocktail_CBTprofiles.zip"))
    rows = []
    for member, study in [("HDBR_BBR2_CBT.csv", "BBR2"), ("HDBR_Cocktail_CBT.csv", "Cocktail")]:
        txt = z.read(member).decode("utf-8", errors="replace")
        d = pd.read_csv(io.StringIO(txt))
        tcol = d.columns[0]
        tmin = pd.to_numeric(d[tcol], errors="coerce").values
        for c in d.columns[1:]:
            m = re.match(r"CBT_(\w+?)\.(w\d+)", c)
            if not m:
                continue
            subj, week = m.groups()
            v = pd.to_numeric(d[c], errors="coerce").values
            ok = ~np.isnan(v)
            # Time [min] 은 하루 중 분(1321분 = 22:01)에서 시작해 32시간 이어진다
            rows.append(pd.DataFrame({
                "dataset": f"hdbr_{study}", "subject_id": f"{study}_{subj}",
                "species": "human",
                # 자기대조: 3주차를 기준, 8주차를 '더 오래 누워 있은 상태' 로 둔다
                "condition": np.where(week == "w3", "baseline", "treatment")[()].repeat(ok.sum())
                if False else (["baseline"] if week == "w3" else ["treatment"]) * int(ok.sum()),
                "gravity": "HDT",
                "t_hours": (tmin[ok] - tmin[0]) / 60,
                "clock_hours": (tmin[ok] / 60) % 24,
                "variable": "tb_core", "value": v[ok]}))
    return pd.concat(rows, ignore_index=True)


# ---------------------------------------------------------------- OSD-595
def load_osd595():
    p = os.path.join(FILES, "LSDS-42_Video_Recording_LSDS-42_video_mhatreTRANSFORMED.csv")
    d = pd.read_csv(p)
    cols = [c for c in d.columns if "integratedpixel" in c]
    rows = []
    for _, r in d.iterrows():
        for c in cols:
            m = re.match(r"day(\d+)_video(\d+)_", c)
            day, vid = int(m.group(1)), int(m.group(2))
            rows.append({
                "dataset": "osd595", "subject_id": str(r["Sample Name"]).replace(" ", ""),
                "species": "drosophila",
                # 궤도상 1G 원심이 대조, uG 가 처리
                "condition": "baseline" if r["condition"] == "1G in centrifuge" else "treatment",
                "gravity": "1G" if r["condition"] == "1G in centrifuge" else "uG",
                "t_hours": day * 24.0 + (vid - 1) * 2.0,   # 명기 12시간을 6등분한 근사
                "clock_hours": np.nan,                      # 실제 시각 미상 -> 위상 계산 금지
                "variable": "activity", "value": float(r[c])})
    return pd.DataFrame(rows)


def main():
    print("=" * 78)
    print("R2  원자료 -> 공통 스키마")
    print("=" * 78)

    # --- manifest
    man = []
    for fn in sorted(os.listdir(FILES)):
        p = os.path.join(FILES, fn)
        if os.path.isfile(p):
            man.append({"file": fn, "bytes": os.path.getsize(p), "sha256": sha256(p)})
    mdf = pd.DataFrame(man)
    mdf.to_csv(os.path.join(OUT, "manifest.csv"), index=False, encoding="utf-8-sig")
    print(f"\nmanifest {len(mdf)}개 파일 고정")
    for _, r in mdf.iterrows():
        print(f"  {r.sha256[:16]}  {r.bytes:>10,}  {r.file}")

    parts = []
    for name, fn in [("helissen", load_helissen), ("hdbr", load_hdbr), ("osd595", load_osd595)]:
        try:
            d = fn()
            parts.append(d)
            print(f"\n[{name}] {len(d):,} 행 / 개체 {d.subject_id.nunique()} / "
                  f"변수 {sorted(d.variable.unique())}")
        except Exception as e:
            print(f"\n[{name}] 실패: {type(e).__name__}: {e}")
            raise

    long = pd.concat(parts, ignore_index=True)
    long.to_csv(os.path.join(OUT, "long.csv"), index=False, encoding="utf-8-sig")

    print("\n" + "-" * 78)
    print("데이터셋 x 조건 x 변수 요약")
    s = (long.groupby(["dataset", "condition", "variable"])
         .agg(n=("value", "size"), subj=("subject_id", "nunique"),
              mean=("value", "mean")).round(3))
    print(s.to_string())

    # --- 군 라벨 유도 검증: HU 시작에서 활동량이 실제로 급락하는가
    print("\n" + "-" * 78)
    print("[검증] Helissen 군 라벨을 날짜에서 유도했다. 데이터로 확인한다.")
    print("       논문: HU 중 평균 활동량 85% 감소. baseline 대비 treatment 비를 본다.")
    for ds in ["helissen2020", "helissen2021", "helissen2022"]:
        a = long[(long.dataset == ds) & (long.variable == "activity")]
        if not len(a):
            continue
        b = a[a.condition == "baseline"].value.mean()
        t = a[a.condition == "treatment"].value.mean()
        r = a[a.condition == "recovery"].value.mean()
        print(f"  {ds}: baseline {b:.4f} -> HU {t:.4f} ({(t/b-1)*100:+.1f}%) -> 회복 {r:.4f}")

    print(f"\n-> {OUT}/long.csv ({len(long):,} 행), {OUT}/manifest.csv")


if __name__ == "__main__":
    main()
