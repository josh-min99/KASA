"""
v2 — 미션 비교를 조직 매칭으로 다시 한다.

문제의식
  표 4(미션별 I4)는 RR-1 8건과 기타 미션 8건을 비교한다. 그런데 두 군의 조직이
  거의 겹치지 않는다.
      RR-1    : 부신·EDL·비복근·신장·사두근·가자미근·전경골근·간  (근육 5 + 장기 3)
      기타 미션: 가자미근×2·비장×2·망막×2·등쪽피부×2            (근육 2 + 비근육 6)
  겹치는 조직은 가자미근 하나뿐이다.
  이것은 v2 에서 v1 결론을 철회할 때 지적한 것과 정확히 같은 교락 구조다
  (REPORT_V2 §1). 같은 실수를 표 4 에서 반복했다.

  게다가 표 3 에서 조직별 I4 를 보면 근육이 망막·피부보다 높다.
  따라서 'RR-1 3.11 vs 기타 0.75' 의 상당 부분이 조직 차이일 수 있다.

이 스크립트가 하는 일
  1) RR-1 8건을 v2 대비 선택 규칙으로 재계산한다.
     (기존 표 4 의 RR-1 열은 v1 파이프라인 산출값이라 기타 미션 열과
      선택 규칙이 다르다. 같은 규칙으로 맞춘다.)
  2) 미션 라벨을 ISA/프로토콜 메타데이터에서 확인한 값으로 붙인다.
  3) 조직을 맞춘 미션 비교를 한다.
  4) 조직 교락의 크기를 정량화한다 (근육 / 비근육 분리).

출력: results/v2/tables/v2_mission_tissue_matched.csv
      results/v2/tables/v2_mission_layers.csv
"""
import os, sys, warnings
import numpy as np
import pandas as pd

warnings.filterwarnings("ignore")
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import v2_tissue_matched as V

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
TAB = os.path.join(ROOT, "results", "v2", "tables")

# RR-1 (Rodent Research-1, SpaceX-4, 2014). Life 2020 이 사용한 8건.
RR1 = [("부신", 98), ("EDL", 99), ("비복근", 101), ("신장", 102),
       ("사두근", 103), ("가자미근", 104), ("전경골근", 105), ("간", 168)]

# 미션 라벨. data/raw 의 ISA·프로토콜 텍스트에서 식별자를 직접 확인한 값이다.
MISSION = {
    "OSD-98": "RR-1", "OSD-99": "RR-1", "OSD-101": "RR-1", "OSD-102": "RR-1",
    "OSD-103": "RR-1", "OSD-104": "RR-1", "OSD-105": "RR-1", "OSD-168": "RR-1",
    "OSD-770": "RR-23", "OSD-506": "RR-23", "OSD-255": "RR-9", "OSD-254": "RR-7",
    "OSD-288": "JAXA MHU-1", "OSD-238": "JAXA MHU-2",
    "OSD-714": "JAXA MHU-1/4/5", "OSD-758": "JAXA MHU-8",
}
# 골격근 여부. 조직 교락을 분리하기 위한 구분이다.
MUSCLE = {"EDL", "비복근", "사두근", "가자미근", "전경골근"}


def main():
    # --- 1) RR-1 을 v2 규칙으로 재계산 -------------------------------------
    print("RR-1 8건을 v2 대비 규칙으로 재계산한다.\n")
    rr = []
    for tis, osd in RR1:
        r = V.analyse(osd, tis, "FLIGHT")
        rr.append(r)
        print(f"  [{tis}] OSD-{osd}: I4={r.get('I4')} DEG={r.get('DEG')} "
              f"채택={r.get('채택')}  대비=({str(r.get('대비_처리'))[:44]})")
    RR = pd.DataFrame(rr)

    # --- 2) 기타 미션 비행 데이터는 조직매칭 실행본에서 가져온다 -------------
    D = pd.read_csv(os.path.join(TAB, "v2_per_dataset.csv"))
    other = D[(D["처리"] == "FLIGHT") & (D["채택"] == True)].copy()
    other = other[~other.OSD_ID.isin(RR.OSD_ID)]

    A = pd.concat([RR[RR["채택"] == True], other], ignore_index=True)
    A["미션"] = A.OSD_ID.map(MISSION)
    A["미션군"] = np.where(A.미션 == "RR-1", "RR-1", "기타 미션")
    A["근육"] = np.where(A.조직.isin(MUSCLE), "근육", "비근육")
    A = A[["미션군", "미션", "조직", "근육", "OSD_ID", "DEG", "I4"]].sort_values(
        ["미션군", "조직", "OSD_ID"])
    A.to_csv(os.path.join(TAB, "v2_mission_tissue_matched.csv"),
             index=False, encoding="utf-8-sig")

    print("\n" + "=" * 92)
    print("비행 데이터셋 전체 (v2 규칙 통일)")
    print("=" * 92)
    print(A.to_string(index=False))

    # --- 3) 조직을 맞춘 미션 비교 ------------------------------------------
    print("\n" + "=" * 92)
    print("조직을 맞춘 미션 비교")
    print("=" * 92)
    rows = []
    for tis, g in A.groupby("조직", sort=False):
        a = g[g.미션군 == "RR-1"]; b = g[g.미션군 == "기타 미션"]
        if len(a) and len(b):
            verdict = ("RR-1 > 기타" if a.I4.min() > b.I4.max()
                       else "기타 > RR-1" if b.I4.min() > a.I4.max() else "겹침")
            print(f"  {tis:<8} RR-1 {list(a.I4)}  vs  기타 {list(b.I4)}   [{verdict}]")
            rows.append({"조직": tis, "n_RR1": len(a), "n_기타": len(b),
                         "RR1_I4": "; ".join(f"{v:.3f}" for v in a.I4),
                         "기타_I4": "; ".join(f"{v:.3f}" for v in b.I4),
                         "판정": verdict})
        else:
            miss = "기타 미션 없음" if not len(b) else "RR-1 없음"
            print(f"  {tis:<8} 비교 불가 ({miss})")
            rows.append({"조직": tis, "n_RR1": len(a), "n_기타": len(b),
                         "RR1_I4": "; ".join(f"{v:.3f}" for v in a.I4) or "-",
                         "기타_I4": "; ".join(f"{v:.3f}" for v in b.I4) or "-",
                         "판정": f"비교 불가({miss})"})
    T = pd.DataFrame(rows)

    # --- 4) 조직 교락의 크기 ----------------------------------------------
    print("\n" + "=" * 92)
    print("조직 교락의 크기 — 층별 중앙값")
    print("=" * 92)
    lay = []

    def add(label, sub):
        if not len(sub):
            return
        lay.append({"층": label, "n": len(sub), "I4_중앙값": round(sub.I4.median(), 3),
                    "최소": round(sub.I4.min(), 3), "최대": round(sub.I4.max(), 3),
                    "조직": ", ".join(sorted(set(sub.조직)))})

    add("RR-1 전체", A[A.미션군 == "RR-1"])
    add("기타 미션 전체", A[A.미션군 == "기타 미션"])
    add("RR-1 · 근육만", A[(A.미션군 == "RR-1") & (A.근육 == "근육")])
    add("기타 미션 · 근육만", A[(A.미션군 == "기타 미션") & (A.근육 == "근육")])
    add("RR-1 · 비근육만", A[(A.미션군 == "RR-1") & (A.근육 == "비근육")])
    add("기타 미션 · 비근육만", A[(A.미션군 == "기타 미션") & (A.근육 == "비근육")])
    L = pd.DataFrame(lay)
    print(L.to_string(index=False))
    L.to_csv(os.path.join(TAB, "v2_mission_layers.csv"), index=False, encoding="utf-8-sig")

    # 비율 요약
    def med(mask):
        s = A[mask]
        return s.I4.median() if len(s) else np.nan
    r_all = med(A.미션군 == "RR-1") / med(A.미션군 == "기타 미션")
    r_mus = (med((A.미션군 == "RR-1") & (A.근육 == "근육")) /
             med((A.미션군 == "기타 미션") & (A.근육 == "근육")))
    print(f"\n  전체 비교의 RR-1/기타 비율      : {r_all:.2f}배")
    print(f"  근육으로 층화한 뒤의 비율       : {r_mus:.2f}배")
    print(f"  -> 차이의 일부는 조직 구성에서 온다. 층화하면 비율이 줄어든다.")

    print("\n" + "=" * 92)
    print(T.to_string(index=False))
    print(f"\n-> {TAB}")


if __name__ == "__main__":
    main()
