"""
v2 보조 — 어제의 'FLIGHT I4 = 3.11' 이 무엇에서 왔는지 확인.

조직 매칭 결과를 보면 비행 데이터셋의 I4 가 크게 갈린다.
어제 FLIGHT-A 로 쓴 8개는 전부 RR-1 미션이었고 I4 가 1.90~5.28 이었다.
그런데 이번에 추가된 다른 미션의 비행 데이터셋은 대부분 1.0 미만이다.
즉 어제 수치가 '우주비행 일반' 이 아니라 'RR-1 특이' 일 가능성이 있다.

RR-1 소속 여부로 나눠 비교한다.
출력: results/v2/tables/v2_mission_effect.csv
"""
import os, warnings
import numpy as np
import pandas as pd

warnings.filterwarnings("ignore")
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
TAB_V2 = os.path.join(ROOT, "results", "v2", "tables")
TAB_V1 = os.path.join(ROOT, "results", "tables")

# 어제 FLIGHT-A 로 쓴 RR-1 8조직
RR1 = {98, 99, 101, 102, 103, 104, 105, 168}


def main():
    v2 = pd.read_csv(os.path.join(TAB_V2, "v2_per_dataset.csv"))
    v1 = pd.read_csv(os.path.join(TAB_V1, "stage4_normalized.csv"))

    # v1 에서 FLIGHT-A(RR-1) 값 가져오기
    a = v1[v1.series == "FLIGHT-A"][["OSD_ID", "tissue", "clock_enrichment", "n_DEG_fdr05"]]
    a = a.rename(columns={"clock_enrichment": "I4", "n_DEG_fdr05": "DEG"})
    a["mission"] = "RR-1"
    a["source"] = "v1 (어제)"

    # v2 에서 채택된 비행 데이터셋
    b = v2[(v2["채택"] == True) & (v2["처리"] == "FLIGHT")][
        ["OSD_ID", "조직", "I4", "DEG"]].rename(columns={"조직": "tissue"})
    b["osd_num"] = b.OSD_ID.str.extract(r"(\d+)").astype(int)
    b["mission"] = np.where(b.osd_num.isin(RR1), "RR-1", "기타 미션")
    b["source"] = "v2 (조직매칭)"

    both = pd.concat([a.assign(osd_num=a.OSD_ID.str.extract(r"(\d+)").astype(int)[0]), b],
                     ignore_index=True)
    both = both.drop_duplicates("OSD_ID", keep="first")
    both = both[["OSD_ID", "tissue", "mission", "I4", "DEG", "source"]].sort_values(
        ["mission", "I4"], ascending=[True, False])

    both.to_csv(os.path.join(TAB_V2, "v2_mission_effect.csv"),
                index=False, encoding="utf-8-sig")

    print("=" * 88)
    print("비행 데이터셋의 I4 — 미션별")
    print("=" * 88)
    print(both.to_string(index=False))

    print("\n" + "=" * 88)
    g = both.groupby("mission").I4.agg(["count", "median", "min", "max"])
    print(g.round(3).to_string())

    rr1 = both[both.mission == "RR-1"].I4
    oth = both[both.mission == "기타 미션"].I4
    print(f"\nRR-1 {len(rr1)}건 중앙값 {rr1.median():.3f} (범위 {rr1.min():.3f}~{rr1.max():.3f})")
    print(f"기타 {len(oth)}건 중앙값 {oth.median():.3f} (범위 {oth.min():.3f}~{oth.max():.3f})")
    print(f"기타 미션 중 I4 > 1.0 인 것: {int((oth > 1.0).sum())}/{len(oth)}")

    # HLU 와 비교
    h = v2[(v2["채택"] == True) & (v2["처리"] == "HLU")].I4
    print(f"\n참고 — HLU {len(h)}건 중앙값 {h.median():.3f} (범위 {h.min():.3f}~{h.max():.3f})")
    print(f"기타 미션 비행 vs HLU: 중앙값 {oth.median():.3f} vs {h.median():.3f}")


if __name__ == "__main__":
    main()
