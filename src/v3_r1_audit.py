"""
R1-c: 리듬 데이터셋 감사표.

R0 에서 고정하고 6건으로 역방향 검증한 규칙(`v3_audit_rules.classify`)을
R1 에서 찾은 후보 전건에 적용한다.

탐색 범위 (전부 실행 기록이 남아 있다)
  - OSDR: 스터디 633건 / 파일 157,801개 전수 (`src/v3_r1_osdr_sweep.py`)
           검색 API 만으로는 411건밖에 안 잡혀 222건을 놓쳤다. ID 전수 순회로 교정했다.
  - NASA LSDA/NLSP: 데이터셋 2,593건. STS-90(Neurolab) 146건 전건 확인
  - Recherche Data Gouv / figshare / GitHub / 논문 부속자료
  - 문헌: PubMed·Google Scholar 교차 검색 (중력 조건 11개 x 리듬 지표 8개)

각 행의 C1~C7 은 전부 **원문 또는 파일 실물에서 확인한 사실**이며 근거를 evidence 에 적었다.
추정으로 채운 칸은 없다. 확인하지 못한 것은 그 사실을 evidence 에 적고 보수적으로 False.

산출: data/rhythm/audit.csv, results/v3/r1_audit_summary.txt
"""
import os
import sys

import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from v3_audit_rules import classify, CRITERIA  # noqa: E402

if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, "data", "rhythm")
RES = os.path.join(ROOT, "results", "v3")
os.makedirs(OUT, exist_ok=True)
os.makedirs(RES, exist_ok=True)


def R(name, arm, species, variables, source, evidence,
      C1, C2, C3a, C3b, C4, C5, C6, C7, phase_contrast, pub_stats):
    return dict(name=name, arm=arm, species=species, variables=variables, source=source,
                evidence=evidence, C1=C1, C2=C2, C3a=C3a, C3b=C3b, C4=C4, C5=C5, C6=C6,
                C7=C7, has_phase_contrast=phase_contrast, has_published_stats=pub_stats)


CANDIDATES = [
    # ---------------------------------------------------------------- 지상 언로딩 (설치류)
    R("Helissen 2023 HLU 마우스 텔레메트리", "지상 언로딩", "마우스",
      "활동량, 피하온도, 심부온도, 심박, 혈압",
      "Life 13:844 / doi:10.57745/QVRW8W",
      "2h bin x 10-11일 연속, 개체별(implant ID). 대조 3일 / HU 5일 / 회복 2일. "
      "LD 12:12(07-19시 점등). Etalab 2.0 오픈. 파일 4개 실제 다운로드 완료. "
      "2020/2021 은 피하온도+활동+혈압(각 n=5), 2022 는 심부온도+활동(n=6)",
      True, True, True, True, True, True, True, True, True, True),

    R("HLU 쥐 5군 분해 (사회격리/구속/언로딩)", "지상 언로딩", "쥐",
      "활동량, 체온, 혈압, 심박",
      "PMC3388052",
      "30분 평균 연속, 3일 baseline + 14일, LD 12:12. 5군(쌍사육/개별대조/언로딩케이지 개별/"
      "수평구속/HDT거상)으로 스트레스 성분을 분해한 유일한 설계. "
      "원자료 공개 진술 없음. 논문에 리듬 기술 수치 있음",
      True, True, True, True, True, True, False, True, True, True),

    R("HLU 마우스 종단 에소그램", "지상 언로딩", "마우스",
      "행동 범주 7종",
      "PMC12842622",
      "관찰이 오후 1-5시에만 이뤄졌고 일주기 분석을 수행하지 않았다고 원문에 명시. "
      "원자료는 요청 시 제공",
      True, False, False, False, True, True, False, True, False, False),

    # ---------------------------------------------------------------- 지상 언로딩 (인간)
    R("60일 6도 HDBR 심부체온 (BBR2-2 + Cocktail)", "지상 언로딩", "인간",
      "심부체온",
      "doi:10.6084/m9.figshare.13633790",
      "3주차·8주차에 각 32시간 연속(6분 간격 320점), 13명(BBR2-2 n=5, Cocktail n=8), "
      "개체별 열, 자기대조(w3 vs w8). CC BY 4.0. 실제 다운로드 완료. "
      "7일 연속 기록은 없어 위상 '추이'는 불가",
      True, True, True, False, True, True, True, True, True, True),

    R("AGBRESA 60일 HDBR + 단완 원심분리 (수면)", "인공중력 되돌림", "인간",
      "PSG (EEG/ECG/EMG/EOG)",
      "npj Microgravity 10 (PMC11621691)",
      "지상 아날로그에 중력을 되돌려 넣은 유일한 통제 실험(대조/연속AG/간헐AG 각 n=8). "
      "그러나 PSG 를 8개 야간에만 측정했고 연속 다일 기록도 actigraphy 도 없다. "
      "원문에 일주기 위상·진폭·acrophase 보고가 없다. 원자료는 요청 시 제공",
      True, False, False, False, True, True, False, True, False, False),

    # ---------------------------------------------------------------- 과중력 (지상)
    R("Martin 2020 2G 전정자극 재동조", "지상 과중력", "쥐",
      "심부체온, 활동량",
      "Sci Rep 10:8646 (PMC7280278)",
      "복강 telemetry 로 6일 baseline + 실험기간 연속 기록. sham n=24 / 양측전정손실 n=24. "
      "2G 1시간 펄스를 매일 7회. LD 6시간 전진 후 재동조 속도 측정. "
      "재동조 16.0+-3.1일(2G) vs 20.4+-0.8일(LD만), p=0.002. 원자료는 요청 시 제공",
      True, True, True, True, True, True, False, True, True, True),

    R("Fuller 만성 원심분리(2G) 쥐 일주기 교란", "지상 과중력", "쥐",
      "심부체온, 활동량",
      "J Appl Physiol (Fuller 1994 계열)",
      "2G 노출 초기 7-10일간 활동 리듬 소실·체온 진폭 감소 후 회복. "
      "1990년대 연구로 원자료 공개 없음. 논문 보고 수치만",
      True, True, True, True, True, True, False, True, True, True),

    # ---------------------------------------------------------------- 실제 비행
    R("OSD-952 / LSDS-168 RR-1 ISS 행동 에소그램", "실제 비행", "마우스",
      "행동 8범주 (DeepEthogram)",
      "OSDR OSD-952",
      "미세중력에서 리듬을 직접 측정한 유일한 공개 데이터. Day 0-32, 고유 프레임 1,175,486. "
      "그러나 Day별 Light/Dark 2점뿐이라 하루 4점 미달, 케이지 단위(개체 식별 불가), "
      "비행 전 baseline 없음. Dark 0.931 vs Light 0.455 로 명암 대비는 존재",
      True, False, False, False, False, False, True, True, True, True),

    R("OSD-595 / LSDS-42 MVP 초파리 (궤도 uG vs 궤도 1G)", "실제 비행 (중력만 분리)", "초파리",
      "활동량 (최대투영 강도)",
      "OSDR OSD-595 / Cell Rep 40:111279",
      "같은 habitat 안 두 원심기(68rpm=1G vs 2.2rpm=0.00095G)라 방사선·격리·발사·온도가 "
      "양쪽 동일 -> 공개 자료 중 중력만 분리되는 유일한 설계. "
      "16일(mission day 13-28) x 하루 6점 x 11모듈(1G 6 / uG 5), 결측 0. LD 12:12. "
      "그러나 원문이 '6 videos per day spanned the 12-hr light period' 이고 "
      "암기 녹화는 적외선 고장으로 미실시 -> 24h 커버리지 없음. 리듬 축 사용 불가",
      True, False, True, True, True, True, True, True, False, True),

    R("OSD-596 / LSDS-43 MVP 초파리 climbing", "실제 비행", "초파리",
      "climbing 성공률",
      "OSDR OSD-596",
      "착륙 후 단일 시점 행동검사. 시계열 아님",
      True, False, False, False, True, True, True, True, False, False),

    R("Shimbo 2021 MHU-1 (궤도 uG vs 인공 1g)", "실제 비행 (중력만 분리)", "마우스",
      "활동비율 (AIS 영상)",
      "Sci Rep 11:2827",
      "궤도상 1g 원심 대조가 있는 마우스 설계(각 n=6, 34일). "
      "그러나 저장용량 제약으로 AG 2개·MG 2개 영상만 남아 24h 연속 없음. "
      "지상대조 주간 32.5% / 야간 63.2% 로 명암 대비는 존재",
      True, True, False, False, True, True, True, True, True, True),

    R("Neurolab STS-90 CNS Control of Rhythms", "실제 비행", "쥐",
      "심부체온, 심박, 활동량",
      "LSDA AINV0000000939 / 실험 9301132",
      "일주기를 목적으로 설계된 유일한 비행 실험. NBS 로 16일 연속 체온·심박 기록. "
      "그러나 LSDA 레코드에 'Availability: Available offline', "
      "'Dataset Description: No data submitted by PI. Refer to publication citations' "
      "로 명시 -> 원자료 존재하지 않음. STS-90 데이터셋 146건 전건 확인",
      True, True, True, True, True, True, False, True, True, True),

    R("Cosmos 2044 / 2229 영장류", "실제 비행", "붉은털원숭이",
      "체온(액와·뇌), 심박, 활동량",
      "J Appl Physiol 81:188 (Fuller 1996)",
      "비행 중 체온 리듬의 위상 지연과 평균 심박 감소를 보고. n=4. "
      "1990년대 연구로 원자료 공개 없음. 논문 보고 수치만",
      True, True, True, True, True, True, False, True, True, True),

    R("우주인 24h Holter 심박변이 일주기", "실제 비행", "인간",
      "심박변이(HRV)",
      "Sci Rep 8:10381 / PMC5192238",
      "6개월 이상 체류 우주인 7명, 발사 전/1개월/2개월/귀환 2주 전/착륙 3개월 후 "
      "각 24시간 Holter. 24·12·8시간 성분이 비행 중에도 유지됨을 보고. "
      "원자료 공개 없음",
      True, True, True, False, True, True, False, False, True, True),

    # ---------------------------------------------------------------- 탈락 (기록용)
    R("OSD-681 / LSDS-83 HDT 쥐 생체원격", "지상 언로딩", "쥐",
      "체온, 두개내압",
      "OSDR OSD-681",
      "cache 파일 실물 확인: Baseline/Day0/30/60/90/120/150/180 의 기간 평균 33행뿐. "
      "일주기 내 해상도 없음",
      True, False, False, False, True, True, True, False, False, False),

    R("OSD-691 ATPase Activity (오탐)", "실제 비행", "마우스",
      "ATPase 활성",
      "OSDR OSD-691",
      "파일명 토큰 'activity' 로 걸렸으나 근육 ATPase 효소활성 측정. 행동·리듬과 무관. "
      "탐색 규칙의 오탐 사례로 기록",
      True, False, False, False, True, True, True, False, False, False),
]


def main():
    rows = []
    for c in CANDIDATES:
        verdict, why = classify(c)
        r = dict(c)
        r["verdict"] = verdict
        r["reason"] = why
        rows.append(r)

    df = pd.DataFrame(rows)
    cols = (["name", "arm", "species", "variables", "verdict", "reason"] +
            CRITERIA + ["has_phase_contrast", "has_published_stats", "source", "evidence"])
    df = df[cols]
    df.to_csv(os.path.join(OUT, "audit.csv"), index=False, encoding="utf-8-sig")

    lines = []
    P = lines.append
    P("=" * 96)
    P("R1  리듬 데이터셋 감사 결과")
    P("=" * 96)
    P("")
    P("탐색 범위")
    P("  OSDR            스터디 633건 / 파일 157,801개 전수")
    P("  NASA LSDA/NLSP  데이터셋 2,593건 (STS-90 Neurolab 146건 전건 확인)")
    P("  외부 저장소     Recherche Data Gouv / figshare / GitHub / 논문 부속자료")
    P("  문헌            중력 조건 11개 x 리듬 지표 8개 교차 검색")
    P("")
    P(f"후보 {len(df)}건 판정")
    for v, n in df.verdict.value_counts().items():
        P(f"  {v:16s} {n}건")
    P("")
    P("-" * 96)
    for v in ["PASS", "CONDITIONAL", "LEVEL_ONLY", "DIRECTION_ONLY", "LITERATURE", "REJECT"]:
        sub = df[df.verdict == v]
        if not len(sub):
            continue
        P(f"\n[{v}]  {len(sub)}건")
        for _, r in sub.iterrows():
            P(f"  - {r['name']}")
            P(f"      arm={r['arm']} / {r['species']} / {r['variables']}")
            P(f"      출처: {r['source']}")
            P(f"      사유: {r['reason']}")
    P("")
    P("=" * 96)
    P("게이트 G1 판정")
    usable = df[df.verdict.isin(["PASS", "CONDITIONAL"])]
    P(f"  원자료로 정량 사용 가능(PASS+CONDITIONAL): {len(usable)}건")
    for _, r in usable.iterrows():
        P(f"    {r['verdict']:12s} {r['name']}")
    if len(usable) >= 3:
        P("  -> 3건 이상. R4 정합성 분석을 계획대로 진행")
    else:
        P("  -> 3건 미만. 계획서에 정한 분기에 따라 R4 는 '요약통계 기반 정합성 표'로 축소하고")
        P("     R5(검정력)·R6(역문제)에 자원을 몰아준다. 감사 결과 자체를 산출물로 삼는다")
    P("")
    P("주목할 부정 결과 (계획서에 그대로 쓸 것)")
    P("  1. 일주기를 목적으로 설계된 유일한 비행 실험(Neurolab STS-90)은")
    P("     NASA 아카이브에 'No data submitted by PI' 로 기록돼 원자료가 존재하지 않는다")
    P("  2. 지상 아날로그에 중력을 되돌려 넣은 유일한 통제 실험(AGBRESA)은")
    P("     일주기 위상을 측정하지 않았다 (PSG 8개 야간만, actigraphy 없음)")
    P("  3. OSDR 전체 157,801개 파일 중 리듬 후보는 52개(5개 스터디)뿐이고,")
    P("     그중 24시간 커버리지를 갖춘 것은 0건이다")

    txt = "\n".join(lines)
    print(txt)
    with open(os.path.join(RES, "r1_audit_summary.txt"), "w", encoding="utf-8") as fh:
        fh.write(txt + "\n")
    print(f"\n-> data/rhythm/audit.csv, results/v3/r1_audit_summary.txt")


if __name__ == "__main__":
    main()
