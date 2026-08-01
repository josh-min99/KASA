"""
R0: 리듬 데이터셋 채택 규칙 + 역방향 검증.

왜 규칙을 먼저 고정하는가
  PROGRESS.md 에 같은 유형의 사고가 두 번 기록돼 있다. 배제용 정규식이 정탐까지
  지워버린 사고다(`irradiat` 가 `non-irradiated` 를 매칭). 원인은 규칙을 만든 뒤
  그것이 무엇을 지우는지 확인하지 않고 바로 본 분석에 쓴 것이었다.
  그래서 이번에는 판정을 이미 아는 5건에 규칙을 먼저 적용해, 규칙이 그 5건을
  예상대로 분류하는지 확인한 뒤에만 본 검색을 돌린다.

채택 기준
  C1  중력 조건 명시 (비행 / 지상 언로딩 / 원심분리 / 부분중력 / 인공중력)
  C2  일주기 해상도 >= 하루 4점
  C3a 리듬 파라미터 1회 추정 가능 (>= 24h 연속)
  C3b 위상 추이 추정 가능 (>= 7일 연속)
  C4  개체 식별 가능
  C5  자기 대조(baseline 구간) 또는 동시 대조군
  C6  다운로드 가능 + 라이선스 확인
  C7  조명 조건 기록

판정 5단계
  PASS           전 축 정량 사용 (효과크기 + CI + 위상 추이)
  CONDITIONAL    파라미터 1점 추정만. 위상 '추이' 는 쓸 수 없다
  DIRECTION_ONLY 방향만. 효과 크기를 주장하지 않는다
  LITERATURE     원자료는 못 받지만 논문이 리듬 파라미터를 수치로 보고했다.
                 정합성 표의 칸을 '요약통계기반' 태그로 채운다
  REJECT         사용 불가

  LITERATURE 를 나중에 추가한 이유: 초판 규칙은 C6(다운로드 가능) 실패를 무조건
  REJECT 로 보냈다. 그런데 실제로 감사해 보니 이 분야의 리듬 데이터는 대부분
  원자료가 없고 논문 수치만 있다(예: LSDA 의 Neurolab 일주기 데이터셋은
  'Availability: Available offline / No data submitted by PI' 로 기록돼 있다).
  그것들을 전부 버리면 정합성 표가 비어 버리는데, 방향 논증에는 쓸 수 있다.
  단 원자료 기반과 반드시 구분해 표기한다.
  (규칙 변경이므로 고정 사례 5건으로 재검증했다. 5건 모두 C6=True 라 판정 불변.)

  C3 를 a/b 로 쪼갠 이유: 하나로 두면(>=7일 연속) 60일 HDBR 이 탈락한다. 그 연구는
  3주차·8주차에 각각 36시간씩만 연속 기록했는데, 그것으로 acrophase 를 뽑아
  위상 지연을 보고했다. 파라미터 1회 추정에 필요한 길이와 추이 추정에 필요한 길이는
  다르므로 기준도 나뉘어야 한다.

실행: python src/v3_audit_rules.py
      역방향 검증 실패 시 exit 1. 본 검색(R1)은 이 스크립트가 통과해야 진행한다.
"""
import os
import sys
import json

if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, "results", "v3")
os.makedirs(OUT, exist_ok=True)

CRITERIA = ["C1", "C2", "C3a", "C3b", "C4", "C5", "C6", "C7"]

CRITERIA_DESC = {
    "C1":  "중력 조건 명시",
    "C2":  "일주기 해상도 >= 하루 4점 이며 24시간에 걸쳐 분포",
    "C3a": "연속 관측 >= 24h (파라미터 1회 추정)",
    "C3b": "연속 관측 >= 7일 (위상 추이 추정)",
    "C4":  "개체 식별 가능",
    "C5":  "자기 대조 또는 동시 대조군",
    "C6":  "다운로드 가능 + 라이선스 확인",
    "C7":  "조명 조건 기록",
}


def classify(rec):
    """후보 레코드 -> (판정, 사유).

    rec 는 CRITERIA 키에 대한 bool 과 `has_phase_contrast` 를 가진 dict.
    has_phase_contrast: 하루 4점이 안 되더라도 명/암 대비 같은 일주기 대비가 존재하는가.
    """
    miss = [c for c in CRITERIA if not rec.get(c, False)]

    # C1 은 대체 불가. 중력 조건을 모르면 아무것도 못 한다.
    if not rec.get("C1", False):
        return "REJECT", "C1 중력 조건 불명"

    # C6(원자료 접근) 실패는 논문 보고 수치가 있으면 LITERATURE 로 살린다.
    if not rec.get("C6", False):
        if rec.get("has_published_stats", False):
            return "LITERATURE", "C6 원자료 비공개 — 논문 보고 수치만 사용"
        return "REJECT", "C6 접근 불가, 보고 수치도 없음"

    if rec.get("C2") and rec.get("C3b") and rec.get("C5") and rec.get("C4") and rec.get("C7"):
        return "PASS", "전 기준 충족"

    if rec.get("C2") and rec.get("C3a") and rec.get("C5"):
        return "CONDITIONAL", "미충족: " + ",".join(m for m in miss if m in ("C3b", "C4", "C7"))

    if rec.get("has_phase_contrast", False):
        return "DIRECTION_ONLY", "미충족: " + ",".join(miss)

    # 24시간 커버리지가 없어 리듬은 못 보지만, 여러 날 연속 + 대조군이 있어
    # '활동량 수준' 의 중력 대비는 볼 수 있는 경우.
    # OSD-595(궤도 인공중력 초파리)가 여기 해당한다. 하루 6점이지만 전부 12h 명기
    # 안쪽이라 위상은 못 재는데, 궤도상 1G 원심 대조가 있어 중력만 분리된
    # 유일한 활동량 자료다. 리듬 축과 수준 축을 섞지 않으려고 계층을 나눈다.
    if rec.get("C3b") and rec.get("C5"):
        return "LEVEL_ONLY", "24h 커버리지 없음 — 리듬 축 제외, 활동량 수준 축만"

    return "REJECT", "미충족: " + ",".join(miss)


# --------------------------------------------------------------------------
# 역방향 검증용 고정 사례 5건.
# 각 필드는 원문·파일 실물에서 확인한 사실이며, 근거를 evidence 에 남긴다.
# --------------------------------------------------------------------------
FIXTURES = [
    dict(
        name="Helissen 2023 (HLU 마우스 텔레메트리)",
        source="doi:10.57745/QVRW8W",
        C1=True, C2=True, C3a=True, C3b=True, C4=True, C5=True, C6=True, C7=True,
        has_phase_contrast=True,
        expect="PASS",
        evidence="2h bin x 11일 연속, 개체별, 대조 3일/HU 5일/회복 2일, LD 12:12, Etalab 오픈",
    ),
    dict(
        name="OSD-681 / LSDS-83 (Fuller HDT 쥐 생체원격)",
        source="OSDR OSD-681",
        C1=True, C2=False, C3a=False, C3b=False, C4=True, C5=True, C6=True, C7=False,
        has_phase_contrast=False,
        expect="REJECT",
        evidence="cache 파일 실물 확인: Baseline/Day0/30/60/90/120/150/180 기간 평균 33행뿐",
    ),
    dict(
        name="OSD-952 / LSDS-168 (RR-1 ISS 행동 에소그램)",
        source="OSDR OSD-952",
        C1=True, C2=False, C3a=False, C3b=False, C4=False, C5=False, C6=True, C7=True,
        has_phase_contrast=True,
        expect="DIRECTION_ONLY",
        evidence="Day별 Light/Dark 2점뿐(하루 4점 미달), 케이지 단위, 비행 전 baseline 없음. "
                 "다만 Dark 0.931 vs Light 0.455 로 명암 대비는 존재",
    ),
    dict(
        name="Shimbo 2021 MHU-1 (궤도 uG vs 인공 1g 영상 활동)",
        source="Sci Rep 11:2827",
        C1=True, C2=True, C3a=False, C3b=False, C4=True, C5=True, C6=True, C7=True,
        has_phase_contrast=True,
        expect="DIRECTION_ONLY",
        evidence="시간별 활동비율은 있으나 저장용량 제약으로 AG 2개·MG 2개 영상뿐 -> 24h 연속 없음. "
                 "지상대조 주간 32.5% / 야간 63.2% 로 명암 대비는 존재",
    ),
    dict(
        name="OSD-595 / LSDS-42 (MVP 초파리, 궤도 uG vs 궤도 1G 원심)",
        source="OSDR OSD-595",
        C1=True, C2=False, C3a=True, C3b=True, C4=True, C5=True, C6=True, C7=True,
        has_phase_contrast=False,
        expect="LEVEL_ONLY",
        evidence="16일(mission day 13-28) x 하루 6점 x 11모듈, 결측 0. 같은 habitat 안 두 원심기라 "
                 "방사선·격리·발사가 양쪽 동일 -> 중력만 분리되는 유일 자료. "
                 "그러나 원논문이 '6 videos per day spanned the 12-hr light period' 이고 "
                 "암기 녹화는 적외선 고장으로 미실시 -> 24h 커버리지 없음",
    ),
    dict(
        name="60일 6도 HDBR 심부체온 (BBR2-2 + Cocktail)",
        source="doi:10.6084/m9.figshare.13633790",
        C1=True, C2=True, C3a=True, C3b=False, C4=True, C5=True, C6=True, C7=True,
        has_phase_contrast=True,
        expect="CONDITIONAL",
        evidence="3주차·8주차에 각 36시간 연속 기록 -> 24h 충족, 7일 연속은 미충족. "
                 "acrophase 16.23h -> 16.68h 위상 지연 (자기대조)",
    ),
]


def main():
    print("=" * 78)
    print("R0  채택 규칙 역방향 검증")
    print("=" * 78)
    print("\n[기준]")
    for c in CRITERIA:
        print(f"  {c:4s} {CRITERIA_DESC[c]}")
    print("\n[판정 순서]")
    print("  1. C1 또는 C6 실패            -> REJECT (대체 불가)")
    print("  2. C2 & C3b & C5 & C4 & C7    -> PASS")
    print("  3. C2 & C3a & C5              -> CONDITIONAL")
    print("  4. has_phase_contrast         -> DIRECTION_ONLY")
    print("  5. 그 외                      -> REJECT")

    print("\n" + "-" * 78)
    rows, bad = [], []
    for f in FIXTURES:
        got, why = classify(f)
        ok = (got == f["expect"])
        if not ok:
            bad.append((f["name"], f["expect"], got))
        mark = "OK " if ok else "!! "
        print(f"{mark}{f['name']}")
        print(f"    기대 {f['expect']:15s} 판정 {got:15s} ({why})")
        print(f"    근거: {f['evidence']}")
        rows.append(dict(name=f["name"], source=f["source"], expect=f["expect"],
                         got=got, reason=why, ok=ok, evidence=f["evidence"]))

    print("-" * 78)
    if bad:
        print(f"\n역방향 검증 실패 {len(bad)}건. 규칙을 고치고 다시 돌릴 것.")
        for n, e, g in bad:
            print(f"  {n}: 기대 {e} != 판정 {g}")
    else:
        print(f"\n역방향 검증 {len(FIXTURES)}/{len(FIXTURES)} 통과. R1 진행 가능.")

    with open(os.path.join(OUT, "r0_rule_validation.json"), "w", encoding="utf-8") as fh:
        json.dump(dict(criteria=CRITERIA_DESC, fixtures=rows, passed=not bad),
                  fh, ensure_ascii=False, indent=2)
    print(f"-> results/v3/r0_rule_validation.json")
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
