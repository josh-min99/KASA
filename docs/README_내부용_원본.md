# 중력을 이용한 Circadian Rhythm 조절 — 드라이랩 재현 패키지

연구계획서 작성을 위한 사전 분석·시뮬레이션. 모든 수치는 아래 스크립트 재실행으로 재현된다.

## 문서 읽는 순서

| 문서 | 내용 |
|---|---|
| **`드라이랩_분석보고.md`** | **주 문서.** 방법과 결과 전체. 절별로 데이터·전처리·검정 방법을 명시 |
| `정리_전체요약.md` | 배경 설명용. 통계·생물 용어를 풀어 씀 + 용어 사전 |
| `TODO.md` | 진행 상황 + 역할별 열린 항목 |
| `연구계획_v2.md` | 개정된 연구계획서 초안 |
| `결과_week1-4.md` | 데이터 감사 · 양성대조 · 행동 리듬 · 진동자 모델 |
| `결과_week5-6.md` | 위상 추정기 · 검정력 분석 · 웻랩 설계 확정안 |
| `계획.txt` | 최초 계획 (개정 전 원본, 보관용) |

## 실행

```bash
pip install -r requirements.txt

python run_all.py            # 전체 파이프라인 재실행
python run_all.py --check    # 산출물·핵심 수치만 검증 (재실행 없음)
```

개별 실행도 가능하다. 스크립트는 서로 독립이며 `cache/`를 공유한다.

```bash
python scripts/build_inventory.py       # Week1  인공중력 스터디 인벤토리
python scripts/clock_dose_response.py   # Week1  clock 용량반응 + 최소검출효과크기
python scripts/positive_control.py      # Week1  양성대조 (파이프라인 검증)
python scripts/ground_and_telemetry.py  # Week1-2 지상 스터디 / 텔레메트리 탐색
python scripts/behavior_rhythm.py       # Week2  ISS 행동 리듬 (OSD-952)
python scripts/oscillator_model.py      # Week3-4 진동자 모델 / PRC / Arnold tongue
python scripts/phase_predictor.py       # Week5  위상 추정기 + 검증
python scripts/power_analysis.py        # Week6  검정력 → 표본수
python scripts/make_figures.py          # Week7-8 그림 5장
```

**최초 실행 시 약 400 MB를 내려받는다** (NASA OSDR + NCBI GEO). 인증 불필요.
이후에는 `cache/`를 재사용한다. `cache/`와 `figures/`는 버전관리에서 제외한다.

## 디렉터리

```
scripts/    분석 스크립트 (독립 실행 가능)
data/       산출 CSV — 계획서 표의 원본
figures/    그림 5장 (PNG)
cache/      내려받은 원본 데이터 (git 제외)
```

## 외부 데이터

| 출처 | 용도 | 접근 |
|---|---|---|
| NASA OSDR | 우주비행 마우스 전사체·행동 데이터 | 인증 불필요, REST API |
| NCBI GEO **GSE54650** | 지상 circadian atlas (12조직 × 24시점) — 위상 추정기 학습 | FTP |

> ⚠️ 중력 조건은 `/osdr/data/osd/meta/{id}` JSON API가 **null로 반환한다.**
> 반드시 `ISA.zip` 안의 `s_*.txt`를 파싱할 것. `GLOpenAPI`는 서비스 중단(404).

## 재현성 확인 항목

`run_all.py --check`가 검증하는 값:

| 항목 | 기대값 |
|---|---|
| 인공중력 마우스 스터디 수 | 11 |
| 중력 단계 3개 이상 (용량반응 가능) | 3 |
| 희생 시각이 기록된 스터디 수 | **0** |
| OSD-758 clock 모듈 permutation p | 0.61 |
| OSD-758 군내 SD | 0.343 log2 |
| 망막 DEG (uG → 1G) | 693 → 56 |
| 위상추정 LOTO 중앙값 오차 | 0.70 h |
| 현실조건 위상 상관 (최고) | 0.16 |
| PRC 최대 지연 / 전진 | −4.13 h / +1.87 h |

## 실무 함정 (전부 실제로 겪은 것)

1. **`Arntl` → `Bmal1`** — MGI 심볼 변경. 옛 이름으로 조회하면 에러 없이 조용히 NaN
2. **중력 조건은 ISA.zip에서만** — meta JSON API는 null
3. **GSE54650은 선형 강도** — log2 변환 없이 쓰면 진폭이 수백 단위로 나옴
4. **cosinor 위상은 `atan2(b2, b1)`** — 부호를 뒤집으면 성능이 무작위 수준으로 붕괴
5. **U+2212(−)는 Malgun Gothic에 없음** — 그림 텍스트에는 ASCII 하이픈
