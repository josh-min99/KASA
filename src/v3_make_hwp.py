"""
진짜 한글 문서(.hwp) 생성 — 설치된 한글을 COM 자동화로 구동한다.

왜 이 방식인가
  HWP 5.0 은 OLE 복합문서 기반 폐쇄 포맷이라 외부에서 유효한 파일을 만들기 어렵다.
  이 PC 에 한글이 설치돼 있고 HWPFrame.HwpObject COM 인터페이스가 등록돼 있으므로,
  한글 자신에게 변환을 시키는 것이 가장 확실하다.
  (src/v3_make_hwp_docx.py 가 만든 .docx 를 열어 .hwp 로 저장한다.)

보안 모듈 안내
  한글 COM 자동화는 기본적으로 파일 접근 시 보안 팝업을 띄운다.
  RegisterModule("FilePathCheckDLL", "FilePathCheckerModule") 로 등록된 모듈이 있으면
  자동 승인되지만, 없으면 팝업이 뜬다. 없을 때는 팝업을 한 번 눌러 주면 된다.

실행: python src/v3_make_hwp.py
"""
import os
import sys
import time

if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DOCS = os.path.join(ROOT, "docs")
SRC = os.path.join(DOCS, "KASA_연구계획서_2연구방법_3예상결과.docx")
DST = os.path.join(DOCS, "KASA_연구계획서_2연구방법_3예상결과.hwp")


def main():
    if not os.path.exists(SRC):
        print(f"원본이 없습니다: {SRC}")
        print("먼저 python src/v3_make_hwp_docx.py 를 실행하세요.")
        return 1

    import win32com.client as win32

    print("한글 실행 중...")
    hwp = win32.gencache.EnsureDispatch("HWPFrame.HwpObject")
    try:
        # 보안 모듈이 등록돼 있으면 파일 접근 팝업이 자동 승인된다
        try:
            hwp.RegisterModule("FilePathCheckDLL", "FilePathCheckerModule")
            print("  보안 모듈 등록됨 (팝업 없이 진행)")
        except Exception:
            print("  보안 모듈 없음 — 접근 허용 팝업이 뜨면 '허용' 을 눌러 주세요")

        # 창을 띄워 둔다. 숨기면 보안 확인 창이 화면 밖에서 뜰 수 있다.
        try:
            hwp.XHwpWindows.Item(0).Visible = True
        except Exception:
            pass

        # 포맷 문자열을 주면 실패한다("MSWord" 등은 이 버전에서 인식되지 않음).
        # 인자 없이 부르면 확장자로 자동 판별되고 정상 동작한다.
        print(f"열기: {os.path.basename(SRC)}")
        ok = hwp.Open(SRC)
        if not ok:
            print("  .docx 열기 실패")
            return 1
        try:
            print(f"  불러온 쪽수: {hwp.PageCount}")
        except Exception:
            pass

        if os.path.exists(DST):
            os.remove(DST)
        print(f"저장: {os.path.basename(DST)}")
        hwp.SaveAs(DST, "HWP", "")
        time.sleep(1.0)
    finally:
        try:
            hwp.Clear(1)          # 저장 여부 묻지 않고 닫기
            hwp.Quit()
        except Exception:
            pass

    if os.path.exists(DST):
        size = os.path.getsize(DST)
        print(f"\n-> {DST}  ({size:,} bytes)")
        # 실제 HWP 5.0 인지 확인 (OLE 복합문서 서명)
        with open(DST, "rb") as fh:
            sig = fh.read(8)
        ole = sig == b"\xd0\xcf\x11\xe0\xa1\xb1\x1a\xe1"
        print(f"   OLE 복합문서 서명: {'정상' if ole else '확인 필요'}")
        return 0
    print("\n저장 실패")
    return 1


if __name__ == "__main__":
    sys.exit(main())
