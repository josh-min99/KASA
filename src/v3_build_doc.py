"""
그림을 파일 안에 박아 넣은 배포용 HTML 을 만든다.

왜 필요한가
  원본 HTML 은 그림을 상대경로로 참조한다. 그 상태로 구글 독스에 붙여넣으면
  그림이 따라오지 않는다. 여기서 PNG 를 base64 로 인코딩해 문서 안에 넣으면,
  브라우저로 열어 전체 선택·복사했을 때 그림이 함께 클립보드에 담긴다.

사용법
  python src/v3_build_doc.py
  -> docs/한글양식_최종_이미지포함.html 생성
  -> 브라우저로 열고 Ctrl+A, Ctrl+C 후 구글 독스나 한글에 붙여넣기

원본은 docs/한글양식_2_연구방법_3_예상결과.html 이다. 내용 수정은 원본에만 한다.
"""
import base64
import os
import re
import sys

if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DOCS = os.path.join(ROOT, "docs")
SRC = os.path.join(DOCS, "한글양식_2_연구방법_3_예상결과.html")
DST = os.path.join(DOCS, "한글양식_최종_이미지포함.html")

WRAP = """<!doctype html>
<html lang="ko"><head><meta charset="utf-8">
<title>KASA 연구계획서 — 연구방법 / 예상되는 결과 (드라이랩)</title>
<style>body{{margin:28px auto; max-width:840px;}}</style>
</head><body>
{body}
</body></html>
"""


def main():
    html = open(SRC, encoding="utf-8").read()
    n, missing = 0, []

    def repl(m):
        nonlocal n
        rel = m.group(1)
        path = os.path.normpath(os.path.join(DOCS, rel))
        if not os.path.exists(path):
            missing.append(rel)
            return m.group(0)
        with open(path, "rb") as fh:
            b64 = base64.b64encode(fh.read()).decode("ascii")
        n += 1
        print(f"  넣음 {os.path.basename(path):28s} {os.path.getsize(path):>9,} B")
        return 'src="data:image/png;base64,%s"' % b64

    out = re.sub(r'src="([^"]+\.png)"', repl, html)

    if missing:
        print("\n그림 파일을 찾지 못했습니다:")
        for m in missing:
            print("  -", m)
        return 1

    with open(DST, "w", encoding="utf-8") as fh:
        fh.write(WRAP.format(body=out))
    print(f"\n그림 {n}개 포함 -> {DST}  ({os.path.getsize(DST) / 1e6:.1f} MB)")
    print("브라우저로 연 뒤 Ctrl+A, Ctrl+C 하여 구글 독스나 한글에 붙여넣으십시오.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
