"""
구글 독스용 HTML 원고를 한글(HWP) 직접 붙여넣기용 텍스트로 변환한다.

왜 자동 변환인가
  같은 내용을 HTML 과 txt 로 따로 손으로 유지하면 반드시 어긋난다.
  HTML 을 원본으로 두고 txt 는 항상 여기서 생성한다.

실행: python src/v3_html_to_text.py
"""
import os
import re
import sys
from html.parser import HTMLParser

if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DOCS = os.path.join(ROOT, "docs")

SRC = os.path.join(DOCS, "한글양식_2_연구방법_3_예상결과.html")
DST = os.path.join(DOCS, "한글양식_2_연구방법_3_예상결과.txt")

HEADER = """\
═══════════════════════════════════════════════════════════════════════
KASA 연구계획서 — 「2 연구의 내용 및 방법」·「3 기대효과 및 활용방안」
드라이랩 담당분 원고 (한글 붙여넣기용)

  ※ 이 파일은 docs/한글양식_2_연구방법_3_예상결과.html 에서 자동 생성됩니다.
     내용을 고칠 때는 HTML 을 고치고 `python src/v3_html_to_text.py` 를 다시 실행하세요.
═══════════════════════════════════════════════════════════════════════

"""


class Conv(HTMLParser):
    def __init__(self):
        super().__init__()
        self.out = []
        self.buf = []
        self.in_table = False
        self.rows = []
        self.row = []
        self.cell = []
        self.in_cell = False

    # ---- 텍스트 수집
    def handle_data(self, d):
        (self.cell if self.in_cell else self.buf).append(d)

    def handle_starttag(self, tag, attrs):
        if tag == "table":
            self.in_table = True
            self.rows = []
        elif tag == "tr" and self.in_table:
            self.row = []
        elif tag in ("td", "th") and self.in_table:
            self.in_cell, self.cell = True, []
        elif tag == "br":
            (self.cell if self.in_cell else self.buf).append(" ")
        elif tag == "hr":
            self.flush()
            self.out.append("\n" + "-" * 71)

    def handle_endtag(self, tag):
        if tag == "table":
            self.render_table()
            self.in_table = False
        elif tag == "tr" and self.in_table:
            if self.row:
                self.rows.append(self.row)
        elif tag in ("td", "th") and self.in_table:
            self.row.append(re.sub(r"\s+", " ", "".join(self.cell)).strip())
            self.in_cell = False
        elif tag in ("p", "div"):
            self.flush()

    # ---- 출력
    def flush(self):
        t = re.sub(r"\s+", " ", "".join(self.buf)).strip()
        self.buf = []
        if not t:
            return
        # 마커로 시작하면 들여쓰기를 준다
        if t.startswith("□"):
            self.out.append("\n" + t)
        elif t.startswith("◦"):
            self.out.append("\n" + t)
        elif t.startswith("-"):
            self.out.append("  " + t)
        elif t.startswith("※"):
            self.out.append("  " + t)
        elif t.startswith("[") and t.endswith("]"):
            self.out.append("  " + t)
        else:
            self.out.append(t)

    def render_table(self):
        rows = [r for r in self.rows if r]
        if not rows:
            return
        n = max(len(r) for r in rows)
        rows = [r + [""] * (n - len(r)) for r in rows]
        w = [max(len(r[i]) for r in rows) for i in range(n)]
        self.out.append("")
        for i, r in enumerate(rows):
            line = "  " + "  ".join(c.ljust(w[j]) for j, c in enumerate(r)).rstrip()
            self.out.append(line)
            if i == 0 or (i == 1 and n > 1):
                self.out.append("  " + "-" * min(sum(w) + 2 * (n - 1), 100))
        self.out.append("")


def main():
    html = open(SRC, encoding="utf-8").read()
    c = Conv()
    c.feed(html)
    c.flush()
    body = "\n".join(c.out)
    body = re.sub(r"\n{3,}", "\n\n", body)
    with open(DST, "w", encoding="utf-8") as fh:
        fh.write(HEADER + body + "\n")
    print(f"-> {DST}  ({len(body):,} 자)")


if __name__ == "__main__":
    main()
