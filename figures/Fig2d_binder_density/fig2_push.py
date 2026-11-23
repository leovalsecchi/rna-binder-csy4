"""Re-export the Figure 2 blocks of captions/Fig2.md into the manuscript doc.

The doc's section 4/5/6 Figure 2 subsections are the 2026-11-16 export and
describe the pre-17-Nov figure (old panel d = GPU cost, old panel b territory
rule, no density panel). Each of the three is replaced wholesale by the current
repo block: a legend rewritten on 11-17 does not diff usefully against its
predecessor paragraph by paragraph.

Three regions in ONE batch, highest index first, because an edit at a higher
index never moves a lower one, so each region's style requests - computed as
insert_index + offset into its own inserted string - are still exact when they
run. The HEADING_2 and the blank paragraph that separates the subsection from
the next heading are left alone; only the body paragraphs between them move.
"""
import json
import re
from pathlib import Path

CAP = Path("/Volumes/wicky/users/lvalsecchi/Leo_Valsecchi/LV_ThesisFigures/captions/Fig2.md")
DOC = "1s5SR9KCyp_DUrpfLzKloAOv8-ml3oWc6J3mWlVE2ttU"
TAB = "t.0"
MD = CAP.read_text()

# (name, delete_from, delete_to, markdown heading, bulleted)
REGIONS = [
    ("6. Open questions",  68595, 69743, "## To confirm", True),
    ("5. Methods notes",   48118, 53683, "## Methods material (not drawn on the figure)", True),
    ("4. Figure descriptions", 30889, 34657, "## Legend", False),
]


def block(heading):
    b = MD.split(heading, 1)[1]
    return re.split(r"\n---\n|\n## ", b, 1)[0].strip()


def build(text_blocks, bulleted):
    """-> (text, bold, italic, bullet_runs) with offsets into text."""
    out, bold, ital, bullets = [], [], [], []

    def emit(s):
        i = 0
        for m in re.finditer(r"\*\*(.+?)\*\*|\*(.+?)\*|`(.+?)`", s, re.S):
            out.append(s[i:m.start()])
            pos = sum(len(x) for x in out)
            inner = (m.group(1) or m.group(2) or m.group(3)).replace("\n", " ")
            out.append(inner)
            (bold if m.group(1) else ital if m.group(2) else []).append(
                (pos, pos + len(inner)))
            i = m.end()
        out.append(s[i:])

    for para in text_blocks:
        start = sum(len(x) for x in out)
        emit(para)
        if bulleted:
            bullets.append([start, sum(len(x) for x in out) + 1])
        out.append("\n")

    runs = []
    for s, e in bullets:
        if runs and runs[-1][1] == s:
            runs[-1][1] = e
        else:
            runs.append([s, e])
    return "".join(out), bold, ital, runs


def main():
    reqs, report = [], []
    for name, lo, hi, heading, bulleted in REGIONS:
        raw = block(heading)
        if bulleted:
            items = [i.strip() for i in re.split(r"\n(?=- )", raw) if i.strip()]
            assert all(i.startswith("- ") for i in items), name
            paras = [" ".join(i[2:].split()) for i in items]
        else:
            paras = [" ".join(p.split()) for p in raw.split("\n\n")]

        text, bold, ital, runs = build(paras, bulleted)
        assert all(ord(c) < 0x10000 for c in text), "astral char breaks offsets"
        assert not text.startswith("\n") and text.endswith("\n")

        reqs.append({"deleteContentRange": {
            "range": {"startIndex": lo, "endIndex": hi, "tabId": TAB}}})
        reqs.append({"insertText": {
            "text": text, "location": {"index": lo, "tabId": TAB}}})
        reqs.append({"updateParagraphStyle": {
            "range": {"startIndex": lo, "endIndex": lo + len(text), "tabId": TAB},
            "paragraphStyle": {"namedStyleType": "NORMAL_TEXT"},
            "fields": "namedStyleType"}})
        for s, e in bold:
            reqs.append({"updateTextStyle": {
                "range": {"startIndex": lo + s, "endIndex": lo + e, "tabId": TAB},
                "textStyle": {"bold": True}, "fields": "bold"}})
        for s, e in ital:
            reqs.append({"updateTextStyle": {
                "range": {"startIndex": lo + s, "endIndex": lo + e, "tabId": TAB},
                "textStyle": {"italic": True}, "fields": "italic"}})
        for s, e in sorted(runs, reverse=True):
            reqs.append({"createParagraphBullets": {
                "range": {"startIndex": lo + s, "endIndex": lo + e, "tabId": TAB},
                "bulletPreset": "BULLET_DISC_CIRCLE_SQUARE"}})
        report.append(
            "%-24s %5d->%5d  was %5d chars, now %5d  %2d paras  "
            "%2d bold %d ital %d bullet runs"
            % (name, lo, hi, hi - lo, len(text), len(paras), len(bold),
               len(ital), len(runs)))

    p = Path(__file__).with_name("fig2_requests.json")
    p.write_text(json.dumps(reqs, ensure_ascii=False))
    print("\n".join(report))
    print("\nrequests %d   payload %.1f KB" % (len(reqs), p.stat().st_size / 1024))
    print("doc", DOC, "\nwrote", p)


main()
