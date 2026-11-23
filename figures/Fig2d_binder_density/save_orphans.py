"""Rescue the Figure 2 notes the 2026-11-17 caption rewrite dropped.

They exist only in the manuscript doc's sections 5 and 6, which are exports of
captions/Fig2.md - so re-exporting the current file would delete them from the
last place they live. Appended verbatim to captions/Fig2.md under a heading the
three builders do not read: they split on "## Legend", "## Methods material
(not drawn on the figure)" and "## To confirm", and each split stops at the
next "\n## ", so a section after them exports nowhere.

Not included: the old "Panel b, the territories ... contoured at 5% of that
generator's own maximum" bullet, which the current file explicitly supersedes
and documents, and the two old panel-d-as-GPU-cost bullets' framing (the
numbers inside them are kept).
"""
import json
from pathlib import Path

SRC = Path("/Users/leo/.claude/projects/-Users-leo/7f89041b-ba35-4810-a02f-7a1a907af824/"
           "tool-results/mcp-1455dd7d-667c-405f-87cd-9fc086771d50-read_doc-1795098943579-yiomjn.txt")
CAP = Path("/Volumes/wicky/users/lvalsecchi/Leo_Valsecchi/LV_ThesisFigures/captions/Fig2.md")
WANT = [(50684, "5"), (51738, "5"), (51983, "5"), (49912, "5"),
        (52952, "5"), (53390, "5"),
        (68595, "6"), (68822, "6"), (68987, "6"),
        (69087, "6"), (69389, "6"), (69587, "6")]

C = json.loads(SRC.read_text())["content"]["tabs"][0]["documentTab"]["body"]["content"]
byidx = {}
for e in C:
    if e.get("paragraph") and e.get("startIndex") is not None:
        byidx[e["startIndex"]] = "".join(
            r["textRun"]["content"] for r in e["paragraph"]["elements"] if "textRun" in r)


def wrap(s, width=88, lead="- ", cont="  "):
    words, lines, cur = s.split(), [], lead
    for w in words:
        if len(cur) + len(w) + 1 > width and cur.strip() not in (lead.strip(),):
            lines.append(cur.rstrip())
            cur = cont + w + " "
        else:
            cur += w + " "
    lines.append(cur.rstrip())
    return "\n".join(lines)


def bolded(s):
    """Put the lead-in sentence back in bold, the way the file writes them."""
    for stop in (". ", ", changed"):
        i = s.find(stop)
        if 0 < i < 120:
            return "**" + s[:i + 1].strip() + "** " + s[i + len(stop):].strip()
    return s


out = ["", "---", "",
       "## Orphaned from the 2026-11-16 export (kept, not exported)", "",
       wrap("These were in sections 5 and 6 of the manuscript doc, exported from this file "
            "on 2026-11-16, and the 2026-11-17 rewrite dropped them. They are still true and "
            "the doc was the only place left holding them, so they are parked here before the "
            "re-export of 2026-11-19 overwrote those blocks. Nothing below is exported - the "
            "three builders read only Legend, Methods material and To confirm, and stop at "
            "this heading. Promote an item into one of those sections to put it back in the "
            "manuscript, or delete it once it is genuinely settled.", lead="", cont=""),
       ""]
for sec in ("5", "6"):
    out.append("### From section %s (%s)" % (
        sec, "methods notes" if sec == "5" else "open questions"))
    out.append("")
    for idx, s in WANT:
        if s != sec:
            continue
        txt = " ".join(byidx[idx].split())
        out.append(wrap(bolded(txt) if sec == "5" else txt))
        out.append("")

add = "\n".join(out)
before = CAP.read_text()
assert "## Orphaned from the 2026-11-16 export" not in before, "already appended"
CAP.write_text(before + add)
print("appended %d chars / %d items to %s" % (len(add), len(WANT), CAP))
print("file now %d lines" % len(CAP.read_text().splitlines()))
