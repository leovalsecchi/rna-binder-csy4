"""Fig 1a draft glyphs: the same seven steps drawn as the things they are.

The shipped glyphs in code/common/glyphs.py are one-colour line icons. On the
assembled page panel a is then the only panel with no object in it, while b, d
and f all draw the campaign's molecules in the campaign's colours.

These redraw the seven in the same 48x48 box, at the same stroke weight, under
COLORS.md rule 6: a glyph that depicts a MOLECULE takes its synbio colour
(binder MAIN, RNA theme 400, mScarlet, SMS kinase teal / RNase amber) with a
neutral-700 outline; a glyph that depicts APPARATUS - funnel, tube, sequencer,
cytometer - stays petrol line art and is given life by the molecule inside it.
Rule 1's single accent governs data series, not molecular identity, which is
what rule 6 exists for.

Kept local to the draft on purpose. If the set is adopted it belongs in
code/common/glyphs.py beside the others, so Fig 2a and readout_workflow can
use it too.
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from common import palette                                   # noqa: E402

SYN = palette.SYNBIO
MAIN = palette.MAIN
RNA = SYN["rna"]
BINDER = SYN["binder"]
SCARLET = SYN["reporter_alt"]
KINASE, RNASE = SYN["enzyme"], SYN["enzyme_alt"]
TARGET = palette.THEME["200"]
OUT = SYN["outline"]
N300, N400, N500 = (palette.NEUTRAL[k] for k in ("300", "400", "500"))

SW = 2.4           # the shared glyph stroke, as common/glyphs.py
THIN = 1.4         # outlines on a filled molecule: heavier and it turns muddy
APP = MAIN         # apparatus - funnel, tube, cytometer - stays petrol line
                   # art. Only a MOLECULE takes a synbio fill with a neutral
                   # outline; drawing the glassware in neutral too drained the
                   # colour out of the whole row.
# The construct, as it appears in panels b and d: binder, linker, SMS kinase,
# SMS RNase, target. Widths are the drawn proportions, not base pairs.
CASSETTE = [(BINDER, 7.0), (N300, 1.6), (KINASE, 6.4), (RNASE, 4.0), (TARGET, 5.0)]


CAP = "stroke-linecap='round'"          # kept out of the f-strings below:
JOIN = "stroke-linejoin='round'"        # a backslash cannot appear in one
CAPJOIN = CAP + " " + JOIN


def _o(extra="", col=None):
    return f'fill="none" stroke="{col or APP}" stroke-width="{SW}" {extra}'


def _f(fill, extra="", col=None, sw=None):
    """Outlined, with a fill. _o() already carries fill="none", so a shape that
    needs one has to come through here or the attribute is written twice."""
    return (f'fill="{fill}" stroke="{col or OUT}" '
            f'stroke-width="{sw or SW}" {extra}')


def _bar(x, y, w, h, segs, sw=THIN):
    """A construct bar: coloured segments under one outline."""
    total = sum(s[1] for s in segs)
    o, cx = [], x
    for col, frac in segs:
        sw_ = w * frac / total
        o.append(f'<rect x="{cx:.2f}" y="{y:.2f}" width="{sw_:.2f}" height="{h:.2f}" '
                 f'fill="{col}"/>')
        cx += sw_
    o.append(f'<rect x="{x:.2f}" y="{y:.2f}" width="{w:.2f}" height="{h:.2f}" '
             f'fill="none" stroke="{OUT}" stroke-width="{sw}" rx="1"/>')
    return "".join(o)


def _hairpin(cx, top, h, w, col=RNA, sw=SW):
    """A stem-loop drawn as panels b and d draw it: two rails with pairing
    rungs, under an open circular loop. A smooth arch instead of the circle
    read as a bent paperclip at this size."""
    half, r = w / 2, w * 0.62
    stem_top = top + 2 * r
    o = [f'<circle cx="{cx:.1f}" cy="{top+r:.1f}" r="{r:.1f}" fill="none" '
         f'stroke="{col}" stroke-width="{sw}"/>',
         f'<path d="M{cx-half:.1f},{stem_top:.1f} L{cx-half:.1f},{top+h:.1f} '
         f'M{cx+half:.1f},{stem_top:.1f} L{cx+half:.1f},{top+h:.1f}" fill="none" '
         f'stroke="{col}" stroke-width="{sw}" stroke-linecap="round"/>']
    y = stem_top + 4.5
    while y < top + h - 2.5:
        o.append(f'<line x1="{cx-half+0.7:.1f}" y1="{y:.1f}" x2="{cx+half-0.7:.1f}" '
                 f'y2="{y:.1f}" stroke="{N300}" stroke-width="{THIN*0.85:.2f}"/>')
        y += 5.6
    return "".join(o)


def _blob(cx, cy, r, col):
    """A compact protein: two overlapping lobes under one outline."""
    return (f'<ellipse cx="{cx:.1f}" cy="{cy:.1f}" rx="{r:.1f}" ry="{r*0.76:.1f}" '
            f'transform="rotate(-18 {cx:.1f} {cy:.1f})" fill="{col}" '
            f'stroke="{OUT}" stroke-width="{THIN}"/>')


def g_binder_rna():
    """Step 1: a designed binder docked on the target hairpin."""
    # the blob overlaps the right rail: set apart they read as two objects
    # side by side rather than a binder docked on its target
    return (_hairpin(16, 6, 38, 9)
            + _blob(29.5, 28, 9.5, BINDER))


def g_funnel():
    """Step 2: many designs in, one out. The rejected ones are neutral."""
    o = [f'<path d="M7,12 L41,12 L28.5,27 L28.5,39 L19.5,43 L19.5,27 Z" '
         f'{_f(palette.NEUTRAL["100"], JOIN, APP)}/>']
    for x, col in ((13, N400), (24, BINDER), (35, N400)):
        o.append(f'<circle cx="{x}" cy="6.5" r="3.1" fill="{col}" stroke="{OUT}" '
                 f'stroke-width="{THIN}"/>')
    o.append(f'<circle cx="24" cy="46" r="3.1" fill="{BINDER}" stroke="{OUT}" '
             f'stroke-width="{THIN}"/>')
    return "".join(o)


def g_oligo_pool():
    """Step 3: the ordered pool, single strands of synthetic DNA."""
    return "".join(
        f'<path d="M7,{10+i*7} q5.5,-4.5 11,0 t11,0 t11,0" fill="none" stroke="{col}" '
        f'stroke-width="{SW}" stroke-linecap="round"/>'
        for i, col in enumerate((MAIN, RNA, MAIN, RNA, MAIN)))


def g_integration():
    """Step 4: the cassette landing in the genome inside a cell."""
    return (f'<path d="M24,2 L24,9" {_o(CAP)}/>'
            f'<path d="M20.5,6 L24,10 L27.5,6" {_o(CAPJOIN)}/>'
            f'<ellipse cx="24" cy="30" rx="20" ry="15" '
            f'{_f(palette.NEUTRAL["100"], "", N400, THIN)}/>'
            f'<line x1="5" y1="30" x2="43" y2="30" stroke="{N400}" '
            f'stroke-width="{THIN}"/>'
            + _bar(9, 25.5, 30, 9, CASSETTE))


def g_tube():
    """Step 5: extracted RNA in a tube."""
    return (f'<path d="M15,8 L33,8 L33,31 L24,43 L15,31 Z" '
            f'{_f(palette.NEUTRAL["100"], JOIN, APP)}/>'
            f'<path d="M12.5,8 L35.5,8" {_o(CAP)}/>'
            + "".join(f'<path d="M17.5,{y} q3.25,-3.2 6.5,0 t6.5,0" fill="none" '
                      f'stroke="{RNA}" stroke-width="{SW}" stroke-linecap="round"/>'
                      for y in (22, 30)))


def g_long_reads():
    """Step 6: long reads tiling the construct they came from."""
    o = [_bar(6, 7, 36, 7, CASSETTE)]
    for x, w, y in ((6, 25, 22), (15, 27, 30), (10, 30, 38)):
        o.append(f'<line x1="{x}" y1="{y}" x2="{x+w}" y2="{y}" stroke="{N500}" '
                 f'stroke-width="{SW}" stroke-linecap="round"/>')
    return "".join(o)


def g_cytometry():
    """Step 7: cells in the stream, scored bright or dark."""
    o = [f'<path d="M14,3 C18.5,16 18.5,32 14,45" {_o(CAP)}/>',
         f'<path d="M32,3 C27.5,16 27.5,32 32,45" {_o(CAP)}/>',
         f'<path d="M36,24 L44,24" {_o(CAP)}/>',
         f'<path d="M40,20.5 L44,24 L40,27.5" {_o(CAPJOIN)}/>']
    for cy, col in ((12, SCARLET), (24, N300), (36, SCARLET)):
        o.append(f'<circle cx="23" cy="{cy}" r="4.6" fill="{col}" stroke="{OUT}" '
                 f'stroke-width="{THIN}"/>')
    return "".join(o)


GLYPHS = {
    "binder-rna": g_binder_rna,
    "funnel": g_funnel,
    "oligo-pool": g_oligo_pool,
    "integration": g_integration,
    "tube": g_tube,
    "long-reads": g_long_reads,
    "cytometry": g_cytometry,
}


def group(name, cx, cy, k=1.34):
    """One glyph placed at (cx, cy), scaled by k. Each path carries its own
    colour, so unlike common/glyphs.icon_group the parent does not set one."""
    h = 24 * k
    return (f'<g transform="translate({cx-h:.2f},{cy-h:.2f}) scale({k})" '
            f'stroke-linecap="round" stroke-linejoin="round">{GLYPHS[name]()}</g>')


def sheet(cols=7, box=60):
    """All seven in a row, for looking at them on their own."""
    w = cols * box
    o = [f'<svg xmlns="http://www.w3.org/2000/svg" width="{w*2}" height="{box*2}" '
         f'viewBox="0 0 {w} {box}"><rect width="{w}" height="{box}" fill="#ffffff"/>']
    for i, name in enumerate(GLYPHS):
        o.append(group(name, box * (i + 0.5), box / 2, 1.0))
    o.append("</svg>")
    return "\n".join(o)


if __name__ == "__main__":
    out = Path(__file__).resolve().parents[2] / "working" / "Fig1a_drafts"
    out.mkdir(parents=True, exist_ok=True)
    (out / "rich_glyph_sheet.svg").write_text(sheet(), encoding="utf-8")
    print("wrote", out / "rich_glyph_sheet.svg")
