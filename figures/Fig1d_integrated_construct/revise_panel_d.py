#!/usr/bin/env python3
"""Fig 1 panel d: Leo's prio-2 revisions, applied to the compact panel SVG.

`working/Fig1d_integrated_construct/Full Panel V4 compact.svg` is what
`Fig1_assembly` places as panel d. `compact_full_panel.py` narrows Leo's
hand-drawn `Full Panel V2.svg` into it; this script corrects that result.

    python code/Fig1d_integrated_construct/revise_panel_d.py

in   working/_prio2_fork/Full Panel V4 compact.source.svg   (pristine; frozen on
     the first run from the live file, then never written again)
out  working/Fig1d_integrated_construct/Full Panel V4 compact.svg

Reading a pristine copy and writing the live path makes this IDEMPOTENT and
leaves the assembly untouched. **If `compact_full_panel.py` is ever re-run it
overwrites the live file, so delete the pristine copy and re-run this after.**

## Coordinates

The file's group transforms are translate-only, so page coordinates are just the
raw numbers plus the accumulated translate. Every target below sits in a group
at translate(-261.5, 0), so page x = raw x - 261.5. The numbers in this script
are RAW (what is written in the d attribute); the comments give page values.
viewBox is 0 0 2520.043 1280.03.

## The four changes (Leo, 2026-11-17)

1. **The straight arrow over the reporter hairpin is removed**, in the mechanism
   inset - the same edit as panel b, and the same drawing. The two curved arrows
   either side stay; only the middle, straight one goes.
2. **The two "Target module" brackets are aligned to their blocks.** Measured off
   a 2520 px render rather than guessed: the bracket spans page x 2274.75-2440.49
   while its light-blue block spans 2254.8-2436.8 (top row) and 2252.8-2438.8
   (lower row). The LEFT edge was inset by 20-22; the right edge was out by
   2-4, which is within the 2.961 stroke and is left alone. Both brackets move
   their left edge 21 left.
3. **A gap where the two leads cross.** The Target-module lead runs horizontally
   at page y 279.41 from x 1512.47 to 2357.62; the SMS lead drops vertically at
   page x 2046.28 from y 196.96 to 354.18. They cross at (2046.28, 279.41) with
   no break. The horizontal is split either side of the crossing - it is the
   longer line, so it reads as the one passing behind. The gap is +/-9, which
   leaves 6.0 units of white either side of the vertical once its 2.961 stroke
   and the round caps are taken off.
4. **An ONT glyph under the splicing-ratio readout**, with an arrow down to it.
   `long-reads` is used, not `nanopore`. `nanopore` is the platform glyph and is
   what `readout_workflow` pairs with "ONT sequencing", but it was tried first
   and fails here: its threading strand is a vertical line through the middle of
   the glyph, which runs straight on from the arrow pointing down at it, so the
   two read as one stroke. `long-reads` is horizontal lines only, so there is
   nothing for the arrow to merge with, and it matches the label's own wording.

   **The panel does not grow** (Leo, 2026-11-17: "keep the panel size and just add
   a small glyph and small text under it"). It goes in the empty pocket inside
   the pale panel below the readout box, which runs from the box foot at page
   y 1112.6 to about y 1278 - 165 units. The label is set on TWO lines because
   "Readout with ONT long read" measures 343.8 units at the body size against an
   inner box 318.7 wide; the two lines are 156.5 and 176.5.
"""
import re
import shutil
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from common import paths, glyphs                                 # noqa: E402

LIVE = paths.ROOT / "working" / "Fig1d_integrated_construct" / "Full Panel V4 compact.svg"
SRC_DIR = paths.ROOT / "working" / "_prio2_fork"
SRC = SRC_DIR / "Full Panel V4 compact.source.svg"

DX = -261.5                      # the translate every target below sits under
page = lambda raw: raw + DX

# 1. the straight arrow in the mechanism inset
ARROW_SHAFT = "M1722.419 408.969 L1722.419 423.706"
ARROW_HEAD = "M1722.419 428.474 L1724.412 423.706 L1720.426 423.706 Z"

# 2. the two Target-module brackets; left edge moves BRACKET_DX
BRACKETS = (
    "M2536.253 131.671 L2536.253 118.669 L2701.995 118.669 L2701.995 131.671",
    "M2536.253 531.715 L2536.253 518.716 L2701.995 518.716 L2701.995 531.715",
)
BRACKET_LEFT_OLD, BRACKET_DX = "2536.253", -21.0

# 3. the horizontal lead and the vertical it crosses (raw x of the crossing)
LEAD_H = "M1773.975 279.412 L2619.123 279.412"
# GAP is the half-width. The caps are round, so the visible white either side
# is GAP - 1.48 (half the crossing stroke) - 1.48 (the cap): 6.0 at GAP 9.0.
# At GAP 5.0 that left 2.0, and the lines read as touching rather than jumping.
CROSS_X, GAP = 2307.784, 9.0

# 4. the ONT readout, in page coordinates (inserted at the top level)
BOX_CX, BOX_FOOT = 2297.16, 1112.61
# The arrow is kept SHORT and the glyph large: at k=0.92 the glyph was 44 units
# (~8.8 pt on the page) and its pore walls closed up. k=1.25 gives 60 units, and
# the arrow loses what the glyph gains so the pocket still holds both.
ARROW_Y0, ARROW_Y1, HEAD_H, HEAD_W = 1118.0, 1138.0, 7.0, 3.1
GLYPH_CY, GLYPH_K = 1180.0, 1.25
LABEL = ("Readout with", "ONT long read")
LABEL_Y0, LABEL_DY = 1236.0, 29.5
LEAD_GREY, INK = "#919699", "#1d2428"
STROKE = 2.961


def bootstrap() -> str:
    if not SRC.exists():
        SRC_DIR.mkdir(parents=True, exist_ok=True)
        shutil.copy2(LIVE, SRC)
        print(f"froze pristine source -> {SRC}")
    return SRC.read_text()


def drop_path(svg: str, d: str, what: str) -> str:
    pat = re.compile(r'\s*<path[^>]*?d="' + re.escape(d) + r'"[^>]*/?>')
    svg, n = pat.subn("", svg)
    assert n == 1, f"{what}: expected 1 path, removed {n}"
    return svg


def sub_d(svg: str, old_d: str, new_d: str, what: str) -> str:
    assert old_d in svg, f"{what}: source d not found"
    return svg.replace(f'd="{old_d}"', f'd="{new_d}"', 1)


def main() -> None:
    svg = bootstrap()

    # 1 ---------------------------------------------------------------------
    svg = drop_path(svg, ARROW_SHAFT, "straight arrow shaft")
    svg = drop_path(svg, ARROW_HEAD, "straight arrow head")

    # 2 ---------------------------------------------------------------------
    new_left = f"{float(BRACKET_LEFT_OLD) + BRACKET_DX:.3f}"
    for d in BRACKETS:
        assert d.count(BRACKET_LEFT_OLD) == 2, "bracket: expected 2 left-edge points"
        svg = sub_d(svg, d, d.replace(BRACKET_LEFT_OLD, new_left), "Target module bracket")

    # 3 ---------------------------------------------------------------------
    assert f'd="{LEAD_H}"' in svg, "horizontal lead not found"
    m = re.search(r'<path([^>]*?)d="' + re.escape(LEAD_H) + r'"([^>]*?)/?>', svg)
    assert m, "horizontal lead element not found"
    attrs = (m.group(1) + m.group(2)).strip()
    split = (f'<path {attrs} d="M1773.975 279.412 L{CROSS_X - GAP:.3f} 279.412"/>'
             f'<path {attrs} d="M{CROSS_X + GAP:.3f} 279.412 L2619.123 279.412"/>')
    svg = svg[:m.start()] + split + svg[m.end():]

    # 4 ---------------------------------------------------------------------
    arrow = (f'<path d="M{BOX_CX} {ARROW_Y0} L{BOX_CX} {ARROW_Y1}" fill="none" '
             f'stroke="{LEAD_GREY}" stroke-width="{STROKE}" stroke-linecap="round"/>'
             f'<path d="M{BOX_CX} {ARROW_Y1 + HEAD_H} L{BOX_CX + HEAD_W} {ARROW_Y1} '
             f'L{BOX_CX - HEAD_W} {ARROW_Y1} Z" fill="{LEAD_GREY}" stroke="{LEAD_GREY}" '
             f'stroke-width="{STROKE}" stroke-linejoin="round"/>')
    glyph = glyphs.icon_group("long-reads", BOX_CX, GLYPH_CY, INK,
                              k=GLYPH_K, sw=STROKE / GLYPH_K)
    text = "".join(
        f'<text transform="matrix(1 0 0 1 {BOX_CX} {LABEL_Y0 + i * LABEL_DY})" x="0" y="0" '
        f'xml:space="preserve" class="cls-10"><tspan x="0" y="0">{line}</tspan></text>'
        for i, line in enumerate(LABEL))
    block = f'<g id="ONT_readout" data-name="ONT readout">{arrow}{glyph}{text}</g>'

    assert svg.rstrip().endswith("</svg>")
    svg = svg.rstrip()[: -len("</svg>")] + block + "\n</svg>\n"

    LIVE.write_text(svg)
    print(f"wrote {LIVE}")
    print(f"  arrow removed at page x {page(1722.419):.1f}")
    print(f"  brackets left edge {page(float(BRACKET_LEFT_OLD)):.1f} -> {page(float(new_left)):.1f}")
    print(f"  lead gap at page ({page(CROSS_X):.1f}, 279.4), +/-{GAP}")
    print(f"  ONT glyph at page ({BOX_CX}, {GLYPH_CY}), k={GLYPH_K}, box foot {BOX_FOOT}")


if __name__ == "__main__":
    main()
