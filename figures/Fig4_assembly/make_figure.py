"""Fig 4, assembled on an A4 page.

    python code/Fig4_assembly/make_figure.py [outdir]   # default working/Fig4_assembly

Layout, from Leo's draft page 4 (2026-11-13):
    a        the screen overview: library cassette and oligo pool, the four
             cell populations, the readout workflow, and the hit gate plot in
             the corner. Already a composite - code/Fig4a_screen_overview
             builds it from three emitted SVGs plus its own gate plot.
    b        the 32 ordered hits, 4 x 8 structure renders
    c | d    the hit analyses, filling the slot the draft marks
             "Analysis of Hits". d was the hairpin footprint until
             2026-11-13; it is now the four-layer epitope funnel from
             Fig2d_epitope_funnel, on Leo's call.

**c is a raster.** Fig4ef_ordered32 emits only a PNG (32 ChimeraX renders), so
it is downsampled to GRID_PX across - 300 dpi at 178 mm - and embedded as
base64 in a wrapper SVG. The full-size original is several MB and would
treble the size of the figure book for resolution no print can use. Since
2026-11-14 the sheet is 7 x 5, not 8 x 4: taller, bigger cells, and it takes
up most of the white the page had under it.

**d is shared with Fig 2 panel e** - the same campaign/gate/library/hits
column panel appears on both pages. Deliberate, but if one has to go it is
the Fig 2 copy, since the hits bar is what earns it a place in Fig 4.

**The hairpin footprint is no longer placed.** Fig4g_hit_epitope/hairpin_footprint
was panel d until 2026-11-13 and is still built; it says the same thing as the
columns in less space but without the campaign and gate layers.

**Not placed, and why.** Fig4g_hit_epitope/epitope_columns, and all three of
Fig4h_hit_geometry (shape_space, density_curves, percentile_strip), are built
and good but do not fit: a and b alone take 174 mm of the 267 mm block, which
leaves one row. Those four are all wide (aspect 2.8 to 3.4) and want a row
each. They need either a second Fig 4 page or a decision about which to drop.

Leo's note on the draft - "background colour + rank by core? remove name" -
is about the grid itself and belongs in Fig4ef_ordered32, not here.
"""
import base64
import re
import subprocess
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "code"))
from common.page_assembly import (GAP_H, MM, PAGE_H, TOP, T_HOUSE, W,   # noqa: E402
                                  Panel, fit, floor, page, pair, report,
                                  variant, write)

OUT = REPO / "output"
A = OUT / "Fig4a_screen_overview/Fig4a_screen_overview.svg"
B = OUT / "Fig4a_screen_overview/Fig4b_screen_gates.svg"
C_PNG = OUT / "Fig4ef_ordered32/Fig4f_grid_5x7_y120_boxed.png"
D = OUT / "Fig4d_stack_hit_rate/Fig4d_stack_hit_rate.svg"
D_SCRIPT = REPO / "code/Fig4d_stack_hit_rate/make_figure.py"
E = OUT / "Fig2d_epitope_funnel/epitope_funnel_columns.svg"
GRID_PX = 2100           # 178 mm at 300 dpi
A_W = 70 * MM            # panel a's width in the top row; b takes the rest

# Panel a is itself a composite: Fig4a_screen_overview places the cassette, the
# populations and the readout workflow at three different internal scales, each
# chosen so its own labels land on 5.5 pt in panel units. So the file holds
# three font-size values that all PRINT the same, and the panel's own unit is
# a point - at scale 1 every label is already 5.5 pt.
#
# That breaks the usual retext(), which rewrites every font-size to one number:
# the workflow's 13.65 sits inside a scale(0.40293), so making it equal to the
# cassette's 5.5 printed it at 40% - the bug this replaced. The panel is
# rescaled PROPORTIONALLY instead (every size times the same factor), which
# preserves each label's printed size whatever transform it sits under.
NEST = {13.65: 13.65 / 5.5, 18.0664: 18.0664 / 5.5}   # not riders: nested scales


def grid_wrapper(workdir):
    """The 4 x 8 grid PNG, downsampled and embedded in an SVG so it can be
    placed like any other panel."""
    workdir.mkdir(parents=True, exist_ok=True)
    small = workdir / "grid_300dpi.png"
    subprocess.run(["magick", C_PNG, "-resize", f"{GRID_PX}x", small], check=True)
    w, h = subprocess.run(["magick", "identify", "-format", "%w %h", small],
                          check=True, capture_output=True, text=True).stdout.split()
    b64 = base64.b64encode(small.read_bytes()).decode()
    svg = (f'<svg xmlns="http://www.w3.org/2000/svg" '
           f'xmlns:xlink="http://www.w3.org/1999/xlink" width="{w}" height="{h}" '
           f'viewBox="0 0 {w} {h}"><image width="{w}" height="{h}" '
           f'xlink:href="data:image/png;base64,{b64}"/></svg>')
    p = workdir / "grid.svg"
    p.write_text(svg, encoding="utf-8")
    print(f"grid: {C_PNG.name} -> {w}x{h}, {len(svg)/1e6:.1f} MB embedded")
    return p


def scaled_text(t, f):
    """Every font-size and letter-spacing times f. Unlike retext this keeps the
    sizes' RATIOS, which is what a panel with nested scales needs."""
    t = re.sub(r'(font-size(?::\s*|="))([0-9.]+)',
               lambda m: f"{m.group(1)}{float(m.group(2)) * f:.4f}", t)
    return re.sub(r'(letter-spacing=")([0-9.]+)',
                  lambda m: f"{m.group(1)}{float(m.group(2)) * f:.4f}", t)


def overview(s):
    """Panel a with its text pre-shrunk by 1/s, so that placing it at s leaves
    every label printing at exactly the size it was drawn to - 5.5 pt."""
    return Panel("a", A, riders=NEST,
                 prep=(lambda u: scaled_text(u, 1 / s)) if s else None)


def rows_v1(grid):
    """Each panel at the size it was drawn."""
    return [[Panel("a", A, riders=NEST).at(1.0), Panel("b", B).at(1.0)],
            [Panel("c", grid).at(1.0)],
            [Panel("d", D).at(1.0), Panel("e", E).at(1.0)]]


def rows_v2(T, grid, workdir):
    # a and b share the top row (2026-11-14): with the gate plot gone, a is
    # about square and full width would make it 175 mm tall, which the page
    # has not got. Paired with the two-frame gate panel it lands at about a
    # third of the block. a keeps the proportional-text trick either way -
    # pair() asks for its panel at the scale it has solved for, so overview()
    # pre-shrinks the text by that same scale and the labels print at T.
    # a and b are both about square once b's two gate frames are stacked, so
    # pair() splits the row near evenly and the two come out ~86 mm tall -
    # 18 mm more page than the side-by-side version, and no white under
    # either. a keeps the proportional-text trick: pair() asks for its panel
    # at the scale it has solved for, so overview() pre-shrinks the text by
    # that same scale and the labels print at T.
    a, b = pair(lambda t, s: overview(s), lambda t, s: Panel("b", B, t, s), T)
    c = Panel("c", grid)                     # a raster: no text to re-set
    c.at(W / c.iw)
    d, e = pair(lambda t, s: Panel("d", D, t, s),
                lambda t, s: Panel("e", E, t, s), T)
    # e is a 4:1 strip, so the pair leaves d about 44 mm wide. Scaled into
    # that slot d shrinks ~0.6x and its text, re-set to T, grows against the
    # bars until "RFD3+LiMPNN" runs off the left edge of the page. Redraw it
    # at the slot size instead, with the type already at T.
    d = variant("d", D_SCRIPT, D.name, d.w, e.h, T, workdir, sys.executable)
    rows = [[a, b], [c], [d, e]]
    for p in sum(rows, []):                  # last, once every scale is final
        p.retouch(lambda u, s=p.s: floor(u, s))
    return rows


def main(outdir):
    outdir.mkdir(parents=True, exist_ok=True)
    grid = grid_wrapper(outdir / "_grid")
    r1 = rows_v1(grid)
    svg1, fit1, bot1 = page(r1, "Fig4 V1 as sketched")
    r2 = rows_v2(T_HOUSE, grid, outdir / "_variant")
    svg2, fit2, bot2 = page(r2, "Fig4 V2 adjusted")
    report("V1 as sketched", r1, fit1, bot1)
    report(f"V2 adjusted, T = {T_HOUSE:.2f} pt", r2, fit2, bot2, T_HOUSE)
    assert fit2 == 1.0, "V2 must fit without shrinking, or its type is no longer T"
    # a's labels are guarded by construction, not by report(): its printed size
    # is (file size / s) * s = the size it was drawn at, whatever nesting.
    a = r2[0][0]
    assert abs(min(a.sizes()) - T_HOUSE) < 0.02, f"a's base label is {min(a.sizes())}"
    print(f"room left below the figure: {(PAGE_H - TOP - bot2) / MM:.1f} mm")
    write(outdir, "Fig4_V1_as_sketched", svg1)
    write(outdir, "Fig4_V2_adjusted", svg2)


if __name__ == "__main__":
    main(Path(sys.argv[1]) if len(sys.argv) > 1 else REPO / "working" / "Fig4_assembly")
