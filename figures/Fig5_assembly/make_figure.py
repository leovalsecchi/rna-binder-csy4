"""Fig 5, assembled on an A4 page.

    python code/Fig5_assembly/make_figure.py [outdir]   # default working/Fig5_assembly

Layout, from Leo's draft page 5 (2026-11-13), with d-f dropped 2026-11-18:
    a | b        a the transfection matrix (RNA target scaffold and control
                 plasmid against the well grid), b the 32 ordered binders
                 ranked by ON/OFF
    c            the four best hits by flow fold change, as surfaces, each
                 with an X-ray window onto the cartoon underneath

The Csy4-chemistry pair that used to sit at d | e (code/Fig5de_csy4_chemistry) is no
longer on this page: the aromatic/Csy4-direction result it showed turned out to be
asparagine in disguise. The folder still builds, it is just not promoted.

**c is code/Fig5c_cutaway_surface** (2026-11-18, Leo: "replace the renderings
in c by the ones in d"): the four best hits by flow fold change - twist 23, 15,
16 and twist 7 - as molecular surfaces, each with a square window cut into it
showing the cartoon through the surface. It began as a second row under
code/Fig5c_best_hits, the cartoon close-ups, and took that row's place: one row
of four, not two rows of the same four. The cameras are still c's own, imported
from Fig5c_best_hits/render.py, so that panel stays the source of the four views
even though it is no longer promoted.
It takes the house W_MM H_MM FONT_PT OUTDIR CLI, so it is placed with
`variant` like d, e and f, drawn straight to W x C_H with its type already at
T. `placeholder` is kept below, unused, for the next panel that needs a slot
held open.

**Which b.** Fig5a_wetlab_on_off emits eighteen candidates. The draft's panel
is the one with "Ordered binder" up the side, "ON/OFF avg Pos / avg Neg"
along the bottom and the 0.48x marker at the foot - that is
Fig5a_PosAvg_vs_NegAvg, matched on the 0.48x. The others are the same plot
against single negatives, correlations between them, or the pass matrix.

Both panels come from code/Fig5de_csy4_chemistry/make_figure.py at 105 mm wide;
`pair` re-scales them to share one height across the 178 mm block and re-sets
their text to the house size at that scale. Neither takes the
W_MM H_MM FONT_PT OUTDIR CLI, so `variant` is not used - `pair` is the right
tool, as it is for Fig 3's c and d.

Panels
  d  output/Fig5de_csy4_chemistry/love_oriented.svg
  e  output/Fig5de_csy4_chemistry/curves.svg
"""
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "code"))
from common import palette                                    # noqa: E402
from common.page_assembly import (ARIAL, GAP_H, GAP_V, MM, PAGE_H, TOP,  # noqa: E402
                                  T_HOUSE, W, Panel, fit, floor, page,
                                  report, variant, write)

ASN = REPO / "code/Fig5def_asn_interface/make_figure.py"
HITS = REPO / "code/Fig5c_best_hits/make_figure.py"
SURF = REPO / "code/Fig5c_cutaway_surface/make_figure.py"
A = REPO / "working/screen_matrix/screen_matrix.svg"
B = REPO / "output/Fig5a_wetlab_on_off/Fig5a_PosAvg_vs_NegAvg.svg"
DE_H = 52 * MM           # d | e, the two list panels side by side
F_H = 46 * MM            # f, the asparagine distributions, full width
C_H = 55 * MM            # the c placeholder. It used to take 95 mm - everything the
                         # page had left - but d/e/f need 103 mm + gaps and the block
                         # is 267 mm, so the slot for the hit structures shrinks to
                         # what is actually left. Re-measure it when they exist.
D_H = 40 * MM            # d, the same four as c but as surfaces. Four across
                         # 178 mm gives ~44 mm each and the renders are ~1.3
                         # wide, so the content needs ~34 mm; the rest is slack
                         # the row gap absorbs.
C_NOTE = "Best hits — structures to come"


def placeholder(workdir, T):
    """Panel c's slot, drawn in points at its final size so it places at scale
    1: a pale container like the others in the figure, and a note at the house
    size. Replace the Panel with the real one when it exists."""
    workdir.mkdir(parents=True, exist_ok=True)
    w, h, sw = W, C_H, 0.75
    n = palette.NEUTRAL
    svg = (f'<svg xmlns="http://www.w3.org/2000/svg" width="{w:.3f}pt" height="{h:.3f}pt" '
           f'viewBox="0 0 {w:.3f} {h:.3f}">'
           f'<rect x="{sw/2}" y="{sw/2}" width="{w-sw:.3f}" height="{h-sw:.3f}" rx="8" '
           f'fill="{n["100"]}" stroke="{n["300"]}" stroke-width="{sw}"/>'
           f'<text x="{w/2:.3f}" y="{h/2 + T*0.35:.3f}" font-family="{ARIAL}" '
           f'font-size="{T}" fill="{n["600"]}" text-anchor="middle">{C_NOTE}</text></svg>')
    f = workdir / "c_placeholder.svg"
    f.write_text(svg, encoding="utf-8")
    return Panel("c", f).at(1.0)


def rows_v1(workdir):
    """Each panel at the size it was drawn."""
    return [[Panel("a", A).at(1.0), Panel("b", B).at(1.0)],
            [Panel("c", REPO / "output/Fig5c_best_hits/fig5c_best_hits.svg").at(1.0)]]


def rows_v2(T, workdir):
    """The row filling the block, both panels REDRAWN by their own script.

    Not `pair`. pair re-sets type in place: it rewrites each font-size to T/scale
    so the label prints at T once the panel is shrunk, but it cannot move the
    labels' x,y, so a right-anchored y label grows leftward out of its own
    viewBox. rsvg will not render outside the viewBox, so ink() cannot see the
    overflow, the panel box comes out too narrow and the page clip eats the start
    of the longest labels ("Aromatic at interface" -> "matic at interface").
    `variant` asks the panel script for the geometry instead, so the type is
    already T at scale 1 and nothing is re-set. Same route as Fig1e and Fig3e.
    """
    half = (W - GAP_H) / 2
    # 2026-11-18 (Leo): "the whole kind of matrix view can be made
    # substantially smaller ... and then, of course, resize figure B as well".
    #
    # This row used `pair`, which solves for ONE shared height and splits the
    # block by the two panels' aspect ratios - so a's width was whatever its
    # own shape asked for (89.9 mm even after the BFP strip came out) and
    # there was no knob to turn. a now takes an explicit width and b the
    # remainder. The two are no longer the same height; `page` places a row's
    # panels at the row's top, so a sits top-aligned with white beneath it,
    # which is the trade this change buys.
    #
    # Neither script takes a size CLI (a is generated SVG, b's takes no args),
    # so both go through `fit`, which rescales the finished panel and re-sets
    # its text to T. Watch b's rotated y label - fit can resize a label, not
    # move it. And see the note under A_W: shrinking a RAISES its text
    # relative to its own glyphs, so the cell captions are what break first.
    A_W = 72 * MM
    a = fit(lambda t, s: Panel("a", A, t, s), T, A_W)
    # 2026-11-18 (Leo): "it's fine to have a gap between the two panels on the
    # same level ... you just had to reduce panel b's size until they were on
    # the same horizontal bottom line." So b is sized by HEIGHT, not width:
    # both panels hang from the row's top (as `page` does by default) and, at
    # one height, they end on one line too. b comes out narrower than the
    # remainder of the block, and the leftover goes into the gap BEFORE b so
    # the row still reaches the right edge like the full-width rows under it.
    b = fit(lambda t, s: Panel("b", B, t, s), T, height=a.h)
    b.dx = W - GAP_H - A_W - b.w
    assert b.dx >= 0, "b is wider than the block's remainder even at a's height"
    # 2026-11-18 (Leo): "probably d-f will be removed. You can remove them for
    # now." The three asparagine panels are off the page; d is now the surface
    # companion to c. To put them back, restore the two rows below - the panels
    # themselves are untouched in code/Fig5def_asn_interface and still build.
    #
    #   mk = lambda L, stem, w, h: variant(L, ASN, stem, w, h, T, workdir,
    #                                      sys.executable)
    #   [mk("d", "metrics.svg", half, DE_H),
    #    mk("e", "aminoacid_volcano.svg", half, DE_H)],
    #   [mk("f", "asndist.svg", W, F_H)]
    #
    # (they would become f, g, h once the surface row keeps d.)
    # One row of four, not two. The cartoon close-ups (HITS) used to be c with
    # the surfaces as d underneath; the surfaces carry the labels now and the
    # cartoon row is off the page. To put it back, add
    #   [variant("c", HITS, "fig5c_best_hits.svg", W, C_H, T, workdir,
    #            sys.executable)],
    # above this row and let the surfaces be d again (LABELS = False in their
    # make_figure.py, which is what that flag is for).
    rows = [[a, b],
            [variant("c", SURF, "fig5c_cutaway_surface.svg", W, C_H, T,
                     workdir, sys.executable)]]
    for p in sum(rows, []):          # last, once every scale is final
        p.retouch(lambda u, s=p.s: floor(u, s))
    return rows


def main(outdir):
    outdir.mkdir(parents=True, exist_ok=True)
    r1 = rows_v1(outdir / "_variant")
    svg1, fit1, bot1 = page(r1, "Fig5 V1 as sketched")
    r2 = rows_v2(T_HOUSE, outdir / "_variant")
    svg2, fit2, bot2 = page(r2, "Fig5 V2 adjusted")
    report("V1 as sketched", r1, fit1, bot1)
    report(f"V2 adjusted, T = {T_HOUSE:.2f} pt", r2, fit2, bot2, T_HOUSE)
    assert fit2 == 1.0, "V2 must fit without shrinking, or its type is no longer T"
    room = (PAGE_H - TOP - bot2) / MM
    print(f"room left below the figure: {room:.0f} mm")
    write(outdir, "Fig5_V1_as_sketched", svg1)
    write(outdir, "Fig5_V2_adjusted", svg2)


if __name__ == "__main__":
    main(Path(sys.argv[1]) if len(sys.argv) > 1 else REPO / "working" / "Fig5_assembly")
