"""Fig 1a drafts: the dry/wet split made structural instead of just labelled.

    python code/Fig1a_screen_pipeline/make_drafts.py [outdir]
    default outdir: working/Fig1a_drafts/

The shipped panel (make_figure.py) is one row of seven steps with the two
phases marked only by a name set into a hairline. On the assembled page it is
the one panel that reads as a different document: a 27 mm strip of single
colour line glyphs over five panels of coloured objects.

These two drafts change ONLY the structure, not the glyphs, so the layout can
be judged on its own. Swapping the abstract glyphs for the campaign's real
cartoons is a separate, additive step.

    rows     two rows, dry lab above, wet lab below, the DNA pool between them
             as the handover. Phases still named by a hairline with the name in
             a gap, one rule per row, so nothing but the arrangement changes.

    grounds  the same two rows, with the dry half on the pale neutral ground
             (NEUTRAL 100) the other panels already use for their containers
             and the wet half on white. Phase names sit top-left in each band,
             so the split is visible before a word is read and no rule is
             needed.

NOTE on 'grounds' vs the 2026-11-10 decision. Fig1a's four pale rounded boxes
were removed then as "the slide-deck tell". That judgement was made on the
panel standing alone. Panels b, d and f now all carry pale NEUTRAL-100 rounded
containers, so one here is consistent with the page rather than against it.
It is still a reversal, and it is a draft for that reason.

Both drafts keep the house rules: no boxes around individual steps, chevrons
not arrows, one accent colour, sizes in points, no panel letter.
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import yaml                                                   # noqa: E402
from common import paths, palette, glyphs                     # noqa: E402

REPO = paths.ROOT
_cfg = yaml.safe_load((REPO / "config" / "style.yaml").read_text())
FONT = f'{_cfg["font_family"]}, Helvetica, sans-serif'
DPI = int(_cfg["dpi"])
WIDTH_MM = float(_cfg["width_double_mm"])
LABEL_PT = float(_cfg["font_size_pt"])
PHASE_PT = LABEL_PT
HAIR_PT = float(_cfg["line_width_pt"])

ACC = palette.MAIN
INK = palette.INK
MUT = palette.NEUTRAL["600"]
RULE = palette.NEUTRAL["300"]
CHEV = palette.NEUTRAL["400"]
GROUND = palette.NEUTRAL["100"]      # the pale container the other panels use
GROUND_EDGE = palette.NEUTRAL["300"]
BG = palette.SURFACE

# Two rows and the handover between them. One label line each: with four steps
# per row instead of seven the slots are wide enough, and one line keeps the
# panel short enough to fit the page.
DRY = [("binder-rna", "RNA-protein binder generation"),
       ("funnel", "Selection of favourable designs")]
HANDOVER = ("oligo-pool", "Synthetic DNA pool ordering")
WET = [("integration", "Bxb1 integration"),
       ("tube", "RNA extraction"),
       ("long-reads", "Long-read NGS and analysis"),
       ("cytometry", "Flow cytometry validation")]
DRY_NAME, WET_NAME = "IN SILICO", "EXPERIMENTAL"

W = 1000
U = W * 0.352778 / WIDTH_MM          # viewBox units per point
X0 = 30
SLOT = (W - 2 * X0) / len(WET)       # the wet row sets the grid
ICON_K = 1.34                        # as the shipped panel: ~11.5 mm printed
CX = W / 2

# One row: name tier, glyph, one label line. Offsets from the row's top.
NAME_DY = 12
ICON_DY = 62
LAB_DY = 108
ROW_H = 116
# The handover tier, from the dry row's bottom.
HAND_CHEV_DY = 16
HAND_ICON_DY = 56
HAND_LAB_DY = 102
HAND_H = 118
BAND_OVER = 8          # how far the pale ground runs past the row, bottom and
                       # top. Any less and the dry labels sit on its edge; any
                       # more and it swallows the handover chevron below.
TOP = 8
CHEV_S = 7.0


DRY_SPREAD = 0.34       # the dry pair's separation as a fraction of the span.
                        # Spread over the full width like the wet row they sit
                        # a third of the panel apart with only a chevron
                        # between, which reads as a gap rather than a step.


def centres(n, spread=None):
    """n steps across the grid: evenly when spread is None (the wet row), or
    as a pair straddling the centre line at that fraction of the span."""
    span = W - 2 * X0
    if spread is None:
        return [X0 + span * (i + 0.5) / n for i in range(n)]
    return [CX + span * spread * (i - (n - 1) / 2) / max(n - 1, 1) for i in range(n)]


def T(x, y, s, pt=LABEL_PT, col=INK, w=400, anch="middle", ls=0):
    return (f'<text x="{x:.1f}" y="{y:.1f}" font-family="{FONT}" font-size="{pt*U:.2f}" '
            f'fill="{col}" font-weight="{w}" text-anchor="{anch}" '
            f'letter-spacing="{ls}">{s}</text>')


def hair(x1, x2, y):
    return (f'<line x1="{x1:.1f}" y1="{y:.1f}" x2="{x2:.1f}" y2="{y:.1f}" '
            f'stroke="{RULE}" stroke-width="{HAIR_PT*U:.2f}"/>')


def chevron(x, y, s=CHEV_S, down=False):
    """A chevron pointing right, or down for the step between the two rows."""
    a, b = (s * 0.55, s) if not down else (s, s * 0.55)
    d = (f'M{x-a:.1f},{y-b:.1f} L{x:.1f},{y+b:.1f} L{x+a:.1f},{y-b:.1f}'
         if down else
         f'M{x-a:.1f},{y-b:.1f} L{x+a:.1f},{y:.1f} L{x-a:.1f},{y+b:.1f}')
    return (f'<path d="{d}" fill="none" stroke="{CHEV}" stroke-width="{1.1*U:.2f}" '
            f'stroke-linecap="round" stroke-linejoin="round"/>')


def phase_rule(x1, x2, y, name):
    """Hairline across the row with the phase name set into a gap in it."""
    sz, ls = PHASE_PT * U, 1.6
    cx = (x1 + x2) / 2
    half = (len(name) * sz * 0.60 + len(name) * ls) / 2 + 9
    return (hair(x1, cx - half, y) + hair(cx + half, x2, y)
            + T(cx, y + sz * 0.36, name, PHASE_PT, MUT, 600, "middle", ls))


def row(steps, top, cs):
    """One row of steps with chevrons between them."""
    o = []
    for (glyph, label), cx in zip(steps, cs):
        o.append(glyphs.icon_group(glyph, cx, top + ICON_DY, ACC, ICON_K))
        o.append(T(cx, top + LAB_DY, label))
    for a, b in zip(cs, cs[1:]):
        o.append(chevron((a + b) / 2, top + ICON_DY))
    return o


def handover(top):
    glyph, label = HANDOVER
    return [chevron(CX, top + HAND_CHEV_DY, down=True),
            glyphs.icon_group(glyph, CX, top + HAND_ICON_DY, ACC, ICON_K),
            T(CX, top + HAND_LAB_DY, label)]


def build(variant):
    dry_top = TOP
    hand_top = dry_top + ROW_H
    wet_top = hand_top + HAND_H
    H = wet_top + ROW_H + 8
    dry_cs, wet_cs = centres(len(DRY), DRY_SPREAD), centres(len(WET))

    o = [f'<svg xmlns="http://www.w3.org/2000/svg" width="{WIDTH_MM}mm" '
         f'height="{WIDTH_MM*H/W:.2f}mm" viewBox="0 0 {W} {H:.0f}">',
         f'<rect width="{W}" height="{H:.0f}" fill="{BG}"/>']

    if variant == "grounds":
        # Full width, so both phase names start on the same left edge; the
        # dry band is the only thing filled, which is the whole point.
        o.append(f'<rect x="{X0-14}" y="{dry_top-6}" width="{W-2*(X0-14)}" '
                 f'height="{ROW_H+BAND_OVER}" rx="10" fill="{GROUND}" '
                 f'stroke="{GROUND_EDGE}" stroke-width="{HAIR_PT*U:.2f}"/>')
        for name, top in ((DRY_NAME, dry_top), (WET_NAME, wet_top)):
            o.append(T(X0, top + NAME_DY + PHASE_PT * U * 0.36, name,
                       PHASE_PT, MUT, 600, "start", 1.6))
    else:
        o.append(phase_rule(dry_cs[0] - SLOT / 2 + 10, dry_cs[-1] + SLOT / 2 - 10,
                            dry_top + NAME_DY, DRY_NAME))
        o.append(phase_rule(wet_cs[0] - SLOT / 2 + 10, wet_cs[-1] + SLOT / 2 - 10,
                            wet_top + NAME_DY, WET_NAME))

    o += row(DRY, dry_top, dry_cs)
    o += handover(hand_top)
    o += row(WET, wet_top, wet_cs)
    o.append("</svg>")
    return "\n".join(o), H


def _check(variant, svg, H):
    assert f'font-size="{LABEL_PT*U:.2f}"' in svg, "labels are not the house size"
    widest = max(len(s) for _, s in DRY + WET + [HANDOVER])
    assert widest * LABEL_PT * U * 0.52 < SLOT, "a label is wider than its slot"
    assert H / W > 0.2, "still a letterbox strip"
    # glyphs may carry rects of their own, so count the ground by its fill
    assert svg.count(f'fill="{GROUND}"') == (1 if variant == "grounds" else 0)
    print(f"  {variant}: {WIDTH_MM:.0f} x {WIDTH_MM*H/W:.1f} mm, labels {LABEL_PT} pt")


def main(outdir):
    outdir.mkdir(parents=True, exist_ok=True)
    for variant in ("rows", "grounds"):
        svg, H = build(variant)
        _check(variant, svg, H)
        p = outdir / f"Fig1a_draft_{variant}.svg"
        p.write_text(svg, encoding="utf-8")
        print("wrote", p)


if __name__ == "__main__":
    main(Path(sys.argv[1]) if len(sys.argv) > 1
         else REPO / "working" / "Fig1a_drafts")
