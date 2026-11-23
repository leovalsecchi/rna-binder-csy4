"""Bxb1 cassette integration: donor plasmid, landing pad, integrated locus.

The three rows of the standard promoterless-integration figure, built from
`common/constructs.py` so every module is the shape and colour it has
everywhere else. The donor carries attB sites and no promoter; the genomic
landing pad carries the promoter and attP; recombination leaves hybrid attL and
attR flanking the cassette, now driven by the genomic promoter.

The crossed strokes between the rows are the recombination itself: attB-GT
pairs with attP-GT and attB-GA with attP-GA, which is what makes the
integration directional.

    python code/cartoon_library/integration.py <outdir>
"""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from common import palette
from common import constructs as K

INK = palette.INK
MUT = palette.NEUTRAL["600"]
RULE = palette.NEUTRAL["300"]
LINE = palette.SYNBIO["outline"]
BG = palette.SURFACE
FONT = K.FONT

M = 26.0           # page margin
ROW = 104.0        # between rows
LABEL = 15.0       # label sits this far above a row's strand
ROWLAB = 40.0      # row titles clear the site labels above the strand

DONOR = ["attB_GT", "at", "binder", "sms", "scaffold", "attB_GA"]
# the landing pad carries mNeonGreen between its att sites, so an
# unintegrated cell is green and the donor lands in its place
PAD = ["promoter", "spacer", "attP_GT", "mneongreen", "attP_GA", "pa"]
INTEGRATED = ["promoter", "attL", "at", "binder", "sms", "f30", "attR", "pa"]


def _t(x, y, s, sz=12.5, col=INK, weight=600, anchor="start"):
    return (f'<text x="{x:.1f}" y="{y:.1f}" font-family="{FONT}" '
            f'font-size="{sz}" fill="{col}" font-weight="{weight}" '
            f'text-anchor="{anchor}">{s}</text>')


def _arrow(x, y0, y1, col=MUT, w=2.0, head=7.0):
    return (f'<path d="M{x:.1f},{y0:.1f} L{x:.1f},{y1-head:.1f}" fill="none" '
            f'stroke="{col}" stroke-width="{w}" stroke-linecap="round"/>'
            f'<path d="M{x:.1f},{y1:.1f} L{x-head*0.62:.1f},{y1-head:.1f} '
            f'L{x+head*0.62:.1f},{y1-head:.1f} Z" fill="{col}"/>')


def _exchange(x_from, y_from, x_to, y_to, col=MUT):
    """One recombination stroke, bowed so the GT and GA pairs cross."""
    my = (y_from + y_to) / 2
    return (f'<path d="M{x_from:.1f},{y_from:.1f} C{x_from:.1f},{my:.1f} '
            f'{x_to:.1f},{my:.1f} {x_to:.1f},{y_to:.1f}" fill="none" '
            f'stroke="{col}" stroke-width="1.6" stroke-dasharray="4 4" '
            f'stroke-linecap="round"/>')


def _site_x(x0, items, name):
    """Centre of a named site inside a bar starting at x0."""
    cursor = x0
    for key in items:
        w = K._w(key)
        if key == name:
            return cursor + w / 2
        if not K._is_tick(key):
            cursor += w + K.GAP
    raise KeyError(name)


def build():
    w_donor, w_pad, w_int = (K.width(b) for b in (DONOR, PAD, INTEGRATED))
    inner = max(w_donor, w_pad, w_int)
    W = inner + 2 * M + 150          # room for the row labels on the right
    x_donor = M + (inner - w_donor) / 2
    x_pad = M + (inner - w_pad) / 2
    x_int = M + (inner - w_int) / 2

    y_donor, y_pad = 86.0, 86.0 + ROW
    y_int = y_pad + ROW + 54
    H = y_int + 70

    o = [f'<svg xmlns="http://www.w3.org/2000/svg" width="{W:.0f}" '
         f'height="{H:.0f}" viewBox="0 0 {W:.0f} {H:.0f}">',
         f'<rect width="{W:.0f}" height="{H:.0f}" fill="{BG}"/>']

    # the donor sits inside its plasmid, drawn as an open rounded outline
    pw, ph = w_donor + 56, 62.0
    o.append(f'<rect x="{x_donor-28:.1f}" y="{y_donor-ph/2:.1f}" '
             f'width="{pw:.1f}" height="{ph:.1f}" rx="26" fill="none" '
             f'stroke="{RULE}" stroke-width="1.6"/>')
    o.append(_t(x_donor - 28, y_donor - ph / 2 - 10, "Promoterless donor plasmid"))
    o.append(K.bar(x_donor, y_donor, DONOR))

    o.append(_t(M, y_pad - ROWLAB, "Genomic landing pad"))
    o.append(K.bar(x_pad, y_pad, PAD, strand=40))

    # recombination: each attB finds the attP of the same dinucleotide
    for a, b in (("attB_GT", "attP_GT"), ("attB_GA", "attP_GA")):
        o.append(_exchange(_site_x(x_donor, DONOR, a), y_donor + K.H / 2 + 4,
                           _site_x(x_pad, PAD, b), y_pad - K.H / 2 - 20))
    o.append(_t(mid_x := M + inner / 2, (y_donor + y_pad) / 2 - 10,
            "Bxb1", 12.5, MUT, 600, "middle"))

    mid = mid_x
    o.append(_arrow(mid, y_pad + K.H / 2 + 16, y_int - K.H / 2 - 48))
    o.append(_t(mid + 12, y_pad + ROW / 2 + 14,
                "Landing-pad integration in HEK-Spartan", 12, MUT, 500))

    o.append(_t(M, y_int - ROWLAB, "Integrated cassette, hybrid sites"))
    o.append(K.bar(x_int, y_int, INTEGRATED, strand=40))
    o.append("</svg>")
    return "".join(o)


if __name__ == "__main__":
    out = Path(sys.argv[1] if len(sys.argv) > 1 else ".")
    (out / "cassette_integration.svg").write_text(build(), encoding="utf-8")
    print("wrote cassette_integration.svg to", out)
