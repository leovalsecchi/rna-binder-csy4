"""The F30 three-way junction carrying two Csy4 hairpins.

Drawn the way the F30 literature draws it, and the way RNA2Drawer calls the
"flat base" style: secondary-structure line art, not filled tubes. Every helix
is two parallel rails with one rung per base pair, a hairpin closes with an open
rounded turn, and the three-way junction is left rung-free where the strands
cross it. Nothing is filled, so colour is carried by the stroke alone and each
module reads as one colour.

Colour follows COLORS.md rule 6's line that structure inside a cartoon is
neutral: the scaffold is grey, so the only saturated thing in the drawing is the
Csy4 hairpin, which is what the figure is about. The host transcript is not
drawn at all - the scaffold stands on its own, and whatever places it decides
what it sits in. Note this is NOT the pale blue the F30 module carries in a construct map
- there the box is a fill with a neutral outline, here it is the stroke itself.

Proportions are taken off the published F30-aptamer figure: a long bottom stem
of seven pairs, a short rung-free junction, arms running out at about 35 degrees
from vertical. The arms here are short because a Csy4 hairpin is 16 nt - a
five-pair stem and a loop - against the ~50 nt aptamers that figure carries.

    python code/cartoon_library/f30.py <outdir>
"""
import math
import sys
from pathlib import Path
_HERE = Path(__file__).resolve()
sys.path.insert(0, str(_HERE.parents[1]))
sys.path.insert(0, str(_HERE.parent))

import hairpin as AT
from common import palette

F30_COL = palette.NEUTRAL["500"]   # the scaffold: structure, so neutral
CSY4_COL = palette.THEME["600"]    # the hairpin we present
INK = palette.INK
MUT = palette.NEUTRAL["600"]
BG = palette.SURFACE
FONT = "Arial, Helvetica, sans-serif"

# A helix here is the AT hairpin's helix, at its weight, so the two drawings
# can sit together without one looking like a sketch of the other. The numbers
# come out of that module rather than being copied, so it stays true if the
# hairpin is retuned. Everything else below was drawn against HW 17, so K
# carries those lengths over.
HW = AT.STEM_HALF * 2 * AT.SCALE   # rail separation: the AT stem, as rendered
SW = AT.SW * AT.SCALE              # rails and turns, likewise
RUNG = AT.RUNG * AT.SCALE          # a base pair
K = HW / 17.0                      # this drawing's original duplex width
ANG = math.radians(35)

# A duplex is as long as the pairs in it, at the AT hairpin's pitch: the two
# drawings are the same molecule family, so a base pair has to be the same size
# in both or one looks stretched beside the other. This drawing used a pitch
# about 1.8x the hairpin's, which is what made its arms look drawn out.
PITCH = AT.SP_RUNG_D * AT.SCALE    # one base pair along a helix, as rendered

J = (170.0 * K, 188.0 * K)         # where the three helices meet
STEM_N = 5                         # bottom duplex
JUNC_N = 3                         # the short duplex out of the junction
ARM_N = 4                          # the Csy4 stem: four pairs
STEM_LEN, JUNC_LEN, ARM_LEN = (n * PITCH for n in (STEM_N, JUNC_N, ARM_N))
# the internal loop before the hairpin: the bow has to shrink with the loop,
# or a shorter loop pinches instead of bowing (it kept the old 5/18 ratio)
BULGE_LEN = 1.5 * PITCH
BULGE_OUT = BULGE_LEN * 5.0 / 18.0
LOOP_FLARE = 8.0 * K               # how far the terminal loop swells past the stem
JUNC_GAP = HW * 1.45               # how far the helices stand off the junction;
                                   # with HW it sets the multiloop's radius
FS, LS = 13.0 * K, 12.0 * K        # label sizes, in the same units

_DL = (-math.sin(ANG), -math.cos(ANG))
_DR = (math.sin(ANG), -math.cos(ANG))
_DD = (0.0, 1.0)

_pts = []          # every coordinate drawn, so the canvas can crop to it


def _p(x, y):
    _pts.append((x, y))
    return x, y


def _perp(d):
    return (-d[1], d[0])


def _along(p, d, n):
    return (p[0] + d[0] * n, p[1] + d[1] * n)


def _seg(a, b, col, w=SW):
    _p(*a), _p(*b)
    return (f'<path d="M{a[0]:.1f},{a[1]:.1f} L{b[0]:.1f},{b[1]:.1f}" '
            f'fill="none" stroke="{col}" stroke-width="{w}" '
            f'stroke-linecap="round"/>')


def helix(start, d, length, n_rungs, col):
    """Two rails plus one rung per pair. Returns (svg, rail ends, axis end)."""
    px, py = _perp(d)
    end = _along(start, d, length)
    a0 = (start[0] + px * HW / 2, start[1] + py * HW / 2)
    b0 = (start[0] - px * HW / 2, start[1] - py * HW / 2)
    a1 = (end[0] + px * HW / 2, end[1] + py * HW / 2)
    b1 = (end[0] - px * HW / 2, end[1] - py * HW / 2)
    o = [_seg(a0, a1, col), _seg(b0, b1, col)]
    for i in range(n_rungs):
        c = _along(start, d, length * (i + 0.5) / n_rungs)
        o.append(_seg((c[0] + px * HW / 2, c[1] + py * HW / 2),
                      (c[0] - px * HW / 2, c[1] - py * HW / 2), col, RUNG))
    return "".join(o), (a0, b0), (a1, b1), end


def bulge(start, d, col, length=BULGE_LEN, out=BULGE_OUT):
    """An internal loop: the two rails bow apart, with no pair between them."""
    px, py = _perp(d)
    end = _along(start, d, length)
    o = []
    for sign in (1, -1):
        a = (start[0] + px * sign * HW / 2, start[1] + py * sign * HW / 2)
        b = (end[0] + px * sign * HW / 2, end[1] + py * sign * HW / 2)
        w = (HW / 2 + out * 1.33) * sign
        c1 = _along((start[0] + px * w, start[1] + py * w), d, length * 0.28)
        c2 = _along((end[0] + px * w, end[1] + py * w), d, -length * 0.28)
        _p(*a), _p(*b), _p(*c1), _p(*c2)
        o.append(f'<path d="M{a[0]:.1f},{a[1]:.1f} C{c1[0]:.1f},{c1[1]:.1f} '
                 f'{c2[0]:.1f},{c2[1]:.1f} {b[0]:.1f},{b[1]:.1f}" fill="none" '
                 f'stroke="{col}" stroke-width="{SW}" stroke-linecap="round"/>')
    return "".join(o), end


def loop_turn(end, d, col, reach=1.2, flare=LOOP_FLARE):
    """The open rounded turn that closes a hairpin.

    `flare` pushes the control points sideways as well as forward, so the loop
    swells a little past the stem instead of capping it flush - a loop holds
    more nucleotides than the helix below it is wide.
    """
    px, py = _perp(d)
    a = (end[0] + px * HW / 2, end[1] + py * HW / 2)
    b = (end[0] - px * HW / 2, end[1] - py * HW / 2)
    ca = _along((a[0] + px * flare, a[1] + py * flare), d, HW * reach)
    cb = _along((b[0] - px * flare, b[1] - py * flare), d, HW * reach)
    _p(*_along(end, d, HW * reach))
    _p(ca[0], ca[1]), _p(cb[0], cb[1])
    return (f'<path d="M{a[0]:.1f},{a[1]:.1f} C{ca[0]:.1f},{ca[1]:.1f} '
            f'{cb[0]:.1f},{cb[1]:.1f} {b[0]:.1f},{b[1]:.1f}" fill="none" '
            f'stroke="{col}" stroke-width="{SW}" stroke-linecap="round"/>')


def _arc(a, b, col):
    """Join two rail ends along the junction circle.

    All six rail ends sit the same distance from J - the three helices share a
    width and a stand-off - so a circular arc of that radius through any two of
    them is centred on J, and the multiloop closes as a true circle.
    """
    r = math.dist(a, J)
    ax, ay = a[0] - J[0], a[1] - J[1]
    bx, by = b[0] - J[0], b[1] - J[1]
    sweep = 1 if ax * by - ay * bx > 0 else 0      # the short way, bulging out
    _p(*a), _p(*b)
    _p(J[0] - r, J[1] - r), _p(J[0] + r, J[1] + r)
    return (f'<path d="M{a[0]:.1f},{a[1]:.1f} A{r:.1f},{r:.1f} 0 0 {sweep} '
            f'{b[0]:.1f},{b[1]:.1f}" fill="none" stroke="{col}" '
            f'stroke-width="{SW}" stroke-linecap="round"/>')


def build(labelled=False):
    _pts.clear()
    o, ends = [], {}

    for side, d in (("L", _DL), ("R", _DR)):
        start = _along(J, d, JUNC_GAP)
        svg, (a0, b0), _, j_end = helix(start, d, JUNC_LEN, JUNC_N, F30_COL)
        o.append(svg)
        svg, j_end = bulge(j_end, d, F30_COL)      # the internal loop
        o.append(svg)
        svg, _, _, tip = helix(j_end, d, ARM_LEN, ARM_N, CSY4_COL)
        o.append(svg)
        o.append(loop_turn(tip, d, CSY4_COL))
        ends[side] = (a0, b0)

    stem_start = _along(J, _DD, JUNC_GAP)
    svg, (sa0, sb0), (sa1, sb1), _ = helix(stem_start, _DD, STEM_LEN,
                                           STEM_N, F30_COL)
    o.append(svg)

    # close the junction: left arm over to right arm, each arm down to the stem
    # close the multiloop as a circle through all six rail ends
    o.append(_arc(ends["L"][0], ends["R"][1], F30_COL))
    o.append(_arc(ends["L"][1], sa0, F30_COL))
    o.append(_arc(sb0, ends["R"][0], F30_COL))

    # No host transcript: the scaffold is drawn on its own. The two rails of
    # the bottom stem end where the stem ends; what it is embedded in belongs
    # to whatever places this drawing.

    if labelled:
        def t(x, y, s, sz=FS, col=INK, anchor="middle"):
            # the crop has to know where the text actually lands, and an
            # anchored label runs to one side of its x, not around it
            w = len(s) * sz * 0.58
            left = {"middle": x - w / 2, "start": x, "end": x - w}[anchor]
            _p(left, y - sz), _p(left + w, y + 4)
            return (f'<text x="{x:.0f}" y="{y:.0f}" font-family="{FONT}" '
                    f'font-size="{sz}" fill="{col}" font-weight="600" '
                    f'text-anchor="{anchor}">{s}</text>')
        tip_l = _along(_along(J, _DL, HW * 0.5 + JUNC_LEN), _DL, ARM_LEN)
        o.append(t(tip_l[0] - 28 * K, tip_l[1] - 16 * K, "Csy4 hairpin", FS,
                   CSY4_COL, "end"))
        o.append(t(J[0] + 78 * K, J[1] + 46 * K, "F30 array", FS, F30_COL,
                   "start"))

    m = 10.0 * K
    xs, ys = [p[0] for p in _pts], [p[1] for p in _pts]
    x0, y0 = min(xs) - m, min(ys) - m
    w, h = max(xs) - min(xs) + 2 * m, max(ys) - min(ys) + 2 * m
    return (f'<svg xmlns="http://www.w3.org/2000/svg" width="{w:.0f}" '
            f'height="{h:.0f}" viewBox="{x0:.1f} {y0:.1f} {w:.1f} {h:.1f}">'
            + "".join(o) + "</svg>")


W, H = 360, 330          # nominal only; build() crops to the drawing


if __name__ == "__main__":
    out = Path(sys.argv[1] if len(sys.argv) > 1 else ".")
    (out / "f30_scaffold.svg").write_text(build(), encoding="utf-8")
    (out / "f30_scaffold_labelled.svg").write_text(build(True), encoding="utf-8")
    print("wrote f30_scaffold{,_labelled}.svg to", out)
