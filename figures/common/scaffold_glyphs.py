"""Scaffold line glyphs (fig04 panel a) drawn in mm on a caller-supplied axes.

Both glyphs assume an inverted y axis (ylim(H, 0)), so "up" is decreasing y.
A three-way junction splits SYMMETRICALLY: both arms the same length, both at
the same angle either side (Leo, 2026-11-17).

`dot_col` is the hairpin tip: the Csy4 site colour for a positive construct,
a neutral grey for a negative. The tip dot is centred ON THE END of its arm, so
it caps the stick. It used to sit a further 1.0 beyond it - and in `array` and
`tree` that offset was applied straight up rather than along the arm, so on a
diagonal arm the dot floated off the end of the stick instead of capping it
(Leo, 2026-11-17).

Every glyph takes a scale `k`. It defaults to 1.0, so the sizes fig04 panel a
was drawn at are unchanged; a caller that needs the same glyph small (the Fig 1
workflow strip draws four of them in a 2 x 2) passes k < 1 instead of redrawing
the shapes at another size and letting the two drift apart.
"""
import math
from matplotlib.patches import Circle



def _seg(ax, p, q, col, lw=1.3):
    ax.plot([p[0], q[0]], [p[1], q[1]], color=col, lw=lw,
            solid_capstyle="round", zorder=2)


def _dot(ax, x, y, col, r=0.8):
    ax.add_patch(Circle((x, y), r, color=col, lw=0, zorder=3))


def array(ax, cx, y1, col, dot_col, x0, x1, n=4, stem=2.4, lw=1.3, k=1.0,
          dot_r=None):
    """mod1 / mod3: n Y units on one strand, two hairpins each.

    `stem` is what separates mod1 from mod3 - the same array with long or short
    arms - so pass it rather than drawing a second glyph."""
    yb = y1 - 4.5 * k
    rr = 0.8 * k if dot_r is None else dot_r
    _seg(ax, (x0, yb), (x1, yb), col, lw)          # the strand runs through
    for j in range(n):
        xc = cx - 2.5 * k * (n - 1) + 5.0 * k * j
        _seg(ax, (xc, yb), (xc, yb - stem * k), col, lw)
        for s in (-1, 1):
            tip = (xc + s * 1.4 * k, yb - (stem + 3.2) * k)
            _seg(ax, (xc, yb - stem * k), tip, col, lw)
            _dot(ax, tip[0], tip[1], dot_col, rr)


def ring(ax, cx, y1, col, dot_col, x0, x1, n=8, r=3.2, lw=1.3, k=1.0,
         dot_r=None, stem=3.3):
    """mod4 / circular scaffold: n hairpins on a ring, stem down to the strand.

    `stem` is the bare helix between the strand and the ring, measured from the
    strand to the BOTTOM OF THE RING - so it stays what it says whatever `r` is.
    It was 5.8, which left the ring floating a ring-diameter above the
    transcript and read as a lollipop; 3.3 closes that up (Leo, 2026-11-17).
    The lowest hairpin spikes hang 3.98 below the ring centre, so the clearance
    left under them is stem + r - 3.98, about 2.5 units - still clear of the
    strand at every size the glyph is drawn."""
    rr = 0.8 * k if dot_r is None else dot_r
    r = r * k
    yb = y1 - 3.5 * k
    yc = yb - stem * k - r
    _seg(ax, (x0, yb), (x1, yb), col, lw)
    _seg(ax, (cx, yc + r), (cx, yb), col, lw)
    ax.add_patch(Circle((cx, yc), r, fill=False, ec=col, lw=lw, zorder=2))
    for j in range(n):
        a = math.radians(90 + 360 * (j + 1) / (n + 1))
        ca, sa = math.cos(a), math.sin(a)
        tip = (cx + (r + 2.0 * k) * ca, yc + (r + 2.0 * k) * sa)
        _seg(ax, (cx + r * ca, yc + r * sa), tip, col, lw)
        _dot(ax, tip[0], tip[1], dot_col, rr)


def tree(ax, cx, y1, col, dot_col, x0, x1, spec, lw=1.3, k=1.0,
         bp=0.42, aper=180.0, decay=0.75, dot_r=None):
    """mod2 / branched scaffold, drawn from its REAL branching.

    `spec` is the molecule's own tree, (helix_bp, [children]) with a leaf
    standing for a Csy4 hairpin - derived from the dot-bracket target, not
    invented here, so the glyph shows the shape the design actually has. mod2
    is NOT a balanced tree: it cascades, each junction throwing off a
    two-hairpin arm and carrying on, and only the last one divides evenly. A
    balanced tree drawn in its place is a different molecule.

    Segment lengths are PROPORTIONAL to base pairs (`bp` units per bp), which is
    why the hairpin arms are long against the connecting helices: 16 bp against
    6-13.

    ANGLES: a cascade cannot be drawn at one fixed junction angle - whatever the
    value, the deep side folds back over the shallow one and the arms cross (128
    degrees, chosen on tip distance alone, crossed four segments). So each
    branch carries an angular WEDGE instead: a junction halves its own wedge and
    sends each child down the middle of one half. That keeps the junction
    SYMMETRICAL - equal angle either side of the incoming helix - while the two
    sub-trees sit in wedges that cannot overlap. `decay` shrinks the wedge a
    child inherits. 180 / 0.75 was searched over aperture x decay against three
    tests at once: no segment may cross another, no arm may dip below the
    transcript strand, and the eight tips must clear each other - they clear by
    0.9 mm at the size Fig 1 draws them.

    The cascade leans, because it always carries on to the same side. The glyph
    is therefore laid out first and drawn second, shifted so its BOUNDING BOX is
    centred on `cx`: the foot meets the strand off-centre, which is what an
    asymmetric molecule does, instead of the drawing hanging off one side.
    """
    yb = y1 - 4.5 * k
    rr = 0.8 * k if dot_r is None else dot_r
    _seg(ax, (x0, yb), (x1, yb), col, lw)          # the strand runs through

    segs, tips = [], []

    def grow(x, y, ang, wedge, node):
        length, kids = node
        nx = x + length * bp * k * math.sin(ang)
        ny = y - length * bp * k * math.cos(ang)
        segs.append(((x, y), (nx, ny)))
        if not kids:                                # a Csy4 hairpin
            tips.append((nx, ny))
            return
        half = wedge / 4.0                          # centre of each half-wedge
        for s_, kid in zip((-1, 1), kids):
            grow(nx, ny, ang + s_ * half, wedge * decay, kid)

    grow(0.0, yb, 0.0, math.radians(aper), spec)
    xs = [q[0] for s_ in segs for q in s_]
    dx = cx - (min(xs) + max(xs)) / 2.0
    for a, b in segs:
        _seg(ax, (a[0] + dx, a[1]), (b[0] + dx, b[1]), col, lw)
    for tx, ty in tips:                             # cap each hairpin
        _dot(ax, tx + dx, ty, dot_col, rr)
