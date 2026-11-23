"""Fig 4b - the four cell populations the library lands in.

One population per cell cluster, laid out 2x2 in the SAME arrangement as the
target call-out in Fig 4a: A and B are the two scaffold geometries across, on
and off down. A reader can put the two panels side by side and map a
population onto the target it carries.

Both drawings are reused whole and unedited:

* the cluster is `output/cartoon_library/svg/cell_cluster.svg` - three cells,
  one large and two small, nuclei and all.
* the target is the Fig1f glyph for that construct, the same four files Fig
  4a's call-out uses.

The target is not drawn loose in the cell. A small circle marks a patch of
the large cell's cytoplasm and a cone opens from it into a zoom circle
alongside, and the scaffold is drawn in there. That way the magnification is
explicit - these are molecules in the cytoplasm, seen closer - instead of one
giant scaffold floating inside a cell it could not fit in.

No landing pad is drawn: the integration is carried by Fig 4a and the caption,
and a locus glyph at this size would be a smudge.

    python code/Fig4b_cell_populations/make_figure.py
"""
import math
import re
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from common import glyphs, palette                                    # noqa: E402
from common.paths import ROOT                                 # noqa: E402

CLUSTER = ROOT / "output/cartoon_library/svg/cell_cluster.svg"
GLYPHS = ROOT / "output/Fig1f_scaffold_splicing/elements"
OUT = ROOT / "output/Fig4b_cell_populations"

INK, OUTLINE = palette.INK, palette.SYNBIO["outline"]

# All of these are in cell_cluster.svg's own units. Its large cell is at
# (116, 80) with r 46 and a nucleus of r 19, leaving a cytoplasm ring only
# about 27 wide; the patch is placed in the widest clear part of it, 7 units
# off both the nucleus and the membrane, on the side the zoom sits.
SRC = (149.9, 83.0, 8.0)                     # cx, cy, r - the patch zoomed
ZOOM = (298.0, 92.0, 74.0)                    # cx, cy, r - where it is shown

# The zoom reads top to bottom as a reaction: the target, a plus, and the one
# binder the cell drew from the pool - shown as three variants, "/"
# between them, since which one a given cell got is exactly what is unknown.
# the two scaffolds have very different aspects (the Y array is 1.5 wide, the
# circular one square), so the box is kept near-square or the square glyph
# scales to its height and comes out half the size of the other
GLYPH_DY, GLYPH_BOX = -30.0, (0.575, 0.412)   # scaffold: offset, then w/h of r
PLUS_DY, PLUS_ARM = 14.0, 5.2
BINDER_DY, BINDER_DX = 42.0, 36.0             # the row of variants, and pitch
# 2026-11-18 (Leo, from his supervisor): "color them for the generation
# methods, backbone design method stacks we used" - three steps of one blue
# said "three members of one library" but not that the library is a MIXTURE
# OF METHODS, which is the point of the panel. They now carry the campaign's
# own stack colours, the same ones Fig 4c's key uses, so the two panels agree.
#
# Three of the five, not all five: the row sits 42 units below the zoom
# circle's centre, where the circle is only +/-60.9 wide, and five helices
# plus their separators need a pitch of ~22 - the glyphs would touch. These
# three are canonical-order and well separated in hue, and they span both
# design modes (two diffusion, one hallucination). They are an illustration
# of "a mixture", not an enumeration - say so in the caption.
BINDER_STACKS = ("RFD3 + LigandMPNN", "BoltzGen + LigandMPNN",
                 "BoltzDesign + LigandMPNN")

# A is the Y array (mod1), B the circular scaffold (mod4). Laid out 2x2 in
# Fig 4a's arrangement: the geometries across, positive and negative down.
POPS = [("positive A", "f30_variation_positive", 0, 0),
        ("positive B", "circular_scaffold_positive", 1, 0),
        ("negative A", "f30_variation_negative", 0, 1),
        ("negative B", "circular_scaffold_negative", 1, 1)]

POP_W = 84.0              # one population's drawn width, in figure units
GAP_X, GAP_Y = 7.0, 4.0
LABEL_DROP = 3.8


def bbox(body):
    """The ink's own extent. The cluster's viewBox is much larger than its
    drawing - about a quarter of it is empty - so laying the grid out on the
    viewBox leaves the labels stranded and the clusters off-centre."""
    xs, ys = [], []
    for d in re.findall(r'<path[^>]*\sd="([^"]+)"', body):
        n = [float(v) for v in re.findall(r"-?[\d.]+", d)]
        xs += n[0::2]
        ys += n[1::2]
    for cx, cy, r in re.findall(
            r'<circle cx="([\d.-]+)" cy="([\d.-]+)" r="([\d.-]+)"', body):
        cx, cy, r = float(cx), float(cy), float(r)
        xs += [cx - r, cx + r]
        ys += [cy - r, cy + r]
    pad = 1.5                                  # half the heaviest stroke
    return min(xs) - pad, min(ys) - pad, max(xs) + pad, max(ys) + pad


def cluster_parts():
    """The cluster exactly as it is, and the box it actually fills."""
    t = CLUSTER.read_text()
    body = t[t.index(">", t.index("<svg")) + 1:t.rindex("</svg>")]
    return bbox(body), body


def cone(src, zoom):
    """The two outer tangents from the source patch to the zoom circle.

    Both lines touch each circle, so the cone reads as one solid magnifier
    rather than two lines aimed near it. For circles a distance d apart the
    tangent points lie at gamma +- acos((r1 - r2) / d) around each centre.
    """
    x1, y1, r1 = src
    x2, y2, r2 = zoom
    d = math.hypot(x2 - x1, y2 - y1)
    gamma = math.atan2(y2 - y1, x2 - x1)
    alpha = math.acos(max(-1.0, min(1.0, (r1 - r2) / d)))
    out = []
    for a in (gamma + alpha, gamma - alpha):
        ca, sa = math.cos(a), math.sin(a)
        out.append(f"M{x1 + r1 * ca:.2f} {y1 + r1 * sa:.2f} "
                   f"L{x2 + r2 * ca:.2f} {y2 + r2 * sa:.2f}")
    return " ".join(out)


def glyph(name, x0, y0, x1, y1):
    """One Fig1f scaffold glyph, fitted to the box by its INK.

    Fitting these by their viewBox is what kept the targets looking small: a
    matplotlib figure is mostly margin, and the Y array's drawing is 62 x 21
    inside an 88 x 59 box - three quarters of what was being scaled was empty
    padding. Measure the paths instead.
    """
    t = (GLYPHS / f"{name}.svg").read_text()
    # matplotlib's last </g> closes figure_1, so slicing to it leaves axes_1 whole
    body = t[t.index('<g id="axes_1">'):t.rindex("</g>")]
    body = re.sub(r'\sclip-path="url\(#[^)]*\)"', "", body)
    body = re.sub(r'\sid="[^"]*"', "", body)
    xs, ys = [], []
    for d in re.findall(r'\sd="([^"]+)"', body):
        n = [float(v) for v in re.findall(r"-?[\d.]+", d)]
        xs += n[0::2]
        ys += n[1::2]
    pad = 0.7                                  # half the heaviest stroke
    ix0, iy0 = min(xs) - pad, min(ys) - pad
    ix1, iy1 = max(xs) + pad, max(ys) + pad
    k = min((x1 - x0) / (ix1 - ix0), (y1 - y0) / (iy1 - iy0))
    tx = (x0 + x1) / 2 - k * (ix0 + ix1) / 2
    ty = (y0 + y1) / 2 - k * (iy0 + iy1) / 2
    return (f'<g transform="translate({tx:.3f} {ty:.3f}) scale({k:.5f})">'
            f"{body}</g>")


def binders(cx, cy, fs):
    """Three binder variants in a row, "/" between them.

    Each is one helix from `common.glyphs.ribbon_helix`, which builds a
    two-tone ribbon off real helix geometry. A cell takes a single insert
    from the oligo pool, so it holds one of these - which is what the "/"
    says. They sit beside the target rather than docked on it: the negative
    scaffolds have no site to dock to.
    """
    o = []
    for i, stack in enumerate(BINDER_STACKS):
        x = cx + (i - 1) * BINDER_DX
        o.append(f'<g color="{palette.stack_color(stack)}">'
                 f'{glyphs.ribbon_helix(x, cy, turns=2.0, R=5.4, pitch=13.2, W=6.9, ang=-78)}'
                 f"</g>")
        if i:
            # 2026-11-18 (Leo): "don't do an or. Just put like a slash in
            # between them." Same role, a third of the width, and it does not
            # read as a word competing with the labels around it.
            o.append(f'<text class="orr" x="{x - BINDER_DX / 2:g}" '
                     f'y="{cy + fs * 0.36:g}" font-size="{fs:g}">/</text>')
    return "".join(o)


def population(body, name, fs):
    """One cluster, its zoomed patch of cytoplasm, and what is in it."""
    sx, sy, sr = SRC
    zx, zy, zr = ZOOM
    gw, gh = zr * GLYPH_BOX[0], zr * GLYPH_BOX[1]
    gy, py = zy + GLYPH_DY, zy + PLUS_DY
    return (body
            + f'<path class="cone" d="{cone(SRC, ZOOM)}"/>'
            + f'<circle class="patch" cx="{sx:g}" cy="{sy:g}" r="{sr:g}"/>'
            + f'<circle class="zoom" cx="{zx:g}" cy="{zy:g}" r="{zr:g}"/>'
            + glyph(name, zx - gw, gy - gh, zx + gw, gy + gh)
            + f'<path class="plus" d="M{zx - PLUS_ARM:g} {py:g} '
              f'L{zx + PLUS_ARM:g} {py:g} M{zx:g} {py - PLUS_ARM:g} '
              f'L{zx:g} {py + PLUS_ARM:g}"/>'
            + binders(zx, zy + BINDER_DY, fs))


def main():
    (cx0, cy0, cx1, cy1), body = cluster_parts()
    zx, zy, zr = ZOOM
    x0, y0 = min(cx0, zx - zr), min(cy0, zy - zr)
    x1, y1 = max(cx1, zx + zr), max(cy1, zy + zr)
    k = POP_W / (x1 - x0)                      # cluster units -> figure units
    ph = (y1 - y0) * k
    # The label hangs under the CELLS, not under the whole unit. The zoom
    # circle reaches far below the cluster, so measuring from the unit's
    # bottom left the label stranded - and closing that gap only pushed it
    # into the circle's lower arc. Under the cluster there is nothing to its
    # right at that height, so it can sit close.
    lx = (cx0 + cx1) / 2 - x0
    ly = cy1 - y0
    col_w, row_h = POP_W + GAP_X, ph + GAP_Y

    o = []
    for label, name, col, row in POPS:
        x, y = col * col_w, row * row_h
        o.append(f'<g transform="translate({x:g} {y:g}) scale({k:.5f}) '
                 f'translate({-x0:g} {-y0:g})">'
                 f'{population(body, name, 4.5 / k)}</g>')
        o.append(f'<text class="lbl" x="{x + lx * k:g}" '
                 f'y="{y + ly * k + LABEL_DROP:g}">{label}</text>')

    w, h = col_w + POP_W, row_h + ph + 2
    svg = (f'<svg xmlns="http://www.w3.org/2000/svg" '
           f'viewBox="-4 -4 {w + 8:g} {h + 8:g}"><style>'
           f'.lbl{{font-family:Arial;font-size:5.5px;fill:{INK};'
           f'text-anchor:middle}}'
           f'.cone{{fill:none;stroke:{palette.NEUTRAL["400"]};stroke-width:1.4}}'
           f'.patch{{fill:none;stroke:{palette.NEUTRAL["400"]};stroke-width:1.4}}'
           f'.zoom{{fill:{palette.SURFACE};stroke:{OUTLINE};stroke-width:2.2}}'
           f'.plus{{fill:none;stroke:{INK};stroke-width:1.7;'
           f'stroke-linecap:round}}'
           f'.orr{{font-family:Arial;fill:{palette.NEUTRAL["500"]};'
           f'text-anchor:middle}}'
           f'</style>' + "".join(o) + "</svg>")

    OUT.mkdir(parents=True, exist_ok=True)
    p = OUT / "cell_populations.svg"
    p.write_text(svg, encoding="utf-8")
    subprocess.run(["rsvg-convert", "-w", "1700", "-b", "white", str(p),
                    "-o", str(p.with_suffix(".png"))], check=False)
    subprocess.run(["inkscape", str(p), "--export-type=pdf",
                    f"--export-filename={p.with_suffix('.pdf')}"], check=False,
                   stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    stray = sorted({c.lower() for c in re.findall(r"#[0-9a-fA-F]{6}", svg)} -
                   {c.lower() for c in (*palette.THEME.values(),
                                        *palette.NEUTRAL.values(),
                                        *palette.SYNBIO.values(),
                                        *palette.DIVERGING, palette.INK,
                                        palette.GRID, palette.SURFACE)})
    print(f"wrote {p}" + (f"  LEFTOVER {stray}" if stray else "  palette ok"))


if __name__ == "__main__":
    main()
