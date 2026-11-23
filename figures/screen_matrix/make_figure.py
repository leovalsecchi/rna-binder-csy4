"""The co-transfection matrix: every target condition against every binder.

Down the left, the RNA target scaffold
that goes with it; across the top, the binder series; and the block they meet
in is one cell per combination.

The target glyphs are NOT redrawn here. `output/Fig1f_scaffold_splicing/
elements/` already holds them, so this panel and Fig 1f show the same objects:

    Target 1   f30_variation_positive      the F30 repeat, Csy4 sites coloured
    Target 2   circular_scaffold_positive  the circularised repeat
    Neg. 1     f30_variation_negative      the same F30 geometry, sites grey
    Neg. 2     drawn here                  a bare strand, crossed out

**One scale for every borrowed glyph.** Each is drawn at stroke-width 1.3 in
its own units, so fitting them individually to the slot gave them different
line weights - the circular one came out both smaller and finer. They now
share a scale, set by the tallest, and the non-target is drawn at that same
1.3 so all four rows match.

**One square.** CELL is the only size in the figure: every cell of every
strip, the header and the grid is CELL x CELL, so the whole thing is one
module repeated. The left strips are one cell wide, the header one cell tall,
and the ellipsis has a full cell like any binder. Every cell boundary, in both
directions and in all three containers, is a near-transparent grey hairline.

**The header cells are empty** - no captions, nothing: each is a square
waiting for a PDB render of its binder, plugged in later.

Colours are not settled. Everything the figure invents takes a named constant
at the top; the borrowed glyphs keep the colours Fig 1f gave them.

    python code/screen_matrix/make_figure.py [outdir]
"""
import base64
import io
import math
import re
import sys
from pathlib import Path

_HERE = Path(__file__).resolve()
REPO = _HERE.parents[2]
sys.path.insert(0, str(REPO / "code"))

from common import palette                                   # noqa: E402

OUT = REPO / "working" / "screen_matrix"
ELEMENTS = REPO / "output" / "Fig1f_scaffold_splicing" / "elements"
CLUSTER = REPO / "output" / "cartoon_library" / "svg" / "cell_cluster.svg"

CELL = 60.0            # ONE square. Every cell in the figure is this, and
                       # only this: the strips, the header, the grid.
SLOT = 34.0            # the glyph inside a cell, with its caption under it
CONTENT_Y = CELL / 2   # EVERYTHING sits on the centre of its square - glyph,
                       # plasmid, well, ellipsis alike. The slot is sized so
                       # the caption still fits below it inside the square
                       # rather than pushing the glyph off centre.
PAD = 7.0
GUTTER = 10.0
FS = 7.0
LINE = 8.5
GLYPH_LW = 1.3         # what scaffold_glyphs.py draws with, in glyph units

ELL_W = 30.0           # the ellipsis column, half a cell: it is a break in
                       # the series, not a condition, so it is not a square
WELL_H = 40.0          # the cell cluster in a well, drawn to this height
BINDERS = REPO / "working" / "screen_matrix" / "binders"   # render_binders.py
BINDER_PX = 300        # the embedded icon's longer side, in px. It prints
                       # about 9.5 mm wide, so this is still past 600 dpi.
def _tint(hex_col, f):
    """`hex_col` mixed f of the way into white."""
    r, g, b = (int(hex_col[i:i + 2], 16) for i in (1, 3, 5))
    return "#%02x%02x%02x" % tuple(round(c + (255 - c) * f) for c in (r, g, b))


POS_FILL = _tint(palette.DIVERGING[4], 0.55)   # a positive well: pale red,
                                               # the ramp's warm end let down
NEG_FILL = palette.NEUTRAL["200"]        # a negative one: grey
INK = palette.INK
BOX_FILL = palette.NEUTRAL["100"]
BOX_EDGE = palette.NEUTRAL["300"]
WALL = palette.NEUTRAL["400"]
WALL_LW = 0.9             # the dividers: a hairline, nearly transparent
WALL_OP = 0.35
STRAND = palette.NEUTRAL["600"]
FONT = "Arial, Helvetica, sans-serif"

# (caption, glyph, is this a real target?) - the flag colours the wells in
# that row: a binder meeting a real target reads out, a negative does not.
ROWS = [("Target 1", "f30_variation_positive", True),
        ("Target 2", "circular_scaffold_positive", True),
        ("Neg. 1", "f30_variation_negative", False),
        ("Neg. 2", None, False)]
COLS = [("Binder", "1"), ("Binder", "2"), None, ("Binder", "32")]
# 2026-11-18 (Leo): "the BFP row can be removed completely". The transfection
# control plasmid was identical in every well, so its strip carried no
# contrast - one constant repeated down four rows beside a matrix whose whole
# job is to vary. The plasmid glyph and its colour went with it; git has both
# if the strip is ever wanted back.
TITLES = {"target": ("RNA target", "array")}


def _text(x, y, s, col=INK, size=FS, anchor="middle"):
    return (f'<text x="{x:.1f}" y="{y:.1f}" font-family="{FONT}" '
            f'font-size="{size}" fill="{col}" text-anchor="{anchor}">{s}</text>')


def _box(x, y, w, h, name):
    return (f'<rect data-name="{name}" x="{x:.1f}" y="{y:.1f}" '
            f'width="{w:.1f}" height="{h:.1f}" fill="{BOX_FILL}" '
            f'stroke="{BOX_EDGE}" stroke-width="0.8"/>')


def _glyph_box(stem):
    t = (ELEMENTS / f"{stem}.svg").read_text()
    return [float(v) for v in
            re.search(r'viewBox="([^"]+)"', t).group(1).split()]


def _scale():
    """One scale for every borrowed glyph, set by the tallest.

    Fitting each glyph to the slot on its own is what made their lines
    disagree: a glyph half the size gets half the stroke.
    """
    return SLOT / max(_glyph_box(s)[3] for _, s, _p in ROWS if s)


def _embed(stem, cx, cy, k, tag):
    """A Fig 1f element at scale `k`, centred on (cx, cy).

    Matplotlib numbers its ids from 1 in every file, so two of them in one
    document would collide; each embed gets its own prefix. Its opaque page
    rect goes too, or it would punch a white hole in the container.
    """
    t = (ELEMENTS / f"{stem}.svg").read_text()
    vb = _glyph_box(stem)
    body = t[t.index(">", t.index("<svg")) + 1:t.rindex("</svg>")]
    body = re.sub(r"<metadata>.*?</metadata>", "", body, flags=re.S)
    body = re.sub(r'<g id="patch_1">.*?</g>', "", body, flags=re.S, count=1)
    body = re.sub(r'id="([^"]+)"', lambda m: f'id="{tag}-{m.group(1)}"', body)
    body = re.sub(r"url\(#([^)]+)\)", lambda m: f"url(#{tag}-{m.group(1)})",
                  body)
    return (f'<g data-name="{stem}" '
            f'transform="translate({cx - vb[2] * k / 2:.2f},'
            f'{cy - vb[3] * k / 2:.2f}) scale({k:.4f}) '
            f'translate({-vb[0]},{-vb[1]})">{body}</g>')


def _non_target(cx, cy, k, w=46.0):
    """A bare strand, crossed out: present, but nothing to bind.

    Its stroke is the glyphs' own 1.3 taken through the same scale, so this
    row's lines weigh the same as the three borrowed ones.
    """
    lw = GLYPH_LW * k
    a = 6.0
    return (f'<g data-name="Non-target">'
            f'<path d="M{cx - w / 2:.1f},{cy + 8:.1f} '
            f'L{cx + w / 2:.1f},{cy + 8:.1f}" fill="none" stroke="{STRAND}" '
            f'stroke-width="{lw:.2f}" stroke-linecap="round"/>'
            f'<path d="M{cx - a:.1f},{cy + 8 - a:.1f} '
            f'L{cx + a:.1f},{cy + 8 + a:.1f} M{cx + a:.1f},{cy + 8 - a:.1f} '
            f'L{cx - a:.1f},{cy + 8 + a:.1f}" fill="none" stroke="{STRAND}" '
            f'stroke-width="{lw:.2f}" stroke-linecap="round"/></g>')


def _content_box(svg):
    """The drawn extent of a hand-written glyph, not its viewBox.

    cell_cluster.svg is drawn off-centre in its own box, so centring it on
    the viewBox would leave every well visibly off its square. The numbers
    come out of the path data and the circles instead.
    """
    xs, ys = [], []
    for d in re.findall(r'\sd="([^"]+)"', svg):
        n = [float(v) for v in re.findall(r"-?\d+\.?\d*", d)]
        xs += n[0::2]
        ys += n[1::2]
    for c in re.finditer(r'<circle cx="([-\d.]+)" cy="([-\d.]+)" r="([-\d.]+)"',
                         svg):
        x, y, r = (float(v) for v in c.groups())
        xs += [x - r, x + r]
        ys += [y - r, y + r]
    return min(xs), min(ys), max(xs) - min(xs), max(ys) - min(ys)


def _cluster(cx, cy, name, fill):
    """A cluster of cells: what actually goes into a well.

    Borrowed whole from output/cartoon_library/svg/cell_cluster.svg, drawn to
    WELL_H and centred on its own content, so it sits on the square's centre
    like everything else. Only the cell bodies take `fill`; the outlines and
    nuclei keep the glyph's own colours.
    """
    svg = CLUSTER.read_text()
    body = svg[svg.index(">", svg.index("<svg")) + 1:svg.rindex("</svg>")]
    body = body.replace('fill="#ffffff"', f'fill="{fill}"')   # the cell bodies
    x0, y0, w, h = _content_box(body)
    k = WELL_H / h
    return (f'<g data-name="{name}" '
            f'transform="translate({cx - w * k / 2:.2f},{cy - h * k / 2:.2f}) '
            f'scale({k:.4f}) translate({-x0:.2f},{-y0:.2f})">{body}</g>')


def _binder(cx, cy, cid):
    """One candidate's cartoon, centred on (cx, cy) and fitted to SLOT.

    Written by `render_binders.py`: binder only, no hairpin, in the Fig 2b
    shared frame at y120. The PNG is cropped to its own ink before embedding,
    so the binder sits on the square's centre rather than on the centre of
    whatever empty frame ChimeraX saved around it, and it is downsampled on the
    way in - the icon prints at about 9.5 mm and the render is 600 px.

    Returns "" when the render is missing, so this script still builds a sheet
    with empty header cells, which is what it did before the renders existed.
    """
    f = BINDERS / f"candidate{cid:02d}.png"
    if not f.exists():
        return ""
    from PIL import Image
    im = Image.open(f).convert("RGBA")
    im = im.crop(im.getchannel("A").getbbox() or im.getbbox())
    im.thumbnail((BINDER_PX, BINDER_PX), Image.LANCZOS)
    buf = io.BytesIO()
    im.save(buf, "PNG", optimize=True)
    uri = base64.b64encode(buf.getvalue()).decode()
    k = SLOT / max(im.size)
    w, h = im.size[0] * k, im.size[1] * k
    return (f'<image data-name="Candidate {cid}" x="{cx - w / 2:.2f}" '
            f'y="{cy - h / 2:.2f}" width="{w:.2f}" height="{h:.2f}" '
            f'xlink:href="data:image/png;base64,{uri}"/>')


def _ellipsis(cx, cy, col=INK):
    d = 6.5
    return "".join(f'<circle cx="{cx + i * d:.1f}" cy="{cy:.1f}" r="1.7" '
                   f'fill="{col}"/>' for i in (-1, 0, 1))


def _wall(x0, y0, x1, y1, name):
    """A divider: one unbroken hairline, nearly transparent.

    It was dotted; at this scale the dots rendered as a ragged dashed line
    rather than a dotted one, so the rule is solid and carries its lightness
    in its opacity instead of in its gaps.
    """
    return (f'<path data-name="{name}" d="M{x0:.1f},{y0:.1f} '
            f'L{x1:.1f},{y1:.1f}" stroke="{WALL}" stroke-width="{WALL_LW}" '
            f'stroke-opacity="{WALL_OP}" fill="none"/>')


def _walls(edges, y, h):
    """Dividers on every column boundary, floor to ceiling."""
    return "".join(_wall(x, y, x, y + h, "Lane wall") for x in edges[1:])


def _row_walls(x, w, y0, n):
    """The same, lying down: one on every cell boundary."""
    return "".join(_wall(x, y0 + i * CELL, x + w, y0 + i * CELL,
                         "Row divider") for i in range(1, n))


def build():
    n_rows = len(ROWS)
    k = _scale()
    widths = [ELL_W if c is None else CELL for c in COLS]
    grid_w = sum(widths)
    body_h = n_rows * CELL

    ax = 0.0
    gx = CELL + GUTTER                        # the grid takes the freed column
    gy = CELL + GUTTER                        # under the header, itself one
    row_y = [gy + i * CELL + CONTENT_Y for i in range(n_rows)]
    edges, x = [], gx                              # column left edges
    for w in widths:
        edges.append(x)
        x += w
    col_x = [e + w / 2 for e, w in zip(edges, widths)]
    cap_dy = SLOT / 2 + FS + 1.5

    o = [_box(ax, gy, CELL, body_h, "Target construct"),
         _box(gx, 0, grid_w, CELL, "Binder series"),
         _box(gx, gy, grid_w, body_h, "Combinations"),
         _walls(edges, 0, CELL),
         _walls(edges, gy, body_h),
         _row_walls(ax, CELL, gy, n_rows),
         _row_walls(gx, grid_w, gy, n_rows)]

    tcx = ax + CELL / 2
    for n, line in enumerate(TITLES["target"]):
        o.append(_text(tcx, gy - 13 + n * LINE, line))

    for i, (label, stem, _pos) in enumerate(ROWS):
        y = row_y[i]
        o.append(_embed(stem, tcx, y, k, f"e{i}") if stem
                 else _non_target(tcx, y, k))
        o.append(_text(tcx, y + cap_dy, label))

    # The header: one candidate cartoon per column, captioned under it exactly
    # as the target rows caption their glyphs, and the ellipsis alone in its
    # own half-square. The cells were empty squares until 2026-11-15, held for
    # these renders.
    for j, col in enumerate(COLS):
        if col is None:
            o.append(_ellipsis(col_x[j], CELL / 2))
        else:
            o.append(_binder(col_x[j], CONTENT_Y, int(col[1])))
            o.append(_text(col_x[j], CONTENT_Y + cap_dy, f"Candidate {col[1]}"))

    # one cell cluster per well
    for i in range(n_rows):
        for j, col in enumerate(COLS):
            if col is not None:
                o.append(_cluster(col_x[j], row_y[i],
                                  f"Well r{i + 1} {col[0]} {col[1]}",
                                  POS_FILL if ROWS[i][2] else NEG_FILL))

    w, h, m = gx + grid_w, gy + body_h, 3.0
    return (f'<svg xmlns="http://www.w3.org/2000/svg" '
            f'xmlns:xlink="http://www.w3.org/1999/xlink" '
            f'width="{(w + 2 * m) * 3:.0f}" height="{(h + 2 * m) * 3:.0f}" '
            f'viewBox="{-m} {-m} {w + 2 * m:.1f} {h + 2 * m:.1f}">'
            + "".join(o) + "</svg>")


def main(outdir):
    outdir.mkdir(parents=True, exist_ok=True)
    p = outdir / "screen_matrix.svg"
    p.write_text(build(), encoding="utf-8")
    print("wrote", p)


if __name__ == "__main__":
    main(Path(sys.argv[1]) if len(sys.argv) > 1 else OUT)
