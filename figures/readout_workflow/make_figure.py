"""The readout workflow: RNA out of the cells, through to filtered reads.

An L. Four steps climb the left column from the bench end, and at the top the
flow turns right into sequencing and analysis. Drawn in the idiom
Fig1a_screen_pipeline established - line glyphs in the thesis petrol, a short
label per step, thin grey connectors, no boxes. A filled box adds an outline
and a fill that say nothing the icon and its label do not, and a column of
them joined by arrows is a slide, not a figure panel.

**Every step is one icon with its label to the RIGHT of it**, in the column
and along the row alike. That is not a style choice, it is what makes the L
possible: with labels UNDER the icons, the upward arrow from one step to the
next runs straight through the label of the step it is pointing at. Putting
the label beside the icon leaves the icon's own centre line free, so the
connectors have a clear lane in both directions.

Lines are finer than Fig1a's: the glyphs carry SW_ICON rather than the 2.4
default, which is why `icon_group` grew a stroke-width argument instead of
this figure redrawing anything.

**No glyph here is used by another panel.** `tube` and `funnel` would have
done for the first and last steps, but both already mean something else in
Fig1a, so `cell-lysis` and `read-filter` were drawn for these. Shared glyphs
live in code/common/glyphs.py, never in a figure script.

Sizes are in POINTS, converted by U. **U is fixed and the canvas is derived
from the layout**, not the other way round, so the row can be lengthened
without any printed size moving. Its arrows are long on purpose: in Fig 4a
the last two steps carry across the hit plot beneath them.

No panel letter is drawn; letters go on at assembly (STYLE_NOTES section 11).

    python code/readout_workflow/make_figure.py [--check]
"""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import yaml
from common import paths, palette, glyphs, style

out = paths.fig_output("readout_workflow")

_cfg = yaml.safe_load((paths.ROOT / "config" / "style.yaml").read_text())
FONT = f'{_cfg["font_family"]}, Helvetica, sans-serif'
# NOT the house size. This workflow is a component: Fig 4a scales it by
# 5.5 / (its label size in units), so 7 here is a design ratio - labels
# against glyphs - not a print size. Following the house size would grow
# every glyph inside Fig 4a. Print size is set where it is placed.
LABEL_PT = 7.0

ACC = palette.MAIN              # the one accent; COLORS.md rule 1
INK = palette.INK
CONN = palette.NEUTRAL["400"]   # connectors, which are furniture not data
BG = palette.SURFACE

# **U is the primitive here, not the panel width.** It used to be the other
# way round - a fixed 85 mm canvas with U derived from it - which meant the
# row could not be lengthened without either running off the canvas or, by
# widening it, changing U and so every printed size in the drawing at once.
# With U fixed, the canvas is whatever the layout turns out to need and
# nothing printed moves. 1.95 is the old value to five figures, so this swap
# left the drawing identical.
U = 1.95                        # viewBox units per point
LINE_H = 9.5 * U
MARGIN = 16                     # canvas edge to the drawing, in units
# W and WIDTH_MM are set by build() from the layout; an L leaves its lower
# right empty by nature, but dead margin on the right as well would read as a
# mistake rather than as room the figure has left for something.
W = WIDTH_MM = None

ICON_K = 1.0                    # 48 units, ~10.8 mm printed
ICON_R = 24 * ICON_K            # its half-size, so connectors can clear it
SW_ICON = 1.70                  # finer than the 2.4 default
SW_CONN = 0.85 * U

COL_X = 74                      # the climbing column
COL_DY = 101                    # between steps up the column, set from the
                                # pitch Leo spaced the column at by hand
COL_TOP = 46
ARROW_LEN = 90                  # the gap a connector spans. Long, by
                                # request: in Fig 4a the last two steps have
                                # to carry across the plot beneath them, and
                                # the row is the only part of the L that can
                                # stretch without moving anything else.
LAB_DX = 16                     # container edge to the start of its label
LAB_DY = 15                     # container edge down to a label underneath

# A pale container behind each icon, which is what gives the steps a shared
# footprint instead of six free-floating drawings (STYLE_NOTES section 7).
# Theme 100-300 are tints for fills and backgrounds only (COLORS.md rule 1),
# so the tile can echo the accent without introducing a colour.
TILE = "circle"                 # square | circle | card | none
TILE_PAD = 13                   # icon edge to the tile edge
TILE_FILL = palette.THEME["100"]
TILE_EDGE = palette.THEME["200"]
# What the labels and connectors measure from: the container if there is one,
# the icon if there is not. Spacing the column on the ICON while a tile sits
# round it leaves the tiles all but touching and the arrows invisible.
EDGE = ICON_R + TILE_PAD if TILE != "none" else ICON_R

# (glyph, label lines) in the order the bench does them. The fourth is the
# corner: it ends the column and begins the row.
STEPS = [
    ("cell-lysis",  ["RNA extraction"]),
    ("rt-cdna",     ["Reverse", "transcription"]),
    ("pcr",         ["PCR amplification"]),
    ("barcode",     ["Barcoding of", "target conditions"]),
    ("nanopore",    ["ONT sequencing"]),
    ("read-filter", ["Min. reads", "filtering"]),
]
CORNER = 3                      # index of the step the L turns at


def label_w(lines):
    """The widest line, MEASURED. See common.style.text_w for why not
    estimated: a per-character rule leaves a different amount of invisible
    slack in every gap, and the row's spacing then reads as uneven."""
    return max(style.text_w(ln, LABEL_PT) for ln in lines) * U


def beside(i):
    """True if step i wears its label at its side rather than underneath.

    Column steps carry it beside, because an upward arrow would otherwise run
    through the label of the step it points at. Row steps carry it below,
    which is what keeps the row narrow - a label beside each one pushes the
    next step a label's width further right and the panel turns into a strip.
    The corner is a column step that the row leaves from, so it keeps beside.
    """
    return i <= CORNER


def layout():
    """Every step's centre. The row is derived, so the arrows are all equal."""
    pos = [(COL_X, COL_TOP + (CORNER - i) * COL_DY) for i in range(CORNER + 1)]
    x = COL_X + EDGE + LAB_DX + label_w(STEPS[CORNER][1])
    for i in range(CORNER + 1, len(STEPS)):
        x += ARROW_LEN + EDGE
        pos.append((x, COL_TOP))
        x += EDGE
    return pos


POS = None          # filled by build(), since it depends on the label widths


def place(i):
    return POS[i]


def tile(cx, cy, lines):
    """The container behind a step. One footprint for every icon, so the
    column reads as a column and not as six drawings that line up."""
    if TILE == "none":
        return ""
    r = EDGE
    if TILE == "circle":
        return (f'<circle cx="{cx:.1f}" cy="{cy:.1f}" r="{r:.1f}" '
                f'fill="{TILE_FILL}" stroke="{TILE_EDGE}" stroke-width="0.8"/>')
    w = 2 * r
    if TILE == "card":
        w = EDGE + LAB_DX + label_w(lines) + TILE_PAD + r
    return (f'<rect x="{cx-r:.1f}" y="{cy-r:.1f}" width="{w:.1f}" '
            f'height="{2*r:.1f}" rx="{9 if TILE == "card" else 11}" '
            f'fill="{TILE_FILL}" stroke="{TILE_EDGE}" stroke-width="0.8"/>')


def T(x, y, s, pt=LABEL_PT, col=INK, anch="start"):
    return (f'<text x="{x:.1f}" y="{y:.1f}" font-family="{FONT}" '
            f'font-size="{pt*U:.2f}" fill="{col}" text-anchor="{anch}">{s}</text>')


def label(cx, cy, lines, side):
    """Beside the icon (column) or centred under it (row)."""
    sz = LABEL_PT * U
    if side:
        y0 = cy - (len(lines) - 1) * LINE_H / 2 + sz * 0.35
        return "".join(T(cx + EDGE + LAB_DX, y0 + j * LINE_H, ln)
                       for j, ln in enumerate(lines))
    y0 = cy + EDGE + LAB_DY
    return "".join(T(cx, y0 + j * LINE_H, ln, anch="middle")
                   for j, ln in enumerate(lines))


def arrow(x1, y1, x2, y2, head=5.0):
    """A straight connector with a solid head at (x2, y2)."""
    import math
    a = math.atan2(y2 - y1, x2 - x1)
    c, s = math.cos(a), math.sin(a)
    tip = (x2, y2)
    back = (x2 - head * c, y2 - head * s)
    wing = head * 0.46
    p1 = (back[0] - wing * s, back[1] + wing * c)
    p2 = (back[0] + wing * s, back[1] - wing * c)
    return (f'<path d="M{x1:.1f},{y1:.1f} L{back[0]:.1f},{back[1]:.1f}" '
            f'fill="none" stroke="{CONN}" stroke-width="{SW_CONN:.2f}" '
            f'stroke-linecap="round"/>'
            f'<path d="M{tip[0]:.1f},{tip[1]:.1f} L{p1[0]:.1f},{p1[1]:.1f} '
            f'L{p2[0]:.1f},{p2[1]:.1f} Z" fill="{CONN}"/>')


def build():
    global POS, H, W, WIDTH_MM
    POS = layout()
    deepest = max(cy for _, cy in POS) + EDGE
    tail = max(POS[i][1] + EDGE + LAB_DY + len(STEPS[i][1]) * LINE_H
               for i in range(len(STEPS)) if not beside(i))
    H = int(max(deepest, tail) + MARGIN)
    # the widest thing on the right is a container, or a beside-label, or a
    # centred row label overhanging its own container
    right = max(max(cx + EDGE for cx, _ in POS),
                max(POS[i][0] + EDGE + LAB_DX + label_w(STEPS[i][1])
                    if beside(i) else POS[i][0] + label_w(STEPS[i][1]) / 2
                    for i in range(len(STEPS))))
    W = int(right + MARGIN)
    WIDTH_MM = W * 0.352778 / U

    o = [f'<svg xmlns="http://www.w3.org/2000/svg" width="{WIDTH_MM:.2f}mm" '
         f'height="{WIDTH_MM*H/W:.2f}mm" viewBox="0 0 {W} {H}">',
         f'<rect width="{W}" height="{H}" fill="{BG}"/>']

    for i, (glyph, lines) in enumerate(STEPS):
        cx, cy = place(i)
        o.append(tile(cx, cy, lines))
        o.append(glyphs.icon_group(glyph, cx, cy, ACC, ICON_K, SW_ICON))
        o.append(label(cx, cy, lines, beside(i)))

        if i == len(STEPS) - 1:
            continue
        nx, ny = place(i + 1)
        if i < CORNER:                      # climbing: up the icon's own axis
            o.append(arrow(cx, cy - EDGE - 4, cx, ny + EDGE + 4))
        else:                               # along the row
            x1 = (cx + EDGE + LAB_DX + label_w(lines) if beside(i)
                  else cx + EDGE) + 6
            o.append(arrow(x1, cy, nx - EDGE - 6, cy))

    o.append('</svg>')
    return "\n".join(o)


def _check():
    """House label size, no glyph invented here, nothing overlapping."""
    svg = build()                       # also fills POS and H
    assert f'font-size="{LABEL_PT*U:.2f}"' in svg, "label size is not font_size_pt"
    assert abs(WIDTH_MM * U / 0.352778 - W) < 1.0, "W and WIDTH_MM disagree"
    missing = [g for g, _ in STEPS if g not in glyphs.ICONS]
    assert not missing, f"no such glyph: {missing}"

    spans = []
    for i, (_, lines) in enumerate(STEPS):
        cx, cy = place(i)
        if beside(i):
            x0 = cx + EDGE + LAB_DX
            x1 = x0 + label_w(lines)
            y1 = cy + (len(lines) - 1) * LINE_H / 2
        else:
            x0, x1 = cx - label_w(lines) / 2, cx + label_w(lines) / 2
            y1 = cy + EDGE + LAB_DY + (len(lines) - 1) * LINE_H
        spans.append((x0, x1))
        assert 0 < x0 and x1 < W, f"label {i} runs off the canvas"
        assert cy - EDGE > 0 and y1 < H, f"step {i} runs off the canvas"
        assert cx + EDGE < W, f"step {i} container runs off the canvas"

    # along the row, neither the labels nor the containers may collide
    row = [i for i in range(len(STEPS)) if not beside(i)]
    for a, b in zip(row, row[1:]):
        assert spans[a][1] < spans[b][0], f"labels {a} and {b} overlap"
        assert place(a)[0] + EDGE < place(b)[0] - EDGE, f"steps {a},{b} touch"
    print(f"checks pass: {len(STEPS)} steps, labels {LABEL_PT} pt, "
          f"{WIDTH_MM:.0f} x {WIDTH_MM*H/W:.1f} mm, aspect {W/H:.2f}")


if __name__ == "__main__":
    if "--check" in sys.argv:
        _check()
        sys.exit(0)
    for a in sys.argv[1:]:
        if a.startswith("--tile="):
            globals()["TILE"] = a.split("=", 1)[1]
    svg_path = out / "readout_workflow.svg"
    svg_path.write_text(build(), encoding="utf-8")
    glyphs.write_icons(out / "icons")
    print("wrote", svg_path)
