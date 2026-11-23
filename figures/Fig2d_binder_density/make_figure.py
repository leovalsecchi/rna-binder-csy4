"""Fig 2 bottom row: binder density on the target, one column per backbone generator.

    python code/Fig2d_binder_density/make_figure.py [outdir] [inputs]

Four columns side by side, each labelled with the backbone generator that made
the designs, in that stack's colour. The inverse folder moves no atoms, so the
density is a property of the backbone generator, not of the stack: BoltzGen+BIF
is not a fifth column, it shares BoltzGen's backbones.

Top of each column: 800 design backbones, placed in the target's frame and
contoured at the level where 30% of that generator's binders reach (see
analyses/binder_density on Euler). The camera is identical in all four - the
hairpin stands on its own helical axis, truncated ends at the bottom, leaning
25 degrees toward the viewer - so the grey RNA is pixel-identical across the
columns and only the density differs.

Bottom of each column: how many heavy-atom contacts that generator's binders
make with each of the 16 nucleotides of the repeat hairpin - protein-RNA heavy
atom pairs within 4.5 A, averaged over the designs, with the GUAUA loop
(nt7-11) as a pale band behind the bars. All four share one y scale, so bar
heights compare directly, and the nucleotide scale is drawn once, on RFD3.

The renders are design backbones; a contact needs side chains, and three of the
four stacks write backbone-only design files. Those designs are read from the
LigandMPNN side-chain retro-pack instead - same backbone, same sequence - so
every structure counted here is full-atom and the four are comparable. Measured
on the design's OWN pose, not on a refold, which is what the render shows.

Emitted at 1 unit = 1 pt at the 178 mm block, so the assembly places it at
s ~ 1 and `retext` lands the labels on the house size with no scaling surprise.

Two things the page assembly needs that are easy to miss:
  * it floors every stroke and asserts the thinnest is >= 0.5 pt, taking a `min`
    over the strokes it finds - a panel with NO stroke-width at all raises
    ValueError on an empty min(), so the labels carry an explicit (unstroked)
    stroke-width.
  * images are embedded as base64, not linked: the panel's body is inlined into
    the page, where a relative href would no longer resolve.
"""
import base64
import csv
import os
import sys
from collections import OrderedDict
from pathlib import Path

# In the repo this resolves to the repo root. LV_REPO overrides it, which is how
# this was developed while the NAS share was mounted read-only.
REPO = Path(os.environ.get("LV_REPO") or Path(__file__).resolve().parents[2])
sys.path.insert(0, str(REPO / "code"))
from common import csy4, palette, style                       # noqa: E402

MM = 72 / 25.4
W = 210 * MM - 2 * 16 * MM          # 504.567 pt, the house block
T = float(style.CFG["font_size_pt"])
ARIAL = "Arial, Helvetica, sans-serif"

# A render is taller than it is wide (1200 x 1499), so four at full width
# overflow the A4 page: the image width is fixed and the slack goes into the
# gaps. The distribution under it is wider than the render, which is what keeps
# the row looking like one block rather than four stamps - the cell is the wider
# of the two, and the render is centred in it.
IMG_W = float(os.environ.get("IMG_W", 57.0))     # pt, each render
# 57 is the largest that still fits. The page has no slack, so every line of
# text under the row is paid for here: 78 -> 60 bought the contact strip, and
# 60 -> 57 bought the y-quantity line under the sequence.
PLOT_W = float(os.environ.get("PLOT_W", 90.0))   # pt, each distribution
PLOT_H = float(os.environ.get("PLOT_H", 17.0))   # pt, its height
LABEL_GAP = 6.0                     # pt, image bottom -> label baseline
PLOT_GAP = 4.0                      # pt, label baseline -> plot top
AXIS_GAP = 2.5                      # pt, plot baseline -> sequence
YLABEL_GAP = 1.2                    # pt, extra leading under the sequence
TICK = 1.5                          # pt, tick length below the baseline
NT = 16                             # nucleotides of the Csy4 repeat core
BAR_GAP = 0.18                      # share of a slot left as air between bars

# The renders sit on a tinted rounded panel, the same one Fig 2a and Fig 1d use:
# neutral.100 with rx 14 on a 1520-unit viewBox, i.e. 14 * W / 1520 printed, so
# the two boxes on this page carry the same corner radius.
BOX_FILL = palette.NEUTRAL["100"]
BOX_RX = 14.0 * W / 1520.0
PAD_X, PAD_Y = 12.0, 6.0            # pt, inside the box
RULE = palette.NEUTRAL["400"]       # hairline axes
BAND = palette.NEUTRAL["200"]       # the pale loop band behind the bars
MUTED = palette.NEUTRAL["600"]      # the one axis label and the stack mark
LINE = 0.5                          # pt, the assembly's stroke floor

# label -> (file stem, colour name in palette.json "stacks")
PANELS = [
    ("RFD3", "RFD3", "RFD3 + LigandMPNN"),
    ("BoltzGen + LigandMPNN", "BoltzGen", "BoltzGen + LigandMPNN"),
    ("BoltzDesign", "BoltzDesign", "BoltzDesign + LigandMPNN"),
    ("HalluDesign", "HalluDesign", "HalluDesign + LigandMPNN"),
]
# BoltzGen's label carries its stack name because it is the one generator whose
# backbones feed TWO stacks (+LigandMPNN and +BIF): the render is the shared
# backbone, but a contact count belongs to a sequence, so it has to say which
# stack it was measured on. The other three have one stack each, so naming it
# there would be ink for nothing. All four groups ARE the +LigandMPNN stack -
# that belongs in the caption.
#
# What the bar heights are, on its own line UNDER the first column's sequence.
# It was briefly inside the plot at top left, which clears RFD3's bars but sits
# on the pale loop band and reads as a clash (2026-11-19, Leo). Under the
# letters it touches nothing; it costs one line of page height, which the
# renders pay for. The caption still carries the full definition and n.
Y_LABEL = "mean contacts"
Y_LABEL_ON = "RFD3"
#
# The x axis is the sequence itself, a letter under each bar, on every column:
# the nucleotide a bar belongs to is the whole read of this strip, and counting
# across to a scale drawn once on RFD3 is exactly what a reader should not have
# to do. It replaces the 1/16 ticks and the "contacts per nt" label, so the y
# quantity now lives only in the caption.


def png_size(path):
    """(width, height) from the IHDR, so PIL is not needed at assembly time."""
    b = path.read_bytes()[16:24]
    return int.from_bytes(b[:4], "big"), int.from_bytes(b[4:], "big")


def contacts(inputs):
    """{group: (n, [mean contacts per nt, 1..16])} from contacts_table.py's CSV."""
    p = inputs / "nt_contacts_by_group.csv"
    if not p.exists():
        raise SystemExit(f"missing {p} - run contacts_table.py first")
    out = OrderedDict()
    for r in csv.DictReader(p.open()):
        n, vals = out.setdefault(r["group"], (int(r["n"]), []))
        vals.append(float(r["mean_contacts"]))
    for g, (n, vals) in out.items():
        if len(vals) != NT:
            raise SystemExit(f"{g}: {len(vals)} nucleotides, expected {NT}")
    return out


def per_nt_bars(x0, ybase, vals, ymax, colour, width):
    """Mean contacts on each of the 16 nucleotides, on a hairline baseline.

    The GUAUA loop (nt7-11, `common.csy4`) is a pale band behind the bars: the
    whole read of this strip is whether a generator works the loop or the stem,
    and a band costs no height, which a bracket under the axis would."""
    slot = width / NT
    lo, hi = min(csy4.LOOP) - 1, max(csy4.LOOP)       # 0-based, half-open
    o = [f'<rect x="{x0 + lo * slot:.2f}" y="{ybase - PLOT_H:.2f}" '
         f'width="{(hi - lo) * slot:.2f}" height="{PLOT_H:.2f}" fill="{BAND}"/>']
    for i, v in enumerate(vals):
        h = v / ymax * PLOT_H
        o.append(f'<rect x="{x0 + (i + BAR_GAP / 2) * slot:.2f}" y="{ybase - h:.2f}" '
                 f'width="{slot * (1 - BAR_GAP):.2f}" height="{h:.2f}" fill="{colour}"/>')
    o.append(f'<line x1="{x0:.2f}" y1="{ybase:.2f}" x2="{x0 + width:.2f}" '
             f'y2="{ybase:.2f}" stroke="{RULE}" stroke-width="{LINE:.2f}"/>')
    return o


def build(inputs):
    dists = contacts(inputs)
    ymax = max(v for n, vals in dists.values() for v in vals)
    n = len(PANELS)
    cell = max(IMG_W, PLOT_W)
    inner = W - 2 * PAD_X
    gap = (inner - n * cell) / (n - 1)
    if gap < 0:
        raise SystemExit(f"cell {cell} too wide for the {inner:.1f} pt inside the box")

    imgs = []
    for label, stem, _ in PANELS:
        p = inputs / f"density_{stem}.png"
        if not p.exists():
            raise SystemExit(f"missing {p}")
        iw, ih = png_size(p)
        imgs.append((label, p, IMG_W, IMG_W * ih / iw))
    tall = max(h for _, _, _, h in imgs)

    y_label = PAD_Y + tall + LABEL_GAP + T * 0.716      # label baseline
    y_plot = PAD_Y + tall + LABEL_GAP + T + PLOT_GAP    # plot top
    y_base = y_plot + PLOT_H
    y_axis = y_base + TICK + AXIS_GAP + T * 0.716       # sequence baseline
    y_ylabel = y_axis + T + YLABEL_GAP                  # the y quantity, below it
    height = y_ylabel + T * 0.284 + PAD_Y

    out = [f'<svg xmlns="http://www.w3.org/2000/svg" '
           f'xmlns:xlink="http://www.w3.org/1999/xlink" '
           f'viewBox="0 0 {W:.3f} {height:.3f}" width="{W:.3f}" height="{height:.3f}">',
           f'<rect x="0" y="0" width="{W:.3f}" height="{height:.3f}" '
           f'rx="{BOX_RX:.3f}" fill="{BOX_FILL}"/>']
    for i, (label, p, w, h) in enumerate(imgs):
        cx = PAD_X + i * (cell + gap)
        mid = cx + cell / 2
        colour = palette.stack_color(PANELS[i][2])

        x = cx + (cell - w) / 2
        y = PAD_Y + tall - h              # bottom-aligned, so the hairpins line up
        href = "data:image/png;base64," + base64.b64encode(p.read_bytes()).decode()
        out.append(f'<image x="{x:.3f}" y="{y:.3f}" width="{w:.3f}" height="{h:.3f}" '
                   f'preserveAspectRatio="xMidYMax meet" xlink:href="{href}"/>')
        out.append(f'<text x="{mid:.3f}" y="{y_label:.3f}" '
                   f'font-family="{ARIAL}" font-size="{T:.4f}" text-anchor="middle" '
                   f'fill="{colour}" stroke="none" stroke-width="1">{label}</text>')

        px0 = cx + (cell - PLOT_W) / 2
        _, vals = dists[PANELS[i][1]]   # keyed by the group, not the drawn label
        out += per_nt_bars(px0, y_base, vals, ymax, colour, PLOT_W)

        slot = PLOT_W / NT
        for k, base in enumerate(csy4.SEQUENCE):
            out.append(f'<text x="{px0 + (k + 0.5) * slot:.3f}" y="{y_axis:.3f}" '
                       f'font-family="{ARIAL}" font-size="{T:.4f}" '
                       f'text-anchor="middle" fill="{MUTED}" stroke="none" '
                       f'stroke-width="1">{base}</text>')

        if PANELS[i][1] == Y_LABEL_ON:
            out.append(f'<text x="{px0:.3f}" y="{y_ylabel:.3f}" '
                       f'font-family="{ARIAL}" font-size="{T:.4f}" fill="{MUTED}" '
                       f'stroke="none" stroke-width="1">{Y_LABEL}</text>')
    out.append("</svg>")
    return "".join(out)


def main():
    outdir = Path(sys.argv[1]) if len(sys.argv) > 1 else \
        REPO / "output" / "Fig2d_binder_density"
    inputs = Path(sys.argv[2]) if len(sys.argv) > 2 else REPO / "inputs" / "binder_density"
    outdir.mkdir(parents=True, exist_ok=True)
    svg = build(inputs)
    p = outdir / "binder_density_row.svg"
    p.write_text(svg, encoding="utf-8")
    # the longest label has to fit its cell, or two columns run together
    wide = max(PANELS, key=lambda t: style.text_w(t[0], T))
    w = style.text_w(wide[0], T)
    cell = max(IMG_W, PLOT_W)
    print(f"wrote {p}  ({len(svg) / 1e6:.2f} MB, {W / MM:.0f} mm wide), "
          f"widest label {wide[0]!r} {w:.1f} pt in a {cell:.0f} pt cell")
    if w > cell:
        print(f"  WARNING: that label is wider than the cell")


if __name__ == "__main__":
    main()
