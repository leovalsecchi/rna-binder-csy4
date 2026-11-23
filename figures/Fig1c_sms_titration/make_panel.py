"""Csy4-H29A-SMS scaffold titration: cartoon, micrograph, mScarlet, per construct.

Three rows on one column grid - the F30 scaffold drawn with its real number of
units, the mScarlet micrograph, and the normalized mScarlet signal.

The scaffold cartoons are `common.scaffold_glyphs.array`, the same glyph fig04
uses for the F30 variation ("n Y units on one strand, two hairpins each"), so
1x / 3x / 8x are n = 1 / 3 / 8 of it rather than new art. Every column is drawn
at one scale: the 1x scaffold is small because it is short. Csy4-H29A-Only has
no scaffold plasmid at all, so its slot is deliberately empty.

Numbers: `inputs/sms_titration_flow.csv`, the 10 ng group of Matteo's
E.MLS11.104.1 (`Csy4mRNAsensingMultipleHairpins`), wells B3-C4. `Norm mScarlet`
is his own normalized column; the construct order is his chart's x axis, and
his FoldChange column (1x at F30-2xNegCtrlLoop, 2.098 / 74.97 / 135.65) is what
ties the two together.

The brackets are the ratio of every column to Neg Ctrl (RATIOS = "control"),
which is how Matteo's own chart reads: 2.1x / 124x / 133x. Set
RATIOS = "adjacent" to bracket each column over the one left of it instead,
which localises the gain (60x between 1x and 3x, 1.1x after).

2026-11-17, Leo: brackets moved from "adjacent" to "control"; the control index
was stale (1, from when Csy4-H29A-Only was column 0) and is now CONTROL_I, read
from COLUMNS. Nested control brackets are stacked one per level instead of
sharing one line. The left spine is cut at the last tick rather than running on
unlabelled to the top of the ylim, the bar baseline is drawn, and the
"Csy4-H29A-SMS, 10 ng" footer is dropped (it is caption material).

Micrographs: `inputs/micrographs/sms_titration/`, cut from his own eight-panel
montage of the same wells, so all five share one display range - see the README.

    python code/Fig1c_sms_titration/make_panel.py

The figure number is a guess, as in fig01b: change `fig_output` if it belongs
elsewhere.
"""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import numpy as np
import pandas as pd
from PIL import Image
import matplotlib.pyplot as plt

from common import style, paths, palette, scaffold_glyphs

style.apply()

TILES = paths.INPUTS / "micrographs" / "sms_titration"
FLOW = paths.INPUTS / "sms_titration_flow.csv"
DOSE = "10 ng"
RATIOS = "control"       # every column over Neg Ctrl, as Matteo's chart does
                         # (Leo, 2026-11-17). "adjacent" brackets each column
                         # over the one left of it, which localises the gain.

# label, construct in the CSV, micrograph, Csy4 sites (vs neutral loops). The
# strip imaged five of his eight conditions; x2, x4 and x5 were not imaged.
# Csy4-H29A-Only (well B3) is measured and in the CSV, but is not shown: it
# and the neutral-loop scaffold differ by 1.1x, so the scaffold-free control
# adds a column without adding a statement.
COLUMNS = [
    ("Neg Ctrl",        "F30-2xNegCtrlLoop",     "B4_negctrl.png", False),
    ("1×",         "F30-2xCsy4Hairpin",     "B5_1x.png",      True),
    ("3×",         "(F30-2xCsy4Hairpin)x3", "C1_3x.png",      True),
    ("8×",         "(F30-2xCsy4Hairpin)x8", "C4_8x.png",      True),
]

PITCH = 30.0                      # column pitch; the panel is as wide as it needs
TILE_W_MM = 28.0                  # micrograph width
BAR_W_MM = 10.0                   # bar width, deliberately narrower than the tile
LABEL_PT, SMALL_PT = 6.0, 5.5     # under the 7 pt base: the panel is 120 mm wide
VALUE_PT = 7.0
# A type variant for page assembly; the canonical panel is untouched:
#   make_panel.py FONT_PT OUTDIR    (every label at FONT_PT, geometry unchanged)
VARIANT = len(sys.argv) == 3
if VARIANT:
    LABEL_PT = SMALL_PT = VALUE_PT = float(sys.argv[1])
CART_MM, BAR_MM = 9.0, 26.0       # cartoon row, bar row
# Clear air between the top of the bar axes and the micrographs. The bracket
# labels are matplotlib text, which does NOT clip to the axes, so the topmost
# one can ride up over the tiles. It does not show at this panel's own size: at
# assembly the panel is drawn at ~74 mm instead of its native 120 while every
# label is re-set to the house 5.5 pt, so the type is ~1.6x taller in data units
# than it is here, and the stack grew into the tiles (Leo, 2026-11-17).
BAR_GAP_MM = 3.0
LABEL_MM, GAP_MM, FOOT_MM = 7.0, 1.6, 9.0

RNA = palette.SYNBIO["rna"]                 # the scaffold strand
SITE = palette.DIVERGING[-1]                # Csy4 hairpin tip, as in fig04
LOOP = palette.NEUTRAL["500"]               # neutral loop, no Csy4 site
BAR = palette.SYNBIO["reporter_alt"]        # mScarlet, the channel imaged
MUT = palette.NEUTRAL["600"]

flow = pd.read_csv(FLOW).query("dose == @DOSE").set_index("construct")
UNITS = [int(flow.loc[c, "f30_units"]) or None for _, c, _, _ in COLUMNS]
VALUES = [float(flow.loc[c, "norm_mscarlet"]) for _, c, _, _ in COLUMNS]
IMAGES = [np.asarray(Image.open(TILES / f).convert("RGB")) for _, _, f, _ in COLUMNS]

shape = {im.shape for im in IMAGES}
assert len(shape) == 1, f"micrographs differ in size: {shape}"
H_PX, W_PX = IMAGES[0].shape[:2]
TILE_H_MM = TILE_W_MM * H_PX / W_PX

W_MM = PITCH * len(COLUMNS)
CX = [PITCH * (i + 0.5) for i in range(len(COLUMNS))]

# ---- figure ------------------------------------------------------------------
FIG_W = W_MM + 16.0                       # room for the y axis on the left
FIG_H = LABEL_MM + CART_MM + GAP_MM + TILE_H_MM + BAR_MM + BAR_GAP_MM + FOOT_MM
fig = plt.figure(figsize=(style.mm(FIG_W), style.mm(FIG_H)))
L = 14.0 / FIG_W                          # left edge of the column grid


def row(y_mm, h_mm):
    ax = fig.add_axes([L, y_mm / FIG_H, W_MM / FIG_W, h_mm / FIG_H])
    ax.set_xlim(0, W_MM)
    return ax


y_bar = FOOT_MM
y_tile = y_bar + BAR_MM + BAR_GAP_MM
y_cart = y_tile + TILE_H_MM + GAP_MM

# micrographs, one axes each so the pitch is exact
for im, cx in zip(IMAGES, CX):
    ax = fig.add_axes([(14.0 + cx - TILE_W_MM / 2) / FIG_W, y_tile / FIG_H,
                       TILE_W_MM / FIG_W, TILE_H_MM / FIG_H])
    ax.imshow(im, interpolation="lanczos", aspect="auto")
    ax.axis("off")

# cartoons: one axes, one scale for every column, so 8x reads as eight of 1x
ax_c = row(y_cart, CART_MM)
ax_c.axis("off")
span = lambda n: 2 * (2.5 * (n - 1) + 3.4)                 # glyph width, units
S = (TILE_W_MM - 1.0) / span(max(n for n in UNITS if n))
ax_c.set_xlim(0, W_MM / S)
ax_c.set_ylim(-8.2 + CART_MM / S / 2, -8.2 - CART_MM / S / 2)   # inverted: up
for (name, _, _, positive), n, cx in zip(COLUMNS, UNITS, CX):
    if n is None:
        continue                                            # no scaffold at all
    u, half = cx / S, span(n) / 2
    scaffold_glyphs.array(ax_c, u, 0.0, RNA, SITE if positive else LOOP,
                          u - half, u + half, n=n, lw=1.3 * S)

# labels, above the cartoons
ax_l = row(y_cart + CART_MM, LABEL_MM)
ax_l.set_ylim(0, 1); ax_l.axis("off")
for (name, _, _, _), cx in zip(COLUMNS, CX):
    ax_l.text(cx, 0.12, name, ha="center", va="bottom",
              fontsize=LABEL_PT, color=palette.INK)

# bars
ax_b = row(y_bar, BAR_MM)
for v, cx in zip(VALUES, CX):
    ax_b.bar(cx, v, width=BAR_W_MM, color=BAR, lw=0)
    ax_b.text(cx, v + max(VALUES) * 0.035, f"{v:.4f}".rstrip("0"),
              ha="center", va="bottom", fontsize=VALUE_PT, color=palette.INK)
# ratio brackets, all on one line above the bars, as on Matteo's own chart
top = max(VALUES)
# The control is found by label, not by position: columns have been added and
# removed here, and a hard-coded index silently brackets the wrong one.
CONTROL_I = next(i for i, c in enumerate(COLUMNS) if c[0] == "Neg Ctrl")
pairs = ([(i, i + 1) for i in range(len(CX) - 1)] if RATIOS == "adjacent"
         else [(CONTROL_I, j) for j in range(len(CX)) if j != CONTROL_I])
tick = top * 0.030
# "control" brackets land on the COLUMN CENTRES (Leo, 2026-11-17): they all
# start on Neg and nest, so each one points at the column it is the ratio of.
# "adjacent" brackets run end to end instead, and there the inset is what
# stops two neighbours reading as one long rule.
inset = 0.0 if RATIOS == "control" else PITCH * 0.10
# "adjacent" brackets are disjoint and share one line; "control" brackets all
# start at the same column and nest, so each gets its own level.
# 0.40, not 0.17: the step has to clear a LABEL, and at assembly the panel is
# drawn at ~0.62 of its native width with the type re-set to the house size, so
# a label is ~0.07 of the ylim there against ~0.04 here. Sized for the assembled
# case; it only looks airy in this panel's own PNG.
STEP = top * 0.40 if RATIOS == "control" else 0.0
y_base = top * (1.22 if RATIOS == "control" else 1.26)
for level, (i, j) in enumerate(pairs):
    r = VALUES[j] / VALUES[i]
    y_br = y_base + level * STEP
    x0, x1 = CX[i] + inset, CX[j] - inset
    ax_b.plot([x0, x0, x1, x1], [y_br - tick, y_br, y_br, y_br - tick],
              color=MUT, lw=0.8, solid_capstyle="butt")
    ax_b.text((x0 + x1) / 2, y_br + tick * 0.4,
              f"{r:.0f}\u00d7" if r >= 10 else f"{r:.2g}\u00d7",
              ha="center", va="bottom", fontsize=SMALL_PT, color=palette.INK)

ax_b.set_ylim(0, y_base + max(len(pairs) - 1, 0) * STEP + top * 0.45)
# Ticked up to 0.4, not stopped at 0.2 (Leo, 2026-11-17). The ylim carries the
# bracket stack and reaches ~0.59, so three ticks left two thirds of the spine
# bare; the spine is still cut at the LAST tick rather than running to the ylim.
YTICKS = [0, 0.1, 0.2, 0.3, 0.4]
ax_b.set_yticks(YTICKS)
ax_b.set_xticks([])
# The ylim carries the bracket stack, so the spine would run on unlabelled well
# above the last tick; bound it to the ticks and let the brackets sit clear.
ax_b.spines["left"].set_bounds(YTICKS[0], YTICKS[-1])
ax_b.spines["bottom"].set_visible(True)
ax_b.spines["bottom"].set_bounds(0, W_MM)
ax_b.tick_params(axis="x", length=0)
ax_b.set_ylabel("mScarlet, norm.", fontsize=LABEL_PT)
ax_b.tick_params(axis="y", labelsize=SMALL_PT)

out = Path(sys.argv[2]) if VARIANT else paths.fig_output("Fig1c_sms_titration")
style.save(fig, out, "sms_titration")
fig.savefig(out / "sms_titration.svg", bbox_inches="tight")
print(f"tiles {W_PX}x{H_PX} px -> {TILE_W_MM} x {TILE_H_MM:.1f} mm, "
      f"pitch {PITCH:.1f} mm, glyph scale {S:.3f} mm/unit")
print("wrote", out / "sms_titration.png")
