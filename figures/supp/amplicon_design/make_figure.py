"""Amplicon design and the primer test that chose it.

Two bare panels. There is no derived table: the lengths come from the PCR
worksheet, the module spans from the construct files, the gel from an image.
All three are named in the README.

Panel (a) draws the transcript with its three modules - recorder, sensor,
target - and BOTH candidate designs under it on the same base-pair axis, so a
reader sees what each amplicon covers. That is the whole reason one design is
~850 bp and the other ~2,650: the short one stops inside the sensor, the long
ones had to run through the sensor and the target.

Panel (b) is the LV01.04 gel with the ladder turned into a size axis and the
expected product positions marked, so the panel is evidence and not a photo.

Module spans are read off the files, never sketched. Their provenance, and the
one boundary that comes from the part rather than from this template, is in the
README under "Where each module span comes from".

**The glyphs are main figure 1's, not SBOL bars** (changed 2026-11-19, Leo:
"take the elements from main figure 1"). The construct in Fig 1 is
`working/Fig1d_integrated_construct/Fig1d_integrated_construct.svg` - it
reaches the page as panel 1e, see `code/Fig1_assembly/make_figure.py` - and it
draws modules as ROUNDED RECTANGLES with brackets above and label ticks below,
not as the CDS pentagons `common/constructs.py` builds. That module's bars are
on Leo's excluded list (2026-11-12: "use Matteo-derived elements only, never
the ones we drew from scratch"), so every proportion and colour below is
carried from the Fig 1 file instead and `_selfcheck` asserts the fills still
agree with `constructs.py`, which is what stops the two drifting apart.
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

import numpy as np
import matplotlib.image as mpimg
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, Rectangle

from common import paths, style, palette
from common import constructs as K

FIG = Path(__file__).resolve().parent.name
GELS = Path("/Volumes/wicky/users/lvalsecchi/Leo_Valsecchi/LV_Gels")

INK = palette.INK
MUTED = palette.NEUTRAL["600"]
FS = 5.5

# LV_PCR/LV1.4.1.xlsx, sheet 1: one forward primer, three reverse primers,
# three expected sizes, three lanes. Journal LV01.04 (14.07.26) is the same
# experiment and names the same pairs.
FWD = "PR10148"

# PR250 anneals 20 bp 3' of PR249. Read off
# LV_Cloning/pLV_NNN/pLV001/pLV001_part158.gb, where PR249_BindingSite is
# 3253..3272 and PR250_BindingSite is 3273..3292 - the two sites ABUT, so the
# offset is exactly 20 bp and PR250 gives the longer product.
#
# That file is a library DONOR, not this template: its AT cassette is 117 bp
# against this template's 159 and its binder slot holds the ccdB placeholder,
# so the donor's absolute coordinates imply ~3,360 bp where the worksheet says
# ~2,650 (README, "Still open"). The ORDER carries across and the POSITIONS do
# not, so the offset is drawn and both products keep the worksheet's ~2,650.
REV_OFFSET = {"PR249": 0, "PR250": 20}
CANDIDATES = [
    #  reverse    expected  lane
    ("PR10137", 850, 3),
    ("PR249", 2650, 5),
    ("PR250", 2650, 7),
]

# --- the transcript, in amplicon coordinates ------------------------------
# Base 1 is the 5' end of PR10148: the first base of the LV1.27.1 ONT
# reference (LV_Sequencing/LV1.27.1_ONT/MLS12_18_1_AllSamples.fa, 865 bp),
# which LV1.27.xlsx says was made with PR10148 + PR10137. Aligning that
# reference to pLV018 gives every boundary below; see the README.
#   1..62     Twist 5' UTR          (== pLV018 1180..1241)
#   63..221   AT splicing cassette  (159 bp absent from pLV018, 3' 30mer in pLV001)
#   223..232  Kozak, ATG at 229     (== pLV018 1241..1250)
#   241..801  Csy4-H29A             (== pLV018 1259..1819)
#   808..843  (GGS-GGT)x2 linker    (== pLV018 1826..1861)
#   847..2328 MgIRE1                (== pLV018 1865..3346)
#   2329..2580 target module        (see TARGET_END)
#   2581..2650 3' UTR + the reverse primer site: BARE STRAND, as Fig 1 draws
#             its 3' lead. This is the stretch the reverse primers anneal in.
#
# Fig 1's own names are kept for the blocks: the binder SLOT is what the map
# draws, and on MLS12.18.1.1 that slot holds Csy4-H29A (said in the caption,
# not on the canvas). "Binder" and "SMS" are therefore the same objects here
# as in Fig 1, which is the point of the revision.
# Where the target module STOPS, and so where the bare 3' strand the reverse
# primers anneal in begins.
#
# Leo, 2026-11-19, zooming in: *"The reverse primers should be after so more
# right than the target module."* He is right, and the block was the fault:
# it ran to 2,650, which is the PRODUCT's end, so it swallowed the 3' UTR and
# the primer binding site and no reverse primer mark could sit clear of it.
# (The two lanes were already on one x axis - both go through bp() with the
# same x0 and scale - so nothing about the axis needed changing.)
#
# Derived, not sketched. The product ends at the reverse primer's 5' end, so
# PR249's 20 bp footprint is 2631..2650. pLV001_part158.gb puts 50 bp between
# the F30 array's 3' end (3202) and PR249's site (3253): Twist 3' UTR 40,
# X6-GTCA overhang 4, 3' UTR 6. So the array ends 50 + 20 bp before 2,650.
#
# ⚠️ That 50 bp is the DONOR's 3' architecture. The README already carried a
# `VERIFY:` on the target's extent and it still stands - pLV018's own Twist
# 3' UTR is 64 bp, not 50, so the two constructs' 3' ends are not identical.
# This is a better-founded estimate than "the remainder up to the primers",
# which is what it replaces, and it is the one that puts the primers where the
# files say they anneal: 3' of the target module.
REV_FOOTPRINT = 20               # PR249_BindingSite 3253..3272, PR250 3273..3292
TARGET_3UTR = 50                 # 3203..3252, between the array and that site
TARGET_END = 2650 - REV_FOOTPRINT - TARGET_3UTR

BLOCKS = [
    # key         start        end   fill                      label
    ("utr",           1,        62, palette.NEUTRAL["200"],   None),
    ("at",           63,       221, palette.SYNBIO["intron"], None),
    ("binder",      241,       801, palette.SYNBIO["binder"], "Binder"),
    ("linker",      808,       843, palette.NEUTRAL["300"],   None),
    ("sms",         847,      2328, palette.SYNBIO["enzyme"], "SMS"),
    ("f30",        2329, TARGET_END, palette.THEME["200"],    None),
]
# bracket spans and their captions, Fig 1's three module names
BRACKETS = [(63, 221, "Recorder"), (229, 2328, "Sensor"),
            (2329, TARGET_END, "Target array")]

# The RNase lobe's share of the SMS block, carried verbatim from the Fig 1
# file: its SMS outline runs x 69.2315..120.6 and the RNase path starts at
# 100.696, so the lobe is 19.904 of 51.3685.
#
# ⚠️ This is Fig 1's SCHEMATIC proportion, NOT a sequence coordinate. Every
# other boundary on this panel is base pairs off the alignment; this one is
# not, because no kinase/RNase boundary for MgIRE1 was read off this template.
# Said in the README and in the caption.
RNASE_FRAC = (120.6 - 100.696) / (120.6 - 69.2315)

# --- Fig 1's construct geometry, in its own units -------------------------
# One unit is one point: the file sets captions at 5.5 and blocks at 6, the
# BLOCK_PER_TEXT = 6/5.5 ratio Fig1d/make_figure.py asserts at build time. All
# offsets below are measured off Fig1d_integrated_construct.svg, whose strand
# sits at y 30.28, and are quoted as units from that strand.
PT = 25.4 / 72.0                 # mm per point
BLOCK_H = 6.0                    # 27.28..33.28
BLOCK_R = 0.8                    # the corner arcs in every block path
LINK_H = 3.0                     # the linker is drawn THIN: 28.78..31.78
BRK_TOP, BRK_ARM = -10.0, -7.4   # bracket spine 20.28, arm feet 22.88
BRK_CAP = -12.78                 # its caption's baseline, 17.5
BRK_INSET = 0.35                 # mm each side. Fig 1 insets its brackets
                                 # from the blocks (sensor 42.2..120.4 over
                                 # blocks 42..120.6); here it also has work to
                                 # do, because the recorder ends at 221 and the
                                 # sensor starts at 229 - 0.18 mm apart at this
                                 # scale, so without it the two feet merge into
                                 # one rule. _selfcheck measures the gap.
TICK_TOP, TICK_BOT = 4.0, 6.1    # label ticks 34.28..36.38
TICK_CAP = 11.92                 # their captions' baseline, 42.2
LW_OUT, LW_STRAND, LW_HAIR = 0.4, 0.4, 0.35

# The primer lane. Fig 1 puts a primer BELOW the strand at 35.3 - "in the lane
# the label ticks used" (Fig1d/make_pcr.py), and strips the ticks to free it.
# This panel keeps the ticks, so the lane drops below their captions instead.
# Either way a primer mark can never touch a module bar, which is what Leo
# asked for: PR10137 anneals INSIDE the sensor, so "before or after the
# elements" is not something the x axis can deliver.
PRIMER_Y = 15.5
PRIMER_CAP = 21.8                # 6.3 below the arrow, as Fig 1's 35.3 -> 41.6
PRIMER_LEN = 2.8                 # mm. SCHEMATIC: a 20-nt primer is 0.44 mm
                                 # here. Fig 1 draws its primers schematic for
                                 # the same reason ("segment widths are
                                 # schematic, NOT base pairs"). Nothing
                                 # constrains this one - it is in its own lane
                                 # on open strand - so it is drawn to read.
                                 # 2.2 -> 2.8 with the bigger head (HEAD_SCALE)
                                 # so the arrow keeps a visible shaft instead
                                 # of becoming all triangle; _selfcheck asserts
                                 # the shaft survives.

# Arrowhead size, shared by the forward primer and the three bar termini so
# every primer mark on the panel is one glyph at one size.
#
# Leo, 2026-11-19: *"these primer arrows could be a little larger. I can barely
# see it. The reverse primer arrow for the lane 3 piece you made."* At
# mutation_scale 4 the head was 1.6 pt = 0.56 mm long and 0.56 mm WIDE, which
# is exactly the amplicon bar's own thickness (1.6 pt) - so on lane 3, whose
# bar is neutral-700 against the head's near-black ink, the triangle had no
# silhouette to read by and merged into the bar end. Doubling it to 8 makes the
# head 1.13 mm, twice the bar's thickness, so it reads by SHAPE and no longer
# depends on the ink-against-neutral-700 contrast.
#
# 🚨 8 is the CEILING, not a preference. matplotlib draws "-|>" with
# head_length = 0.4 * mutation_scale points, so the head is 0.4*8 pt = 1.13 mm
# and it has to fit inside REV_HEAD (1.2 mm), which in turn cannot exceed
# 1.32 mm or the head stops clearing the target module's right edge. See
# REV_HEAD. _selfcheck asserts the fit, so raising this alone aborts the build.
HEAD_SCALE = 8
HEAD_LEN = 0.4 * HEAD_SCALE * PT  # the drawn head's length, mm

# The REVERSE primers are not drawn as this glyph, and deliberately so. They
# are the amplicon bars' own termini: a bar ends at the product's 3' end,
# which IS where its reverse primer anneals, and a leftward head there says
# "this end, read leftward" - the ordinary way an amplicon is drawn under a
# map. So the head belongs AT the bar's end, not a shaft's length inside it.
#
# 🚨 That shaft was the whole fault. It was 2.6 mm = 118 bp, so both heads
# landed ~2.5 mm deep inside the target block. Just the head now.
#
# Two things were tried in between and are recorded so they are not retried:
# a full-length arrow shrunk to fit the 70 bp of bare 3' strand is ~1.2 mm,
# which is all head and no shaft and reads as a stray tick; and moving that
# arrow off the bar to sit below it reads as a smudge hanging under the bar
# end. The original glyph was right; only its position was wrong.
#
# 1.2 mm (was 0.8) to carry the bigger head. It is within 0.12 mm of the HARD
# ceiling and cannot grow further: the head's tip must clear the target
# module's right edge by 0.2 mm, the bar ends 69 bp = 1.52 mm past that edge,
# so REV_HEAD < 1.32. The bare 3' strand the mark sits on is only 70 bp
# (1.54 mm) and that is the second ceiling. Both are asserted.
REV_HEAD = 1.2

# --- the gel ---------------------------------------------------------------
# Lane geometry measured off the raw image (see README): well centres fit
# x = 239.7 + 109.39*(lane-1) to within 4 px over lanes 1-7, and lane 3's
# product band sits at x = 458.7 against a predicted 458.5. The well row is
# level (-0.016 deg), so the crop is not rotated.
GEL = dict(path=GELS / "gel_2026-07-15_primer-test_inverted.png",
           crop=(182, 150, 948, 744),          # x0, y0, x1, y1
           lane_x0=239.7, lane_dx=109.39, lane_w=100.0, box_dx=-4.0)
# `lane_w` is the expected-band BOX's width, not a measurement of the lane. It
# was 76 and the one band on this gel is 90 px wide (see BAND3_BOX), so the box
# cut 11 px into the band's left end - Leo, 2026-11-19: *"as you can see with
# the leftmost and lower one, it's not even going around the band fully."*
#
# 🚨 100 is the CEILING, and it is a tight one. Two walls close on this number
# from opposite sides and a 106 px box was tried first and rejected by looking
# at the render:
#   * lane 3's band runs 409..499, so the box must reach at least 408.5 on the
#     left, i.e. 100 px about the fitted lane centre 458.5;
#   * lane 7's centre is 896 and the crop ENDS at 948, because lane 8's well
#     rim starts at ~948 and the panel shows lanes 1-7. At 106 px that box's
#     right edge landed at 949 - over the edge of the gel, hanging in the
#     panel's margin, which is what the 106 px render showed.
# A 100 px box centred on the fitted well centre satisfies both arithmetically
# - and the render showed it still grazing the band, because the clearance it
# leaves on the left is 0.5 px = 0.05 mm and the dashed rule is 0.5 pt = 0.18
# mm wide, so the RULE sits on the band's left tail even though the box edge
# technically clears it. Hence `box_dx`.
#
# `box_dx` offsets the box from the fitted well centre. The band is NOT
# centred on its lane by extent: it reaches 49.5 px left of the fitted centre
# and 40.5 px right (its darkness centroid is 458.7, so the ink is centred;
# the left SMEAR is not). -4 px puts a 100 px box 4.5 px clear on the left and
# 5.5 on the right, and still leaves lane 7's box 6 px inside the crop.
#
# ⚠️ It is a drawing offset, not a correction to the lane model. The well fit
# is good - re-measured here, well 3's bright interior centres on 460.5 against
# the model's 458.5 - and 4 px is inside the +/-4 px residual the fit already
# carries over lanes 1-7, so the box is placed at one end of that residual
# rather than against the fit. The lane LABELS and everything else still use
# the fitted centre. _selfcheck asserts both ends of every box. A box also
# stays inside the 109.39 px lane pitch, so it can never reach a neighbour.

# GeneRuler 1 kb DNA Ladder (Thermo), size -> darkness centroid in lane 1 of
# the raw image. The protocol records "1 kb DNA ladder" with no supplier, so
# the identification is from the pattern: twelve bands, eight bunched at the
# well end, then 2000/1500/1000/750, and - the clincher - the three darkest
# within their neighbourhoods are 1000, 3000 and 6000 (0.485, 0.316, 0.278
# against 0.20-0.29 for their neighbours), which are exactly GeneRuler 1 kb's
# three reference bands. Its 500 and 250 bands are not visible on this run.
LADDER = [(10000, 524.1), (8000, 510.4), (6000, 489.4), (5000, 474.9),
          (4000, 453.9), (3500, 439.3), (3000, 419.9), (2500, 397.8),
          (2000, 365.1), (1500, 318.9), (1000, 246.0), (750, 194.2)]
LABELLED = (10000, 2000, 1000)      # sparse, and clear of the 850 / 2,650 ink
BAND3_Y = 208.4                     # the one product, darkness centroid, lane 3

# The same band's EXTENT, raw-image pixels (x0, x1, y0, y1), which is what an
# expected-band box has to enclose and what the old one did not. Measured as
# the 15%-of-depth contour against a local background taken from 20 rows above
# and below: that is the first contour that is STABLE (x is 409..499 at 15, 20
# and 25% and 410..498 at 35%, so it is the band's edge and not a threshold
# artefact) and it is also where the band stops being visible on screen. Below
# 15% the contour runs away into the lane's smear - 127 px at 10%, 160 at 5%.
BAND3_BOX = (409, 499, 196, 221)

# Where an expected product would run, drawn as a BOX rather than a rule
# (2026-11-19, Leo: "could you put a box instead of the lines? box doesnt have
# to be precise and can be bigger than the exact location of the band"). The
# box spans the expected size +/- BAND_TOL, symmetric in log bp because that is
# how the size axis interpolates. A hairline implies a precision this panel
# does not have - ~850 and ~2,650 are the worksheet's nominal figures, and the
# one product that did appear reads 812 - so a band of sizes is the honest
# mark. The window is a drawing choice, stated in the caption, not a measurement.
#
# 1.15 -> 1.25, with Leo's *"make them literally a bit taller and wider"*. The
# height is NOT just this number: see bp_to_y, whose np.interp used to CLAMP
# below the ladder's smallest band (750 bp at y 194.2), so every tolerance past
# 1.133 drew the same lower edge and the box could not be made taller at all.
# That clamp, not the tolerance, is why the box cut the band's lower edge off.
BAND_TOL = 1.25
LANE_LABEL = "expected band position"

W_PANEL = 85.0
LGUT, RGUT, TOP = 11.0, 0.8, 4.2    # size axis, slack, lane labels


def bp_to_y(bp):
    """Raw-image y for a size, interpolated along the ladder in log bp.

    🚨 Off the ends it EXTRAPOLATES on the terminal interval's slope instead of
    clamping, which `np.interp` does. The clamp was a real defect, not a
    nicety: the ladder's smallest visible band is 750 bp, the expected-band box
    runs to 850/BAND_TOL, and every tolerance past 1.133 therefore drew its
    lower edge at the 750 bp row - y 194.2, 2 px off the band's own lower edge
    at 196 - so the box could not be made taller no matter what BAND_TOL said.
    The extrapolation is 0.04 of a decade past the last band (680 against 750)
    and is used ONLY to place a drawn box edge. `y_to_bp` deliberately keeps
    the clamp: reading a SIZE off this gel outside the ladder is not something
    the ladder supports, and nothing here asks it to.
    """
    lo = np.log10([b for b, _ in LADDER][::-1])
    ys = np.array([y for _, y in LADDER][::-1], float)
    x = np.log10(bp)
    i = int(np.clip(np.searchsorted(lo, x) - 1, 0, len(lo) - 2))
    return float(ys[i] + (x - lo[i]) * (ys[i + 1] - ys[i]) / (lo[i + 1] - lo[i]))


def y_to_bp(y):
    ys = np.array([y for _, y in LADDER][::-1], float)
    lo = np.log10([b for b, _ in LADDER][::-1])
    return float(10 ** np.interp(y, ys, lo))


def mm_fig(w, h):
    """A bare canvas measured in mm, y growing downward."""
    fig = plt.figure(figsize=(style.mm(w), style.mm(h)))
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(0, w)
    ax.set_ylim(h, 0)
    ax.axis("off")
    return fig, ax


def tw(s, pt=FS):
    """Width of `s` at `pt`, in mm. Measured, never estimated."""
    return style.text_w(s, pt) / 72.0 * 25.4


# --- Fig 1's module glyph --------------------------------------------------
# A rounded rectangle, outlined in synbio.outline at 0.4, with the corner
# radius the Fig 1 paths use. NOT constructs.py's CDS pentagon: Fig 1's own
# construct map does not use one, and its SBOL bars are on the excluded list.
def _rounded(x, w, y, h, fill, edge=None, lw=LW_OUT, z=3):
    r = min(BLOCK_R * PT, w * 0.4, h * 0.4)
    return FancyBboxPatch((x + r, y - h / 2 + r), w - 2 * r, h - 2 * r,
                          boxstyle=f"round,pad={r}", facecolor=fill,
                          edgecolor="none" if edge is None else edge,
                          linewidth=lw, zorder=z)


def _bracket(ax, xa, xb, y, name):
    """Fig 1's module bracket: a spine with two feet, caption centred above."""
    top, foot = y + BRK_TOP * PT, y + BRK_ARM * PT
    ax.plot([xa, xa, xb, xb], [foot, top, top, foot], color=INK,
            lw=LW_HAIR, solid_joinstyle="miter", zorder=4)
    ax.text((xa + xb) / 2, y + BRK_CAP * PT, name, ha="center", va="baseline",
            fontsize=FS, color=INK)


def _tick(ax, x, y, name):
    """Fig 1's label tick: a hairline under the block, caption below it."""
    ax.plot([x, x], [y + TICK_TOP * PT, y + TICK_BOT * PT],
            color=palette.NEUTRAL["500"], lw=LW_HAIR, zorder=4)
    ax.text(x, y + TICK_CAP * PT, name, ha="center", va="baseline",
            fontsize=FS, color=INK)


def _axis():
    """Panel (a)'s base-pair axis: left edge, right edge, mm per base pair.

    The gutters are set by MEASURED type (common/style.text_w), not by an
    em-per-character guess, so the drawing is as wide as the labels leave it.
    Shared with _selfcheck, which measures gaps in millimetres and so needs
    the same scale the panel is drawn at.
    """
    lg = max(tw(f"lane {l}") for _, _, l in CANDIDATES) + 1.6
    rg = max(tw(f"{r}, ~{n:,} bp") for r, n, _ in CANDIDATES) + 1.6
    x0, x1 = lg + 1.0, W_PANEL - rg - 0.8
    # The axis spans the LONGEST product drawn, PR250's, which is 20 bp past
    # the worksheet's 2,650. Everything on the panel - blocks, strand, bars,
    # primer marks - is placed through bp() on this one scale, so a mark at a
    # coordinate lands where that coordinate is and nowhere else.
    span = max(n + REV_OFFSET.get(r, 0) for r, n, _ in CANDIDATES)
    return x0, x1, (x1 - x0) / span


def panel_a(out):
    W, H = W_PANEL, 26.5
    fig, ax = mm_fig(W, H)

    x0, x1, sc = _axis()
    ty, h = 6.6, BLOCK_H * PT                           # transcript baseline
    bp = lambda b: x0 + b * sc

    # the strand. Fig 1 runs a hairline through the blocks and caps it "5'" at
    # the 5' end; there is no 3' caption here, because the transcript does NOT
    # end where the amplicon does and saying so would be false.
    ax.plot([x0, x1], [ty, ty], color=INK, lw=LW_STRAND,
            solid_capstyle="round", zorder=2)
    ax.text(x0 - 0.9, ty, "5′", ha="right", va="center", fontsize=FS, color=INK)

    # the blocks, to scale, in Fig 1's shapes and Fig 1's colours
    for key, a, b, fill, label in BLOCKS:
        x, w = bp(a), (b - a + 1) * sc
        if key == "linker":                 # Fig 1 draws it thin, not full height
            ax.add_patch(_rounded(x, w, ty, LINK_H * PT, fill,
                                  palette.SYNBIO["outline"]))
            continue
        if key == "sms":                    # kinase + RNase lobe, ONE outline
            outer = _rounded(x, w, ty, h, fill, palette.SYNBIO["outline"])
            ax.add_patch(outer)
            cut = x + w * (1 - RNASE_FRAC)
            lobe = Rectangle((cut, ty - h / 2), w * RNASE_FRAC, h,
                             facecolor=palette.SYNBIO["enzyme_alt"], lw=0,
                             zorder=4)
            ax.add_patch(lobe)
            lobe.set_clip_path(outer)
            ax.plot([cut, cut], [ty - h / 2, ty + h / 2],
                    color=palette.SYNBIO["outline"], lw=0.3, zorder=5)
            ax.add_patch(_rounded(x, w, ty, h, "none",
                                  palette.SYNBIO["outline"], z=6))
        else:
            ax.add_patch(_rounded(x, w, ty, h, fill, palette.SYNBIO["outline"]))
        if label:
            _tick(ax, x + w / 2, ty, label)

    for a, b, name in BRACKETS:
        _bracket(ax, bp(a) + BRK_INSET, bp(b) - BRK_INSET, ty, name)

    # The forward primer, in its own lane BELOW the strand - Fig 1's
    # convention, and what keeps every primer mark clear of the module bars.
    py = ty + PRIMER_Y * PT
    ax.annotate("", xy=(x0 + PRIMER_LEN, py), xytext=(x0, py),
                arrowprops=dict(arrowstyle="-|>", lw=0.7,
                                color=palette.NEUTRAL["700"],
                                mutation_scale=HEAD_SCALE,
                                shrinkA=0, shrinkB=0),
                zorder=4)
    ax.text(x0 + PRIMER_LEN / 2, ty + PRIMER_CAP * PT, FWD, ha="center",
            va="baseline", fontsize=FS, color=INK)

    # the three candidate amplicons, on the same axis. Adopted against
    # abandoned is LIGHTNESS, not hue: the hues on this panel are molecules.
    for i, (rev, nt, lane) in enumerate(CANDIDATES):
        y = 17.2 + i * 3.6
        xe = bp(nt + REV_OFFSET.get(rev, 0))
        ax.plot([x0, xe], [y, y],
                color=palette.NEUTRAL["700"] if nt < 1000
                else palette.NEUTRAL["400"], lw=1.6, solid_capstyle="butt",
                zorder=3)
        # the bar's terminus: a leftward head AT the product's 3' end, which
        # is where the reverse primer anneals. Its tip is its leftmost point,
        # so asserting the tip clears the target block clears the whole mark.
        ax.annotate("", xy=(xe - REV_HEAD, y), xytext=(xe, y),
                    arrowprops=dict(arrowstyle="-|>", lw=0.7, color=INK,
                                    mutation_scale=HEAD_SCALE,
                                    shrinkA=0, shrinkB=0))
        ax.text(xe + 1.0, y, f"{rev}, ~{nt:,} bp", ha="left", va="center",
                fontsize=FS, color=INK)
        ax.text(x0 - 1.0, y, f"lane {lane}", ha="right", va="center",
                fontsize=FS, color=MUTED)

    style.save(fig, out, "panel_a_amplicon_design", tight=False)
    plt.close(fig)


def panel_b(out):
    """The gel, with the ladder as a size axis and the expectations marked."""
    img = mpimg.imread(GEL["path"])
    cx0, cy0, cx1, cy1 = GEL["crop"]
    sub = img[cy0:cy1, cx0:cx1][::-1]            # wells to the top
    ph, pw = sub.shape[:2]

    W = W_PANEL
    iw = W - LGUT - RGUT
    k = iw / pw                                  # mm per raw pixel
    ih = ph * k
    H = TOP + ih
    fig, ax = mm_fig(W, H)

    def px(x):                                   # raw x -> panel mm
        return LGUT + (x - cx0) * k

    def py(y):                                   # raw y -> panel mm (flipped)
        return TOP + (cy1 - y) * k

    ax.imshow(sub, cmap="gray", aspect="auto",
              extent=(LGUT, LGUT + iw, TOP + ih, TOP), interpolation="none",
              zorder=1)

    # size axis: a tick at every ladder band, a label on four of them
    for bp, y in LADDER:
        yy = py(y)
        ax.plot([LGUT - 1.1, LGUT - 0.3], [yy, yy], color=MUTED, lw=0.4,
                zorder=4)
        if bp in LABELLED:
            s = f"{bp:,} bp" if bp == max(LABELLED) else f"{bp:,}"
            ax.text(LGUT - 1.5, yy, s, ha="right", va="center", fontsize=FS,
                    color=MUTED)

    # the expected products: a tick in ink on the axis, a box in the lane
    for size in sorted({n for _, n, _ in CANDIDATES}):
        yy = py(bp_to_y(size))
        ax.plot([LGUT - 1.6, LGUT - 0.3], [yy, yy], color=INK, lw=0.6, zorder=4)
        ax.text(LGUT - 2.0, yy, f"{size:,}", ha="right", va="center",
                fontsize=FS, color=INK)
    for rev, nt, lane in CANDIDATES:
        c = GEL["lane_x0"] + (lane - 1) * GEL["lane_dx"] + GEL["box_dx"]
        xa, xb = px(c - GEL["lane_w"] / 2), px(c + GEL["lane_w"] / 2)
        ya, yb = py(bp_to_y(nt * BAND_TOL)), py(bp_to_y(nt / BAND_TOL))
        ax.add_patch(Rectangle((xa, ya), xb - xa, yb - ya, facecolor="none",
                               edgecolor=INK, lw=0.5,
                               linestyle=(0, (1.4, 1.2)), zorder=4))
        # ONE label for the set of three, not one per box (Leo, 2026-11-19:
        # "maybe you can label the right ones like expected location of a
        # band"). It sits on the box that has a band under it, which is the
        # one that shows what the mark means; the other two are the same mark.
        # CENTRED on that box and 0.5 mm above its top rule (Leo,
        # 2026-11-20: "Center align the expected band position band to the
        # box"), so it labels it by position and needs no leader. Centring
        # puts half the label either side of the box centre, so `_selfcheck`
        # measures the width with style.text_w and asserts BOTH halves fit -
        # the right half must stop before lane 5's box, the left half must
        # stay on the cropped gel.
        if rev == "PR10137":
            ax.text((xa + xb) / 2, ya - 0.5, LANE_LABEL, ha="center",
                    va="baseline", fontsize=FS, color=INK, zorder=5)

    # lane labels, from the measured well centres
    for lane, label in [(1, "ladder"), (3, "PR10137"), (5, "PR249"),
                        (7, "PR250")]:
        c = px(GEL["lane_x0"] + (lane - 1) * GEL["lane_dx"])
        ax.text(c, TOP - 0.6, label, ha="center", va="bottom", fontsize=FS,
                color=INK)

    style.save(fig, out, "panel_b_primer_test_gel", tight=False)
    plt.close(fig)


def main():
    style.apply()
    out = paths.fig_output(FIG)
    panel_a(out)
    panel_b(out)
    for rev, nt, lane in CANDIDATES:
        print(f"{FWD} + {rev}: ~{nt:,} bp expected, lane {lane}, "
              f"ladder y {bp_to_y(nt):.0f}, drawn to "
              f"{nt + REV_OFFSET.get(rev, 0):,}")
    print(f"observed band, lane 3: raw y {BAND3_Y} -> {y_to_bp(BAND3_Y):.0f} "
          f"bp against ~850 expected")
    print(f"expected-band boxes span +/-{(BAND_TOL - 1) * 100:.0f}% in bp: "
          + ", ".join(f"{n / BAND_TOL:.0f}-{n * BAND_TOL:.0f}"
                      for n in sorted({n for _, n, _ in CANDIDATES})))
    print("LV_PCR/LV1.4.1.xlsx: KAPA HiFi, 50 uL, 24 cycles, "
          "template MLS12.18.1.1 cDNA")
    for key, a, b, _, _ in BLOCKS:
        print(f"  {key:9s} amplicon coords {a:>5}..{b:<5} = {b - a + 1:>5} bp")
    print("wrote", out)


def _selfcheck():
    """What the size axis asserts, checked before anything is drawn."""
    # The ladder must be monotone: a bigger fragment never runs further. A
    # transposed or mistyped entry shows up here and nowhere else, because
    # bp_to_y interpolates the table and so always reproduces it.
    assert LADDER == sorted(LADDER, reverse=True), "ladder is out of order"
    assert all(y1 > y2 for (_, y1), (_, y2) in zip(LADDER, LADDER[1:])), \
        "a smaller fragment must sit further from the well"
    # Mobility in agarose is smooth in log size: no two adjacent intervals
    # may differ by more than 2x in px per decade, which is what a band
    # assigned to the wrong size looks like.
    import math
    rate = [(y1 - y2) / (math.log10(b1) - math.log10(b2))
            for (b1, y1), (b2, y2) in zip(LADDER, LADDER[1:])]
    assert max(r2 / r1 for r1, r2 in zip(rate, rate[1:])) < 2.0, \
        f"ladder spacing jumps: {[round(r) for r in rate]} px per decade"
    # and it must put the one product where it was measured
    assert 800 < y_to_bp(BAND3_Y) < 830, "lane 3 product should read ~815 bp"
    # the box is the expectation; the one product that appeared has to be
    # inside its own, or the window is narrower than the panel's own evidence
    assert 850 / BAND_TOL < y_to_bp(BAND3_Y) < 850 * BAND_TOL, \
        "lane 3's box must contain the band that actually ran"
    assert tw("10,000 bp") + 2.2 < LGUT, "size labels must clear the gutter"

    # 🚨 And it must contain the WHOLE band, not just its centroid. The box
    # used to be 76 px wide against a 90 px band and its lower edge sat 2 px
    # inside the band's, so it cut through both - which is what Leo saw. In raw
    # image pixels, where BAND3_BOX was measured, with 2 px of clearance so a
    # box edge never grazes an edge of the ink it is meant to enclose.
    def _box_x(lane):
        c = GEL["lane_x0"] + (lane - 1) * GEL["lane_dx"] + GEL["box_dx"]
        return c - GEL["lane_w"] / 2, c + GEL["lane_w"] / 2

    bx0, bx1, by0, by1 = BAND3_BOX
    la, lb = _box_x(3)
    assert la < bx0 - 2 and lb > bx1 + 2, \
        (f"lane 3's box spans x {la:.0f}..{lb:.0f} and the band spans "
         f"{bx0}..{bx1} - the box has to go AROUND the band")
    # and every box has to stay on the GEL. Lane 7's centre sits 52 px from the
    # crop's right edge, where lane 8's well rim starts, so the width that
    # clears lane 3's band on the left is also the width that fits inside the
    # image on the right - this caught a 106 px box hanging off the panel edge.
    for _, _, lane in CANDIDATES:
        a, b = _box_x(lane)
        assert GEL["crop"][0] < a and b < GEL["crop"][2], \
            (f"lane {lane}'s box runs off the cropped gel ({a:.0f}..{b:.0f} "
             f"against the crop's {GEL['crop'][0]}..{GEL['crop'][2]})")
    assert bp_to_y(850 / BAND_TOL) < by0 - 2 and bp_to_y(850 * BAND_TOL) > by1 + 2, \
        (f"lane 3's box spans y {bp_to_y(850 / BAND_TOL):.1f}..."
         f"{bp_to_y(850 * BAND_TOL):.1f} and the band spans {by0}..{by1} - "
         f"the box has to go AROUND the band")
    # The width is shared by all three boxes, so that half of the check carries
    # to lanes 5 and 7. The height cannot: neither lane produced a band, so
    # there is nothing there to contain and no honest check to write.
    # What IS checkable for every lane is that a box stays in its own lane.
    assert GEL["lane_w"] < GEL["lane_dx"], \
        "a box wider than the lane pitch would run into the next lane"

    # the one label for the set is CENTRED on lane 3's box, so half of it
    # sits either side of that box's centre. Both halves must fit: the right
    # half over empty lane-4 gel and stopping before lane 5's box, the left
    # half staying on the cropped gel. Measured, not estimated (STYLE_NOTES).
    mmpx = (W_PANEL - LGUT - RGUT) / (GEL["crop"][2] - GEL["crop"][0])
    c3 = sum(_box_x(3)) / 2
    half = tw(LANE_LABEL) / 2
    run_r = (_box_x(5)[0] - c3) * mmpx
    run_l = (c3 - GEL["crop"][0]) * mmpx
    assert half < run_r, \
        (f"'{LANE_LABEL}' half-width is {half:.1f} mm and only {run_r:.1f} mm "
         f"runs from lane 3's box centre to lane 5's box")
    assert half < run_l, \
        (f"'{LANE_LABEL}' half-width is {half:.1f} mm and only {run_l:.1f} mm "
         f"of gel runs left from lane 3's box centre")

    # --- panel (a): the blocks tile the amplicon, in order, to 2,650 -------
    assert [a for _, a, _, _, _ in BLOCKS] == sorted(a for _, a, _, _, _
                                                     in BLOCKS), \
        "blocks must be in 5' to 3' order"
    for (_, _, b, _, _), (k, a, _, _, _) in zip(BLOCKS, BLOCKS[1:]):
        assert a > b, f"{k} overlaps the block before it"
    assert BLOCKS[0][1] == 1, "the blocks must start at base 1"

    # 🚨 The fault Leo found by zooming in, 2026-11-19, and the one thing on
    # this panel that was silently wrong: BOTH reverse primers anneal 3' of
    # the target module - that is what pLV001_part158.gb shows and what the
    # caption says - but their arrowheads were landing 2.6 and 2.2 mm INSIDE
    # the pale blue block, because the block ran to the product's end.
    #
    # Measured in millimetres on the panel's own axis, so it tests the drawing
    # and not the intention behind it.
    x0, _, sc = _axis()
    f30 = next(b for b in BLOCKS if b[0] == "f30")
    right = x0 + f30[1] * sc + (f30[2] - f30[1] + 1) * sc
    for rev, nt, _ in CANDIDATES:
        if rev == "PR10137":        # anneals at ~850, INSIDE the sensor: its
            continue                # head belongs under a block
        head = x0 + (nt + REV_OFFSET.get(rev, 0)) * sc - REV_HEAD
        assert head > right + 0.2, (
            f"{rev}'s arrowhead is {right - head:+.2f} mm from the target "
            f"module's right edge - it must sit AFTER the block, not on it")
    assert REV_HEAD < (REV_FOOTPRINT + TARGET_3UTR) * sc, \
        "the bar's head is longer than the bare 3' strand it sits on"
    # 🚨 The glyph matplotlib actually DRAWS has to fit in the segment those
    # two assertions measure. mutation_scale sets the head's length in points
    # (0.4 x scale for "-|>") independently of REV_HEAD, so without this the
    # head could be enlarged straight back over the target block while both
    # checks above still passed - they measure the tip, not the triangle.
    assert HEAD_LEN <= REV_HEAD, (
        f"the drawn arrowhead is {HEAD_LEN:.2f} mm (mutation_scale "
        f"{HEAD_SCALE}) and the segment it sits in is {REV_HEAD:.2f} mm - "
        f"it would overrun the bar's end")
    # and the forward primer must keep a shaft, or it stops reading as an arrow
    assert PRIMER_LEN - HEAD_LEN > 1.0, \
        f"{FWD}'s arrow is {PRIMER_LEN - HEAD_LEN:.2f} mm of shaft: all head"
    for a, b, name in BRACKETS:
        assert any(a <= s and e <= b for _, s, e, _, _ in BLOCKS), \
            f"the {name} bracket spans no block"
    # Adjacent brackets must read as two, not one. The recorder and the sensor
    # are 8 bp apart on this template, which is 0.18 mm here, so without the
    # inset their feet merge into a single rule with a notch in it - which is
    # what the first render of this revision did.
    _, _, sc = _axis()
    for (_, b, n1), (a, _, n2) in zip(BRACKETS, BRACKETS[1:]):
        gap = (a - b) * sc + 2 * BRK_INSET
        assert gap > 0.6, f"{n1} and {n2} brackets are {gap:.2f} mm apart"
    # and the primer lane must sit clear of the deepest block ink, which is
    # the label captions hanging under the blocks
    assert PRIMER_Y > TICK_CAP, "the primer lane overlaps the label captions"

    # Fig 1's colours, not a second opinion about them: constructs.py is the
    # module table the main figures agree on, so a palette edit that moved one
    # of these would have to move both. The SHAPES are Fig 1's and deliberately
    # not constructs.py's (see the docstring); only the fills are shared.
    fills = {k: f for k, _, _, f, _ in BLOCKS}
    for key in ("utr", "at", "binder", "f30"):
        assert fills[key] == K.MODULES[key][1], \
            f"{key} has drifted from constructs.py"
    assert fills["sms"] == K.MODULES["sms"][1], "SMS fill drifted"
    assert palette.SYNBIO["enzyme_alt"] == K.MODULES["sms"][4], \
        "the RNase lobe drifted from constructs.py"
    assert 0.3 < RNASE_FRAC < 0.45, \
        f"the RNase lobe is {RNASE_FRAC:.3f} of SMS - re-read the Fig 1 file"

    # PR250's product must be drawn LONGER than PR249's, by the 20 bp its
    # binding site sits 3' of PR249's in pLV001_part158.gb. If that offset is
    # ever zeroed, the panel silently stops answering the question Leo asked.
    assert REV_OFFSET["PR250"] - REV_OFFSET["PR249"] == 20, \
        "PR250 anneals 20 bp 3' of PR249 (3273..3292 against 3253..3272)"

    # the two marked expectations must not collide with a labelled ladder tick
    mm = (W_PANEL - LGUT - RGUT) / (GEL["crop"][2] - GEL["crop"][0])
    rows = sorted([bp_to_y(b) for b in LABELLED]
                  + [bp_to_y(n) for _, n, _ in CANDIDATES])
    gaps = [(b - a) * mm for a, b in zip(rows, rows[1:]) if b - a > 1e-6]
    assert min(gaps) > 1.95, f"size-axis labels overlap: {min(gaps):.2f} mm"


if __name__ == "__main__":
    _selfcheck()
    main()
