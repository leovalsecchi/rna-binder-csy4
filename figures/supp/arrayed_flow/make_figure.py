"""arrayed_flow - the LV2.33/LV2.34.1 arrayed re-test of the 32 ordered candidates.

Three bare panels over one variant set and ONE variant order (32 candidates ranked by
ON/OFF, then the two Csy4-H29A control doses past a rule):

  (a) every well: normalized mScarlet per variant x target, as log2 of the variant's
      own mean over its four targets. Absolute level is a transfection property and
      spans ~120x across variants, so the raw number would show nothing but that.
      HORIZONTAL since 2026-11-19 - the 34 variants run left to right and the four
      targets are four rows.
  (b) ranked ON/OFF with the spread it is actually entitled to - there is ONE
      biological replicate, so the band is an empirical null measured on the same
      wells, not a replicate SD, and the per-variant line is the range over the four
      positive x negative pairings.
  (c) those four pairings, one row per variant. This is the panel the pLV15/pLV16
      substitution makes necessary: only PosA/NegA pairs a target with a negative of
      its own scaffold family, so the four ratios are not four measurements of one
      quantity and must not be averaged silently.

LAYOUT, 2026-11-19 (Leo): "the arrayed flow, you can make the A horizontal, and then B
and C, you can really expand the x-axis so that they're almost like a square size plot,
each of them." So (a) is transposed and (b)/(c) are widened from 66/62 mm to BC_W_MM.
That costs the one thing the old block bought - (a)'s rows no longer line up with
(b)/(c)'s, a reader reads (a) across and (b)/(c) down - and buys the x range the ranked
panels were short of. The ORDER is untouched, so a variant is at the same place along
(a) as it is down (b). The three sizes are not free: make_supp_book.py packs panels into
rows by their EMITTED size, so BOOK_* below is what decides that (a) takes its own row
and (b) and (c) sit side by side. main() asserts it, and asserts the A4 page still fits.

(b) and (c) are built like the main text's Fig5a_PosAvg_vs_NegAvg and share ONE x
range, which THE 32 CANDIDATES SET: the Csy4-H29A controls (10.8x, 96.7x) and the two
failed wells (0.06) are reported but never allowed into the limits, because between
them they used to stretch the axis over 0.05-200 and crush the designs into a fifth of
it. See xlimits(), clamp() and ctrl_text().

Condition colours are theme steps and a neutral, never a route colour (COLORS.md
rules 2 and 3); route colour appears only in (b), where the mark is a design.

scipy is broken in this venv (_propack dlopen), so the two tests in report() are
hand-rolled: an exact two-sided binomial by enumeration and an exact Wilcoxon
signed-rank by DP over the rank set. Both are exact, not normal approximations.

Run: ~/.venvs/binder_map/bin/python3 code/supp/arrayed_flow/make_figure.py
"""
import math
import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import FuncFormatter

from common import palette, paths, style

FIG = Path(__file__).resolve().parent.name

CUT = 1.2                 # proposed binders >= 1.2; below, inactive designs
MIN_KEEP = 0.5            # under this a ratio is a dead WELL, not a weak design
TARGETS = ["PosA", "PosB", "NegA", "FullNeg"]
# Short on the panel, expanded in the caption (STYLE_NOTES section 3).
TLAB = {"PosA": "posA", "PosB": "posB", "NegA": "NegA", "FullNeg": "no scaf"}
# Pairing -> (which positive, which negative). Only the first is scaffold-matched.
PAIRS = [("PosA", "NegA"), ("PosA", "FullNeg"), ("PosB", "NegA"), ("PosB", "FullNeg")]
PAIR_LAB = {("PosA", "NegA"): "PosA / NegA  (matched)",
            ("PosA", "FullNeg"): "PosA / no array",
            ("PosB", "NegA"): "PosB / NegA  (unmatched)",
            ("PosB", "FullNeg"): "PosB / no array"}
POS_STEP = {"PosA": "700", "PosB": "400"}     # colour = positive target
FILLED = {"NegA": True, "FullNeg": False}     # fill = which negative
CTRL = {"P6_C4": "Csy4-H29A 100 ng", "P6_C5": "Csy4-H29A 10 ng"}
# Horizontal (a) hangs the control label under a ~3.5 mm column, and
# "Csy4-H29A 100 ng" measures 16.1 mm - rotated, that is 16 mm of empty panel
# below 32 two-digit numbers. Split it over two tiers instead: the dose on the
# column, the name once under the pair. Same two labels, 6 mm instead of 16.
# (b) and (c) keep CTRL whole, where it is a horizontal row label with room.
CTRL_NAME = "Csy4-H29A"
CTRL_DOSE = {"P6_C4": "100 ng", "P6_C5": "10 ng"}
assert all(CTRL[k] == f"{CTRL_NAME} {d}" for k, d in CTRL_DOSE.items())
PITCH_MM = 2.4            # 5.5 pt type is ~1.94 mm tall; this is the floor
PT_MM = 25.4 / 72.0       # style.text_w returns points when size is in points

# --- authored sizes, and the book budget that sets them ---------------------
# make_supp_book.py places panels left to right, in filename order, into rows no
# wider than the A4 text block (210 - 2 x 14 mm margin) with a 6 mm gutter, and
# it packs by the EMITTED size, i.e. after bbox_inches="tight" has cropped each
# panel. So these two numbers, not the figsizes, decide the page: (b) and (c)
# are side by side only while each stays inside (182 - 6) / 2, and (a) keeps its
# own row above them only while it is too wide to follow (a) + 6 + (b) <= 182.
# main() asserts both on the files it just wrote.
BOOK_AVAIL, BOOK_GUTTER, BOOK_MARGIN, BOOK_TITLE = 182.0, 6.0, 14.0, 10.0
A4_H = 297.0
BC_BUDGET = (BOOK_AVAIL - BOOK_GUTTER) / 2        # 88.0 mm each
# (a) horizontal: 35 columns (32 candidates, a blank, 2 controls) along 152 mm,
# so a column is ~4 mm and the 1-2 digit binder numbers fit upright under it
# (measured: 0.51-2.02 mm). Wide enough to hold its own book row, well inside
# the 182 mm text block so the book never has to scale the page down.
A_W_MM, A_H_MM = 152.0, 38.0
A_GRID_TOP = 0.78        # tight_layout rect: the top 22% is the colourbar's
# (b)/(c): as wide as the budget allows, which is what "almost square" costs.
# Was 66 / 62 mm against an 82 mm tall row block - tall and narrow.
BC_W_MM = 86.0
BC_H_MM = PITCH_MM * 34 + 22
# Heatmap scale. Every one of the 128 candidate wells lands inside ±0.6 log2 of
# its variant's mean EXCEPT the two dead positive wells (binder 4 PosB 0.0061,
# binder 6 PosA 0.0029), which reach 2.9 and 2.2 and would flatten everything
# else if the scale followed them. They are clipped and outlined instead;
# panel_a asserts that no third well joins them.
VMAX = 0.6
# Fig5a_wetlab_on_off's tick set, which is the main-text panel this one is built
# like. It was sparse here (decades only) because the axis used to span 0.05-200
# and 0.7 / 1.5 collided at 5.5 pt; on the clipped axis five ticks fit with ~9 mm
# between them, and a ranked panel over half a decade needs them.
TICKS = [0.05, 0.1, 0.2, 0.3, 0.5, 0.7, 1, 1.5, 2, 3, 5, 10]


# --------------------------------------------------------------- load + check
def load():
    d = pd.read_csv(paths.INPUTS / "lv2_34_1_flow_on_off.csv", index_col=0)
    idx = pd.read_csv(paths.INPUTS / "candidates32_index.tsv", sep="\t")

    assert len(d) == 34, f"{len(d)} rows: the kept plate columns are 35 - 1 (P6_C3)"
    assert list(d.index[d.design_id.isna()]) == list(CTRL), "controls moved"
    # The stored ratios are re-derived from the four raw condition columns rather
    # than trusted: this is the one check that fails if the table is ever rebuilt
    # with a different averaging convention.
    on = d[["PosA", "PosB"]].mean(axis=1)
    off = d[["NegA", "FullNeg"]].mean(axis=1)
    # rtol, not atol: the stored columns are rounded to 5 dp and the smallest
    # wells are ~0.003, so the ratio of two of them carries ~5e-4 relative error.
    assert np.allclose(on / off, d.ONOFF_avg, rtol=1e-3), "ONOFF_avg does not re-derive"
    assert np.allclose(d.PosA / d.NegA, d.ONOFF_PosA_NegA, rtol=1e-3)

    cand = d[d.design_id.notna()].merge(
        idx.set_index("flow_well")[["twist_id", "generator", "inverse_folding",
                                    "hts_comparison"]],
        left_index=True, right_index=True, how="left")
    assert cand.twist_id.notna().all() and len(cand) == 32
    cand["route"] = cand.generator + "+" + cand.inverse_folding
    cand["label"] = cand.twist_id.astype(int).astype(str)
    cand["color"] = cand["stack"].map(palette.stack_color)
    assert cand.color.notna().all(), "a route has no palette colour"

    ctrl = d.loc[list(CTRL)].copy()
    ctrl["label"] = [CTRL[i] for i in ctrl.index]
    ctrl["color"] = palette.NEUTRAL["600"]
    # One order for all three panels: strongest ON/OFF at the top.
    cand = cand.sort_values("ONOFF_avg")
    ctrl = ctrl.sort_values("ONOFF_avg")
    return cand, ctrl


def slots(cand, ctrl):
    """Variant positions: candidates 0..31, a blank, then the two controls.

    One set of positions for all three panels, read down the y axis in (b)/(c)
    and along the x axis in the transposed (a). The blank slot is what the rule
    sits in, so no panel has to special-case the gap.
    """
    y_c = np.arange(len(cand), dtype=float)
    y_k = len(cand) + 1.0 + np.arange(len(ctrl))
    return y_c, y_k, len(cand) + 0.5


def ax_mm(fig, ax):
    """The DRAWN axes box, in mm - what "almost square" is a claim about.

    tight_layout decides it from the labels, so it is not the figsize and it is
    not knowable before a draw.
    """
    fig.canvas.draw()
    bb = ax.get_window_extent().transformed(fig.dpi_scale_trans.inverted())
    return bb.width * 25.4, bb.height * 25.4


def emitted_mm(out, stem):
    """Size of the panel just saved, in mm, read back off its SVG header.

    bbox_inches="tight" crops to the drawn content, so the size the book packs
    by is only knowable from the file. matplotlib always writes pt; be loud if
    it ever does not (same check make_supp_book.svg_size makes).
    """
    head = (out / f"{stem}.svg").read_text()[:600]
    got = [re.search(rf'{k}="([\d.]+)pt"', head) for k in ("width", "height")]
    assert all(got), f"{stem}.svg: no pt width/height in the header"
    return tuple(float(m.group(1)) * PT_MM for m in got)


def null_sd(cand):
    """Spread of the SAME statistic on a contrast with no binding in it.

    Both negatives are present for every variant, so log2(FullNeg/NegA) is one
    ON/OFF measured on two plasmids that differ by a scaffold rather than by a
    hairpin. Its SD is this experiment's well-to-well noise, which a single
    replicate cannot otherwise give. Its CENTRE is a plasmid offset, not noise
    (the two negatives are different constructs), so only the spread is used and
    the band is centred on the no-effect line. After lv2_33_stats.py, which does
    the same thing on the per-cell mScarlet-positive fraction.
    """
    return float(np.log2(cand.FullNeg / cand.NegA).std(ddof=1))


def pairings(frame):
    """The four positive x negative ratios, one column each."""
    return pd.DataFrame({f"{p}/{n}": frame[p] / frame[n] for p, n in PAIRS})


def xlimits(ratios, extra):
    """The x range for (b) and (c). THE CANDIDATES SET THE SCALE, nothing else.

    2026-11-19 (Leo): "cut the controls and super large deviations on the x
    axis". Two things were setting the extent and neither is a candidate
    effect: the Csy4-H29A doses at 10.8x and 96.7x, and the two failed wells
    (binder 4 posB, binder 6 posA) which drag a ratio down to 0.06. Between
    them the axis ran 0.05-200 and the 32 designs - 0.48 to 1.56, a factor of
    3 - were a vertical smear in a fifth of it.

    Both are reported, neither is allowed into the limits: the controls as
    their own number printed in their own row (what panel (a) already does),
    the off-scale points as a clamp mark on the edge carrying the true value.
    Rule and the 1.35 pad are Fig5a_wetlab_on_off's, so the supplementary
    panels and the main-text one read at one scale.

    `ratios` and `extra` must already be the candidates alone.
    """
    v = pd.concat([ratios.stack(), extra])
    kept = v[v >= MIN_KEEP]
    assert len(kept) > len(v) - 8, "more off-scale points than the known dead wells"
    return float(kept.min() / 1.35), float(kept.max() * 1.35)


def clamp(ax, x, y, lo, hi, color, label=None, ms=3.5):
    """An off-scale point, drawn on the edge it ran off, carrying its value.

    Fig5a_wetlab_on_off's idiom: a clipped axis must never hide that something
    sits outside it. `clip_on=False` so the mark survives on the spine.
    """
    over = x > hi
    at = hi if over else lo
    ax.plot(at, y, marker=">" if over else "<", ms=ms, c=color, clip_on=False,
            zorder=4)
    if label:
        # Lifted half a row: unlike the main-text panel this one draws a range
        # line through every row, and a centred label sits on it as a
        # strikethrough. Half a row clears the line and still clears the row
        # above (~6 pt of row pitch against 5.5 pt type).
        ax.annotate(label, (at, y), textcoords="offset points", xytext=(5, 2.4),
                    ha="left", va="center", color=palette.INK,
                    annotation_clip=False)


def ctrl_name(fig, texts, below):
    """"Csy4-H29A" once, off the MEASURED edge of the two dose labels.

    The name is split from the dose (CTRL_NAME / CTRL_DOSE) because whole it is
    16.1 mm, and it was being spent twice: 16 mm of hanging rotated label under
    (a)'s 3.9 mm columns, and 16 mm of (b)/(c)'s left margin, where it was the
    single largest reason those panels were narrow. Split, it costs 6 mm and a
    word, and says the same thing.

    Placed off the rendered boxes of the dose labels rather than a guessed
    offset: they are rotated in (a) and upright in (b)/(c), so their extent is
    not the same number twice, and in (b)/(c) it has to clear the tick labels
    without landing on the y label. Call it AFTER the layout is settled.
    """
    bb = [t.get_window_extent(fig.canvas.get_renderer()) for t in texts]
    px = fig.dpi / 25.4
    W, H = fig.get_figwidth() * 25.4, fig.get_figheight() * 25.4
    mid = lambda f: sum(f(b) for b in bb) / len(bb)
    if below:
        fig.text(mid(lambda b: (b.x0 + b.x1) / 2) / (px * W),
                 (min(b.y0 for b in bb) - 0.7 * px) / (px * H),
                 CTRL_NAME, ha="center", va="top")
    else:
        # ha="right", not "center": rotation_mode "default" rotates first and
        # aligns the ROTATED box, so a centred anchor puts half the label's ink
        # height back over the gap - which is how the first render had the "4"
        # of Csy4 sitting on the "1" of "100 ng".
        fig.text((min(b.x0 for b in bb) - 0.7 * px) / (px * W),
                 mid(lambda b: (b.y0 + b.y1) / 2) / (px * H),
                 CTRL_NAME, ha="right", va="center", rotation=90)


def ctrl_text(ax, ctrl, y_k, lo, value):
    """Controls as a printed number in their own row, not a dot on the axis.

    They are the scale reference and have to stay findable, but at 10.8x and
    96.7x they cannot also set the limits. `value` maps the control's index to
    the string for its row. Panel (a) already reports them this way.
    """
    for yy, i in zip(y_k, ctrl.index):
        ax.annotate(value(i), (lo * 1.06, yy), ha="left", va="center",
                    color=palette.NEUTRAL["600"])


def yaxis(ax, y_c, y_k, split, cand, ctrl, labels=True):
    # Dose only; ctrl_name() adds "Csy4-H29A" once beside the pair, after the
    # layout. Whole, the two names were 16.1 mm of a panel 86 mm wide.
    ax.set_yticks(np.r_[y_c, y_k])
    ax.set_yticklabels(list(cand.label) + [CTRL_DOSE[i] for i in ctrl.index]
                       if labels else [])
    ax.set_ylim(-0.9, y_k[-1] + 0.9)
    ax.axhline(split, color=palette.GRID, lw=0.6, zorder=1)
    ax.tick_params(axis="y", length=0)


def logx(ax, lo, hi, label):
    ax.set_xscale("log")
    ax.set_xlim(lo, hi)
    ax.set_xticks([t for t in TICKS if lo <= t <= hi])
    ax.set_xticks([], minor=True)
    ax.xaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{v:g}"))
    ax.set_xlabel(label)


def band(ax, sd, lo, hi):
    """+-2 SD of the empirical null, centred on no effect.

    Pale fill plus a hairline at each edge, not a mid grey: once the axis was
    clipped to the candidates (2026-11-19) the band covers two thirds of the
    frame, and at NEUTRAL 200 it read as the panel's background and pulled the
    contrast out of every dot sitting on it. The edges are what carry the
    +-2 SD claim, so they are drawn and the fill is only a tint.
    """
    f = 2 ** (2 * sd)
    ax.axvspan(max(lo, 1 / f), min(hi, f), color=palette.NEUTRAL["100"], lw=0,
               zorder=0)
    for e in (1 / f, f):
        if lo < e < hi:
            ax.axvline(e, color=palette.NEUTRAL["300"], lw=0.5, zorder=1)
    ax.axvline(1, color=palette.NEUTRAL["500"], lw=0.75, ls="--", zorder=2)
    return f


def check_columns(fig, ax, span, labels, n_up):
    """Every column label fits its column, MEASURED, never estimated.

    style.text_w, not 0.52 em per character: that stand-in is out by 10-33% and
    by a different amount for every string, which is exactly how a laid-out row
    ends up unevenly spaced for no visible reason (style.py, "measuring type").

    Two different fits, because the 34 labels are not one kind of label:
      * the 32 binder numbers stay UPRIGHT, which is what a reader wants on a
        number, so their INK WIDTH has to clear the column pitch;
      * "Csy4-H29A 100 ng" is 16.1 mm and clears no pitch this panel could have,
        so the two control names turn 90 deg. Rotated, the width stops being the
        constraint and the ink HEIGHT becomes it.
    The pitch is only known after a draw (tight_layout sets the axes box from
    the labels), so this runs on the laid-out figure, and the last assert is the
    catch-all: the real rendered boxes of all 34, in order, must not overlap.
    """
    w_ax, _ = ax_mm(fig, ax)
    pitch = w_ax / span
    up = {s: style.text_w(s) * PT_MM for s in labels[:n_up]}
    assert max(up.values()) + 0.6 <= pitch, (
        f"upright column labels want {max(up.values()):.2f} + 0.6 mm of a "
        f"{pitch:.2f} mm column (widest {max(up, key=up.get)!r})")
    rot = {s: style.text_w(s) * PT_MM for s in labels[n_up:]}
    boxes = [t.get_window_extent(fig.canvas.get_renderer())
             for t in ax.get_xticklabels()]
    px_mm = fig.dpi / 25.4
    gaps = [(b.x0 - a.x1) / px_mm for a, b in zip(boxes, boxes[1:])]
    assert min(gaps) > 0, (
        f"column labels overlap by {-min(gaps):.2f} mm at "
        f"{labels[int(np.argmin(gaps))]!r}/{labels[int(np.argmin(gaps)) + 1]!r}")
    return pitch, min(gaps), max(rot.values())


# ------------------------------------------- (a) every well, per-variant scale
def panel_a(cand, ctrl, out):
    """Horizontal: variants along x, the four targets as four rows.

    Transposed 2026-11-19 on Leo's "you can make the A horizontal". What moves
    is only the orientation: same wells, same order (strongest ON/OFF at the
    RIGHT, where it used to be at the top), same +-0.6 clip with the same
    assert, same two control doses past the rule with their own number printed
    rather than drawn as a colour the scale cannot carry. The four target names
    stop being rotated 45 deg - as rows they are horizontal and need no trick -
    and the 34 variant labels take over the rotation question, which
    check_columns() settles by measuring rather than guessing.
    """
    all_ = pd.concat([cand, ctrl])
    v = np.log2(all_[TARGETS].div(
        np.exp(np.log(all_[TARGETS]).mean(axis=1)), axis=0))
    vmax = VMAX
    over = v.loc[cand.index].abs().max(axis=1) > vmax
    assert set(cand.index[over]) == {"P1_C4", "P1_C6"}, (
        f"new wells outside the ±{vmax} scale: {list(cand.index[over])}")

    x_c, x_k, split = slots(cand, ctrl)
    x = np.r_[x_c, x_k]
    # One grid, with the gap COLUMN left blank, so the cells are contiguous and
    # the rule has somewhere to sit.
    grid = np.full((len(TARGETS), int(x_k[-1]) + 1), np.nan)
    grid[:, x.astype(int)] = v[TARGETS].values.T

    fig, ax = plt.subplots(figsize=(style.mm(A_W_MM), style.mm(A_H_MM)))
    mesh = ax.pcolormesh(np.arange(grid.shape[1] + 1) - 0.5,
                         np.arange(len(TARGETS) + 1) - 0.5,
                         np.ma.masked_invalid(grid), cmap="thesis_div",
                         vmin=-vmax, vmax=vmax)
    # Clipped cells carry a hairline so the scale cannot be misread as the value.
    for j, t in enumerate(TARGETS):
        for xx in x[(v[t].abs() > vmax).values]:
            ax.add_patch(plt.Rectangle((xx - 0.5, j - 0.5), 1, 1, fill=False,
                                       ec=palette.INK, lw=0.35, zorder=3))

    # The control columns run off the candidates' scale by design; their own
    # numbers are what matters there, so they are printed rather than implied by
    # colour - above the column now that the row they used to sit beside is a
    # column. "11x"/"97x" measure 2.76/3.03 mm against a ~4 mm pitch, so they
    # are the one annotation here that can stay upright.
    for xx, i in zip(x_k, ctrl.index):
        ax.annotate(f"{ctrl.ONOFF_avg[i]:.0f}x", (xx, -0.62), ha="center",
                    va="bottom", color=palette.NEUTRAL["600"],
                    annotation_clip=False)
    ax.set_yticks(range(len(TARGETS)), [TLAB[t] for t in TARGETS])
    # Inverted, so posA is the TOP row: the target order reads the way it did
    # when it was the column order, top to bottom.
    ax.set_ylim(len(TARGETS) - 0.5, -0.5)
    labels = list(cand.label) + [CTRL_DOSE[i] for i in ctrl.index]
    ax.set_xticks(x, labels)
    for t in ax.get_xticklabels()[len(cand):]:
        t.set_rotation(90)
        t.set_va("top")
    ax.set_xlim(-0.9, x_k[-1] + 0.9)
    # Down the MIDDLE of the blank column, not on `split`. (b)/(c) put the rule
    # on split because nothing else is at that y; here split is the boundary the
    # first control cell is drawn on, and all eight control cells are outlined
    # (off-scale by construction), so the rule landed under an INK edge and
    # vanished - the gap read as a gap with no rule in it.
    ax.axvline(split - 0.5, color=palette.GRID, lw=0.6, zorder=1)
    ax.set_xlabel("Ordered binder")
    ax.tick_params(length=0)
    for s in ax.spines.values():
        s.set_visible(False)

    # `rect` keeps the top band of the canvas clear for the colourbar, and
    # tight_layout runs BEFORE the bar exists. An inset placed above the axes
    # and then tight_laid-out is a feedback loop: the inset is positioned in
    # axes fractions, so shrinking the axes to fit it moves it down and the
    # layout shrinks again - which flattened this grid to 4 mm the first time.
    fig.tight_layout(rect=[0, 0, 1, A_GRID_TOP])
    fit = check_columns(fig, ax, x_k[-1] + 1.8, labels, len(cand))

    # Colourbar above the grid, at the left end. The free space on a 140 x 17 mm
    # panel is along its LONG edge, and the right end of that band is taken by
    # the two control numbers. Beside the bar it would have to be rotated, and
    # the label (17.7 mm) is longer than a bar bounded by the grid height can be,
    # so it came out clipped.
    cax = ax.inset_axes([0.0, 1.15, 0.19, 0.08])
    cb = fig.colorbar(mesh, cax=cax, orientation="horizontal")
    cax.xaxis.set_ticks_position("top")
    cax.xaxis.set_label_position("top")
    cb.set_label("log2 vs design mean")
    cb.outline.set_visible(False)
    cb.ax.tick_params(length=2, width=0.5)
    ctrl_name(fig, ax.get_xticklabels()[len(cand):], below=True)
    style.save(fig, out, "panel_a_wells_by_target")
    plt.close(fig)
    return vmax, v, fit


# ------------------------------------ (b) ranked ON/OFF, with honest spread
def panel_b(cand, ctrl, ratios, sd, xlim, out):
    """Built like the main text's Fig5a_PosAvg_vs_NegAvg - same 32 designs,
    same CSV, same ranking, same clipped scale, same clamp idiom, same colour
    by stack - and carrying two things that panel does not: the empirical null
    band and the per-design range over the four pairings. Those are the
    evidence the supplement exists for, so they stay (see README, "Open").

    Widened to BC_W_MM on 2026-11-19 ("really expand the x-axis so that they're
    almost like a square size plot"): the 34-row block fixes the height at
    PITCH_MM, so the only lever on the aspect is the width, and the book's
    two-panel budget is how much of it there is.
    """
    y_c, y_k, split = slots(cand, ctrl)
    lo, hi = xlim
    r = ratios.loc[cand.index]

    fig, ax = plt.subplots(figsize=(style.mm(BC_W_MM), style.mm(BC_H_MM)))
    f = band(ax, sd, lo, hi)
    # Range over the four pairings, cut at the frame rather than allowed to set it.
    ax.hlines(y_c, r.min(axis=1).clip(lower=lo), r.max(axis=1).clip(upper=hi),
              color=palette.NEUTRAL["400"], lw=0.6, zorder=2)
    inside = cand.ONOFF_avg.between(lo, hi).values
    ax.scatter(cand.ONOFF_avg[inside], y_c[inside], s=14, c=cand.color[inside],
               edgecolors=palette.SURFACE, linewidths=0.4, zorder=3)
    for yy, i in zip(y_c[~inside], cand.index[~inside]):
        clamp(ax, cand.ONOFF_avg[i], yy, lo, hi, cand.color[i],
              f"{cand.ONOFF_avg[i]:.2f}x")
    # A row whose RANGE runs off the frame gets a caret in the line's own
    # colour - except where the dot already clamped there, since one mark on
    # one row already says "this continues off-scale".
    for yy, i, ins in zip(y_c, cand.index, inside):
        if ins and (r.loc[i].min() < lo or r.loc[i].max() > hi):
            clamp(ax, r.loc[i].min() if r.loc[i].min() < lo else r.loc[i].max(),
                  yy, lo, hi, palette.NEUTRAL["400"], ms=2.5)
    ctrl_text(ax, ctrl, y_k, lo, lambda i: f"{ctrl.ONOFF_avg[i]:.0f}x")
    ax.axvline(CUT, color=palette.THEME["600"], lw=0.6, ls=":", zorder=2)
    # Both labels sit inside the frame, above the top row: the foot of a ranked
    # panel is where the tick labels are.
    ax.annotate(f"≥ {CUT:g}", (CUT * 1.08, split + 0.3), color=palette.THEME["600"],
                ha="left", va="bottom")
    ax.annotate(f"null ±2 SD ({f:.2f}x)", (1, y_k[-1] + 0.55),
                color=palette.NEUTRAL["600"], ha="center", va="bottom")
    logx(ax, lo, hi, "ON/OFF   avg Pos / avg Neg")
    yaxis(ax, y_c, y_k, split, cand, ctrl)
    ax.set_ylabel("Ordered binder")
    handles = [Line2D([], [], ls="", marker="o", ms=3, mfc=palette.stack_color(k),
                      mec=palette.SURFACE, label=palette.SHORT[k])
               for k in palette.STACK_ORDER if (cand["stack"] == k).any()]
    ax.legend(handles=handles, loc="lower center", bbox_to_anchor=(0.5, 1.0),
              ncol=2, frameon=False, handletextpad=0.3, columnspacing=1.0)
    fig.tight_layout()
    box = ax_mm(fig, ax)
    ctrl_name(fig, ax.get_yticklabels()[-len(ctrl):], below=False)
    style.save(fig, out, "panel_b_onoff_ranked")
    plt.close(fig)
    return box


# ------------------------------------------- (c) the four pairings, per variant
def panel_c(cand, ctrl, ratios, sd, xlim, out):
    """Same scale, same width and the same two exclusions as (b) - controls
    printed, dead wells clamped - so a row reads straight across from (b) to
    (c), which is what the book's side-by-side placement is for."""
    y_c, y_k, split = slots(cand, ctrl)
    lo, hi = xlim
    r = ratios.loc[cand.index]

    fig, ax = plt.subplots(figsize=(style.mm(BC_W_MM), style.mm(BC_H_MM)))
    band(ax, sd, lo, hi)
    for p, n in PAIRS:
        col = palette.THEME[POS_STEP[p]]
        kw = dict(facecolors=col, edgecolors=col) if FILLED[n] else \
            dict(facecolors="none", edgecolors=col)
        v = r[f"{p}/{n}"]
        ok = v.between(lo, hi).values
        ax.scatter(v[ok], y_c[ok], s=11, linewidths=0.5, zorder=3, **kw)
    # Both off-scale pairings in a row share their positive (binder 4 posB,
    # binder 6 posA - the two failed wells), so ONE caret and one range label
    # per row says it; two triangles would land on the same point.
    off = r < lo
    for yy, i in zip(y_c, cand.index):
        if not off.loc[i].any():
            continue
        v = r.loc[i][off.loc[i].values]
        pos = sorted({c.split("/")[0] for c in v.index})
        # One caret per row only holds while both off-scale pairings share a
        # positive; if a third dead well ever appears this has to be revisited.
        assert len(pos) == 1, f"{i}: off-scale pairings span two positives"
        # Caret only, no value: a label long enough to carry both numbers runs
        # into the same row's in-range dots (row 4's nearest sits at 0.97), and
        # these are dead-well readings, which the caption already names. The
        # averaged ON/OFF of the same two rows is labelled in (b), where the
        # left of the row is empty.
        clamp(ax, v.min(), yy, lo, hi, palette.THEME[POS_STEP[pos[0]]], ms=3)
    ctrl_text(ax, ctrl, y_k, lo, lambda i: "\u2013".join(
        f"{v:.0f}" for v in (ratios.loc[i].min(), ratios.loc[i].max())) + "x")
    logx(ax, lo, hi, "ON/OFF   by pairing")
    yaxis(ax, y_c, y_k, split, cand, ctrl)
    ax.set_ylabel("Ordered binder")
    handles = [Line2D([], [], ls="", marker="o", ms=3,
                      mfc=palette.THEME[POS_STEP[p]] if FILLED[n] else "none",
                      mec=palette.THEME[POS_STEP[p]], label=PAIR_LAB[(p, n)])
               for p, n in PAIRS]
    ax.legend(handles=handles, loc="lower center", bbox_to_anchor=(0.5, 1.0),
              ncol=1, frameon=False, handletextpad=0.3, labelspacing=0.25)
    fig.tight_layout()
    box = ax_mm(fig, ax)
    ctrl_name(fig, ax.get_yticklabels()[-len(ctrl):], below=False)
    style.save(fig, out, "panel_c_pairings")
    plt.close(fig)
    return box


# ----------------------------------------------------------- numbers + tests
def binom_two_sided(k, n, p=0.5):
    """Exact two-sided binomial by enumeration (method of small p). No scipy."""
    pm = [math.comb(n, i) * p ** i * (1 - p) ** (n - i) for i in range(n + 1)]
    return sum(v for v in pm if v <= pm[k] * (1 + 1e-12))


def wilcoxon_exact(diff):
    """Exact two-sided Wilcoxon signed-rank by DP over the rank set. No scipy."""
    d = np.asarray([v for v in diff if v != 0], float)
    rk = pd.Series(np.abs(d)).rank().values
    assert np.allclose(rk, np.round(rk)), "ties in |d|: the exact test is invalid"
    rk = rk.round().astype(int)
    w, tot = int(rk[d > 0].sum()), int(rk.sum())
    dp = np.zeros(tot + 1)
    dp[0] = 1
    for r in rk:                      # count rank subsets by their sum
        dp[r:] += dp[:tot + 1 - r]
    dp /= dp.sum()
    lo = min(w, tot - w)
    return len(d), w, float(min(dp[:lo + 1].sum() + dp[tot - lo:].sum(), 1.0))


def report(cand, ctrl, ratios, sd, vmax, xlim):
    r = cand.ONOFF_avg
    n_up = int((r > 1).sum())
    print(f"candidates n={len(r)}  >1.0 {n_up}/{len(r)} ({100*n_up/len(r):.0f}%)  "
          f">={CUT} {int((r >= CUT).sum())}/{len(r)}  max {r.max():.4f} ({r.idxmax()})  "
          f"median {r.median():.4f}  IQR {r.quantile(.25):.4f}-{r.quantile(.75):.4f}")
    print(f"  exact two-sided binomial {n_up}/{len(r)} vs 0.5: p = "
          f"{binom_two_sided(n_up, len(r)):.3g}")
    n_, w_, p_ = wilcoxon_exact(cand.ON_avg - cand.OFF_avg)
    print(f"  exact Wilcoxon signed-rank ON_avg vs OFF_avg: n={n_} W+={w_} "
          f"p = {p_:.3g}")
    print(f"empirical null log2(FullNeg/NegA): SD {sd:.3f} = {2**sd:.2f}x, "
          f"±2 SD = {2**(2*sd):.2f}x; candidates beyond +2 SD: "
          f"{int((np.log2(r) > 2*sd).sum())}/{len(r)}, beyond -2 SD: "
          f"{int((np.log2(r) < -2*sd).sum())}/{len(r)}")
    negr = cand.NegA / cand.FullNeg
    print(f"the two negatives disagree by median {negr.median():.3f}x "
          f"(range {negr.min():.3f}-{negr.max():.3f}), i.e. the choice of negative "
          f"moves ON/OFF as much as the biology does")
    for p, n in PAIRS:
        v = ratios.loc[cand.index, f"{p}/{n}"]
        print(f"  {p}/{n:8s} >1.0 {int((v>1).sum()):2d}/32  >={CUT} "
              f"{int((v>=CUT).sum()):2d}/32  median {v.median():.3f}")
    a, b = np.log2(cand.PosA / cand.NegA), np.log2(cand.PosB / cand.NegA)
    rp = float(np.corrcoef(a, b)[0, 1])
    rng = np.random.default_rng(0)
    perm = np.array([np.corrcoef(a, rng.permutation(b.values))[0, 1]
                     for _ in range(20000)])
    print(f"do the two positives agree? Pearson {rp:+.3f} over NegA, "
          f"permutation p = {(np.abs(perm) >= abs(rp)).mean():.3f} (20,000 draws)")
    print(f"controls: " + ", ".join(
        f"{CTRL[i]} {ctrl.ONOFF_avg[i]:.2f}x (max candidate = "
        f"{100*r.max()/ctrl.ONOFF_avg[i]:.1f}% of it, "
        f"{100*(r.max()-1)/(ctrl.ONOFF_avg[i]-1):.1f}% of its lift)"
        for i in ctrl.index))
    print(f"routes: {cand.route.value_counts().to_dict()}")
    print(f"pooled-screen comparison the design was picked on: "
          f"{cand.hts_comparison.value_counts().to_dict()}")
    for c in sorted(cand.hts_comparison.unique()):
        s = cand[cand.hts_comparison == c]
        print(f"  {c}: PosA/NegA median {(s.PosA/s.NegA).median():.3f}  "
              f"PosB/NegA median {(s.PosB/s.NegA).median():.3f}  (n={len(s)})")
    lo, hi = xlim
    r4 = ratios.loc[cand.index]
    n_off = int((r4 < lo).sum().sum() + (r4 > hi).sum().sum())
    print(f"(b)/(c) x axis {lo:.3f}-{hi:.3f}, set by the 32 candidates alone: "
          f"controls ({', '.join(f'{ctrl.ONOFF_avg[i]:.2f}x' for i in ctrl.index)}) "
          f"printed in their rows, not plotted; {n_off}/{r4.size} pairing points "
          f"and {int((~cand.ONOFF_avg.between(lo, hi)).sum())}/{len(cand)} ON/OFF "
          f"dots off-scale, each drawn clamped on the edge with its value "
          f"(rows {', '.join(sorted(cand.label[((r4 < lo) | (r4 > hi)).any(axis=1)]))})")
    w = pd.concat([cand, ctrl])[TARGETS]
    print(f"heatmap clip ±{vmax:g} log2; absolute normalized mScarlet over the "
          f"128 candidate wells {w.loc[cand.index].min().min():.5f}-"
          f"{w.loc[cand.index].max().max():.5f} "
          f"({w.loc[cand.index].max().max()/w.loc[cand.index].min().min():.0f}x)")


def layout_report(out, fit, box_b, box_c):
    """The layout claims, checked on the three files just written.

    Everything here is a number the README quotes, so it is printed rather than
    remembered, and the three asserts are the ones that break the page: (b) and
    (c) stop sitting side by side, (a) stops holding its own row, or the stack
    runs off A4.
    """
    (wa, ha), (wb, hb), (wc, hc) = (emitted_mm(out, s) for s in (
        "panel_a_wells_by_target", "panel_b_onoff_ranked", "panel_c_pairings"))
    pitch, gap, rot = fit
    assert max(wb, wc) <= BC_BUDGET, (
        f"(b)/(c) emit {wb:.1f}/{wc:.1f} mm against a {BC_BUDGET:.0f} mm "
        f"budget - the book would stack them and the page would overflow")
    assert wa + BOOK_GUTTER + wb > BOOK_AVAIL, (
        f"(a) at {wa:.1f} mm would share its book row with (b) at {wb:.1f} mm")
    assert wa <= BOOK_AVAIL, f"(a) at {wa:.1f} mm forces the book to scale the page"
    page = BOOK_MARGIN + BOOK_TITLE + ha + BOOK_GUTTER + max(hb, hc) + BOOK_GUTTER
    assert page <= A4_H - BOOK_MARGIN, (
        f"the two rows run {page:.0f} mm, past the {A4_H - BOOK_MARGIN:.0f} mm "
        f"A4 page body")
    print(f"layout: (a) {wa:.1f} x {ha:.1f} mm horizontal, own book row; "
          f"(b) {wb:.1f} x {hb:.1f} and (c) {wc:.1f} x {hc:.1f} side by side "
          f"({wb + BOOK_GUTTER + wc:.1f} of {BOOK_AVAIL:.0f} mm); page body "
          f"{page:.0f} of {A4_H - BOOK_MARGIN:.0f} mm")
    print(f"  (b)/(c) plot box {box_b[0]:.1f} x {box_b[1]:.1f} and "
          f"{box_c[0]:.1f} x {box_c[1]:.1f} mm, w/h {box_b[0]/box_b[1]:.2f} / "
          f"{box_c[0]/box_c[1]:.2f} (was 0.63 / 0.57 at 66 / 62 mm wide)")
    print(f"  (a) column pitch {pitch:.2f} mm, tightest label gap {gap:.2f} mm; "
          f"32 binder numbers upright (widest ink 2.02 mm, measured), the two "
          f"control doses rotated 90 deg ({rot:.1f} mm, over the pitch even "
          f"split from \"{CTRL_NAME}\", which alone is "
          f"{style.text_w(CTRL_NAME) * PT_MM:.1f} mm)")


def main():
    style.apply()
    plt.rcParams.update({"axes.labelsize": 6.0})
    out = paths.fig_output(FIG)
    cand, ctrl = load()
    sd = null_sd(cand)
    ratios = pairings(pd.concat([cand, ctrl]))
    # ONE x range for (b) and (c): they plot the same quantity on rows that
    # line up, so a reader reads straight across the two panels.
    xlim = xlimits(ratios.loc[cand.index], cand.ONOFF_avg)
    vmax, _, fit = panel_a(cand, ctrl, out)
    box_b = panel_b(cand, ctrl, ratios, sd, xlim, out)
    box_c = panel_c(cand, ctrl, ratios, sd, xlim, out)
    report(cand, ctrl, ratios, sd, vmax, xlim)
    layout_report(out, fit, box_b, box_c)


if __name__ == "__main__":
    main()
