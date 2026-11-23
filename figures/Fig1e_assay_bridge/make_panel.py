#!/usr/bin/env python3
"""Fig 1 FINAL panel: the target-presentation readouts, both assays on one axis pair.

Six bars: flow on the left axis, ONT on the right, so the fluorescence screen and
the splicing screen can be read against one another.

2026-11-17 (Leo): "the panel e to be the final panel showing both the flow
cytometry of the f30 and the ont splicing ratio readouts showed against it", and
the scaffold bars from the old panel f merge in here.

That merge REMOVES A DUPLICATION rather than just relocating ink. This panel's
two ONT bars were LV1.27.B.005 (pART242, "12 hairpins", 0/267) and LV1.27.B.001
(pART158, F30, 117/300). Those are the same two samples, from the same LV1.27.1
run, that the scaffold panel plotted as its Neg and F30 bars - the figure drew
them twice. The ONT side now carries all four scaffolds once, and the assertion
below fails the build if the two source tables ever stop agreeing.

  flow  normalized mScarlet, column P6_C5 of LV2_34_1 (the csy4 control),
        FullNeg vs PosB -> 91x
  ONT   splicing ratio (removed/(intact+removed)) from LV1.27.1, all four
        scaffolds, ascending: Neg (pART242, 12 sites) -> F30 variation (mod1,
        8 sites) -> Circular scaffold (mod4, 8 sites) -> F30 (pART158, 16 sites).
        The fold change is still read Neg -> F30, the two extremes.

No micrographs. An earlier version carried a row of image slots over the flow
pair; they are gone on purpose (STYLE_NOTES section 0 - a placeholder is ink
that carries nothing). If the flow images are wanted later they belong beside
this panel at assembly, not inside it.

NOTE ON THE TWO AXES: a second y-scale normally lets the author manufacture
agreement by choosing the relative scaling. That is avoided here by giving both
axes the SAME numeric range (0 to YTOP). The two quantities are different
things - a fluorescence ratio and a fraction of transcripts cleaved - so they
keep separate axes, but no stretching is applied.

2026-11-17 (Leo): the AXES are now neutral ink, not assay-coloured. The colour
that keys each axis to its pair is carried by the bars and by the fold-change
bracket over each pair, which is enough; two coloured spines on a four-bar
panel read as decoration. Only AXIS furniture lost its colour - bars and
brackets keep FLOW / ONT.

NOTE ON THE ONT FOLD CHANGE: the negative gave zero cleaved reads, so the ratio
is undefined, not large. A +1 pseudocount gives the negative a floor of 1/(n+1),
and the positive over that floor is 105x - but because that floor is an upper
bound on a rate never actually observed, the panel writes "> 105x", not "105x"
(Leo, 2026-11-17). The flow pair keeps a plain "91x": its negative is a real,
non-zero measurement.
"""
import sys
from pathlib import Path

import csv

import pandas as pd
import matplotlib.pyplot as plt

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from common import palette, paths, style                        # noqa: E402

YTICKS = [0.0, 0.1, 0.2, 0.3, 0.4]   # the spine is cut at the last one
# The bracket stack is spaced in MILLIMETRES and converted, not taken as a
# fraction of the data range. YTOP is then solved so the stack fits, once the
# bars are known - see "bracket geometry" below.
BR_GAP_MM = 3.2      # clear air between one bracket level and the next
BR_LIFT_MM = 2.2     # top of a bar's VALUE LABEL to the bracket's descending tick
BR_PAD_MM = 1.6      # above the topmost label
BR_TICK_MM = 0.7     # how far the bracket's end ticks hang below its line
VAL_GAP_MM = 1.0     # bar top (and its dot marker) to the value label
# Every bar carries a value label, so a bracket lifted from the BAR lands on
# that label - the tick came down through "0.393" and the dot sat on the digits
# (Leo, 2026-11-18: "the bar is clashing"). Lift brackets from the top of the
# label instead, and give the dot its own clearance.
# 130 not 100: the ONT side went from 2 bars to 4 (2026-11-17).
FIG_W_MM, FIG_H_MM = 130.0, 79.0
# A size variant for page assembly; the canonical panel is untouched:
#   make_panel.py W_MM H_MM FONT_PT OUTDIR
VARIANT = len(sys.argv) == 5
if VARIANT:
    FIG_W_MM, FIG_H_MM = float(sys.argv[1]), float(sys.argv[2])
BAR_RECT = (0.105, 0.175, 0.80, 0.720)      # left, bottom, w, h in figure coords
FLOW, FLOW_LT = palette.THEME["700"], palette.THEME["300"]
ONT, ONT_LT = palette.DIVERGING[-1], palette.DIVERGING[-3]
AXIS = palette.INK            # axis furniture is neutral; the bars carry the key

style.apply()
if VARIANT:
    plt.rcParams["font.size"] = float(sys.argv[3])
out = Path(sys.argv[4]) if VARIANT else paths.fig_output("Fig1e_assay_bridge")
d = pd.read_csv(paths.INPUTS / "csy4_control_both_assays.csv")
v = {(r.assay, r.role): r for r in d.itertuples()}

flow_fc = v[("flow", "positive")].value / v[("flow", "negative")].value
ont_n = int(v[("ont", "negative")].n)
# 2026-11-16 (Leo): the negative gave ZERO cleaved reads, so positive/negative
# is undefined and used to be reported as ">=35x" from the rule of three (a 3/n
# upper bound on the unseen rate). Now a POINT ESTIMATE from a +1 pseudocount:
# one cleaved read is added to the negative, which adds one to its total too, so
# its floor is 1/(n+1) instead of 0. That is a convention for dividing by zero,
# not a measurement - the negative is still 0/267 - but it gives one honest
# number rather than a bound, and it is the fraction the +1 actually implies.
ont_neg = 1.0 / (ont_n + 1)
ont_fc = v[("ont", "positive")].value / ont_neg

# The four ONT scaffolds, ascending by ratio, from the panel that used to own them.
DISPLAY = {"mod1": "F30 variation", "mod4": "Ring array"}
EXCLUDE = {"pART215"}                       # the 4x negative, never shown
_sc = [r for r in csv.DictReader((paths.INPUTS / "scaffold_splicing.csv").open())
       if r["construct"] not in EXCLUDE]
for r in _sc:
    r["ratio"] = int(r["spliced"]) / (int(r["spliced"]) + int(r["unspliced"]))
_sc.sort(key=lambda r: r["ratio"])          # ascending; Neg = 0 lands first

# The merge is only honest if the two tables describe the same experiment. The
# bridge CSV's ONT pair and the scaffold CSV's Neg / F30 rows are the same two
# samples, so they must agree to the rounding in the bridge table.
_by_con = {r["construct"]: r for r in _sc}
for _role, _con in (("negative", "pART242"), ("positive", "pART158")):
    _a, _b = v[("ont", _role)].value, _by_con[_con]["ratio"]
    assert abs(_a - _b) < 5e-3, (
        f"ONT {_role} disagrees between the two tables: "
        f"csy4_control_both_assays={_a} vs scaffold_splicing={_b} ({_con})")

ONT_X0, ONT_DX = 2.3, 1.0
SLOTS = [("flow", "negative", 0.0, FLOW_LT), ("flow", "positive", 1.0, FLOW)]
SC_X = [ONT_X0 + ONT_DX * i for i in range(len(_sc))]
XLIM = (-0.65, SC_X[-1] + 0.65)


fig = plt.figure(figsize=(style.mm(FIG_W_MM), style.mm(FIG_H_MM)))
ax = fig.add_axes(BAR_RECT)    # no gridlines: STYLE_NOTES section 7
X = [s[2] for s in SLOTS] + SC_X
vals = [v[(s[0], s[1])].value for s in SLOTS] + [r["ratio"] for r in _sc]
# ONT bars all take the ONT axis colour: the axis is colour-matched to its bars
# (STYLE_NOTES, dual axes), and the four scaffolds are told apart by their tick
# labels, not by hue. The floor bar is the light tint, the scaffolds the full one.
cols = [s[3] for s in SLOTS] + [ONT_LT if r["ratio"] == 0 else ONT for r in _sc]
# ---- bracket geometry -------------------------------------------------------
# The stack has to clear a LINE OF TYPE. Tying its step to the data range is
# what made every label strike through the line above it once the panel was
# drawn at the assembly's height instead of its own (Leo, 2026-11-17: "the
# brackets in panel f look bad"). Solve YTOP so the whole stack fits:
#   YTOP = bar / (1 - head_mm / axes_mm)
AX_H_MM = FIG_H_MM * BAR_RECT[3]
LINE_MM = plt.rcParams["font.size"] * 25.4 / 72 * 1.25
N_BR = len(_sc) - 1                       # one ONT bracket per scaffold
# A bracket's horizontal line sits this far above the bar it spans: clear the
# bar's own value label (VAL_GAP + one line of type), then the air, then the
# tick that hangs back down off the line.
LABEL_MM = VAL_GAP_MM + LINE_MM           # bar top -> top of its value label
BR_BASE_MM = LABEL_MM + BR_LIFT_MM + BR_TICK_MM
HEAD_MM = BR_BASE_MM + (N_BR - 1) * BR_GAP_MM + LINE_MM + BR_PAD_MM
assert HEAD_MM < AX_H_MM, "the bracket stack does not fit the axes"
YTOP = max(vals) / (1 - HEAD_MM / AX_H_MM)
U = YTOP / AX_H_MM                        # data units per mm

ax.bar(X, vals, width=0.52, color=cols, zorder=3)
for x, val in zip(X, vals):
    ax.text(x, val + VAL_GAP_MM * U, f"{val:.3f}", ha="center", va="bottom", color=palette.INK)

# One ONT bracket PER SCAFFOLD over the negative, stacked, not a single bracket
# over the two extremes (Leo, prio-2: "add FC brackets to the columns"). Every
# scaffold has its own fold change over the floor; the extremes bracket reported
# only F30's.
#
# All of them are written "> N x", not "N x": the negative gave zero cleaved
# reads, so 1/(n+1) is a FLOOR on its rate and every ratio built on it is a
# lower bound, not a point estimate (Leo, 2026-11-17). Flow keeps a plain "91x"
# - its negative is a real, non-zero measurement.
BR_STEP = BR_GAP_MM * U
_ont_top = max(r["ratio"] for r in _sc)
_brackets = [((X[0], X[1]), f"{flow_fc:.0f}×", FLOW,
              max(vals[0], vals[1]) + BR_BASE_MM * U)]
for _i in range(1, len(_sc)):
    _brackets.append(((SC_X[0], SC_X[_i]),
                      f"> {_sc[_i]['ratio'] / ont_neg:.0f}×", ONT,
                      _ont_top + BR_BASE_MM * U + (_i - 1) * BR_STEP))
# every bracket's lowest ink must clear the value label of the tallest bar it spans
for (a, b), _t, _c, _y in _brackets:
    _tallest = max(v_ for x_, v_ in zip(X, vals) if a - 1e-9 <= x_ <= b + 1e-9)
    assert _y - BR_TICK_MM * U >= _tallest + LABEL_MM * U - 1e-9, \
        f"bracket '{_t}' tick lands on the {_tallest:.3f} bar's value label"
for (a, b), txt, col, y in _brackets:
    ax.plot([a, a, b, b], [y - BR_TICK_MM * U, y, y, y - BR_TICK_MM * U], lw=0.75, color=col)
    ax.text((a + b) / 2, y + 0.35 * U, txt, ha="center", va="bottom", color=col)

# One dot per bar, on the bar's own value: every bar here is a single
# measurement - one ONT run per scaffold, one flow column per control - not a
# triplicate mean (Leo, prio-2; flow pair dotted on his say-so 2026-11-17).
ax.plot(X, vals, linestyle="none", marker="o",
        markersize=2.2, color=palette.INK, zorder=4)

ax.set_xticks(X, ["negative\nno array", "positive\nPosB"] +
              ["\n".join(DISPLAY.get(r["label"], r["label"]).split(" ")
                         + [f'{r["hairpin_sites"]} hairpins']) for r in _sc])
ax.set_ylim(0, YTOP)
ax.set_yticks(YTICKS)
ax.spines["left"].set_bounds(YTICKS[0], YTICKS[-1])
ax.set_xlim(*XLIM)
ax.set_ylabel("Normalized mScarlet", color=AXIS)
ax.tick_params(axis="y", colors=AXIS)
ax.spines["left"].set_color(AXIS)

r = ax.twinx()                      # same range on both: no relative stretching
r.set_ylim(0, YTOP)
r.set_yticks(YTICKS)
r.spines["right"].set_bounds(YTICKS[0], YTICKS[-1])
r.set_ylabel("Splicing ratio (fraction cleaved)", color=AXIS, rotation=270,
              labelpad=14 * plt.rcParams["font.size"] / 7)   # 14 at the 7 pt default
r.tick_params(axis="y", colors=AXIS)
r.spines["right"].set_color(AXIS)
r.spines["right"].set_visible(True)
r.spines["top"].set_visible(False)

# method labels head their pair, above the panel, not under the tick labels
for x, lab in ((0.5, "Flow cytometry"), ((SC_X[0] + SC_X[-1]) / 2, "ONT splicing")):
    ax.text(x, 1.075, lab, transform=ax.get_xaxis_transform(), ha="center",
            va="bottom", color=palette.NEUTRAL["700"], fontweight="bold")

style.save(fig, out, "Fig1e_assay_bridge")
# SVG as well: this panel is assembled with the construct schematic in the deck
fig.savefig(out / "Fig1e_assay_bridge.svg", bbox_inches="tight")
plt.close(fig)
print(f"flow {flow_fc:.1f}x   ONT {ont_fc:.0f}x "
      f"(0/{ont_n} in the negative, floor 1/{ont_n + 1} = {ont_neg:.5f})")
print("  ONT bars:", [(DISPLAY.get(r["label"], r["label"]), round(r["ratio"], 3))
                      for r in _sc])
