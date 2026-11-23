"""Fig3d: how the two selection rules differ in what they pick.

  b  fold diversity: how many of the pool's fold clusters each rule holds,
     against an equal-sized random draw from the same pool.
  c  what each rule does to the stack mix: the share of each set coming from
     each generator, campaign against the two selections.

Cluster count against a random-draw null is the diversity convention in the
structure-generation literature (Foldseek clusters at TM 0.5). It replaced two
attempts at the same question through secondary structure, both of which died
on the same fact: 74% of the library and 84% of ipTM have no strand at all.
Helix content could not separate beta-rich from floppy, and an alpha/beta
balance axis turned that 74/84% into a point mass the smoother drew as a wall.
Fold clusters sidestep secondary structure altogether. The class split is
still printed, since binning is the honest way to show a point mass.

Cut on 2026-11-11, in git history if they are wanted back:
  a  the map of which rule holds each backbone, over generator territories
  c  the log2 abundance ratio per fold cluster, by dominant generator
and with them the difference-in-density strip under the violin. The helix
violin that stood where b is went the same way on 2026-11-12, replaced by the
balance one; the column chart took the lower-right slot the old c had. The
left column is still deliberately empty: something new goes there.

Warm is library, cool is ipTM.

NO P-VALUES. Both arms are deterministic selections from one pool, so there is
no sampling error to test. Every number here is descriptive.

Data: inputs/counterfactual_k30.csv, inputs/stackmix_k30.csv,
      inputs/funnel_real_gate_v2.csv, inputs/sse_arms.csv.
Both k30 files are built by figures/k30_fig3c_inputs.py on Euler, from
library_selection/hier_k30_max_selections.parquet and v6k30_pool.parquet.
Writes output/Fig3d_selection_difference/selection_difference.{svg,pdf,png}.
Nothing is written to final/ - assembly exports there.
"""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import numpy as np
import pandas as pd
import yaml
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from common import paths, palette

FIG = Path(__file__).resolve().parent.name
out = paths.fig_output(FIG)
_cfg = yaml.safe_load((paths.ROOT / "config" / "style.yaml").read_text())
plt.rcParams.update({
    "font.family": _cfg["font_family"], "pdf.fonttype": 42, "svg.fonttype": "none",
    "savefig.dpi": 600, "axes.grid": False,
})
def mm(x): return x / 25.4
palette.register_colormaps()

CMAP = matplotlib.colormaps["thesis_div"]
INK = palette.INK
LIB_C, TOP_C = CMAP(0.88), CMAP(0.12)

POOL_C = palette.NEUTRAL["400"]

# ------------------------------------------------------- fold coverage
# How many of the pool's fold clusters each rule ends up holding, against what
# an equal-sized random draw from the same pool would hold. Cluster count is
# the diversity metric the structure-generation literature uses (Foldseek
# clusters at TM 0.5 in Proteina and Protein FID); the random draw is the null
# that turns a bare count into a statement.
#
# PARTITION, rebuilt 2026-11-16. This panel is now on the k30 partition
# (foldseek_map/results_sw_k30_f50_mut, TM >= 0.50 mutual k=30, Leiden r8.0),
# which is the clustering that actually drove the library selection: 55,706
# gated designs on 31,763 backbones in 917 folds. It replaces the campaign-wide
# MAP partition (139 folds, counterfactual3.csv) which could not answer the
# question, because it was built over all 6,992,515 designs and carried only
# 71% of the ipTM arm. Here both arms are exactly 12,000 designs drawn from
# exactly that pool by library_selection/hier_k30_max_selections.parquet, so
# the comparison is like for like. Do not reintroduce counterfactual3.csv.
#
# THE SPLIT IS THE PANEL. The 917 folds are 256 multi-member clusters plus 661
# single-backbone ones, and the two arms behave completely differently on the
# two. Aggregated, ipTM holds 326/917 and looks far below chance; split, it is
# AT chance on populated folds (206 against a null of 210 +/- 5) and below it
# only on isolated backbones (120 against 194 +/- 11). Never draw the aggregate
# 917 on its own: it says something the data does not support.
#
# The null is DESIGN-matched (300 draws of 12,000 designs from the pool's
# designs, seed 0, computed in figures/k30_fig3c_inputs.py), not
# backbone-matched. That is what lets the backbone row carry a null bar at all:
# a backbone-matched draw of N backbones holds N backbones and only restates
# its own size.
_k = pd.read_csv(paths.INPUTS / "counterfactual_k30.csv").set_index("arm")
POOL_BB = int(_k.loc["pool", "backbones"])
POOL_MULTI = int(_k.loc["pool", "multi"])
POOL_SINGLE = int(_k.loc["pool", "single"])
POOL_FOLDS = int(_k.loc["pool", "folds"])
assert POOL_MULTI + POOL_SINGLE == POOL_FOLDS

_ROW = {"ipTM": "ipTM", "lib": "library", "rand": "random"}
METRICS = ("designs", "backbones", "multi", "single", "folds")
cov = {}
for _a, _name in _ROW.items():
    _r = _k.loc[_name]
    cov[_a] = {m: float(_r[m]) for m in METRICS}
    cov[_a].update({"sd_" + m: float(_r["sd_" + m])
                    for m in ("backbones", "multi", "single")})

for _a, _lab in (("lib", "library"), ("ipTM", "ipTM-ranked"), ("rand", "random")):
    v = cov[_a]
    print(f"{_lab:<12} {v['designs']:>6,.0f} designs on {v['backbones']:>8,.0f} of "
          f"{POOL_BB:,} backbones ({100*v['backbones']/POOL_BB:.0f}%); folds "
          f"{v['multi']:>5.0f}/{POOL_MULTI} populated "
          f"({100*v['multi']/POOL_MULTI:.0f}%), {v['single']:>5.0f}/{POOL_SINGLE} "
          f"one-member ({100*v['single']/POOL_SINGLE:.0f}%), "
          f"{v['folds']:>5.0f}/{POOL_FOLDS} total")
print("library coverage of both fold kinds is 100% BY CONSTRUCTION: the selection "
      "is a fold round-robin and 12,000 picks over 917 folds guarantees it.")

# ------------------------------------------------- per-stack composition
# Canonical stack order (COLORS.md rule 2), short tick labels as the rest of
# the deck writes them.
STACKS = [palette.SHORT[k] for k in palette.STACK_ORDER]
TICK = {"RFD3+LiMPNN": "RFD3", "BG+LiMPNN": "BG+LiM", "BG+BIF": "BG+BIF",
        "BD1+LiMPNN": "BD1", "HD+LiMPNN": "HD"}
_f = pd.read_csv(paths.INPUTS / "funnel_real_gate_v2.csv").set_index("step")
_s = pd.read_csv(paths.INPUTS / "stackmix_k30.csv").set_index("stack")
cnt = pd.DataFrame({
    "all": _f.loc["All designs", STACKS].astype(float),
    "gate": _s["pool"].reindex(STACKS).astype(float),
    "ipTM": _s["ipTM"].reindex(STACKS).astype(float),
    "lib": _s["library"].reindex(STACKS).astype(float),
}).astype(float)
# The gate column of stackmix_k30 IS the funnel's last step, checked here so the
# two files cannot drift apart silently.
assert (cnt["gate"].values == _f.iloc[-1][STACKS].astype(float).values).all(), \
    "k30 pool and the funnel gate row disagree"
assert not cnt.isna().any().any(), "a stack is missing from one of the sets"
TOT = cnt.sum().astype(int).to_dict()
share = {k: (cnt[k] / cnt[k].sum()).values for k in cnt}
print("per-stack share of each set (%):")
print((100 * cnt / cnt.sum()).round(1).to_string())

# ------------------------------------------------- SSE class, not drawn
# The composition companion to b, kept as a printed check because it is the
# obvious alternative panel and the README quotes it. Classes are the CATH
# split: no strand at all, no helix at all, or both.
#
# This is also the honest form of the alpha/beta balance violin that stood
# here until 2026-11-12. That violin was unusable: 74% of the library and 84%
# of ipTM have zero strand, so on a (helix-strand)/(helix+strand) axis they
# are one point mass at +1 and the smoother turns it into a wall. Binning
# into classes puts that population in a bar, where a point mass belongs.
sse = pd.read_csv(paths.INPUTS / "sse_arms.csv")
ss = sse[sse.sel.isin(["library", "top_iptm"])]
klass = np.where(ss.strand < 0.05, "all-alpha",
                 np.where(ss.helix < 0.05, "all-beta", "alpha+beta"))
mix = pd.crosstab(ss.sel, klass, normalize="index").mul(100).round(1)
print("SSE class, % of each arm (not drawn):")
print(mix.to_string())


# =========================================================== figure
fig = plt.figure(figsize=(mm(178), mm(100)))

# --------------------------------------- b what each rule reduces
# Two reductions side by side, as a share of the pool each rule drew from:
# how many backbones survive, and how many distinct fold clusters they still
# cover. The point is the gap between the two rows - a rule can keep a quarter
# of the backbones and most of the folds, or a quarter and barely half.
bx = fig.add_axes([0.545, 0.545, 0.430, 0.365])
# Three rows, not two. The fold row is split into populated clusters and
# one-member clusters because the two arms separate on exactly that line, and
# an aggregate fold count hides it. See the PARTITION note above.
GROUPS = [("unique backbones", "backbones", POOL_BB, "{:,.0f}"),
          ("populated fold clusters (\u22652 backbones)", "multi", POOL_MULTI, "{:,.0f}"),
          ("one-member fold clusters", "single", POOL_SINGLE, "{:,.0f}")]
# Arm names, used by both panels. Deliberately parallel: both arms are
# selections of the same size from the same pool, differing only in the rule.
# "(ordered)" marks the one that was actually synthesised, since otherwise
# "diversity-optimised" reads as a hypothetical arm rather than the library.
ARMS = [("ipTM-ranked 12,000", "ipTM", TOP_C),
        ("diversity-optimised 12,000 (ordered)", "lib", LIB_C)]
RAND = ("random draw", "rand", POOL_C)
ALL_ARMS = ARMS + [RAND]
BH = 0.26
# Group labels sit above their own pair, inside the axes. As y tick labels
# they overhung the axes to the left and would collide with whatever
# eventually fills the left column.
_bases, _y = [], 0.0                     # built bottom-up
for _ in GROUPS:
    _bases.insert(0, _y)
    _y += len(ALL_ARMS) * BH + 0.55
for gi, (gname, key, pool_n, fmt) in enumerate(GROUPS):
    base = _bases[gi]
    bx.text(0, base + len(ALL_ARMS) * BH + 0.02,
            f"{gname}, {fmt.format(pool_n)} in the pool",
            fontsize=5.6, color=INK, ha="left", va="bottom", weight="bold")
    for ai, (alab, akey, col) in enumerate(ALL_ARMS):
        v = cov[akey][key]
        y = base + (len(ALL_ARMS) - 1 - ai) * BH
        bx.barh(y, 100 * v / pool_n, BH * 0.86, color=col, linewidth=0, zorder=2)
        if akey == "rand":               # the null: its own label, and a +/- 1 SD whisker
            sd = cov["rand"]["sd_" + key]
            bx.errorbar(100 * v / pool_n, y, xerr=100 * sd / pool_n, fmt="none",
                        ecolor=palette.NEUTRAL["600"], elinewidth=0.6, capsize=1.3,
                        capthick=0.6, zorder=3)
            txt = f"{v:,.0f} \u00b1 {sd:.0f} at random"
            pad = 1.2 + 100 * sd / pool_n
        else:
            txt = f"{fmt.format(v)}  ({100*v/pool_n:.0f}%)"
            pad = 1.2
        bx.text(100 * v / pool_n + pad, y, txt,
                fontsize=5.4, color=col, ha="left", va="center", weight="bold")
for side in ("top", "right", "left"):
    bx.spines[side].set_visible(False)
bx.spines["bottom"].set_color(palette.NEUTRAL["400"])
bx.spines["bottom"].set_linewidth(0.6)
bx.set_xlim(0, 128)   # room for the two "(100%)" labels in the fold rows
bx.set_ylim(-0.22, _bases[0] + len(ALL_ARMS) * BH + 0.42)
bx.set_yticks([])
bx.set_xticks([0, 50, 100])
bx.set_xticklabels(["0", "50", "100%"], fontsize=5.6,
                   color=palette.NEUTRAL["600"])
bx.tick_params(axis="x", length=2, width=0.6, pad=1.5,
               color=palette.NEUTRAL["400"])
bx.tick_params(axis="y", length=0, pad=2)
bx.set_xlabel("kept, as a share of the pool", fontsize=5.8,
              color=palette.NEUTRAL["600"], labelpad=1.5)
# Key on one line ABOVE the axes: stacked inside, its second row landed on
# the first group label.
_kx = 0.0
for alab, akey, col in ARMS:
    bx.text(_kx, 1.10, "■", transform=bx.transAxes, fontsize=5.5, color=col,
            ha="left", va="center")
    bx.text(_kx + 0.035, 1.10, alab, transform=bx.transAxes, fontsize=5.2,
            color=INK, ha="left", va="center")
    _kx += 0.035 + 0.0125 * len(alab) + 0.045

# ------------------------------------------- c what each rule does to the mix
# Share of each set that comes from each stack. Shares, not counts: the three
# sets are 55,706 / 12,000 / 12,000 designs, so counts would only show that
# the campaign is big. The campaign is near flat by construction - the five
# stacks were run to roughly equal size - so any departure from flat here is
# the selection rule talking, not the campaign.
#
# Stripped on 2026-11-12: no section title, no campaign note under the axis,
# no n= in the key. The bars, their numbers and the three names are the whole
# panel; the counts and the campaign baseline live in the README and the
# caption. The axes grew into the space all of that used to take.
cx = fig.add_axes([0.545, 0.105, 0.430, 0.320])
x = np.arange(len(STACKS))
W = 0.27
for i, (lab, vals, col) in enumerate([
        ("passed the gate", 100 * share["gate"], POOL_C),
        (ARMS[0][0], 100 * share["ipTM"], TOP_C),
        (ARMS[1][0], 100 * share["lib"], LIB_C)]):
    xs = x + (i - 1) * W
    cx.bar(xs, vals, W * 0.92, color=col, linewidth=0, zorder=2)
    for xx, v in zip(xs, vals):
        cx.text(xx, v + 1.0, f"{v:.0f}", fontsize=4.8, color=palette.NEUTRAL["600"],
                ha="center", va="bottom")
    cx.text(0.44, 0.97 - i * 0.085, "■", transform=cx.transAxes, fontsize=5.5,
            color=col, ha="left", va="center")
    cx.text(0.475, 0.97 - i * 0.085, lab, transform=cx.transAxes, fontsize=5.2,
            color=INK, ha="left", va="center")
for side in ("top", "right", "left"):
    cx.spines[side].set_visible(False)
cx.spines["bottom"].set_color(palette.NEUTRAL["400"])
cx.spines["bottom"].set_linewidth(0.6)
cx.set_yticks([])
cx.set_xticks(x)
cx.set_xticklabels([TICK[s] for s in STACKS], fontsize=5.6, color=INK)
cx.tick_params(length=0, pad=2)
cx.set_xlim(-0.6, len(STACKS) - 0.4)
cx.set_ylim(0, 56)
cx.set_ylabel("% of the designs in each set", fontsize=5.8,
              color=palette.NEUTRAL["600"], labelpad=2)

for ext in ("svg", "pdf", "png"):
    fig.savefig(out / f"selection_difference.{ext}", transparent=False,
                dpi=600 if ext == "png" else None)
plt.close(fig)
print("->", out / "selection_difference.png")
