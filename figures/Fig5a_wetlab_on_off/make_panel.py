#!/usr/bin/env python3
"""Fig 5 panel B: arrayed-flow ON/OFF of the 32 ordered candidates.

Five ranked dotplots, one per positive/negative pair:
  Fig5a_PosA_vs_NegA       Fig5a_PosA_vs_FullNeg
  Fig5a_PosB_vs_NegA       Fig5a_PosB_vs_FullNeg
  Fig5a_PosAvg_vs_NegAvg    <- avg(PosA, PosB) / avg(NegA, FullNeg)
plus three PosA-vs-PosB scatters, one per negative:
  Fig5a_corr_PosA_PosB_{NegA,FullNeg,NegAvg}
and one scatter of the two negatives against the same positive, with the
count of binders each negative calls ON beside it:
  Fig5a_corr_PosA_NegA_vs_FullNeg
  Fig5a_above1_by_negative_gate{1,1.2,1.4}
and the full pass/fail matrix over both positives x both negatives:
  Fig5a_pass_matrix_gate{1,1.2,1.4}

Values are normalized mScarlet (mScarlet mean x %mScarlet+ / BFP mean x %BFP+),
one well per condition, from experiment LV2_34_1; PosA/PosB/NegA/FullNeg are
plate rows A-D and plate column 1-6 is ordered binder 1-6, so rows carry the
Twist/ID_Number. Dot colour is the generator stack.

The two csy4-binder control columns (P6_C4, P6_C5) are EXCLUDED throughout:
they reach 10-120x and are not ordered designs. Scales therefore come from the
32 candidates alone. The one near-zero ratio in each pair is still clamped to
the edge carrying its true value, because it comes from a failed Pos well
(binder 4 PosB = 0.006 vs NegA 0.088) rather than from the binder.
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import FuncFormatter

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from common import palette, paths, style                        # noqa: E402

MIN_KEEP = 0.5      # below this a ratio is a dead well, not an inactive design
# The binder call. avg Pos / avg Neg >= 1.2 is the flow ground truth for the
# ordered 32 - 15 proposed binders, 17 inactive designs - and the cut Fig 5 d-f are built
# at. It is defined on the averaged pair only, so the bands below are drawn
# on that panel alone and not on the single-condition ones.
# Naming, both halves hedged on purpose (Leo, 2026-11-15): "proposed binders",
# because one well per condition and a cut nothing in this experiment fixes is
# a proposal, not a confirmed binder; and "inactive designs", never "non-binders",
# because 1.0-1.19 is weak or absent signal, not evidence of no binding.
CUT = 1.2
BANDS = ("Inactive designs", "Proposed binders")
TICKS = (0.05, 0.1, 0.2, 0.3, 0.5, 0.7, 1, 1.5, 2, 3, 5, 10)
PAIRS = [("PosA", "NegA"), ("PosA", "FullNeg"),
         ("PosB", "NegA"), ("PosB", "FullNeg"),
         ("PosAvg", "NegAvg")]
NICE = {"PosAvg": "avg Pos", "NegAvg": "avg Neg"}

style.apply()
# smaller than the 7 pt base: these panels carry 32 tick labels in one column
LABEL_PT, SMALL_PT = 6.0, 5.5
plt.rcParams.update({"font.size": SMALL_PT, "axes.labelsize": LABEL_PT,
                     "legend.fontsize": SMALL_PT,
                     "xtick.labelsize": SMALL_PT, "ytick.labelsize": SMALL_PT,
                     "axes.titlesize": LABEL_PT})
out = paths.fig_output("Fig5a_wetlab_on_off")


def save(fig, stem):
    """PDF and PNG as everywhere else, plus SVG for hand assembly."""
    style.save(fig, out, stem)
    fig.savefig(out / f"{stem}.svg", bbox_inches="tight")

d = pd.read_csv(paths.INPUTS / "lv2_34_1_flow_on_off.csv", index_col=0)
pc = d.index.str.extract(r"P(\d+)_C(\d+)").astype(int)
d["ID"] = ((pc[0] - 1) * 6 + pc[1]).values
d = d[d.ID <= 32].copy()            # drop the two csy4 control columns
assert len(d) == 32 and d["stack"].notna().all()
d["label"] = d.ID.astype(str)
d["PosAvg"] = d[["PosA", "PosB"]].mean(axis=1)
d["NegAvg"] = d[["NegA", "FullNeg"]].mean(axis=1)

# dot colour means generator stack, as everywhere else in the thesis
d["color"] = d["stack"].map(palette.STACKS)
HANDLES = [Line2D([], [], ls="", marker="o", ms=3.5, mfc=palette.STACKS[k],
                  mec=palette.SURFACE, label=palette.SHORT[k])
           for k in palette.STACK_ORDER if (d["stack"] == k).any()]   # never sort


def clamp_marks(ax, vals, colors, lo, hi, at, horizontal):
    """Draw out-of-range points on the edge, carrying their true value."""
    for pos, v, col in zip(at, vals, colors):
        over = v > hi
        edge = hi if over else lo
        xy = (edge, pos) if horizontal else (pos, edge)
        ax.plot(*xy, marker=(">" if over else "<") if horizontal else ("^" if over else "v"),
                ms=4, c=col, clip_on=False, zorder=4)
        ax.annotate(f"{v:.0f}x" if v >= 10 else f"{v:.2f}x", xy,
                    textcoords="offset points", xytext=(6, 0 if horizontal else 4),
                    ha="left", va="center", color=palette.INK, annotation_clip=False)


for pos, neg in PAIRS:
    s = d.assign(ratio=d[pos] / d[neg]).sort_values("ratio").reset_index(drop=True)
    y = np.arange(len(s))
    kept = s.ratio[s.ratio >= MIN_KEEP]
    lo, hi = kept.min() / 1.35, kept.max() * 1.35

    # 2.4 mm per binder, not 2.9: the row is what sets the a|b height on the
    # page, and 32 labels at 5.5 pt still clear each other at this pitch.
    # 2026-11-18 (Leo): "make the figure b smaller vertically". The PITCH is
    # already at its floor - a 5.5 pt line is about 1.94 mm tall, so 2.4 mm
    # leaves 0.46 mm between the 32 y labels and anything less collides. The
    # height came out of the CONSTANT instead, 26 -> 16 mm, which is the
    # margin the one-row legend above no longer needs.
    # 2026-11-18 (Leo): "resize figure B as well". Fig 5's a|b row gives this
    # panel the width a's matrix no longer needs. Widening it HERE rather than
    # scaling it on the page is the point: scaling a finished SVG to 100 mm
    # would have scaled its height to match, and the row's height is what sets
    # the whole figure against the page - at 100 mm scaled it overflowed and
    # the assembler shrank every panel to 0.959, dropping the type to 5.27 pt.
    # Redrawn at the new width the height formula is untouched.
    w_mm = 100 if (pos, neg) == ("PosAvg", "NegAvg") else 85
    fig, ax = plt.subplots(figsize=(style.mm(w_mm),
                                    style.mm(2.4 * len(s) + 16)))
    ax.grid(axis="x", color=palette.GRID, lw=0.5)
    ax.set_axisbelow(True)
    if (pos, neg) == ("PosAvg", "NegAvg"):
        # Behind the grid (zorder 0). The boundary between the two fills IS the
        # 1.2 cut, so no extra rule is drawn for it. Conditions, not design
        # stacks, so a theme step and a neutral, never a stack colour
        # (COLORS.md rules 2 and 3). The labels straddle the cut along the FOOT
        # of the frame: the dots are ranked, so everything near 1.2 sits
        # mid-frame and the bottom strip is empty on both sides of the line.
        # The outer top corner is not - "Proposed binders" runs into the
        # leading dot there.
        ax.axvspan(lo, CUT, color=palette.NEUTRAL["200"], lw=0, zorder=0)
        ax.axvspan(CUT, hi, color=palette.THEME["100"], lw=0, zorder=0)
        # Each label carries the cut on the side that touches the line, so the
        # threshold is read where it is drawn and needs no tick of its own.
        for x, ha, name in ((CUT / 1.04, "right", f"{BANDS[0]} < {CUT:g}"),
                            (CUT * 1.04, "left", f"\u2265 {CUT:g} {BANDS[1]}")):
            ax.text(x, 0.008, name, transform=ax.get_xaxis_transform(), ha=ha,
                    va="bottom", color=palette.NEUTRAL["600"], zorder=5)
    ax.axvline(1, color=palette.NEUTRAL["500"], lw=0.75, ls="--", zorder=2)
    inside = s.ratio.between(lo, hi)
    ax.scatter(s.ratio[inside], y[inside.values], s=16, c=s.color[inside],
               edgecolors=palette.SURFACE, linewidths=0.5, zorder=3)
    clamp_marks(ax, s.ratio[~inside], s.color[~inside], lo, hi,
                y[~inside.values], horizontal=True)
    ax.set_yticks(y, s.label)
    ax.set_ylim(-0.9, len(s) - 0.1)
    ax.set_xscale("log")
    ax.set_xlim(lo, hi)
    ax.set_xticks([t for t in TICKS if lo <= t <= hi])
    ax.set_xticks([], minor=True)
    ax.xaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{v:g}"))
    ax.set_xlabel(f"ON/OFF   {NICE.get(pos, pos)} / {NICE.get(neg, neg)}")
    ax.set_ylabel("Ordered binder")
    ax.tick_params(axis="y", length=0)
    # 2026-11-18 (Leo): "make the design method labels horizontal". One row
    # of five, not a 2 x 3 block. It also buys back two lines of height, which
    # is half of what the shortening below needs. The five labels measure
    # about 73 mm of the panel's 100 mm, so they fit on one row at 5.5 pt.
    ax.legend(handles=HANDLES, loc="lower center", bbox_to_anchor=(0.5, 1.0),
              ncol=len(HANDLES), frameon=False, handletextpad=0.3,
              columnspacing=1.2)
    fig.tight_layout()
    save(fig, f"Fig5a_{pos}_vs_{neg}")
    plt.close(fig)
    print(f"{NICE.get(pos,pos)+'/'+NICE.get(neg,neg):18s} median {s.ratio.median():.2f}  "
          f"{(s.ratio > 1).sum()}/32 above 1  clamped {(~inside).sum()}")

# --------------------------------------------------------------- scatters --
def scatter(x, y, xlabel, ylabel, stem, subtitle=""):
    """One log-log ON/OFF scatter, 1:1 line, out-of-range points on the edge."""
    keep = (x >= MIN_KEEP) & (y >= MIN_KEEP)
    lo = min(x[keep].min(), y[keep].min()) / 1.35
    hi = max(x[keep].max(), y[keep].max()) * 1.35

    fig, ax = plt.subplots(figsize=(style.mm(85), style.mm(95)))
    ax.grid(color=palette.GRID, lw=0.5)
    ax.set_axisbelow(True)
    # dotted lines at ON/OFF = 1 on both axes: a real binder sits upper right
    ax.axvline(1, color=palette.NEUTRAL["500"], lw=0.75, ls=":", zorder=2)
    ax.axhline(1, color=palette.NEUTRAL["500"], lw=0.75, ls=":", zorder=2)
    ax.plot([lo, hi], [lo, hi], color=palette.GRID, lw=0.75, ls="--", zorder=1)
    ax.scatter(x[keep], y[keep], s=18, c=d.color[keep], edgecolors=palette.SURFACE,
               linewidths=0.5, zorder=3)
    for i in d.index[~keep]:                      # clamp on whichever axis is out
        cx, cy = min(max(x[i], lo), hi), min(max(y[i], lo), hi)
        ax.plot(cx, cy, marker="<" if x[i] < lo else "v", ms=4, c=d.color[i],
                clip_on=False, zorder=4)
        ax.annotate(f"b{d.ID[i]} ({x[i]:.2f}, {y[i]:.2f})", (cx, cy),
                    textcoords="offset points", xytext=(6, 4), fontsize=6,
                    color=palette.INK, annotation_clip=False)
    rho = x.rank().corr(y.rank())
    ax.set_xscale("log"); ax.set_yscale("log")
    ticks = [t for t in TICKS if lo <= t <= hi]
    for axis, setter in ((ax.xaxis, ax.set_xticks), (ax.yaxis, ax.set_yticks)):
        setter(ticks)
        setter([], minor=True)
        axis.set_major_formatter(FuncFormatter(lambda v, _: f"{v:g}"))
    ax.set_xlim(lo, hi); ax.set_ylim(lo, hi)
    ax.set_aspect("equal")
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(f"Spearman {rho:+.2f}   n = {len(d)}{subtitle}",
                 loc="left", color=palette.INK)
    ax.legend(handles=HANDLES, loc="lower center", bbox_to_anchor=(0.5, 1.06),
              ncol=2, frameon=False, handletextpad=0.3, columnspacing=1.2)
    fig.tight_layout()
    save(fig, stem)
    plt.close(fig)
    return rho, (~keep).sum()


def shared_null(num, den_x, den_y, n=4000, seed=0):
    """Correlation left when only the shared term is real.

    Both axes divide (or are divided by) one common column, which spans ~100x
    across binders, so they correlate even with the pairing destroyed. Permute
    the other two together and take the median rho: that is the floor the real
    number has to beat, not zero.
    """
    rng = np.random.default_rng(seed)
    a, b, c = num.to_numpy(), den_x.to_numpy(), den_y.to_numpy()
    out = np.empty(n)
    for k in range(n):
        j = rng.permutation(len(a))
        out[k] = pd.Series(a / b[j]).rank().corr(pd.Series(a / c[j]).rank())
    return float(np.median(out))


# ------------------------------------------------- PosA vs PosB, same negative
# Do the two positive conditions agree on a binder's ON/OFF strength? NOTE both
# axes share that denominator, which spans ~100x across binders and inflates the
# correlation on its own: shuffling the positives while keeping each row's
# negative still gives rho ~ +0.44. Read the printed rho against that, not zero.
for neg in ("NegA", "FullNeg", "NegAvg"):
    rho, clamped = scatter(d.PosA / d[neg], d.PosB / d[neg],
                           f"PosA / {NICE.get(neg, neg)}",
                           f"PosB / {NICE.get(neg, neg)}",
                           f"Fig5a_corr_PosA_PosB_{neg}")
    print(f"corr PosA vs PosB over {NICE.get(neg,neg):8s}: spearman {rho:+.2f}  "
          f"clamped {clamped}")

# ------------------------------------------- the two negatives, same positive
# Does the choice of negative change a binder's ON/OFF? Same numerator on both
# axes this time, so the inflation runs the other way: PosA is shared. The null
# permutes the two negatives together, keeping each row's PosA.
null = shared_null(d.PosA, d.NegA, d.FullNeg)
rho, clamped = scatter(d.PosA / d.NegA, d.PosA / d.FullNeg,
                       "ON/OFF   PosA / NegA", "ON/OFF   PosA / no array",
                       "Fig5a_corr_PosA_NegA_vs_FullNeg",
                       subtitle=f"   shared-PosA null {null:+.2f}")
ratio = (d.PosA / d.FullNeg) / (d.PosA / d.NegA)      # = NegA / FullNeg
print(f"corr NegA vs FullNeg under PosA : spearman {rho:+.2f} "
      f"(shared-PosA null {null:+.2f})  clamped {clamped}")
print(f"  FullNeg-based ON/OFF is {ratio.median():.2f}x the NegA-based one "
      f"(range {ratio.min():.2f}-{ratio.max():.2f}); "
      f"above the 1:1 line for {(ratio > 1).sum()}/32")


# ---------------------------------------------- pass/fail at a chosen gate --
# "Passing" is ON/OFF above a threshold. 1.0 is the no-effect line; 1.2 and 1.4
# ask for a 20% or 40% lift before a binder counts. All three are drawn because
# which one is defensible depends on the replicate spread, and there is one well
# per condition here - nothing in this experiment fixes the threshold.
GATES = (1.0, 1.2, 1.4)
COMP = [("PosA", "NegA"), ("PosA", "FullNeg"),
        ("PosB", "NegA"), ("PosB", "FullNeg")]
RATIO = pd.DataFrame({f"{q}/{n}": d[q] / d[n] for q, n in COMP}, index=d.index)

for GATE in GATES:
    tag = f"gate{GATE:g}"

    # --------------------------- which negative calls a binder ON, PosA only
    # Categories are conditions, not design stacks, so they take theme steps
    # and a neutral for the null case (COLORS.md rules 2 and 3).
    a, b = RATIO["PosA/NegA"] > GATE, RATIO["PosA/FullNeg"] > GATE
    CATS = [("Both", a & b, palette.THEME["700"]),
            ("NegA\nonly", a & ~b, palette.THEME["400"]),
            ("FullNeg\nonly", ~a & b, palette.THEME["200"]),
            ("Neither", ~a & ~b, palette.NEUTRAL["300"])]
    counts = [int(m_.sum()) for _, m_, _ in CATS]
    assert sum(counts) == len(d)

    fig, ax = plt.subplots(figsize=(style.mm(68), style.mm(58)))
    x = np.arange(len(CATS))
    ax.bar(x, counts, width=0.62, color=[c for _, _, c in CATS], lw=0)
    for xi, n_ in zip(x, counts):
        ax.text(xi, n_ + len(d) * 0.03, str(n_), ha="center", va="bottom",
                color=palette.INK)
    ax.set_xticks(x, [name for name, _, _ in CATS])
    ax.set_ylim(0, len(d) * 1.05)
    ax.set_ylabel("Ordered binders")
    ax.set_xlabel(f"ON/OFF > {GATE:g} under")
    ax.grid(axis="y", color=palette.GRID, lw=0.5)
    ax.set_axisbelow(True)
    ax.tick_params(axis="x", length=0)
    fig.tight_layout()
    save(fig, f"Fig5a_above1_by_negative_{tag}")
    plt.close(fig)
    print(f"\nabove {GATE:g} (PosA): " + "  ".join(
        f"{name.replace(chr(10), ' ')} {n_}"
        for (name, _, _), n_ in zip(CATS, counts)))
    for pos in ("PosB", "PosAvg"):
        p_, q_ = d[pos] / d.NegA > GATE, d[pos] / d.FullNeg > GATE
        print(f"above {GATE:g} ({pos}): Both {int((p_ & q_).sum())}  NegA only "
              f"{int((p_ & ~q_).sum())}  FullNeg only {int((~p_ & q_).sum())}  "
              f"Neither {int((~p_ & ~q_).sum())}")

    # ------------------------- every positive x negative combination, PosA+B
    # A filled dot clears the gate, in the generator colour; an open ring does
    # not. Rows are sorted by how many of the four clear it, so the groups read
    # top to bottom and the rules between them are the group edges.
    m = RATIO.copy()
    m["npass"] = (m > GATE).sum(axis=1)
    m["ID"], m["color"] = d.ID, d.color
    m["mean"] = RATIO.apply(np.log).mean(axis=1)
    m = m.sort_values(["npass", "mean"]).reset_index(drop=True)

    fig, ax = plt.subplots(figsize=(style.mm(62), style.mm(2.9 * len(m) + 28)))
    for j, (q, n) in enumerate(COMP):
        ok = m[f"{q}/{n}"] > GATE
        ax.scatter(np.full(ok.sum(), j), np.arange(len(m))[ok.values], s=26,
                   c=m.color[ok], edgecolors=palette.SURFACE, linewidths=0.5,
                   zorder=3)
        ax.scatter(np.full((~ok).sum(), j), np.arange(len(m))[(~ok).values],
                   s=26, facecolors="none", edgecolors=palette.NEUTRAL["400"],
                   linewidths=0.8, zorder=3)
    for k in np.flatnonzero(np.diff(m.npass.to_numpy())) + 1:    # group edges
        ax.axhline(k - 0.5, color=palette.GRID, lw=0.6, zorder=1)
    ax.set_xticks(range(len(COMP)), [f"{q}\n{n}" for q, n in COMP])
    ax.set_yticks(range(len(m)), m.ID.astype(str))
    ax.set_xlim(-0.6, len(COMP) - 0.4)
    ax.set_ylim(-0.7, len(m) - 0.3)
    ax.set_ylabel("Ordered binder")
    ax.tick_params(length=0)
    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.xaxis.set_ticks_position("top")
    ax.xaxis.set_label_position("top")
    ax.legend(handles=HANDLES, loc="lower center", bbox_to_anchor=(0.5, 1.06),
              ncol=2, frameon=False, handletextpad=0.3, columnspacing=1.2)
    fig.tight_layout()
    # below the plot, not a title: the legend already occupies the top strip
    fig.text(0.02, 0.004, f"filled = ON/OFF > {GATE:g}", fontsize=7,
             color=palette.NEUTRAL["600"], va="bottom")
    save(fig, f"Fig5a_pass_matrix_{tag}")
    plt.close(fig)

    pat = (m[[f"{q}/{n}" for q, n in COMP]] > GATE).astype(int).astype(str) \
        .agg("".join, axis=1)
    print(f"pass matrix at {GATE:g}, "
          f"pattern = PosA/NegA PosA/FullNeg PosB/NegA PosB/FullNeg")
    for k, grp in sorted(m.groupby(pat.values), key=lambda t: -int(t[0], 2)):
        print(f"  {k}  n={len(grp):2d}  {sorted(grp.ID.tolist())}")
