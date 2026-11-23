"""Fig2d, summary form: what the campaign -> gate -> library funnel did to the
epitope, as one panel.

Five shares, three stages, drawn as a slopegraph because the x axis is three
ordered stages and the message is direction and magnitude, not level at any one
stage. Two families, kept apart by colour temperature AND line style:

  where on the hairpin (warm, solid)     U10 . the GUAUA loop, nt7-11
  what it touches (petrol, dashed)       base . sugar . phosphate

Both families are shares of the same denominator - a design's protein-RNA
heavy-atom contact pairs within 4.5 A - so they belong on one axis, but they
are not mutually exclusive (U10 is inside the loop; base/sugar/phosphate sum
to 1 across the whole hairpin). The caption has to say so.

The point of the panel: the gate carries no epitope term at all (ipTM,
ligand_RMSD, DockQ, seq_run, foot-clash) yet it concentrated contacts onto the
loop and onto bases, and it took every bit of that out of the sugar while
leaving phosphate untouched.

Population and measurement caveats are the same as make_columns.py in this
folder - backbones not designs, campaign is the 499,999-node sample, design
coordinates not a refold. See the README.

Data: inputs/epitope_layers_union.csv.
Writes output/Fig2d_epitope_funnel/epitope_funnel_summary.{svg,pdf,png}.
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
    "font.family": _cfg["font_family"], "font.size": _cfg["font_size_pt"],
    "pdf.fonttype": 42, "svg.fonttype": "none", "savefig.dpi": 600,
    "axes.grid": False, "text.color": palette.INK,
    "axes.labelcolor": palette.INK, "axes.edgecolor": palette.INK,
    "xtick.color": palette.INK, "ytick.color": palette.INK,
    "axes.linewidth": _cfg["line_width_pt"],
    "axes.spines.top": False, "axes.spines.right": False,
})
def mm(x): return x / 25.4

INK = palette.INK
WARM_A, WARM_B = palette.DIVERGING[-1], palette.DIVERGING[-2]
LAYERS = ["campaign", "gate", "library"]
TICKS = ["Campaign", "Gate", "Library"]

# label, colour, linestyle, how to pull it out of the row
SERIES = [
    ("U10",        WARM_A, "-",  lambda r: nt_share(r, [10])),
    ("Loop 7–11", WARM_B, "-", lambda r: nt_share(r, range(7, 12))),
    ("Base",       palette.THEME["800"], (0, (4, 1.6)), lambda r: moiety(r, "base")),
    ("Sugar",      palette.THEME["500"], (0, (4, 1.6)), lambda r: moiety(r, "sugar")),
    ("Phosphate",  palette.THEME["300"], (0, (4, 1.6)), lambda r: moiety(r, "phosphate")),
]
MOIETIES = ("phosphate", "sugar", "base")


def nt_share(r, nts):
    return float(sum(r[f"nt{i:02d}_{m}"] for i in nts for m in MOIETIES))


def moiety(r, m):
    return float(sum(r[f"nt{i:02d}_{m}"] for i in range(1, 17)))


def nudged(values, gap):
    """Label positions: the values, pushed apart to `gap` where they collide.

    Base and phosphate start 0.005 apart and their labels sit on top of each
    other otherwise. Returns positions in the original order.
    """
    out = {}
    last = -1e9
    for i in sorted(range(len(values)), key=lambda k: values[k]):
        out[i] = max(values[i], last + gap)
        last = out[i]
    return [out[i] for i in range(len(values))]


def main():
    df = pd.read_csv(paths.INPUTS / "epitope_layers_union.csv").set_index("layer")
    x = np.arange(3)

    fig, ax = plt.subplots(figsize=(mm(98), mm(58)))
    ends = []
    for label, colour, ls, get in SERIES:
        y = np.array([get(df.loc[l]) for l in LAYERS])
        ax.plot(x, y, ls=ls, color=colour, lw=1.1, zorder=2)
        ax.plot(x, y, "o", ms=3.4, mfc=colour, mec="white", mew=0.6, zorder=3)
        ends.append([label, colour, y])

    # start value on the left, name + end value on the right; both sets nudged
    # apart so near-equal series do not print on top of each other
    left = nudged([e[2][0] for e in ends], 0.028)
    right = nudged([e[2][-1] for e in ends], 0.030)
    for i, (label, colour, y) in enumerate(ends):
        ax.annotate(f"{y[0]:.2f}", (x[0], left[i]), textcoords="offset points",
                    xytext=(-7, 0), ha="right", va="center", fontsize=6,
                    color=colour, annotation_clip=False)
        ax.annotate(f"{label}  {y[-1]:.2f}", (x[-1], right[i]),
                    textcoords="offset points", xytext=(9, 0), ha="left",
                    va="center", fontsize=6.5, color=colour,
                    annotation_clip=False)

    ax.set_xticks(x)
    ax.set_xticklabels(TICKS, fontsize=7)
    ax.set_xlim(-0.33, 2.33)
    ax.set_ylim(0, 0.56)
    ax.set_yticks([0, 0.1, 0.2, 0.3, 0.4, 0.5])
    ax.set_ylabel("Share of contacts")
    ax.tick_params(length=2, width=0.6, pad=1.5)
    ax.spines["bottom"].set_visible(False)
    ax.tick_params(axis="x", length=0)

    # n under each stage: the funnel is the x axis, so its sizes belong on it
    for xi, lay in zip(x, LAYERS):
        ax.annotate(f"n = {int(df.loc[lay].n_nodes):,}", (xi, 0),
                    textcoords="offset points", xytext=(0, -20), ha="center",
                    va="top", fontsize=6, color=palette.NEUTRAL["600"],
                    annotation_clip=False)

    fig.subplots_adjust(left=0.11, right=0.74, top=0.97, bottom=0.20)
    for ext in ("svg", "pdf", "png"):
        fig.savefig(out / f"epitope_funnel_summary.{ext}", bbox_inches="tight")
    print(f"wrote {out}/epitope_funnel_summary.svg|pdf|png")
    for label, _, y in ends:
        print(f"  {label:11s} " + "  ".join(f"{v:.3f}" for v in y))


if __name__ == "__main__":
    main()
