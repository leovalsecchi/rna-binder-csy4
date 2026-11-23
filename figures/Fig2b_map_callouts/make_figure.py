"""Fig02b (thesis Figure 2b): campaign fold space, one map, dots coloured by design stack.

A single panel in the same shared UMAP as fig05_fold_clustering, so a position
means the same backbone in both figures.

The two BoltzGen stacks are COLLAPSED into one colour. BoltzIF is an inverse
folding step on BoltzGen backbones, so BG+LiMPNN and BG+BIF are the same
backbone source and differ only in how the sequence was drawn. Note that they
are nonetheless disjoint by coordinate hash: the 87 plotted BG+BIF backbones
are not byte-identical to any BG+LiMPNN backbone (repacking moves the
coordinates), which is why they appear as separate rows upstream. Merging is a
statement about provenance, not about coordinate identity.

Nodes are unique backbones (coordinate-hash collapsed), NOT designs: 6,992,515
designs collapse to 551,516 distinct backbones because LigandMPNN draws many
sequences per backbone. Verified upstream (export_fig.py): every coord_hash is
produced by exactly ONE stack, so the colour is unambiguous.

The 52,452 degree-0 isolates (9.5%) are DROPPED. They have no mutual k-NN edge
surviving the TM 0.5 floor, so UMAP's random init throws them onto the
perimeter: they form a ring whose angular position carries no meaning, and at
this scale it is a near-solid band that crowds the real structure inward. They
are real backbones, so the dropped count is stated on the panel.

Draw order is RANDOMISED across stacks. With half a million overlapping marks,
whichever stack is drawn last would otherwise dominate the visible colour,
which is a plotting artifact, not a result.

Data: inputs/large/fold_map_scatter.parquet (see inputs/MANIFEST.md).
Writes output/Fig2b_map_callouts/{dots,labels}.{svg,pdf,png}
as two exactly registered layers.
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

# COLORS.md rule 2: fixed colour per stack, fixed order, never cycled or sorted
# by value. Order follows palette.STACK_ORDER with the two BoltzGen entries
# collapsed in place. The merged group keeps the BG+LiMPNN hex, which is the
# palette colour of 149,846 of its 149,933 backbones (99.94%); no new hex is
# invented here.
MERGE = {"BG+BIF": "BoltzGen"}          # folded into the BG+LiMPNN colour
RENAME = {"BG+LiMPNN": "BoltzGen"}
COL = {}
ORDER = []
for long in palette.STACK_ORDER:
    short = palette.SHORT[long]
    name = RENAME.get(short, MERGE.get(short, short))
    if name not in COL:
        COL[name] = palette.stack_color(long if short not in MERGE else
                                        "BoltzGen + LigandMPNN")
        ORDER.append(name)
RNG = np.random.default_rng(0)          # fixed seed: the figure is reproducible

d = pd.read_parquet(paths.INPUTS / "large" / "fold_map_scatter.parquet")
d["stack"] = d["stack"].astype(str)
NISO = int(d["deg0"].sum())
d = d[~d["deg0"]].copy()                # perimeter ring is a UMAP artifact
d["group"] = d["stack"].replace({**RENAME, **MERGE})
N = d.groupby("group").size().to_dict()
assert set(N) == set(ORDER), f"group mismatch: {set(N) ^ set(ORDER)}"
assert N["BoltzGen"] == int((d["stack"].isin(["BG+LiMPNN", "BG+BIF"])).sum())
print(f"{len(d):,} backbones plotted ({NISO:,} isolates dropped) | "
      + " ".join(f"{g}={N[g]:,}" for g in ORDER))

X, Y, G = d["x"].values, d["y"].values, d["group"].values
X0, X1, Y0, Y1 = X.min(), X.max(), Y.min(), Y.max()
padx, pady = 0.03 * (X1 - X0), 0.03 * (Y1 - Y0)

# Two registered layers, so the dots can be placed and recoloured in
# Illustrator without the text riding along:
#   dots.*    the point cloud, nothing else
#   labels.*  axis cues, swatch key and caption, transparent
#
# Registration is the whole point, so the axes rect is FIXED with add_axes and
# both layers are saved with bbox_inches=None. tight_layout + bbox_inches="tight"
# would crop each layer to its own content, and the layers would not line up.
RECT = [0.055, 0.055, 0.925, 0.925]
idx = np.arange(len(d))
RNG.shuffle(idx)                        # interleave stacks: none sits on top


def panel():
    """A figure and axes identical in both layers."""
    fig = plt.figure(figsize=(mm(178), mm(178)))
    ax = fig.add_axes(RECT)
    ax.set_xlim(X0 - padx, X1 + padx)
    ax.set_ylim(Y0 - pady, Y1 + pady)
    ax.set_aspect("equal")
    for sp in ax.spines.values():
        sp.set_visible(False)
    ax.set_xticks([]); ax.set_yticks([])
    return fig, ax


def save(fig, stem, transparent):
    for ext in ("svg", "pdf"):
        fig.savefig(out / f"{stem}.{ext}", transparent=transparent)
    fig.savefig(out / f"{stem}.png", dpi=1200, transparent=transparent)
    plt.close(fig)
    print("wrote", stem, "svg/pdf/png")


# ------------------------------------------------------------- dots only
fig, ax = panel()
ax.scatter(X[idx], Y[idx], s=0.5, c=[COL[g] for g in G[idx]],
           alpha=0.55, linewidths=0, rasterized=True, zorder=2)
save(fig, "dots", False)

# ----------------------------------------------------------- labels only
fig, ax = panel()
ax.patch.set_alpha(0)
# UMAP axes carry no quantity, so they get a direction cue, not a scale.
ax.set_xlabel("UMAP 1", fontsize=7, color=palette.NEUTRAL["600"], labelpad=2)
ax.set_ylabel("UMAP 2", fontsize=7, color=palette.NEUTRAL["600"], labelpad=2)
# Direct swatch labels, not a colour-only key: COLORS.md rule 5 (gold, rose and
# light purple are not separable in grayscale).
for i, g in enumerate(ORDER):
    yy = 0.020 + (len(ORDER) - 1 - i) * 0.028
    ax.text(0.015, yy, "\u25cf", transform=ax.transAxes,
            fontsize=8, color=COL[g], ha="left", va="center")
    ax.text(0.042, yy, f"{g}   {N[g]:,}", transform=ax.transAxes,
            fontsize=7, color=palette.INK, ha="left", va="center")
ax.text(0.015, 0.985, f"Campaign fold space\n{len(d):,} unique backbones"
        f"\n{NISO:,} unclustered isolates not shown",
        transform=ax.transAxes, fontsize=7, color=palette.NEUTRAL["600"],
        ha="left", va="top", linespacing=1.45)
save(fig, "labels", True)

print("->", out)
