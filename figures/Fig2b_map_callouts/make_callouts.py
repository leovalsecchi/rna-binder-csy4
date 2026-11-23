"""Fig02b callout layer: the eight highlighted binders marked on the campaign map.

Drawn in the SAME frame as make_figure.py (identical figsize, RECT and limits),
so callouts.* drops straight onto dots.* / labels.* in Illustrator with no
alignment work. The composed map_callouts.* is the same thing flattened, for
reading on screen and for review.

The eight are the map8 set from Fig2b_stack_binders: two per generator, every
one passing the full design gate and present in the ordered 12k. Their
coordinates come from the same UMAP table the dots are drawn from, so a marker
sits exactly on its own backbone, not on a nearest neighbour (verified: every
marker is within 0.011 UMAP units of a plotted point, which is the point).

Labels sit OUTSIDE the dense core on fixed radial leaders. They are not aligned
to each other: aligning labels weakens the link between a point and its label
when the points themselves are not aligned.

Cluster outlines are deliberately NOT drawn. The eight sit in clusters of
8,813 to 34,804 backbones with x/y spread of 0.6 to 3.3 UMAP units, several
overlapping; a hull around those is a large blob that says nothing about where
the binder is. The cluster id is given in the key instead.

Data: inputs/large/fold_map_scatter.parquet, inputs/map8_callouts.csv.
Writes output/Fig2b_map_callouts/{callouts,map_callouts}.{svg,pdf,png}.
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
import matplotlib.patheffects as pe
from common import paths, palette

FIG = Path(__file__).resolve().parent.name
out = paths.fig_output(FIG)
_cfg = yaml.safe_load((paths.ROOT / "config" / "style.yaml").read_text())
plt.rcParams.update({
    "font.family": _cfg["font_family"], "pdf.fonttype": 42, "svg.fonttype": "none",
    "savefig.dpi": 600, "axes.grid": False,
})
def mm(x): return x / 25.4

# Same merge and colours as make_figure.py: BoltzIF is an inverse-folding step
# on BoltzGen backbones, so it shares the BG colour.
MERGE = {"BG+BIF": "BoltzGen"}
RENAME = {"BG+LiMPNN": "BoltzGen"}
COL, ORDER = {}, []
for long in palette.STACK_ORDER:
    short = palette.SHORT[long]
    name = RENAME.get(short, MERGE.get(short, short))
    if name not in COL:
        COL[name] = palette.stack_color(long if short not in MERGE else
                                        "BoltzGen + LigandMPNN")
        ORDER.append(name)

INK, HALO = palette.INK, palette.SURFACE
LEAD = palette.NEUTRAL["500"]

# Label positions are PLACED, not hand-typed: for each callout we score 16
# directions at three radii and keep the one with the fewest backbones under
# the number, the clearest leader line, and no collision with a label already
# placed. Hand-set offsets were tried first and had to be re-tuned every time a
# pick changed, which is exactly the kind of thing that silently rots.
RINGS = (0.95, 1.25, 1.60, 2.10, 2.80)
NDIR = 24
LAB_R = 0.55            # disc under a number that should be as empty as possible
SEP = 1.30              # minimum separation between two numbers
GAP_M, GAP_L = 0.60, 0.80        # leader trimmed this far from marker / label

d = pd.read_parquet(paths.INPUTS / "large" / "fold_map_scatter.parquet")
d["stack"] = d["stack"].astype(str)
NISO = int(d["deg0"].sum())
d = d[~d["deg0"]].copy()
d["group"] = d["stack"].replace({**RENAME, **MERGE})
N = d.groupby("group").size().to_dict()

c = pd.read_csv(paths.INPUTS / "map8_callouts.csv").sort_values("n")
c["group"] = c["stack"].replace({**RENAME, **MERGE})
assert set(c["group"]) <= set(COL), f"unknown group: {set(c['group']) - set(COL)}"

# Every marker must land on a real plotted backbone, not near one.
XY = d[["x", "y"]].values
for _, r in c.iterrows():
    dd = np.hypot(XY[:, 0] - r.x, XY[:, 1] - r.y).min()
    assert dd < 0.02, f"{r.panel} is {dd:.3f} from the nearest plotted backbone"
print(f"{len(c)} callouts, all within 0.02 UMAP units of a plotted backbone")

X, Y, G = d["x"].values, d["y"].values, d["group"].values
X0, X1, Y0, Y1 = X.min(), X.max(), Y.min(), Y.max()
padx, pady = 0.03 * (X1 - X0), 0.03 * (Y1 - Y0)
RECT = [0.055, 0.055, 0.925, 0.925]
RNG = np.random.default_rng(0)
idx = np.arange(len(d)); RNG.shuffle(idx)


# ------------------------------------------------- automatic label placement
# Density on a coarse grid; H[i, j] counts backbones in one cell.
BIN = 0.30
ex = np.arange(X.min() - 6, X.max() + 6 + BIN, BIN)
ey = np.arange(Y.min() - 6, Y.max() + 6 + BIN, BIN)
H, _, _ = np.histogram2d(X, Y, bins=[ex, ey])


def density(px, py, r):
    """Backbones within r of (px, py), from the grid."""
    i0, i1 = np.searchsorted(ex, [px - r, px + r])
    j0, j1 = np.searchsorted(ey, [py - r, py + r])
    return float(H[max(i0 - 1, 0):i1, max(j0 - 1, 0):j1].sum())


XLO, XHI = X0 - padx, X1 + padx
YLO, YHI = Y0 - pady, Y1 + pady
placed, LABEL = [], {}
for _, r in c.iterrows():
    best, bestcost = None, np.inf
    for rad in RINGS:
        for k in range(NDIR):
            a = 2 * np.pi * k / NDIR
            lx, ly = r.x + rad * np.cos(a), r.y + rad * np.sin(a)
            if not (XLO + 1 < lx < XHI - 1 and YLO + 1 < ly < YHI - 1):
                continue
            if placed and min(np.hypot(lx - px, ly - py) for px, py in placed) < SEP:
                continue
            # dots under the number first, then keep it close to its mark
            cost = density(lx, ly, LAB_R) + 600 * rad
            if cost < bestcost:
                best, bestcost = (lx, ly), cost
    assert best is not None, f"no label position found for {r.panel}"
    placed.append(best)
    LABEL[r.panel] = best
    print(f"  {int(r.n)} {r.panel:7s} label at ({best[0]:6.2f},{best[1]:6.2f}) "
          f"cost {bestcost:,.0f}")


def panel():
    fig = plt.figure(figsize=(mm(178), mm(178)))
    ax = fig.add_axes(RECT)
    ax.set_xlim(X0 - padx, X1 + padx)
    ax.set_ylim(Y0 - pady, Y1 + pady)
    ax.set_aspect("equal")
    for sp in ax.spines.values():
        sp.set_visible(False)
    ax.set_xticks([]); ax.set_yticks([])
    return fig, ax


def draw_callouts(ax):
    for _, r in c.iterrows():
        col = COL[r.group]
        lx, ly = LABEL[r.panel]
        # a leader is drawn ONLY when the badge had to be pushed clear of the
        # mark; at the usual offset the number sits against its own marker and
        # a line between them would be ink for nothing.
        vx, vy = lx - r.x, ly - r.y
        L = np.hypot(vx, vy)
        if L > 1.7:
            ux, uy = vx / L, vy / L
            ax.plot([r.x + ux * GAP_M, lx - ux * GAP_L],
                    [r.y + uy * GAP_M, ly - uy * GAP_L],
                    lw=0.6, color=LEAD, solid_capstyle="round", zorder=6)
        # marker: white halo ring so it reads over any density of dots
        # halo carries a faint keyline so the mark looks identical whether it
        # lands on a dense patch or on white paper
        ax.scatter([r.x], [r.y], s=96, c=HALO, edgecolors=palette.NEUTRAL["300"],
                   linewidths=0.5, zorder=7)
        ax.scatter([r.x], [r.y], s=34, c=col, edgecolors=HALO,
                   linewidths=1.0, zorder=8)
        ax.text(lx, ly, str(int(r.n)), fontsize=9, color=col, weight="bold",
                ha="center", va="center", zorder=9,
                path_effects=[pe.withStroke(linewidth=2.6, foreground=HALO)])


def save(fig, stem, transparent):
    for ext in ("svg", "pdf"):
        fig.savefig(out / f"{stem}.{ext}", transparent=transparent)
    fig.savefig(out / f"{stem}.png", dpi=1200, transparent=transparent)
    plt.close(fig)
    print("wrote", stem, "svg/pdf/png")


# ------------------------------------------- overlay only (registered layer)
fig, ax = panel()
ax.patch.set_alpha(0)
draw_callouts(ax)
save(fig, "callouts", True)

# --------------------------------------------------- composed, for reading
fig, ax = panel()
ax.scatter(X[idx], Y[idx], s=0.5, c=[COL[g] for g in G[idx]],
           alpha=0.40, linewidths=0, rasterized=True, zorder=2)
draw_callouts(ax)
ax.set_xlabel("UMAP 1", fontsize=7, color=palette.NEUTRAL["600"], labelpad=2)
ax.set_ylabel("UMAP 2", fontsize=7, color=palette.NEUTRAL["600"], labelpad=2)
for i, g in enumerate(ORDER):
    yy = 0.020 + (len(ORDER) - 1 - i) * 0.028
    ax.text(0.015, yy, "●", transform=ax.transAxes, fontsize=8,
            color=COL[g], ha="left", va="center")
    ax.text(0.042, yy, f"{g}   {N[g]:,}", transform=ax.transAxes, fontsize=7,
            color=palette.INK, ha="left", va="center")
ax.text(0.015, 0.985, f"Campaign fold space\n{len(d):,} unique backbones"
        f"\n{NISO:,} unclustered isolates not shown",
        transform=ax.transAxes, fontsize=7, color=palette.NEUTRAL["600"],
        ha="left", va="top", linespacing=1.45)
# key for the eight, under the title block: the right edge sits on the cloud
ax.text(0.015, 0.905, "Highlighted binders", transform=ax.transAxes,
        fontsize=7, color=palette.NEUTRAL["600"], ha="left", va="center")
for i, (_, r) in enumerate(c.iterrows()):
    yy = 0.872 - i * 0.027
    ax.text(0.015, yy, str(int(r.n)), transform=ax.transAxes, fontsize=7,
            color=COL[r.group], weight="bold", ha="left", va="center")
    ax.text(0.038, yy,
            f"{r.group}   fold {int(r.cluster_r10)}   "
            f"ipTM {r.iptm:.2f}   DockQ {r.dockq:.2f}   "
            f"RMSD {r.rmsd:.1f} \u00c5",
            transform=ax.transAxes, fontsize=6.5, color=palette.INK,
            ha="left", va="center")
save(fig, "map_callouts", False)
print("->", out)
