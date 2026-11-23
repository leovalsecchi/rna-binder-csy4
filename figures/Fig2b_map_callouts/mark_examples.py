"""Fig02b overlay: where the rendered interaction examples sit on the fold map.

A third layer on the SAME canvas and axes rect as dots.* and labels.*, so it
drops straight onto them in Illustrator and the rings land on the right points.

Reads inputs/large/fig02b_examples.csv, produced by pick_examples.py on Euler.
Writes output/Fig2b_map_callouts/callouts.{svg,pdf,png}: numbered rings only,
transparent, no point cloud. The numbers key to the shortlist table.
"""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

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

RECT = [0.055, 0.055, 0.925, 0.925]        # MUST match make_figure.py

d = pd.read_parquet(paths.INPUTS / "large" / "fold_map_scatter.parquet")
d = d[~d["deg0"]]
X0, X1, Y0, Y1 = d.x.min(), d.x.max(), d.y.min(), d.y.max()
padx, pady = 0.03 * (X1 - X0), 0.03 * (Y1 - Y0)

ex = pd.read_csv(paths.INPUTS / "large" / "fig02b_examples.csv")
ex = ex.reset_index(drop=True)
ex["n"] = ex.index + 1

fig = plt.figure(figsize=(mm(178), mm(178)))
ax = fig.add_axes(RECT)
ax.set_xlim(X0 - padx, X1 + padx)
ax.set_ylim(Y0 - pady, Y1 + pady)
ax.set_aspect("equal")
for sp in ax.spines.values():
    sp.set_visible(False)
ax.set_xticks([]); ax.set_yticks([])
ax.patch.set_alpha(0)

for _, r in ex.iterrows():
    ax.plot(r.x, r.y, marker="o", markersize=7, markerfacecolor="none",
            markeredgecolor=palette.INK, markeredgewidth=0.9, zorder=5)
    ax.annotate(str(int(r.n)), (r.x, r.y), textcoords="offset points",
                xytext=(7, 5), fontsize=7, color=palette.INK, zorder=6)

for ext in ("svg", "pdf"):
    fig.savefig(out / f"callouts.{ext}", transparent=True)
fig.savefig(out / "callouts.png", dpi=1200, transparent=True)
plt.close(fig)
print(f"wrote callouts svg/pdf/png with {len(ex)} rings ->", out)
print(ex[["n", "mode", "design_id", "stack", "x", "y"]].to_string(index=False))
