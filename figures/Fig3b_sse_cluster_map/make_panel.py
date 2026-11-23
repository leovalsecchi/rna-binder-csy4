#!/usr/bin/env python3
"""Static thesis panel: helix content vs strand content, one mark per fold cluster.

Reuses the thesis repo's own style and palette (code/common) rather than
restating colours here, so this panel matches the rest of the figures. Falls
back to inline values only if the NAS share is not mounted.
"""
import csv, pathlib, sys
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

HERE = pathlib.Path(__file__).resolve().parent
# Adopted into LV_ThesisFigures 2026-11-09: inputs come from inputs/, outputs go
# to output/<fig>/, like every other figure here. Falls back to the script's own
# directory so the scripts still run standalone if copied out again.
ROOT = HERE.parents[1]
IN = ROOT / "inputs" if (ROOT / "inputs" / "cluster_sse.tsv").exists() else HERE
OUT = ROOT / "output" / HERE.name
OUT.mkdir(parents=True, exist_ok=True)
REPO = pathlib.Path("/Volumes/wicky/users/lvalsecchi/Leo_Valsecchi/LV_ThesisFigures")

if (REPO / "code" / "common" / "style.py").exists():
    sys.path.insert(0, str(REPO / "code"))
    from common import style, palette
    style.apply()
    INK, GRID = palette.INK, palette.GRID
    mm = style.mm
else:                                            # NAS down: same numbers, inline
    print("NAS not mounted, using inline style")
    matplotlib.rcParams.update({"font.family": "Arial", "font.size": 7,
                                "axes.linewidth": 0.75, "pdf.fonttype": 42,
                                "svg.fonttype": "none", "axes.spines.top": False,
                                "axes.spines.right": False})
    INK, GRID = "#171b1e", "#e3e9ec"
    mm = lambda x: x / 25.4

COLOR = {"RFD3": "#a83902", "BG": "#a376e9", "BD1": "#bd8703", "HD": "#d35182"}
# One shape for everything: the four stack colours carry identity on their own
# (CVD-validated, worst adjacent pair dE 12.8 deutan / 20.2 normal vision).
MARK  = {"RFD3": "o", "BG": "o", "BD1": "o", "HD": "o"}
LABEL = {"RFD3": "RFD3 + LigandMPNN", "BG": "BoltzGen",
         "BD1": "BoltzDesign + LigandMPNN", "HD": "HalluDesign + LigandMPNN"}
ORDER = ["RFD3", "BG", "BD1", "HD"]              # canonical stack order

rows = list(csv.DictReader(open(IN / "cluster_sse.tsv"), delimiter="\t"))
for r in rows:
    for k in ("n_bb", "cluster"):
        r[k] = int(r[k])
    for k in ("helix", "strand"):
        r[k] = float(r[k])
maxbb = max(r["n_bb"] for r in rows)
area = lambda n: 4 + 260 * (n / maxbb)           # area proportional to backbones

fig, ax = plt.subplots(figsize=(mm(120), mm(96)))
# Tinted region per generator, drawn behind the marks. A convex hull is wrong here:
# BD1's two modes sit at opposite ends of the x axis and a hull between them would
# cover most of the panel. A Gaussian density contour splits into separate lobes on
# its own, which is the honest shape. The level is set so every plotted cluster of
# that generator falls inside its own region.
# The region covers the clusters holding CORE of that generator's backbones, biggest
# first. Shading every cluster instead let one stray singleton drag a lobe across the
# panel and left four regions overlapping into mud in the low-strand band.
import numpy as np
gx, gy = np.meshgrid(np.linspace(-6, 102, 360), np.linspace(-5, 48, 220))
BW = 4.5                                          # data units; blob smoothness
CORE = 0.80
for g in ORDER:
    sub = sorted([r for r in rows if r["gen"] == g], key=lambda r: -r["n_bb"])
    need, run, keep = CORE * sum(r["n_bb"] for r in sub), 0, []
    for r in sub:
        keep.append(r); run += r["n_bb"]
        if run >= need:
            break
    if len(keep) < 2:
        keep = sub[:2]
    px = np.array([r["helix"] for r in keep]); py = np.array([r["strand"] for r in keep])
    dens = np.exp(-(((gx[..., None] - px) ** 2 + (gy[..., None] - py) ** 2)
                    / (2 * BW ** 2))).sum(-1)
    at_pts = np.exp(-(((px[:, None] - px) ** 2 + (py[:, None] - py) ** 2)
                      / (2 * BW ** 2))).sum(-1)
    lvl = at_pts.min() * 0.92
    ax.contourf(gx, gy, dens, levels=[lvl, dens.max() * 10],
                colors=[COLOR[g]], alpha=0.11, zorder=1)
    ax.contour(gx, gy, dens, levels=[lvl], colors=[COLOR[g]],
               linewidths=0.5, alpha=0.45, zorder=1)

for g in ORDER:                                   # draw in canonical order
    sub = sorted([r for r in rows if r["gen"] == g], key=lambda r: -r["n_bb"])
    ax.scatter([r["helix"] for r in sub], [r["strand"] for r in sub],
               s=[area(r["n_bb"]) for r in sub], marker=MARK[g],
               facecolor=COLOR[g], alpha=0.62, linewidth=0.5,
               edgecolor="white", zorder=3, label=LABEL[g])

# Direct labels: each generator's LARGEST cluster, one per generator. Labelling the
# six biggest overall left BoltzGen unlabelled entirely -- its clusters are many and
# individually small, so it never reached the top six despite holding 50 of the 78.
# Offsets are hand-placed with leader lines; the big marks sit in a crowded band.
OFFSET = {("BD1", 0): (14, 16), ("RFD3", 3): (16, 24),
          ("HD", 1): (-30, 20), ("BG", 9): (30, 20)}
top_per_gen = {}
for r in rows:
    if r["gen"] not in top_per_gen or r["n_bb"] > top_per_gen[r["gen"]]["n_bb"]:
        top_per_gen[r["gen"]] = r
for r in top_per_gen.values():
    dx, dy = OFFSET[(r["gen"], r["cluster"])]
    ax.annotate(f"{r['gen']} cl{r['cluster']}", (r["helix"], r["strand"]),
                textcoords="offset points", xytext=(dx, dy),
                ha="center", fontsize=5.5, color=INK, zorder=5,
                arrowprops=dict(arrowstyle="-", lw=0.4, color="#919699",
                                shrinkA=1, shrinkB=2))

ax.set_xlabel("Helix content (%)")
ax.set_ylabel("Strand content (%)")
ax.set_xlim(-4, 100); ax.set_ylim(-3, 46)
ax.set_xticks(range(0, 101, 20)); ax.set_yticks(range(0, 41, 10))
ax.set_axisbelow(True)                            # no gridlines; the tinted regions carry the structure

leg = ax.legend(loc="upper right", frameon=False, fontsize=6,
                handletextpad=0.3, borderaxespad=0.2, labelspacing=0.45)
for h in leg.legend_handles:
    h.set_sizes([18])

# size key, drawn as its own legend so it cannot collide with the marks
sizes = [1000, 10000, 50000]
key = [Line2D([], [], ls="", marker="o", markerfacecolor="none",
              markeredgecolor="#919699", markeredgewidth=0.5,
              markersize=(area(s) / 3.14159) ** 0.5 * 2 / 1.6,
              label=f"{s//1000}k" if s >= 1000 else str(s)) for s in sizes]
k2 = ax.legend(handles=key, loc="lower left", frameon=False, fontsize=5.5,
               title="backbones", title_fontsize=5.5, ncol=3,
               handletextpad=0.35, columnspacing=1.0, borderaxespad=0.2)
k2.get_title().set_color("#919699")
ax.add_artist(leg)

fig.tight_layout(pad=0.4)
for ext, kw in (("png", {"dpi": 600}), ("pdf", {}), ("svg", {})):
    fig.savefig(OUT / f"sse_cluster_map.{ext}", bbox_inches="tight", **kw)
print(f"wrote sse_cluster_map.{{png,pdf,svg}} with {len(rows)} marks")
