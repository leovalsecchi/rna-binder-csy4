#!/usr/bin/env python3
"""Two-panel comparison: campaign backbones vs the gated 55k, same cluster frame.

Both panels use CAMPAIGN cluster ids. The passing stage has its own map, but its
ids come from clustering the induced subgraph and do not correspond to campaign
ids, so a side-by-side of the two clusterings would compare nothing. Backbones
nest 100%, so holding the frame fixed and recomputing over survivors is well
defined, and is the only way to read "how the clusters changed".

Area is each cluster's SHARE OF ITS GENERATOR, not an absolute count, so one size
key serves both panels. Sizing by count would shrink panel b by ~20x and hide the
composition shift that is the point.
"""
import csv, pathlib, sys
import numpy as np
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
REPO = ROOT   # was a hardcoded /Volumes/wicky path: wrong under the /Volumes/bsse mount
if (REPO / "code" / "common" / "style.py").exists():
    sys.path.insert(0, str(REPO / "code"))
    from common import style, palette
    style.apply(); INK, GRID, mm = palette.INK, palette.GRID, style.mm
else:
    print("NAS not mounted, using inline style")
    matplotlib.rcParams.update({"font.family": "Arial", "font.size": 7,
                                "axes.linewidth": 0.75, "pdf.fonttype": 42,
                                "svg.fonttype": "none", "axes.spines.top": False,
                                "axes.spines.right": False})
    INK, GRID = "#171b1e", "#e3e9ec"
    mm = lambda x: x / 25.4

COLOR = {"RFD3": "#a83902", "BG": "#a376e9", "BD1": "#bd8703", "HD": "#d35182"}
LABEL = {"RFD3": "RFD3 + LigandMPNN", "BG": "BoltzGen",
         "BD1": "BoltzDesign + LigandMPNN", "HD": "HalluDesign + LigandMPNN"}
ORDER = ["RFD3", "BG", "BD1", "HD"]
GEN_TOT_CAMP = {"RFD3": 152512, "BG": 151996, "BD1": 144468, "HD": 102540}
GEN_TOT_PASS = {"RFD3": 14455, "BG": 9982, "BD1": 1774, "HD": 447}

def load(path, nkey, hkey, ekey, tot):
    out = []
    for r in csv.DictReader(open(IN / path), delimiter="\t"):
        n = int(r[nkey])
        out.append({"gen": r["gen"], "cluster": int(r["cluster"]), "n": n,
                    "helix": float(r[hkey]), "strand": float(r[ekey]),
                    "share": 100 * n / tot[r["gen"]]})
    return out

camp = load("cluster_sse.tsv", "n_bb", "helix", "strand", GEN_TOT_CAMP)
pas  = load("cluster_sse_pass.tsv", "n_bb", "helix", "strand", GEN_TOT_PASS)
maxshare = max(r["share"] for r in camp + pas)
area = lambda s: 4 + 300 * (s / maxshare)

fig, axes = plt.subplots(1, 2, figsize=(mm(178), mm(78)), sharex=True, sharey=True)
gx, gy = np.meshgrid(np.linspace(-6, 102, 320), np.linspace(-5, 48, 200))
BW, CORE = 4.5, 0.80

for ax, rows, ttl in ((axes[0], camp, "campaign"), (axes[1], pas, "gated")):
    for g in ORDER:
        sub = sorted([r for r in rows if r["gen"] == g], key=lambda r: -r["n"])
        if len(sub) < 2:
            continue
        need, run, keep = CORE * sum(r["n"] for r in sub), 0, []
        for r in sub:
            keep.append(r); run += r["n"]
            if run >= need:
                break
        if len(keep) < 2:
            keep = sub[:2]
        px = np.array([r["helix"] for r in keep]); py = np.array([r["strand"] for r in keep])
        dens = np.exp(-(((gx[..., None]-px)**2 + (gy[..., None]-py)**2)/(2*BW**2))).sum(-1)
        at = np.exp(-(((px[:, None]-px)**2 + (py[:, None]-py)**2)/(2*BW**2))).sum(-1)
        lvl = at.min()*0.92
        ax.contourf(gx, gy, dens, levels=[lvl, dens.max()*10], colors=[COLOR[g]],
                    alpha=0.11, zorder=1)
        ax.contour(gx, gy, dens, levels=[lvl], colors=[COLOR[g]], linewidths=0.5,
                   alpha=0.45, zorder=1)
    for g in ORDER:
        sub = sorted([r for r in rows if r["gen"] == g], key=lambda r: -r["n"])
        ax.scatter([r["helix"] for r in sub], [r["strand"] for r in sub],
                   s=[area(r["share"]) for r in sub], marker="o", facecolor=COLOR[g],
                   alpha=0.62, linewidth=0.5, edgecolor="white", zorder=3,
                   label=LABEL[g] if ax is axes[0] else None)
    top = {}
    for r in rows:
        if r["gen"] not in top or r["n"] > top[r["gen"]]["n"]:
            top[r["gen"]] = r
    # Offsets per panel: the gated panel's big clusters all pile into the bottom
    # right, so they need leaders at staggered heights rather than the campaign set.
    OFF = ({"BD1": (30, 6), "RFD3": (16, 24), "HD": (-30, 20), "BG": (28, 18)}
           if ax is axes[0] else
           {"HD": (-20, 20), "BD1": (-6, 42), "BG": (4, 62), "RFD3": (14, 34)})
    for g, r in top.items():
        ax.annotate(f"{g} cl{r['cluster']}", (r["helix"], r["strand"]),
                    textcoords="offset points", xytext=OFF[g], ha="center",
                    fontsize=5.5, color=INK, zorder=5,
                    arrowprops=dict(arrowstyle="-", lw=0.4, color="#919699",
                                    shrinkA=1, shrinkB=2))
    ax.set_xlim(-4, 100); ax.set_ylim(-3, 46)
    ax.set_xticks(range(0, 101, 20)); ax.set_yticks(range(0, 41, 10))
    ax.set_axisbelow(True)                        # no gridlines; the tinted regions carry the structure
    ax.set_xlabel("Helix content (%)")
    n = sum(r["n"] for r in rows)
    tot = sum(GEN_TOT_CAMP.values()) if ax is axes[0] else sum(GEN_TOT_PASS.values())
    ax.text(0.015, 0.03, f"{ttl}\n{n:,d} of {tot:,d} backbones\nin {len(rows)} clusters",
            transform=ax.transAxes, ha="left", va="bottom", fontsize=6,
            color="#6d7276", linespacing=1.5)
axes[0].set_ylabel("Strand content (%)")

leg = axes[0].legend(loc="upper right", frameon=False, fontsize=5.5, handletextpad=0.3,
                     borderaxespad=0.2, labelspacing=0.4)
for h in leg.legend_handles:
    h.set_sizes([14])
key = [Line2D([], [], ls="", marker="o", markerfacecolor="none", markeredgecolor="#919699",
              markeredgewidth=0.5, markersize=(area(s)/3.14159)**0.5*2/1.6, label=f"{s}%")
       for s in (5, 20, 50)]
# Each handle box must hold its circle, or the 50% circle spills over the "20%"
# label beside it and into the title above (2026-11-13). Box sizes are in units
# of the font size, so take them off the biggest marker.
em = max(h.get_markersize() for h in key) / 5.5
k2 = axes[1].legend(handles=key, loc="upper right", frameon=False, fontsize=5.5,
                    title="share of its generator", title_fontsize=5.5, ncol=3,
                    handlelength=em, handleheight=em, handletextpad=0.4,
                    columnspacing=1.2, borderaxespad=0.2)
k2.get_title().set_color("#919699")

fig.tight_layout(pad=0.4, w_pad=1.2)
for ext, kw in (("png", {"dpi": 600}), ("pdf", {}), ("svg", {})):
    fig.savefig(OUT / f"sse_cluster_map_compare.{ext}", bbox_inches="tight", **kw)
print(f"wrote sse_cluster_map_compare.{{png,pdf,svg}}  campaign {len(camp)} / gated {len(pas)} marks")
