"""Fig3e: GPU cost against design quality, one point per design route.

x  designs per GPU-hour = 1 / (GPU-h per design), from the gpubench_v1 benchmark
   (Euler, 2026-11-12): campaign-identical jobs on pinned eu-g6 RTX 4090s,
   covering generation plus all three refold backends.
y  pass rate = the share of a route's OWN designs that reach the 55,706 passing
   set (0.80% of 6,992,515 campaign-wide).

No shading and no reference line: both were removed on 2026-11-17 (see the note
beside the marks below). The x axis is linear from zero, not log, and the
campaign-wide pass rate is not drawn on the panel at all.

Two measurement decisions sit behind the x axis, both in inputs/ header:

* Refolding is measured on an equal 600-design POOL per design route, not per uuid.
  AF3 and ESM load the model once over --input_dir, so per-design cost depends
  on batch size, and HalluDesign only makes 150 designs per uuid -- refolding
  per uuid would have handed it a quarter-size batch and inflated its cost for
  a reason that has nothing to do with HalluDesign. It came out route-
  independent (0.003638-0.003696, a 1.6% spread), which is the check that the
  pooling worked.
* The two BoltzGen routes come out of ONE allocation. Each is charged the full
  SHARED backbone stage plus its OWN inverse folding -- what that route would
  cost run alone, which is the question this panel asks. Measured on a
  timestamped run: backbone 16.7%, BIF 72.6%, LigandMPNN 10.7%, i.e. BIF costs
  6.8x LigandMPNN because it runs a Boltz prediction per design while
  LigandMPNN only assigns sequence. Splitting the job 50/50 instead would be
  wrong twice: it credits a backbone saving that only exists if you want both
  products, and it halves a 6.8:1 stage.
"""
import csv
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import matplotlib.pyplot as plt
import numpy as np
from common import palette, paths, style

FIG = Path(__file__).resolve().parent.name
style.apply()
out = paths.fig_output(FIG)
# A size variant for page assembly; the canonical panel is untouched:
#   make_figure.py W_MM H_MM FONT_PT OUTDIR
VARIANT = len(sys.argv) == 5
if VARIANT:
    plt.rcParams["font.size"] = float(sys.argv[3])
    out = Path(sys.argv[4])

# --- data ------------------------------------------------------------------
rows = []
with (paths.INPUTS / "Fig3e_gpu_efficiency.csv").open() as fh:
    for r in csv.DictReader(l for l in fh if not l.startswith("#")):
        gpuh = float(r["gpuh_gen"]) + float(r["gpuh_refold"])
        rows.append({
            "stack": r["stack"],
            "x": 1.0 / gpuh,
            "y": 100.0 * int(r["passing"]) / int(r["generated"]),
            "gpuh_per_pass": gpuh * int(r["generated"]) / int(r["passing"]),
        })
# canonical stack order, never sorted by value (COLORS.md rule 2)
rows.sort(key=lambda d: palette.STACK_ORDER.index(d["stack"]))

XLIM, YLIM = (0, 300), (0, 2.1)   # linear x from zero (Leo 2026-11-17: no log scale)
X_EDGE, Y_EDGE = 200.0, 100.0 * 55_706 / 6_992_515      # campaign-wide pass rate

# Label offset in points from its mark; the sign picks the side the leader
# leaves from. Tuned per point, NOT a uniform diagonal: a single shared offset
# puts every label in the same relative spot, which crowds the dense corner and
# strands the sparse one. Check the render after moving any mark.
LABEL_POS = {
    "RFD3 + LigandMPNN":        (-11, 9),
    "BoltzGen + LigandMPNN":    (11, 11),
    "BoltzGen + BoltzIF":       (11, 11),
    "BoltzDesign + LigandMPNN": (11, -9),
    "HalluDesign + LigandMPNN": (6, 24),
}

# Double column. A scatter with five direct labels needs the width: at single
# column the 7 pt labels are proportionally huge and collide, and the reference
# figure this follows is a wide panel too.
fig, ax = plt.subplots(figsize=(style.mm(float(sys.argv[1]) if VARIANT else style.CFG["width_double_mm"]),
                                style.mm(float(sys.argv[2]) if VARIANT else 78)))

# No shaded quadrants (Leo, 2026-11-17: "remove the boxes here too"). They were
# two tinted rectangles split at X_EDGE and the campaign-wide pass rate; the
# five marks and their direct labels carry the comparison without them, and a
# background band is the same ink-with-no-data the no-box rule is about.
# NOTE: the campaign-wide pass rate (Y_EDGE) was only ever shown by that shading,
# so it is no longer on the panel. A single hairline at Y_EDGE would put it back
# in idiom (STYLE_NOTES 7: one axis, hairline, palette.GRID, set_axisbelow).

# --- marks + direct labels (STYLE_NOTES 6: no legend box) ------------------
for d in rows:
    c = palette.stack_color(d["stack"])
    ax.scatter([d["x"]], [d["y"]], s=22, color=c, zorder=4,
               edgecolors=palette.SURFACE, linewidths=0.8)
    dx, dy = LABEL_POS[d["stack"]]
    # Label in the SERIES colour, matching Fig3a (and STYLE_NOTES 10, which
    # notes the reference does this against general dataviz advice). Consistency
    # within figure 3 wins; Fig3a labels its stack ribbons the same way.
    ax.annotate(palette.SHORT[d["stack"]], xy=(d["x"], d["y"]), color=c,
                xytext=(dx, dy), textcoords="offset points",
                ha="left" if dx > 0 else "right", va="center", zorder=6,
                arrowprops=dict(arrowstyle="-", color=palette.NEUTRAL["500"],
                                lw=0.5, shrinkA=1, shrinkB=3))

# --- axes ------------------------------------------------------------------
ax.set_xlim(*XLIM)
ax.set_ylim(*YLIM)
# No direction arrows here. STYLE_NOTES 2 is for metrics whose good direction
# is not self-evident (RMSD, H-bond energy). On a throughput and a pass rate,
# "higher is better" is not news and the glyph is just one more mark to read.
ax.set_xlabel("Designs per GPU-hour")
ax.set_ylabel("Pass rate")
ax.set_xticks([0, 100, 200, 300])
ax.xaxis.set_major_formatter(plt.FuncFormatter(lambda v, _: f"{v:,.0f}"))
ax.set_yticks(np.arange(0, 2.1, 0.5))
ax.yaxis.set_major_formatter(plt.FuncFormatter(lambda v, _: f"{v:g}%"))
ax.tick_params(length=2)
ax.set_axisbelow(True)   # no gridlines: STYLE_NOTES 7

style.save(fig, out, FIG)
fig.savefig(out / f"{FIG}.svg", bbox_inches="tight")   # for page assembly
print("wrote", out / f"{FIG}.pdf")

# The derived quantity the panel is a picture of; for the caption, not the canvas.
print(f"\n{'stack':<12}{'designs/GPU-h':>14}{'pass%':>8}{'GPU-h per passing':>19}")
for d in rows:
    print(f"{palette.SHORT[d['stack']]:<12}{d['x']:>14,.0f}{d['y']:>8.3f}"
          f"{d['gpuh_per_pass']:>19.2f}")
