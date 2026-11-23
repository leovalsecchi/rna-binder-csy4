"""Fig 4c: enrichment against splicing for the PoolCraft screen (overlaid).

NOT a volcano, and it was wrongly called one until 2026-11-16. A volcano puts
effect size against SIGNIFICANCE; this puts log2 fold change against the
splicing ratio - two measurements, no p-value anywhere - and its job is to
show the gate and which designs clear it. (Fig 5e IS a volcano; that name
stays.)

Both comparisons overlaid in one panel: x = log2 fold change of ON/OFF splicing
(POS vs matched negative control, UNCAPPED; FC <= 1 stacks on the bottom spine at 0),
x = splicing ratio (fraction of reads spliced = `removed_fraction`), clipped at
0.25 with the overflow stacked on the right spine. Only designs with reads >= 75 in both POS and matched NEG are shown; of those,
points passing the rest of the gate (splicing >= 0.07, FC >= 10) are solid and
coloured by comparison (teal = POS_A/NEG_A, amber = POS_B/NEG_B), the rest faded.
No gridlines, full frame, no per-point labels. The gate lines are drawn but not
annotated and there is no footnote: reads >= 75, splicing >= 0.07, FC >= 10 and
the log2 FC = 0 floor belong in the caption, so the panel carries only the axes and the
legend.

Data: inputs/screen_points.csv (see inputs/MANIFEST.md and extract_from_report.py
for how it was derived and validated). Writes output/Fig4c_fc_vs_splicing/fc_vs_splicing.*
and copies to final/Fig03_Volcano.*.
"""
import sys, csv, math
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import numpy as np
import yaml
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from common import paths, palette

FIG = Path(__file__).resolve().parent.name
out = paths.fig_output(FIG)
_cfg = yaml.safe_load((paths.ROOT / "config" / "style.yaml").read_text())

plt.rcParams.update({
    "font.family": _cfg["font_family"], "pdf.fonttype": 42, "svg.fonttype": "none",
    "savefig.dpi": _cfg["dpi"], "axes.grid": False,
})
def mm(x): return x / 25.4

R_GATE, S_GATE, FC_GATE = 75, 0.07, 10.0
# y = UNCAPPED log2 fold change (fold_change_raw: POS splicing / matched NEG
# splicing, NEG = 1/NEG reads when no NEG read spliced). Leo asked 2026-11-12 to
# drop the old display clips (floor at log2 FC 2, report cap at 50 = 5.64), so
# the top runs to the real max (7.5). The bottom is floored at log2 FC = 0 (Leo,
# 2026-11-12: origin at 0,0): every design with FC <= 1, i.e. no enrichment over
# its control, including FC = 0 (no ON splicing), stacks on the bottom spine.
# x is still clipped: the 12 points above SPL_CAP stack on the right spine
# (SPL_CAP 0.25 keeps 32 of the 34 hits). Display only; the gate is unchanged.
Y_FLOOR = 0.0       # bottom spine, where FC <= 1 designs stack
SPL_CAP = 0.25      # on the x axis (splicing ratio)
# POS_A / POS_B are two screen comparisons, not design stacks, so they must NOT
# take stack colours (rule 2: stack colours mean stack identity and nothing
# else), and the theme ramp is the palette's only other data-mark hue. Two theme
# steps are a lightness pair, which never passes a categorical check. Leo asked
# for round marks on both (2026-11-12), so the pair rests on lightness alone and
# the legend count beside each label is the only other separation - the
# STYLE_NOTES 4 shape key was tried and rejected.
COL = {"POS_A_vs_NEG_A": palette.THEME["700"], "POS_B_vs_NEG_B": palette.THEME["400"]}
MARK = {"POS_A_vs_NEG_A": "o", "POS_B_vs_NEG_B": "o"}
# 2026-11-23: the labels follow the figure convention (PosA/NegA, as Figs S12-13
# use), not the raw column values, which stay the dict keys.
LAB = {"POS_A_vs_NEG_A": "PosA / NegA", "POS_B_vs_NEG_B": "PosB / NegB"}
ORDER = ["POS_A_vs_NEG_A", "POS_B_vs_NEG_B"]
GREY = palette.NEUTRAL["400"]     # the background is neutral, not a faded series

data = {}
with open(paths.INPUTS / "screen_points.csv") as f:
    for r in csv.DictReader(f):
        fc = r["fold_change_raw"]
        if fc in ("", None): continue
        # show only designs with >=75 reads in BOTH the POS arm and its matched
        # NEG control: the NEG cut is what made the ordered set 32, not 34
        # (two POS_B hits had 64 and 71 NEG reads and were not ordered).
        if float(r["reads"]) < R_GATE or float(r["neg_reads"]) < R_GATE: continue
        fc = float(fc); spl = float(r["splicing_ratio"])
        lfc = max(math.log2(fc), Y_FLOOR) if fc > 0 else Y_FLOOR
        ingate = spl >= S_GATE and fc >= FC_GATE      # gate uses the TRUE values
        d = data.setdefault(r["comparison"], {"x": [], "y": [], "g": []})
        d["x"].append(min(spl, SPL_CAP))     # x = splicing ratio
        d["y"].append(lfc)                   # y = log2 fold change
        d["g"].append(ingate)

fig, ax = plt.subplots(figsize=(mm(118), mm(128)))
# Out-of-gate points in ONE neutral grey layer. They used to be the two series
# colours at alpha 0.16, which at this density read as a solid blue mass and
# disagreed with the grey "outside gate" swatch in the legend.
bg = [(np.array(data[c]["x"]), np.array(data[c]["y"]), np.array(data[c]["g"])) for c in ORDER]
ax.scatter(np.concatenate([x[~g] for x, _, g in bg]),
           np.concatenate([y[~g] for _, y, g in bg]),
           s=6, c=GREY, alpha=0.45, linewidths=0, rasterized=True, zorder=1)
for comp in ORDER:          # in-gate hits
    d = data[comp]; x=np.array(d["x"]); y=np.array(d["y"]); g=np.array(d["g"])
    ax.scatter(x[g], y[g], s=11, c=COL[comp], marker=MARK[comp], alpha=1.0,
               linewidths=0, zorder=3)

ax.set_xlim(0, SPL_CAP + 0.008)
ax.set_xticks([0, 0.05, 0.10, 0.15, 0.20, 0.25])
ax.set_ylim(Y_FLOOR, 8)
ax.set_yticks(range(0, 9, 2))
ax.axhline(math.log2(FC_GATE), ls=":", lw=1.0, color=palette.INK, zorder=2)
ax.axvline(S_GATE, ls=":", lw=1.0, color=palette.INK, zorder=2)
for sp in ax.spines.values():
    sp.set_visible(True); sp.set_linewidth(1.1); sp.set_color(palette.INK)
ax.tick_params(labelsize=9, width=1.0, length=4)
ax.set_xlabel("Splicing ratio (fraction spliced)", fontsize=15, labelpad=6)
ax.set_ylabel(r"log$_2$(Fold Change)", fontsize=15, labelpad=6)
# gate thresholds and the log2 FC = 0 floor are stated in the caption, not on the axes

counts = {c: sum(data[c]["g"]) for c in ORDER}
handles = [Line2D([],[],marker=MARK[c],ls="",ms=8,mfc=COL[c],mec="white",mew=0.6,
           label=f"{LAB[c]}  ({counts[c]})") for c in ORDER]
handles.append(Line2D([],[],marker="o",ls="",ms=7,mfc=GREY,mec="none",
               alpha=0.7,label="outside gate"))
ax.legend(handles=handles, loc="lower center", bbox_to_anchor=(0.5, 1.01),
          ncol=3, frameon=False, fontsize=9.5, handletextpad=0.3,
          columnspacing=1.2, borderaxespad=0)
fig.tight_layout(rect=[0, 0, 1, 0.96])

made = []
for ext in ("svg", "pdf", "png"):
    p = out / f"fc_vs_splicing.{ext}"; fig.savefig(p, bbox_inches="tight"); made.append(p)
# Nothing is copied to final/ automatically. final/ holds only the figures Leo
# has named as going into the paper as they are, placed there by hand.
print("wrote", out)
