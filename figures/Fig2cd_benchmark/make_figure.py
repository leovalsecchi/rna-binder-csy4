"""Fig03: campaign benchmark - per-stack quality metrics and compute cost.

Two separate panels, as the claim needs them kept apart:

    metrics.*             five small bar charts in one row, one per metric
    gpu_per_design.*      GPU-hours per million designs, per stack
    gpu_per_design_log.*  the same cost on a log axis, drawn as dots (see below)
    gpu_hours_bar.*       raw GPU-hours per stack, kept as an alternative

Styling follows the RFD4-Proteina paper: no chart titles (the metric name and its
unit are the axis label, panel letters are added at assembly), class shown by
hatch rather than by colour so colour is free to identify the stack, and bars
labelled with the raw counts rather than with the plotted ratio.

Data comes from inputs/campaign_summary.csv (means over all designs, transcribed
from the campaign summary slide). Nothing here is recomputed from raw scores; if
the numbers in Methods change, change the CSV, not this script.

Caveats that belong in the caption, not in the figure - see README.md for the
suggested wording. In short: BG + BIF is scored against the matching BG + LiMPNN
backbone; GPU-hours exclude the 6,671 h of BoltzGen generation shared by both BG
stacks, which gpu_per_design charges to both; sc_DockQ and ipTM are higher-better
while ligand_RMSD is lower-better; designs_m is rounded as on the slide, so the
per-design figures are good to about +-1%.

Writes PDF + PNG + SVG to output/Fig2cd_benchmark/. Nothing is written to
final/: these are panels of one figure, so they go through the normal assembly
step (Inkscape/Illustrator source in working/vector_sources/).
"""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
from matplotlib.patches import Patch
from matplotlib import font_manager

from common import paths, palette, style

FIG = Path(__file__).resolve().parent.name
out = paths.fig_output(FIG)

style.apply()

# config/style.yaml asks for Arial. It is not installed everywhere (Linux boxes,
# CI), and matplotlib would silently fall back to DejaVu Sans with a warning, so
# pick the first metric-compatible sans that is actually present.
_have = {f.name for f in font_manager.fontManager.ttflist}
for _f in (style.CFG["font_family"], "Arial", "Helvetica", "Liberation Sans",
           "Nimbus Sans", "DejaVu Sans"):
    if _f in _have:
        matplotlib.rcParams["font.family"] = _f
        break

# Colour identifies the STACK; the mechanistic class is carried by hatch instead
# (RFD4-Proteina marks method type by marker shape for the same reason). Palette
# as supplied: Shadow Grey, Twilight Indigo, Rich Cerulean, Powder Blue, Frozen
# Water. It is a lightness ramp rather than a set of distinct hues, so it does not
# meet the categorical chroma check; what was validated and must hold is
# separation - worst adjacent pair dE 17.6 in the bar order below - plus a label
# on every bar and a hairline edge, since Frozen Water is 1.13:1 against white.
# Re-run scripts/validate_palette.js before reordering or reassigning a stack.
# COLORS.md rules 1 and 2: colours come from config/palette.json, one fixed
# colour per stack, never hand-picked. The old Shadow/Indigo/Cerulean/Powder/
# Frozen ramp predates palette.json (v0.1, 2026-11-05) and is gone; the hatch
# still carries the mechanistic class, so colour is free for identity alone.

# campaign_summary.csv labels stacks more tersely than palette.SHORT, so map
# the long names onto the CSV's labels; colour and ORDER then both follow
# palette.STACK_ORDER, which is the only place the canonical order is written.
CSV_LABEL = {"RFD3 + LigandMPNN": "RFD3",
             "BoltzGen + LigandMPNN": "BG+LiM",
             "BoltzGen + BoltzIF": "BG+BIF",
             "BoltzDesign + LigandMPNN": "BD1",
             "HalluDesign + LigandMPNN": "HD"}
STACK_COLOR = {CSV_LABEL[k]: palette.STACKS[k] for k in palette.STACK_ORDER}
ORDER = list(STACK_COLOR)
CLASS_OF = {"RFD3": "diffusion", "BG+LiM": "diffusion", "BG+BIF": "diffusion",
            "BD1": "hallucination", "HD": "hallucination"}
CLASS_HATCH = {"diffusion": "", "hallucination": "///"}

EDGE = palette.NEUTRAL["700"]
ACC = palette.MAIN
BAR = palette.NEUTRAL["300"]
INK = palette.INK
MUT = palette.NEUTRAL["600"]

# How the 6,671 h BoltzGen generation job, shared by both BG stacks, is charged:
#   "each"  - the full job to each BG stack. "What did this stack cost to run?"
#             Double counts 6,671 h, so per-stack costs do not sum to the real
#             total. Shares of a total must therefore say so.
#   "split" - half to each. Totals reconcile, but neither BG stack shows what it
#             would really cost on its own.
SHARED_MODE = "each"

FS = style.CFG["font_size_pt"]   # the one house label size. Ticks, values and
                                 # the key were FS / FS: at 5.5 that is
                                 # 4.5 / 3.5 pt, under the 5 pt floor. They keep
                                 # their muted colour; the size is one.
LW = style.CFG["line_width_pt"]

# column, axis label (unit in parentheses). Labels must stay SHORT - the
# panels are ~19 mm tall, so a rotated y label longer than about 13 characters
# at 7 pt runs past the axis and is clipped. "M" (millions) and "Rg" (radius
# of gyration) are therefore expanded in the caption, not on the panel, as are
# the higher-better / lower-better directions.
METRICS = [
    ("designs_m", "Designs (M)"),
    ("iptm", "ipTM"),
    ("sc_dockq", "sc_DockQ"),
    ("sc_rmsd", "ligand_RMSD (Å)"),
    ("rg", "Rg (Å)"),
]


def load():
    return pd.read_csv(paths.INPUTS / "campaign_summary.csv", comment="#")


def charged(all_df):
    """The five stacks with the shared BoltzGen generation job charged on, plus
    GPU-hours per million designs."""
    d = all_df[all_df["kind"] == "stack"].copy()
    shared = float(all_df.loc[all_df["label"] == "shared gen", "gpu_hours"].iloc[0])
    bg = d["label"].isin(["BG+LiM", "BG+BIF"])
    d.loc[bg, "gpu_hours"] += shared if SHARED_MODE == "each" else shared / 2
    d["per_m"] = d["gpu_hours"] / d["designs_m"]     # GPU-hours per million designs
    # COLORS.md rule 2: canonical stack order in every figure, never sorted by
    # value. Enforced here so no caller has to remember it.
    d["label"] = pd.Categorical(d["label"], ORDER, ordered=True)
    return d.sort_values("label").reset_index(drop=True)


def _frame(ax):
    ax.spines["left"].set_color(INK)
    ax.spines["bottom"].set_color(INK)
    ax.tick_params(width=LW, length=2, colors=INK, labelsize=FS, pad=1.5)


def _class_key(ax, **kw):
    """Frameless two-item key for the hatch encoding, in neutral grey so it reads
    as a class marker and not as one of the stack colours."""
    handles = [Patch(facecolor=palette.NEUTRAL["200"], edgecolor=EDGE, linewidth=0.4,
                     hatch=CLASS_HATCH[c], label=c)
               for c in ("diffusion", "hallucination")]
    ax.legend(handles=handles, frameon=False, fontsize=FS, labelcolor=MUT,
              handlelength=0.9, handleheight=0.9, handletextpad=0.4,
              columnspacing=1.0, borderpad=0, borderaxespad=0.2, **kw)


def metrics_figure(all_df):
    """Five small bar charts in one row, one per metric, double column.

    Height and wspace are tuned together, and both must move at once: the total
    width is fixed at 178 mm, so widening the gaps narrows the boxes. At wspace
    0.55 the five boxes are ~75 px wide, and 25 mm of height leaves ~76 px for the
    box once the axis label and the rotated tick labels have taken theirs
    (measured h/w 0.99). Changing the font size or the label rotation changes that
    overhead - re-measure with ax.get_window_extent() rather than guessing.
    """
    df = all_df[all_df["kind"] == "stack"].reset_index(drop=True)
    w = style.mm(style.CFG["width_double_mm"])
    fig, axes = plt.subplots(1, len(METRICS), figsize=(w, style.mm(25)))
    colors = [ACC if lb == "RFD3" else BAR for lb in df["label"]]
    x = range(len(df))

    for ax, (col, label) in zip(axes, METRICS):
        ax.bar(x, df[col], color=colors, width=0.72, linewidth=0)
        ax.set_xticks(list(x))
        ax.set_xticklabels(df["label"], rotation=45, ha="right",
                           rotation_mode="anchor")
        ax.set_ylim(0, df[col].max() * 1.08)
        ax.set_ylabel(label, fontsize=FS, color=INK, labelpad=2)
        ax.margins(x=0.08)
        _frame(ax)

    fig.subplots_adjust(wspace=0.55)
    # SUPERSEDED by make_violins.py, which owns metrics.* now. Kept for
    # reference but writing to a distinct stem so a stray call cannot
    # overwrite the violins with means.
    style.save(fig, out, "metrics_bars_superseded")
    fig.savefig(out / "metrics_bars_superseded.svg", bbox_inches="tight")
    plt.close(fig)


def gpu_per_design(all_df, size_mm=None, dest=None):
    """Compute cost per design generated: GPU-hours per million designs.

    The comparable cost figure, since the stacks did not all generate the same
    number of designs. HalluDesign costs 5.5x what RFdiffusion3 costs per design,
    which the raw totals understate because HalluDesign also returned the fewest
    designs. Bars are labelled with the raw inputs (GPU-hours / designs) rather
    than with the plotted ratio, which the axis already carries.

    Linear axis, starting at zero: a bar encodes magnitude by its length, so a log
    axis would break the encoding. The reference paper uses log only on its scatter
    panels; its bar panel is linear from zero. See gpu_per_design_log for the log
    view in the form that supports it.
    """
    # barh/hlines put row 0 at the bottom, so reverse to read canonical
    # order top to bottom (COLORS.md rule 2). Never re-sort by value.
    d = charged(all_df).iloc[::-1]

    w_mm, h_mm = size_mm or (style.CFG["width_single_mm"], 40)
    fig, ax = plt.subplots(figsize=(style.mm(w_mm), style.mm(h_mm)))
    fig.subplots_adjust(left=0.40, right=0.97, top=0.88, bottom=0.17)

    y = range(len(d))
    bars = ax.barh(list(y), d["per_m"], height=0.68,
                   color=[STACK_COLOR[lb] for lb in d["label"]],
                   edgecolor=EDGE, linewidth=0.4)
    for b, lb in zip(bars, d["label"]):
        b.set_hatch(CLASS_HATCH[CLASS_OF[lb]])

    for yi, (_, r) in zip(y, d.iterrows()):
        ax.text(r["per_m"] + d["per_m"].max() * 0.02, yi,
                f'{r["gpu_hours"]/1000:.1f}k / {r["designs_m"]:.2f}M',
                va="center", ha="left", fontsize=FS, color=MUT)

    ax.set_yticks(list(y))
    ax.set_yticklabels(d["stack"])
    ax.set_xlim(0, d["per_m"].max() * 1.6)       # room for the count labels: "84.6k /
    # 1.02M" needs ~1/3 of the axis at 63 mm wide (2026-11-13; was 1.34, it overran)
    ax.xaxis.set_major_locator(mticker.MultipleLocator(40000))   # 20k steps crowded
    # thousands separators (STYLE_NOTES 7: "1,000 / 10,000"): per million, the
    # ticks run into five figures and "20000" reads as a code, not a quantity
    ax.xaxis.set_major_formatter(mticker.FuncFormatter(lambda v, _: f"{v:,.0f}"))
    ax.set_xlabel("GPU-hours per million designs", fontsize=FS, color=INK,
                  labelpad=2)
    _frame(ax)
    ax.tick_params(axis="y", length=0)
    _class_key(ax, loc="lower right", bbox_to_anchor=(1.0, 1.0), ncol=2)

    if dest is None:
        style.save(fig, out, "gpu_per_design")
    fig.savefig((dest or out) / "gpu_per_design.svg", bbox_inches="tight")
    plt.close(fig)


def gpu_per_design_log(all_df):
    """The same cost on a log axis - drawn as dots, not bars.

    A log axis is legitimate here only because the mark is a point: position, not
    length, carries the value, so there is no baseline to misread. This is how the
    reference paper uses log for compute. The 5.5x spread compresses visibly, which
    is the trade for showing it on a log scale at all.
    """
    # barh/hlines put row 0 at the bottom, so reverse to read canonical
    # order top to bottom (COLORS.md rule 2). Never re-sort by value.
    d = charged(all_df).iloc[::-1]

    w = style.mm(style.CFG["width_single_mm"])
    fig, ax = plt.subplots(figsize=(w, style.mm(40)))
    fig.subplots_adjust(left=0.40, right=0.97, top=0.88, bottom=0.17)

    y = list(range(len(d)))
    ax.hlines(y, 10, d["per_m"], color=palette.NEUTRAL["300"], linewidth=0.7, zorder=1)
    for yi, (_, r) in zip(y, d.iterrows()):
        ax.plot(r["per_m"], yi, marker="o", markersize=3.4, zorder=2,
                color=STACK_COLOR[r["label"]], markeredgecolor=EDGE,
                markeredgewidth=0.35)
        ax.text(r["per_m"] * 1.10, yi,
                f'{r["gpu_hours"]/1000:.1f}k / {r["designs_m"]:.2f}M',
                va="center", ha="left", fontsize=FS, color=MUT)

    ax.set_xscale("log")
    ax.set_xlim(10, 300)
    ax.set_xticks([10, 20, 50, 100, 200])
    ax.set_xticklabels(["10", "20", "50", "100", "200"])
    ax.set_yticks(y)
    ax.set_yticklabels(d["stack"])
    # half a row of air top and bottom: matplotlib's default 5% margin leaves the
    # bottom row sitting on the x-axis
    ax.set_ylim(-0.65, len(d) - 0.35)
    ax.set_xlabel("GPU-hours per million designs (log)", fontsize=FS, color=INK,
                  labelpad=2)
    _frame(ax)
    ax.tick_params(axis="y", length=0)

    style.save(fig, out, "gpu_per_design_log")
    fig.savefig(out / "gpu_per_design_log.svg", bbox_inches="tight")
    plt.close(fig)


def gpu_bar(all_df):
    """Raw GPU-hours per stack. Alternative to the per-design panel; plots the
    hours as accounted, so it does NOT apply SHARED_MODE."""
    df = all_df[all_df["kind"] == "stack"]
    w = style.mm(style.CFG["width_single_mm"])
    fig, ax = plt.subplots(figsize=(w, style.mm(38)))
    fig.subplots_adjust(left=0.40, right=0.97, top=0.88, bottom=0.17)

    d = df.iloc[::-1]                      # first stack at the top
    k = d["gpu_hours"] / 1000.0
    y = range(len(d))
    bars = ax.barh(list(y), k, height=0.68,
                   color=[STACK_COLOR[lb] for lb in d["label"]],
                   edgecolor=EDGE, linewidth=0.4)
    for b, lb in zip(bars, d["label"]):
        b.set_hatch(CLASS_HATCH[CLASS_OF[lb]])

    ax.set_yticks(list(y))
    ax.set_yticklabels(d["stack"])
    ax.set_xlim(0, k.max() * 1.08)
    ax.set_xlabel("GPU-hours (thousands)", fontsize=FS, color=INK, labelpad=2)
    _frame(ax)
    ax.tick_params(axis="y", length=0)
    _class_key(ax, loc="lower right", bbox_to_anchor=(1.0, 1.0), ncol=2)

    style.save(fig, out, "gpu_hours_bar")
    fig.savefig(out / "gpu_hours_bar.svg", bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__" and len(sys.argv) == 5:
    # W_MM H_MM FONT_PT OUTDIR: 2d alone, drawn at the size it lands on the
    # page with its text at the page size (code/Fig2_assembly, _d_variant)
    FS = float(sys.argv[3])
    gpu_per_design(load(), (float(sys.argv[1]), float(sys.argv[2])), Path(sys.argv[4]))
elif __name__ == "__main__":
    df = load()
    # metrics_figure(df)  <- superseded by make_violins.py (violins over the
    # full per-design distributions, not transcribed means)
    gpu_per_design(df)
    gpu_per_design_log(df)
    gpu_bar(df)
    print("wrote", out)
    print("Panels only. Assemble them into the figure in Inkscape/Illustrator "
          "(source in working/vector_sources/) and export that to final/.")
