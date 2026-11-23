"""Square side panel: U10's share of all contacts against the other 15
nucleotides, across campaign -> gate -> library.

The columns panel beside it has the number in it, but spread over 48 bars; a
reader has to hunt for U10 in three places and hold them in their head. This
says the one thing: **one nucleotide of sixteen takes a fifth of all contacts
in the campaign and over a third by the library.**

Drawn as two bands stacked to 1, because U10 and the rest ARE the whole - a
pair of free-floating lines would invite the reader to wonder what the gap
between them is. The U10 band thickens; that thickening is the panel.

Per-nucleotide the gap is wider than the bands suggest, and the caption should
carry it: the other fifteen average 0.053 -> 0.043 each, so U10 goes from
3.9x to 8.5x the average nucleotide.

Same source, caveats and topology as make_columns.py in this folder -
backbones not designs, campaign is the 499,999-node sample, design
coordinates not a refold. See the README.

Data: inputs/epitope_layers_union.csv.
Writes output/Fig2d_epitope_funnel/u10_rise.{svg,pdf,png}.
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
LAYERS = ["campaign", "gate", "library", "hits"]
# "Hits" not "Hits (32)": the x axis has four labels on 46 mm.
TICKS = ["Campaign", "Gate", "Library", "Hits"]
MOIETIES = ("phosphate", "sugar", "base")
U10_C, REST_C = palette.THEME["700"], palette.NEUTRAL["200"]


def main():
    df = pd.read_csv(paths.INPUTS / "epitope_layers_union.csv").set_index("layer")
    u10 = np.array([sum(df.loc[l, f"nt10_{m}"] for m in MOIETIES) for l in LAYERS])
    rest = 1.0 - u10
    x = np.arange(len(LAYERS))

    fig, ax = plt.subplots(figsize=(mm(46), mm(46)))
    ax.fill_between(x, 0, u10, color=U10_C, lw=0)
    ax.fill_between(x, u10, 1.0, color=REST_C, lw=0)
    ax.plot(x, u10, color="white", lw=0.8)
    for xi, v in zip(x, u10):
        ax.plot([xi], [v], "o", ms=3, mfc="white", mec=U10_C, mew=0.8,
                zorder=3, clip_on=False)
        # first and last sit at the axis edges, so tuck them inward
        last = len(LAYERS) - 1
        ha = "left" if xi == 0 else "right" if xi == last else "center"
        dx = 3 if xi == 0 else -3 if xi == last else 0
        ax.annotate(f"{v:.2f}", (xi, v), textcoords="offset points",
                    xytext=(dx, 6), ha=ha, fontsize=6, color=INK)

    # band labels inside the bands: no legend, and the colours stay unexplained
    # nowhere else on the page
    ax.text(0.08, u10[0] / 2, "U10", ha="left", va="center", fontsize=7,
            color="white", fontweight="bold")
    ax.text(0.08, u10[0] + (1 - u10[0]) / 2, "Other 15 nt", ha="left",
            va="center", fontsize=6.5, color=palette.NEUTRAL["700"])

    ax.set_xticks(x)
    # 5.5, the house T the assembler retexts to anyway: with four stages on
    # 36.6 mm of axis, "Campaign" (8.6 mm at 5.5 pt, 9.4 at 6) clears "Gate"
    # by 1.6 mm at 5.5 and only 0.6 at 6. Measured, not estimated -
    # see common.style.text_w.
    ax.set_xticklabels(TICKS, fontsize=5.5)
    ax.set_xlim(0, len(LAYERS) - 1)
    ax.set_ylim(0, 1)
    ax.set_yticks([0, 0.25, 0.5, 0.75, 1.0])
    ax.set_yticklabels(["0", "", "0.5", "", "1"], fontsize=6)
    ax.set_ylabel("Share of contacts")
    ax.tick_params(length=2, width=0.6, pad=1.5)
    ax.spines["bottom"].set_visible(False)
    ax.tick_params(axis="x", length=0)
    # the end stages sit ON the axis edges; left-align the first label and
    # right-align the last so the panel stays square instead of growing two
    # half-labels of overhang
    labs = ax.get_xticklabels()
    labs[0].set_ha("left")
    labs[-1].set_ha("right")

    fig.subplots_adjust(left=0.19, right=0.985, top=0.94, bottom=0.10)
    for ext in ("svg", "pdf", "png"):
        fig.savefig(out / f"u10_rise.{ext}", bbox_inches="tight")
    print(f"wrote {out}/u10_rise.svg|pdf|png")
    other = rest / 15.0
    for l, a, b in zip(LAYERS, u10, other):
        print(f"  {l:9s} U10 {a:.3f}   other 15 mean {b:.4f}   "
              f"U10 = {a / b:.1f}x the average nucleotide")


if __name__ == "__main__":
    main()
