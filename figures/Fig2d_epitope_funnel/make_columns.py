"""Campaign epitope: one column per nucleotide, the whole design campaign.

Same idiom as make_columns.py but a different population and a different
measurement, and the two must not be read as one series.

  make_columns.py        12k library + the 32 hits, contacts measured on the
                         AF3 refold, base vs backbone.
  this script            campaign / gate / library, contacts measured on the
                         DESIGN model's own coordinates, split three ways into
                         phosphate / sugar / base.

Default output is the campaign alone (`--layers campaign`); `--layers
campaign gate library` draws the funnel as three grouped columns.

## Why this source and not the campaign store

`~/binder_filter_store/` carries all 6,992,515 designs but has NO per-
nucleotide column - only `fbase` and the `foot*` footprint summaries. The
per-nucleotide numbers exist in the binder_map union run on Euler,
`foldseek_map/results_union/design/fp/canonical.parquet`: 526,467 nodes x the
48-dim fingerprint (16 nt x phosphate/sugar/base, heavy atoms within 4.5 A),
with `meta/layers.parquet` flagging each node's population.

Two consequences, both of which belong in the caption:

1. **Nodes are backbones, not designs.** 526,467 coordinate-hash nodes, and
   the campaign layer is the 499,999-node SAMPLE of the campaign, not all of
   it. That sample is unbiased at node level (verified per stack to within
   1 SD) but NOT at fold level, so this panel is a per-backbone statement.
2. **Design layer, not a refold.** These are the design model's own
   coordinates. The 12k numbers here therefore differ from the AF3 ones in
   make_columns.py (U10 0.362 against 0.231) - different structures, not a
   contradiction. All three layers here come from the one source, so the
   comparison inside this panel is internally consistent.

The join is design_id -> coord_hash (hashes/*.tsv) -> layers.parquet, checked
to land 526,467 of 526,467; `node_id` is a row index and does NOT equal
coord_hash, so joining on it positionally would silently mismatch.

Data: inputs/epitope_layers_union.csv (3 rows, means + bootstrap CI, built by
the Euler aggregation recorded in the README).
Writes output/Fig4g_hit_epitope/campaign_epitope_columns.{svg,pdf,png}.
Nothing is written to final/ - assembly exports there.
"""
import sys
import argparse
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import numpy as np
import pandas as pd
import yaml
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from common import paths, palette
from common.csy4 import SEQUENCE, BRACKETS

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
# One hue, three tones: the split is ordered (backbone -> base), not categorical.
MOIETY = [("phosphate", palette.THEME["300"]),
          ("sugar", palette.THEME["500"]),
          ("base", palette.THEME["800"])]
LAYER_LABEL = {"campaign": "Campaign", "gate": "Gate", "library": "Library",
               "hits": "Candidates (32)"}


def main():
    ap = argparse.ArgumentParser()
    # The published Fig2d panel is the three-stage funnel; campaign alone is
    # kept because it is the plain "what does the campaign bind" statement.
    ap.add_argument("--layers", nargs="+",
                    default=["campaign", "gate", "library", "hits"],
                    choices=["campaign", "gate", "library", "hits"])
    a = ap.parse_args()

    df = pd.read_csv(paths.INPUTS / "epitope_layers_union.csv").set_index("layer")
    x = np.arange(16)
    n = len(a.layers)
    w = 0.78 / n

    # 45 mm, not 52: with the stem/loop brackets below the axis the ink comes
    # out ~1.2x the authored height, and at 52 the assembled Fig 2 page
    # overflows by 3 mm (a + b + cd + e against 267 mm of block).
    fig, ax = plt.subplots(figsize=(mm(_cfg["width_double_mm"]), mm(45)))
    tops = []
    for k, lay in enumerate(a.layers):
        r = df.loc[lay]
        off = (k - (n - 1) / 2) * w
        bottom = np.zeros(16)
        for mi, (moiety, colour) in enumerate(MOIETY):
            v = np.array([r[f"nt{i+1:02d}_{moiety}"] for i in range(16)])
            # No fading between layers: the three tones ARE the moiety split,
            # and dimming one group turns its three tones into one grey. Layer
            # identity rides on position plus the direct label instead.
            ax.bar(x + off, v, w * 0.92, bottom=bottom, color=colour,
                   edgecolor="none")
            bottom += v
        # 2026-11-18 (Leo): "dont need the bars, just the mean". The bootstrap
        # 95% interval over backbones is no longer drawn. The `nt<NN>_lo` /
        # `nt<NN>_hi` columns stay in inputs/epitope_layers_union.csv, so this
        # is two lines to restore - and the interval is still what licenses
        # the cautious reading of the hits column (n = 32, about +/-0.024 per
        # nucleotide, overlapping the library almost everywhere). With the
        # bars gone that caveat has to be carried by the caption instead.
        if n > 1:
            tops.append((lay, int(np.argmax(bottom)), off, bottom.max()))

    # Labels sit in a right-aligned column clear of the bars, each joined to
    # its own bar by a hairline leader. Stacked above the bars the shorter
    # layers' labels run straight through the taller layers' bars beside them;
    # the leader lets the text move away without losing which bar it names.
    # matplotlib draws the connector from the text's own bounding box, so the
    # gap comes out right without measuring the string.
    #
    # Library and Hits differ by 0.025 at U10 - under one line height on a
    # 45 mm panel - so the LABEL y is nudged apart while the leader still
    # lands on the true bar top. Four layers is where this starts to bite.
    GAP = 0.035
    ys, last = {}, -1e9
    for lay, _, _, h in sorted(tops, key=lambda t: t[3]):
        ys[lay] = last = max(h, last + GAP)
    for lay, top, off, h in tops:
        ax.annotate(LAYER_LABEL[lay],
                    xy=(top + off - w * 0.46, h),
                    xytext=(top - 0.70, ys[lay]), textcoords="data",
                    ha="right", va="center", fontsize=6.5, color=INK,
                    arrowprops=dict(arrowstyle="-", lw=0.5,
                                    color=palette.NEUTRAL["500"],
                                    shrinkA=2.5, shrinkB=0.5))

    ax.set_xticks(x)
    ax.set_xticklabels([f"{SEQUENCE[i]}{i+1}" for i in range(16)], fontsize=6)
    ax.set_ylabel("Share of contacts")
    ax.set_xlim(-0.75, 15.75)
    ax.set_ylim(0, None)
    ax.tick_params(length=2, width=0.6, pad=1.5)

    ymax = ax.get_ylim()[1]
    for k, (moiety, colour) in enumerate(reversed(MOIETY)):
        yk = ymax * (0.93 - 0.10 * k)
        ax.add_patch(plt.Rectangle((12.6, yk), 0.30, ymax * 0.042,
                                   facecolor=colour, edgecolor="none"))
        ax.text(13.05, yk + ymax * 0.021, moiety, fontsize=6, va="center",
                ha="left", color=INK)

    # as a FRACTION of the axis, not an absolute offset: this panel's y range
    # is smaller than the hits panel's and a fixed offset lands on the ticks
    y = -0.118 * ymax
    for lo_, hi_, lab in BRACKETS:
        ax.plot([lo_, hi_], [y, y], color=palette.NEUTRAL["500"], lw=0.8,
                clip_on=False, solid_capstyle="butt")
        ax.text((lo_ + hi_) / 2, y * 1.9, lab, ha="center", va="top",
                fontsize=6, color=palette.NEUTRAL["600"], clip_on=False)

    stem = ("epitope_funnel_columns"
            if a.layers == ["campaign", "gate", "library", "hits"]
            else "epitope_columns_" + "_".join(a.layers))
    for ext in ("svg", "pdf", "png"):
        fig.savefig(out / f"{stem}.{ext}", bbox_inches="tight")
    print(f"wrote {out}/{stem}.svg|pdf|png")

    for lay in a.layers:
        r = df.loc[lay]
        nt = np.array([sum(r[f"nt{i+1:02d}_{m}"] for m, _ in MOIETY)
                       for i in range(16)])
        tot = {m: sum(r[f"nt{i+1:02d}_{m}"] for i in range(16)) for m, _ in MOIETY}
        print(f"{lay:9s} n={int(r.n_nodes):7,d}  U10={nt[9]:.3f}  "
              f"loop(7-11)={nt[6:11].sum():.3f}  "
              f"phos/sugar/base={tot['phosphate']:.3f}/{tot['sugar']:.3f}/{tot['base']:.3f}")


if __name__ == "__main__":
    main()
