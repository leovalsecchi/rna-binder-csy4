"""The two results of the flow-FC analysis that are well powered and defensible.

  h1null     the ten metrics the field ranks designed binders on, scored on BOTH
             refolders, as P(binder > non-binder) with bootstrap intervals.
             All twenty intervals cover 0.5.
  refolder   how reproducibly each interface measurement comes out when the same
             sequence is refolded by AF3 and by ESMFold2, over all 12,000 designs.
             The measurement that produced the one "significant" structural result
             is the least reproducible of the set.

Together: nothing the field measures separates binders from non-binders here, and
the structural measurements are not stable enough to have supported the one result
that looked like an exception.

Why both backends in panel a, rather than one or the "best" one. Picking the refold
that best reproduces each design and then testing a fidelity metric on it is
circular - it selects on the variable under test, and it manufactured a spurious
ligand-RMSD effect (P = 0.043 selected, versus 0.85 on AF3 alone and 0.26 on ESM alone).
Every mark here is one fixed backend, so no selection enters, and showing both also
makes panel b's point without arguing it.

Why P(binder > non) and not a difference of means: n is 15 and 17, the metrics have
incomparable units, several are bounded or zero-inflated. The probability of
superiority is the rank statistic the permutation tests already use, bounded, with
a fixed null at 0.5, and it reads as "how often does a binder beat a non-binder".

Why Spearman over the 12,000 in panel b rather than the 32: the claim is about the
measurements, not about the binders, so it should use every design available. n is
11,999 per feature, so these correlations are precise; the 32 contribute nothing to
them.

Data: inputs/h1null/, written by export_data.py in this folder, from
~/hits32_browse_flowFC (FLOW_FC_ANALYSIS.md appendices 4-7).

    python code/Fig5fg_h1null_refolder/make_figure.py [outdir]
    python code/Fig5fg_h1null_refolder/make_figure.py W_MM H_MM FONT_PT OUTDIR
"""
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "code"))
from common import palette, style  # noqa: E402
from common.style import mm  # noqa: E402

IN = REPO / "inputs" / "h1null"
PAD_IN = 0.04
_VARIANT = ((float(sys.argv[1]), float(sys.argv[2]), float(sys.argv[3]))
            if len(sys.argv) >= 5 else None)

style.apply()
FS = plt.rcParams["font.size"]
if _VARIANT:
    FS = _VARIANT[2]
    plt.rcParams["font.size"] = FS
VAR_SIZE = (lambda w_mm, h_mm: (mm(_VARIANT[0]), mm(_VARIANT[1])) if _VARIANT
            else (mm(w_mm), mm(h_mm)))

N = palette.NEUTRAL
TH = palette.THEME
MAIN = palette.MAIN


def _bare(ax):
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_linewidth(0.5)
    ax.tick_params(width=0.5, length=1.8, pad=1.2, labelsize=FS)


# ------------------------------------------------------------- panel: h1null
def panel_h1null(out):
    A = pd.read_csv(IN / "h1_both_backends.csv")
    order = list(dict.fromkeys(A.label))
    fig, ax = plt.subplots(figsize=VAR_SIZE(86, 52))
    BE = {"AF3": (MAIN, "o", +0.17), "ESM": (N["600"], "s", -0.17)}
    for i, lab in enumerate(order):
        y0 = len(order) - 1 - i
        for be, (c, mk, dy) in BE.items():
            r = A[(A.label == lab) & (A.backend == be)]
            if r.empty:
                continue
            r = r.iloc[0]
            y = y0 + dy
            ax.plot([r.lo, r.hi], [y, y], "-", lw=0.55, color=c, solid_capstyle="butt",
                    zorder=2)
            for xx in (r.lo, r.hi):
                ax.plot([xx, xx], [y - 0.1, y + 0.1], "-", lw=0.55, color=c, zorder=2)
            ax.plot([r.auc], [y], mk, ms=2.4, mfc=c, mec=c, mew=0, zorder=3)
    ax.axvline(0.5, lw=0.5, color=N["500"], zorder=1)
    ax.set_yticks(range(len(order)))
    ax.set_yticklabels(order[::-1], fontsize=FS)
    ax.set_ylim(-0.75, len(order) - 0.05)
    ax.set_xlim(0, 1)
    ax.set_xticks([0, 0.25, 0.5, 0.75, 1])
    ax.set_xticklabels(["0", "0.25", "0.5", "0.75", "1"], fontsize=FS)
    ax.set_xlabel("P(proposed binder > inactive design)")
    ax.annotate("no information", xy=(0.5, len(order) - 0.28), xytext=(2, 0),
                textcoords="offset points", color=N["600"], fontsize=FS, va="center")
    # the key doubles as the statement, so no legend box
    for be, (c, mk, _) in BE.items():
        ax.plot([], [], mk, ms=2.4, mfc=c, mec=c, mew=0, label=f"{be} refold")
    ax.legend(loc="upper left", bbox_to_anchor=(0.0, -0.17), frameon=False, fontsize=FS,
              handletextpad=0.3, borderpad=0, ncol=2, columnspacing=1.2,
              handlelength=0.6)
    nsig = int((A.p < 0.05).sum())
    # one line, in the empty top-left corner; the rest of the statement is caption
    ax.annotate(f"{nsig} of {len(A)} tests at $P$<0.05", xy=(0.012, len(order) - 0.28),
                xytext=(0, 0), textcoords="offset points", color=N["700"],
                fontsize=FS, va="center", ha="left")
    _bare(ax)
    ax.spines["left"].set_visible(False)
    ax.tick_params(axis="y", length=0)
    fig.savefig(out / "h1null.svg", bbox_inches="tight", pad_inches=PAD_IN)
    fig.savefig(out / "h1null.png", bbox_inches="tight", pad_inches=PAD_IN, dpi=400)
    plt.close(fig)


# ----------------------------------------------------------- panel: refolder
def panel_refolder(out):
    B = pd.read_csv(IN / "backend_agreement.csv").sort_values("rho")
    meta = json.load(open(IN / "acidphos_concord.json"))
    fig, ax = plt.subplots(figsize=VAR_SIZE(86, 52))
    y = np.arange(len(B))
    worst = B.rho.idxmin()
    for yi, r in zip(y, B.itertuples()):
        c = MAIN if r.Index == worst else N["500"]
        ax.plot([0, r.rho], [yi, yi], "-", lw=2.2, color=c, solid_capstyle="butt",
                alpha=0.95, zorder=2)
        ax.annotate(f"{r.rho:.2f}", xy=(r.rho, yi), xytext=(2, 0),
                    textcoords="offset points", color=c, fontsize=FS, va="center")
    ax.set_yticks(y)
    ax.set_yticklabels(B.label, fontsize=FS)
    ax.set_ylim(-1.9, len(B) - 0.4)
    ax.set_xlim(0, 1)
    ax.set_xticks([0, 0.25, 0.5, 0.75, 1])
    ax.set_xticklabels(["0", "0.25", "0.5", "0.75", "1"], fontsize=FS)
    ax.set_xlabel("AF3 vs ESMFold2 agreement (Spearman, 11,999 designs)")
    ax.annotate(f"Acid–phosphate contact flag: AF3 marks "
                f"{100*meta['af3_rate']:.1f}% of designs, ESMFold2 "
                f"{100*meta['esm_rate']:.1f}%,\nbut of the {meta['either']:,} marked by "
                f"either they agree on {100*meta['concord']:.0f}% "
                f"({meta['af3_only']:,} AF3 only, {meta['esm_only']:,} ESM only)",
                xy=(0.0, -1.35), xytext=(0, 0), textcoords="offset points",
                color=N["700"], fontsize=FS, va="center", ha="left", linespacing=1.35)
    _bare(ax)
    ax.spines["left"].set_visible(False)
    ax.tick_params(axis="y", length=0)
    fig.savefig(out / "refolder.svg", bbox_inches="tight", pad_inches=PAD_IN)
    fig.savefig(out / "refolder.png", bbox_inches="tight", pad_inches=PAD_IN, dpi=400)
    plt.close(fig)


def main(out):
    out.mkdir(parents=True, exist_ok=True)
    panel_h1null(out)
    panel_refolder(out)
    print("wrote", out)


if __name__ == "__main__":
    if _VARIANT:
        main(Path(sys.argv[4]))
    else:
        main(Path(sys.argv[1]) if len(sys.argv) > 1
             else REPO / "output" / "Fig5fg_h1null_refolder")
