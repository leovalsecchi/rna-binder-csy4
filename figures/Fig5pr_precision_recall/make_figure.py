"""Panels: precision-recall of the Fig 5 metrics (proposed binders vs all HTS non-hits with > 75 reads)
and AF3-vs-ESM agreement (15 / 32 / 12k). Reads inputs/fig5pr_*.csv from export_data.py.

    ~/.venvs/binder_map/bin/python code/Fig5pr_precision_recall/make_figure.py
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "code"))
from common import palette, paths, style  # noqa: E402
from common.style import mm  # noqa: E402

IN = REPO / "inputs"
OUT = paths.fig_output("Fig5pr_precision_recall")
style.apply()
FS = plt.rcParams["font.size"]
N, TH = palette.NEUTRAL, palette.THEME

# metric identity = colour + line style (never colour alone); backend = marker shape
METRICS = ["ipTM", "ligand-RMSD", "DockQ", "Interface ΔSASA", "Interface H-bonds", "Heavy contacts",
           "Shape compl.", "Packstat", "BUNS", "π-stacking", "Asn (sequence)"]
_LS = ["-", "--", ":", "-.", (0, (5, 1, 1, 1, 1, 1))]
COL = dict(zip(METRICS, [TH["400"], TH["600"], TH["800"], N["500"], N["700"],
                         TH["500"], TH["700"], TH["900"], N["400"], N["800"], palette.MAIN]))
LS = {m: _LS[i % 5] for i, m in enumerate(METRICS)}
LS["Asn (sequence)"] = "-"
MK = {"AF3": "o", "ESM": "s", "seq": "D"}


def bare(ax):
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_linewidth(0.5)
    ax.tick_params(width=0.5, length=1.8, pad=1.2, labelsize=FS)


def save(fig, stem):
    style.save(fig, OUT, stem)
    fig.savefig(OUT / f"{stem}.svg", bbox_inches="tight")  # live text (svg.fonttype none)
    plt.close(fig)


def panel_pr():
    C = pd.read_csv(IN / "fig5pr_v2_pr_curves.csv")
    S = pd.read_csv(IN / "fig5pr_v2_ap_summary.csv")
    S = S[S.primary == 1]
    base = S.baseline.iloc[0] * 100
    npos, nneg = S.n_pos.iloc[0], S.n_neg.iloc[0]
    fig = plt.figure(figsize=(mm(190), mm(78)))
    gs = fig.add_gridspec(1, 4, width_ratios=[1, 1, 0.3, 1.1], wspace=0.35)
    axs = [fig.add_subplot(gs[0, 0]), fig.add_subplot(gs[0, 1]), fig.add_subplot(gs[0, 3])]
    for ax, be in zip(axs[:2], ("AF3", "ESM")):
        for m in METRICS[:-1] + (["Asn (sequence)"] if be == "AF3" else []):
            d = C[(C.label == m) & (C.backend.isin([be, "seq"]))]
            ax.step(d.recall, d.precision * 100, where="post", lw=0.8, color=COL[m], ls=LS[m])
        ax.axhline(base, lw=0.6, color=N["900"], ls=(0, (1, 1.5)), zorder=1)
        ax.set_yscale("log")
        ax.set_ylim(0.08, 120)
        ax.set_xlim(0, 1)
        ax.set_yticks([0.1, 1, 10, 100])
        ax.set_yticklabels(["0.1", "1", "10", "100"])
        ax.set_xticks([0, 0.5, 1])
        ax.set_xlabel("Recall")
        ax.set_ylabel("Precision (%)")
        ax.annotate(f"baseline {base:.2f}%", xy=(0.8, base), xytext=(0.8, 4), ha="center", fontsize=FS,
                    color=N["900"], arrowprops=dict(arrowstyle="-", lw=0.4, color=N["900"], shrinkA=1, shrinkB=0))
        ax.text(0.03, 95, be if be == "ESM" else "AF3 (+ Asn, sequence)", fontsize=FS, va="top", color=N["800"])
        bare(ax)
    ax = axs[2]
    rows = [(m, be) for m in METRICS[:-1] for be in ("AF3", "ESM")] + [("Asn (sequence)", "seq")]
    y = np.arange(len(rows))[::-1]
    for yi, (m, be) in zip(y, rows):
        r = S[(S.label == m) & (S.backend == be)].iloc[0]
        ax.plot([r.ap_lo * 100, r.ap_hi * 100], [yi, yi], "-", lw=0.6, color=COL[m], solid_capstyle="butt")
        ax.plot([r.ap * 100], [yi], MK[be], ms=3, mfc=COL[m], mec=COL[m], mew=0)
    ax.axvline(base, lw=0.6, color=N["900"], ls=(0, (1, 1.5)), zorder=1)
    ax.set_xscale("log")
    ax.set_xlim(0.08, 3)
    ax.set_xticks([0.1, 0.3, 1, 3])
    ax.set_xticklabels(["0.1", "0.3", "1", "3"])
    ax.minorticks_off()
    ax.set_yticks(y)
    ax.set_yticklabels([f"{m}, {be}" if be != "seq" else m for m, be in rows], fontsize=FS)
    ax.set_xlabel("Average precision (%)")
    bare(ax)
    save(fig, "pr_curves")
    print(f"emitted pr_curves (n_pos={npos}, n_neg={nneg}, baseline={base:.3f}%)")


def panel_agreement():
    A = pd.read_csv(IN / "fig5pr_v2_backend_agreement.csv")
    groups = [("12k library", N["600"], "D"), ("32 ordered", TH["400"], "s"), ("15 proposed binders", palette.MAIN, "o")]
    mets = list(dict.fromkeys(A.metric))
    fig, ax = plt.subplots(figsize=(mm(86), mm(78)))
    for i, m in enumerate(mets):
        for j, (g, c, mk) in enumerate(groups):
            r = A[(A.metric == m) & (A.group == g)].iloc[0]
            yy = (len(mets) - 1 - i) + (1 - j) * 0.24
            ax.plot([r.lo, r.hi], [yy, yy], "-", lw=0.6, color=c, solid_capstyle="butt")
            ax.plot([r.spearman], [yy], mk, ms=3, mfc=c, mec=c, mew=0)
    ax.axvline(0, lw=0.5, color=N["400"], zorder=1)
    ax.set_yticks(range(len(mets)))
    ax.set_yticklabels(mets[::-1], fontsize=FS)
    ax.set_xlim(-1, 1)
    ax.set_xticks([-1, -0.5, 0, 0.5, 1])
    ax.set_xticklabels(["−1", "−0.5", "0", "0.5", "1"])
    ax.set_xlabel("AF3 vs ESM, Spearman ρ")
    # direct labels beside the top metric's three marks
    top = len(mets) - 1
    for j, (g, c, mk) in enumerate(groups):
        r = A[(A.metric == mets[0]) & (A.group == g)].iloc[0]
        ax.annotate(g, xy=(r.hi, top + (1 - j) * 0.24), xytext=(3, 0), textcoords="offset points",
                    color=c, fontsize=FS, va="center")
    bare(ax)
    save(fig, "backend_agreement")
    print("emitted backend_agreement")


if __name__ == "__main__":
    panel_pr()
    panel_agreement()
