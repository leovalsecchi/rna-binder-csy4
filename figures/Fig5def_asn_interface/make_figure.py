"""Three steps: the standard metrics are flat, only asparagine moves, and it moves on
both measures.

  metrics   the twelve metrics the field ranks designed binders on, as P(binder > non)
            with bootstrap intervals. All twelve cover 0.5.
  aminoacid all nineteen amino acids present in these designs (Cys is excluded by the
            design protocol), as the same statistic, for sequence composition and for
            interface composition. Asn is the only one that survives correction.
  asndist   asparagine in the 12k library, the flow inactive designs, the flow binders and
            59 natural RBPs - as a % of the sequence and as a % of the interface.

Why P(binder > non) throughout: n is 15 and 17, the quantities have incomparable units
and several are bounded or zero-inflated. It is the rank statistic the permutation tests
use, bounded, with a fixed null at 0.5.

**ligand-RMSD and DockQ appear PER BACKEND, never as the better of the two**, because the
refold used elsewhere is picked by lowest ligand-RMSD and testing the selected value would
test the variable that did the selecting.

**Panel b is 19 tests and says so.** Asn reaches q = 0.046 (sequence) and 0.064
(interface); His is nominally significant (P 0.049 and 0.021) and does NOT survive
correction (q 0.20). Two of nineteen at P < 0.05 where one is expected by chance.

**Panel c is percentages, not counts**: the naturals run 61-1301 aa with ~41-residue
interfaces against the designs' 80 aa and ~21, so counts are not comparable. For a design,
1.25% of sequence is one residue. Large groups are drawn as a density, small ones as their
actual points - with n of 15 a density would invent shape that is not in the data. Both
distributions are zero-inflated, so each row carries its "% with none" as well as a mean.

**The honest shape, which the caption must carry**: the binders are not unusually
Asn-rich. At 5.7% of interface residues they are close to the 600-design library sample at
4.6%. The inactive designs at 1.5% are the anomaly, and 71% of them have no interface
asparagine at all against 41% of the library. This predicts failure from a missing
contact, not success from a clever one.

Data: inputs/asn/, written by export_data.py and export_aa20.py here. Analysis:
~/hits32_browse_flowFC/FLOW_FC_ANALYSIS.md appendices 8-9e. P values are 20,000 label
permutations stratified within design stack; intervals are bootstraps.

    python code/Fig5def_asn_interface/make_figure.py [outdir]
    python code/Fig5def_asn_interface/make_figure.py W_MM H_MM FONT_PT OUTDIR
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

IN = REPO / "inputs" / "asn"
PAD_IN = 0.04
_VARIANT = ((float(sys.argv[1]), float(sys.argv[2]), float(sys.argv[3]))
            if len(sys.argv) >= 5 else None)

style.apply()
FS = plt.rcParams["font.size"]
if _VARIANT:
    FS = _VARIANT[2]
    plt.rcParams["font.size"] = FS
VAR_SIZE = (lambda w, h: (mm(_VARIANT[0]), mm(_VARIANT[1])) if _VARIANT
            else (mm(w), mm(h)))

N = palette.NEUTRAL
TH = palette.THEME
MAIN = palette.MAIN
BIND, NON = MAIN, N["600"]
LIB, NAT = N["400"], TH["400"]


def _bare(ax, left=True):
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_linewidth(0.5)
    ax.tick_params(width=0.5, length=1.8, pad=1.2, labelsize=FS)
    if not left:
        ax.spines["left"].set_visible(False)
        ax.tick_params(axis="y", length=0)


# ------------------------------------------------------------- panel: metrics
def panel_metrics(out):
    B = pd.read_csv(IN / "h1_metrics.csv")
    fig, ax = plt.subplots(figsize=VAR_SIZE(86, 52))
    y = np.arange(len(B))[::-1]
    for yi, r in zip(y, B.itertuples()):
        covers = r.lo <= 0.5 <= r.hi
        c = N["600"] if covers else MAIN
        ax.plot([r.lo, r.hi], [yi, yi], "-", lw=0.6, color=c, solid_capstyle="butt", zorder=2)
        for xx in (r.lo, r.hi):
            ax.plot([xx, xx], [yi - 0.22, yi + 0.22], "-", lw=0.6, color=c, zorder=2)
        ax.plot([r.auc], [yi], "o" if covers else "D", ms=2.6 if covers else 2.6,
                mfc=c, mec=c, mew=0, zorder=3)
        if not covers:
            ax.annotate(f"$P$ = {r.p:.3f}", xy=(r.hi, yi), xytext=(3, 0),
                        textcoords="offset points", color=c, fontsize=FS, va="center")
    ax.axvline(0.5, lw=0.5, color=N["400"], zorder=1)
    ax.set_yticks(y)
    ax.set_yticklabels(B.label, fontsize=FS)
    ax.set_ylim(-0.8, len(B) - 0.2)
    ax.set_xlim(0, 1.0)
    ax.set_xticks([0, 0.25, 0.5, 0.75, 1])
    ax.set_xticklabels(["0", "0.25", "0.5", "0.75", "1"], fontsize=FS)
    ax.set_xlabel("P(proposed binder > inactive design)")
    _bare(ax, left=False)
    fig.savefig(out / "metrics.svg", bbox_inches="tight", pad_inches=PAD_IN)
    fig.savefig(out / "metrics.png", bbox_inches="tight", pad_inches=PAD_IN, dpi=400)
    plt.close(fig)


# ----------------------------------------------------------- panel: aminoacid
def panel_aminoacid(out):
    A = pd.read_csv(IN / "aa20.csv").sort_values("auc_seq")
    fig, ax = plt.subplots(figsize=VAR_SIZE(86, 62))
    y = np.arange(len(A))[::-1]
    for yi, r in zip(y, A.itertuples()):
        hit = r.q_seq < 0.1 or r.q_if < 0.1
        c = MAIN if hit else N["500"]
        ax.plot([r.auc_seq, r.auc_if], [yi, yi], "-", lw=0.5, color=c, alpha=0.5, zorder=2)
        ax.plot([r.auc_seq], [yi], "o", ms=2.6, mfc=c, mec=c, mew=0, zorder=3)
        ax.plot([r.auc_if], [yi], "s", ms=2.4, mfc="none", mec=c, mew=0.6, zorder=3)
        if r.p_seq < 0.05 or r.p_if < 0.05:
            ax.annotate(f"$q$ = {min(r.q_seq, r.q_if):.2f}",
                        xy=(max(r.auc_seq, r.auc_if), yi), xytext=(3.5, 0),
                        textcoords="offset points", color=c, fontsize=FS, va="center")
    ax.axvline(0.5, lw=0.5, color=N["400"], zorder=1)
    ax.set_yticks(y)
    THREE = {"A": "Ala", "R": "Arg", "N": "Asn", "D": "Asp", "Q": "Gln", "E": "Glu",
             "G": "Gly", "H": "His", "I": "Ile", "L": "Leu", "K": "Lys", "M": "Met",
             "F": "Phe", "P": "Pro", "S": "Ser", "T": "Thr", "W": "Trp", "Y": "Tyr",
             "V": "Val"}
    ax.set_yticklabels([THREE.get(a, a) for a in A.aa], fontsize=FS)
    ax.set_ylim(-1.0, len(A) - 0.3)
    ax.set_xlim(0.18, 1.0)
    ax.set_xticks([0.25, 0.5, 0.75, 1])
    ax.set_xticklabels(["0.25", "0.5", "0.75", "1"], fontsize=FS)
    ax.set_xlabel("P(proposed binder > inactive design)")
    ax.plot([], [], "o", ms=2.6, mfc=N["600"], mec=N["600"], mew=0, label="whole sequence")
    ax.plot([], [], "s", ms=2.4, mfc="none", mec=N["600"], mew=0.6, label="at the interface")
    ax.legend(loc="upper left", bbox_to_anchor=(-0.02, -0.11), frameon=False, fontsize=FS,
              handletextpad=0.25, labelspacing=0.3, borderpad=0, ncol=2,
              columnspacing=0.9, handlelength=0.6)
    _bare(ax, left=False)
    fig.savefig(out / "aminoacid.svg", bbox_inches="tight", pad_inches=PAD_IN)
    fig.savefig(out / "aminoacid.png", bbox_inches="tight", pad_inches=PAD_IN, dpi=400)
    plt.close(fig)


# ------------------------------------------------------- panel: aminoacid, volcano
def panel_aminoacid_volcano(out):
    """The same 19 amino acids as a volcano: effect on x, significance on y.

    Volcano is the field's convention for "many things tested, which ones moved"
    (effect size against -log10 P, hits labelled, a threshold rule). Here the effect
    axis is P(binder > inactive design) rather than a log fold change, because the
    quantities are fractions with a natural null at 0.5 and no meaningful ratio.

    Two marks per amino acid, joined: whole-sequence composition and composition among
    interface residues. The dashed rule is the Benjamini-Hochberg threshold across the
    19 tests, not P = 0.05 - the point of the panel is that one feature clears
    correction and the rest do not.
    """
    A = pd.read_csv(IN / "aa20.csv")
    fig, ax = plt.subplots(figsize=VAR_SIZE(86, 58))
    # the BH line: the largest P whose q still clears 0.1, taken over both measures
    allp = np.sort(np.r_[A.p_seq.to_numpy(), A.p_if.to_numpy()])
    mm = len(allp)
    ok = allp <= 0.1 * np.arange(1, mm + 1) / mm
    thr = allp[ok].max() if ok.any() else 0.05
    for r in A.itertuples():
        hit = min(r.q_seq, r.q_if) < 0.1
        c = MAIN if hit else N["500"]
        xs_, ys_ = [r.auc_seq, r.auc_if], [-np.log10(r.p_seq), -np.log10(r.p_if)]
        ax.plot(xs_, ys_, "-", lw=0.5, color=c, alpha=0.45, zorder=2)
        ax.plot([xs_[0]], [ys_[0]], "o", ms=2.8, mfc=c, mec=c, mew=0, zorder=3)
        ax.plot([xs_[1]], [ys_[1]], "s", ms=2.5, mfc="none", mec=c, mew=0.6, zorder=3)
        if min(r.p_seq, r.p_if) < 0.05:
            THREE = {"N": "Asn", "H": "His"}
            ax.annotate(THREE.get(r.aa, r.aa), xy=(max(xs_), max(ys_)), xytext=(3.5, 1.5),
                        textcoords="offset points", color=c, fontsize=FS,
                        ha="left", va="bottom")
    ax.axvline(0.5, lw=0.5, color=N["400"], zorder=1)
    ax.axhline(-np.log10(thr), lw=0.5, ls=(0, (2.5, 2)), color=N["500"], zorder=1)
    ax.annotate("BH $q$ = 0.1", xy=(0.005, -np.log10(thr)), xycoords=("axes fraction", "data"),
                xytext=(0, 1.5), textcoords="offset points", color=N["500"],
                fontsize=FS, ha="left", va="bottom")
    ax.set_xlim(0.26, 0.92)
    ax.set_ylim(0, 3.0)
    ax.set_xticks([0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9])
    ax.set_xlabel("P(proposed binder > inactive design)")
    ax.set_ylabel("\u2212log10 $P$")
    ax.plot([], [], "o", ms=2.8, mfc=N["600"], mec=N["600"], mew=0, label="whole sequence")
    ax.plot([], [], "s", ms=2.5, mfc="none", mec=N["600"], mew=0.6, label="at the interface")
    ax.legend(loc="upper left", bbox_to_anchor=(-0.02, -0.17), frameon=False, fontsize=FS,
              handletextpad=0.25, labelspacing=0.3, borderpad=0, ncol=2,
              columnspacing=0.9, handlelength=0.6)
    _bare(ax)
    fig.savefig(out / "aminoacid_volcano.svg", bbox_inches="tight", pad_inches=PAD_IN)
    fig.savefig(out / "aminoacid_volcano.png", bbox_inches="tight", pad_inches=PAD_IN, dpi=400)
    plt.close(fig)


# ------------------------------------------------------------- panel: asndist
# ------------------------------------------------------------- panel: asndist
def panel_asndist(out):
    """One row per group, twice: as a % of the sequence and as a % of the interface.

    Rows rather than overlaid curves. A kernel density over an integer count needs a wide
    bandwidth to stop the integer grid reading as spikes, and that same width smooths away
    the thing that actually separates the groups - the pile of designs with no asparagine
    at all. In rows the small groups are drawn as the designs themselves, so the pile is
    visible and n is countable.

    Large samples (12,000 sequences, 600 interfaces) are a density; the 15, 17 and 59 are
    their own points, because at that size a density invents shape that is not in the data.
    The short rule on each row is the group mean.
    """
    lib12 = pd.read_csv(IN / "lib12k_aaN.csv.gz")
    L = pd.read_csv(IN / "asn_lib600.csv")
    T = pd.read_csv(IN / "asn32.csv")
    NA = pd.read_csv(IN / "asn_nat59.csv")
    b, n = T[T.binder == 1], T[T.binder == 0]
    fig, axes = plt.subplots(1, 2, figsize=VAR_SIZE(178, 50))
    rng = np.random.default_rng(4)
    blocks = [(axes[0], "Asparagine (% of sequence)", 12.0,
               [("Natural RBPs", NA.all_N * 100, NAT, False),
                ("Proposed binders", b.sel_all_N * 100, BIND, False),
                ("12k library", lib12.aa_N * 100, LIB, True),
                ("Inactive designs", n.sel_all_N * 100, NON, False)]),
              (axes[1], "Asparagine (% of interface residues)", 22.0,
               [("Natural RBPs", NA.if_N * 100, NAT, False),
                ("Proposed binders", b.sel_if_N * 100, BIND, False),
                ("12k library", L.if_N * 100, LIB, True),
                ("Inactive designs", n.sel_if_N * 100, NON, False)])]
    for ax, xlab, xmax, rows in blocks:
        xs = np.linspace(0, xmax, 400)
        for i, (lab, v, c, big) in enumerate(rows):
            yy = len(rows) - 1 - i
            v = np.asarray(v, float); v = v[np.isfinite(v)]
            if big:
                bw = max(xmax / 60, 0.25)
                d = np.exp(-0.5 * ((xs[None, :] - v[:, None]) / bw) ** 2).sum(0)
                d = d / d.max() * 0.62
                ax.fill_between(xs, yy - 0.02, yy - 0.02 + d, color=c, lw=0,
                                alpha=0.5, zorder=2)
                ax.plot(xs, yy - 0.02 + d, lw=0.6, color=c, zorder=3)
            else:
                ax.plot(np.clip(v, 0, xmax), yy + 0.16 + rng.uniform(-0.1, 0.1, len(v)),
                        "o", ms=2.1, mfc=c, mec=c, mew=0, alpha=0.85, zorder=3)
            ax.plot([v.mean(), v.mean()], [yy - 0.1, yy + 0.52], "-", lw=1.0,
                    color=c, zorder=4)
            ax.annotate(f"{lab}   $n$ = {len(v):,}", xy=(xmax * 0.995, yy + 0.60),
                        xytext=(0, 0), textcoords="offset points", color=c,
                        fontsize=FS, ha="right", va="bottom")
        ax.set_xlim(-xmax * 0.012, xmax)
        ax.set_ylim(-0.25, len(rows) - 0.1)
        ax.set_yticks([])
        ax.set_xticks(np.arange(0, xmax + 0.1, 2 if xmax < 15 else 5))
        ax.set_xlabel(xlab)
        for sp_ in ("top", "right", "left"):
            ax.spines[sp_].set_visible(False)
        ax.spines["bottom"].set_linewidth(0.5)
        ax.tick_params(width=0.5, length=1.8, pad=1.2, labelsize=FS)
        ax.tick_params(axis="y", length=0)
    fig.subplots_adjust(wspace=0.10)
    fig.savefig(out / "asndist.svg", bbox_inches="tight", pad_inches=PAD_IN)
    fig.savefig(out / "asndist.png", bbox_inches="tight", pad_inches=PAD_IN, dpi=400)
    plt.close(fig)


def main(out):
    out.mkdir(parents=True, exist_ok=True)
    panel_metrics(out)
    panel_aminoacid(out)
    panel_aminoacid_volcano(out)
    panel_asndist(out)
    print("wrote", out)


if __name__ == "__main__":
    if _VARIANT:
        main(Path(sys.argv[4]))
    else:
        main(Path(sys.argv[1]) if len(sys.argv) > 1
             else REPO / "output" / "Fig5def_asn_interface")
