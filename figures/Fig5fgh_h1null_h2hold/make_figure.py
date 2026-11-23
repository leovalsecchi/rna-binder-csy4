"""H1 is null; of H2, only the sequence part holds.

Four bare panels, in reading order:

  scan      every feature tested binder-vs-non-binder, as a uniform QQ plot. If no
            feature carries signal the P values ARE uniform, so the series lies on
            the diagonal inside the order-statistic envelope. Both scorings are
            drawn: AF3-for-every-design (409 features) and the best-representing
            refold (247).
  h1metrics the ten metrics the field actually ranks RNA-binder designs on, as
            P(binder > non-binder) with a bootstrap interval. 0.5 is "no
            information". Nine of ten intervals cover it.
  h2legs    each leg of H2 under three stresses: AF3 always, the best refold, and
            the best refold on the 23 designs whose design structure is
            unambiguous. The structural legs collapse; the sequence legs cannot
            move, because no structure enters them.
  arom      the one thing left standing: aromatic fraction of the sequence, for
            the library, the 32 split by flow outcome, and the 59 natural RBPs.

Why P(binder > non) and not a difference of means: n is 15 and 17, the features
have incomparable units, and several are bounded or zero-inflated. The probability
of superiority is the rank statistic the permutation tests already use, it is
bounded in [0,1] with a fixed null at 0.5, and it reads directly as "how often does
a binder beat a non-binder on this metric".

Why a QQ plot for panel 1: the claim is a NEGATIVE one, and the honest way to show
that a batch of tests found nothing is to show the P values are indistinguishable
from uniform. A volcano or a ranked bar would show the smallest few P values out of
context and invite the reader to believe them. The envelope is what makes it a test
rather than an assertion: 2,000 draws of n uniforms, the 2.5th and 97.5th
percentile of each order statistic, so a series inside it is consistent with pure
noise at every point along its range.

Both axes are -log10 so the small-P tail, the only part anyone argues about, is
legible. The diagonal is y = x, not a fit.

Panel 3's three stresses are not three independent experiments and the panel does
not claim they are. They are nested robustness checks on ONE pre-specified test
(AF3, cut 1.2), and the figure's point is which legs survive them, not which cell
has the smallest P.

Data: inputs/h1h2/, written by export_data.py in this folder from
~/hits32_browse_flowFC (FLOW_FC_ANALYSIS.md appendices 4-6). P values are 20,000
label permutations, stratified within design stack where the text says so.

    python code/Fig5fgh_h1null_h2hold/make_figure.py [outdir]
    python code/Fig5fgh_h1null_h2hold/make_figure.py W_MM H_MM FONT_PT OUTDIR   # one size
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

IN = REPO / "inputs" / "h1h2"
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

INK = palette.INK
N = palette.NEUTRAL
TH = palette.THEME
MAIN = palette.MAIN
SEQ_C = MAIN                 # a leg that needs no structure
STRUCT_C = N["500"]          # a leg measured on a predicted structure


def _bare(ax):
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_linewidth(0.5)
    ax.tick_params(width=0.5, length=1.8, pad=1.2, labelsize=FS)


# ---------------------------------------------------------------- panel: scan
def panel_scan(out):
    fig, ax = plt.subplots(figsize=VAR_SIZE(86, 46))
    rng = np.random.default_rng(0)
    series = [("Best-representing refold", np.load(IN / "p_best.npy"), MAIN, "o", "-"),
              ("AF3 for every design", np.load(IN / "p_af3.npy"), N["600"], "^", "--")]
    # order-statistic envelope for the larger n, drawn once
    nmax = max(len(p) for _, p, *_ in series)
    draws = np.sort(rng.random((2000, nmax)), axis=1)
    lo, hi = np.percentile(draws, [2.5, 97.5], axis=0)
    exp = (np.arange(1, nmax + 1) - 0.5) / nmax
    ax.fill_between(-np.log10(exp), -np.log10(hi), -np.log10(lo),
                    color=N["200"], lw=0, zorder=0)
    top = 0
    for i, (label, p, c, mk, ls) in enumerate(series):
        n = len(p)
        e = (np.arange(1, n + 1) - 0.5) / n
        x, y = -np.log10(e), -np.log10(np.clip(p, 1 / 20001, 1))
        ax.plot(x, y, mk, ms=1.8, mfc=c, mec=c, mew=0, alpha=0.85, zorder=3)
        top = max(top, y.max())
        nsig, chance = int((p < 0.05).sum()), 0.05 * n
        # frameless key in the empty upper-left, marker then two lines of text
        yk = 0.99 - i * 0.40
        ax.plot([1.09], [yk - 0.035], mk, ms=1.8, mfc=c, mec=c, mew=0,
                transform=ax.transAxes, clip_on=False)
        ax.text(1.13, yk, f"{label}\n{n} features\n{nsig} at $P$<0.05\n"
                          f"({chance:.0f} expected by chance)",
                transform=ax.transAxes, color=c, fontsize=FS, va="top", ha="left",
                linespacing=1.35)
    lim = max(top, -np.log10(0.5 / nmax)) * 1.04
    ax.plot([0, lim], [0, lim], "-", lw=0.5, color=N["400"], zorder=1)
    ax.set_xlim(0, lim * 1.02)
    ax.set_ylim(0, lim * 1.02)
    ax.set_xlabel("Expected \u2212log10 P under the null")
    ax.set_ylabel("Observed \u2212log10 P")
    ax.set_aspect("equal")
    _bare(ax)
    fig.savefig(out / "scan.svg", bbox_inches="tight", pad_inches=PAD_IN)
    fig.savefig(out / "scan.png", bbox_inches="tight", pad_inches=PAD_IN, dpi=400)
    plt.close(fig)


# ----------------------------------------------------------- panel: h1metrics
def panel_h1metrics(out):
    B = pd.read_csv(IN / "h1_metrics.csv")
    NICE = {"ligand_RMSD (A)": "ligand-RMSD (\u00c5)", "Interface dSASA": "Interface \u0394SASA",
            "pi-stacking": "\u03c0-stacking", "Shape compl.": "Shape complementarity"}
    B["label"] = B.label.map(lambda x: NICE.get(x, x))
    fig, ax = plt.subplots(figsize=VAR_SIZE(86, 46))
    y = np.arange(len(B))[::-1]
    for yi, r in zip(y, B.itertuples()):
        covers = r.lo <= 0.5 <= r.hi
        c = N["600"] if covers else MAIN
        ax.plot([r.lo, r.hi], [yi, yi], "-", lw=0.6, color=c, solid_capstyle="butt",
                zorder=2)
        for xx in (r.lo, r.hi):                      # thin caps
            ax.plot([xx, xx], [yi - 0.22, yi + 0.22], "-", lw=0.6, color=c, zorder=2)
        ax.plot([r.auc], [yi], "o" if covers else "D", ms=2.6 if covers else 2.4,
                mfc=c, mec=c, mew=0, zorder=3)
        if not covers:                               # selective annotation, section 5
            ax.annotate(f"$P$ = {r.p:.3f}", xy=(r.hi, yi), xytext=(3, 0),
                        textcoords="offset points", color=c, fontsize=FS,
                        va="center", ha="left")
    ax.axvline(0.5, lw=0.5, color=N["400"], zorder=1)
    ax.set_yticks(y)
    ax.set_yticklabels(B.label, fontsize=FS)
    ax.set_ylim(-0.8, len(B) - 0.2)
    ax.set_xlim(0, 1)
    ax.set_xticks([0, 0.25, 0.5, 0.75, 1])
    ax.set_xticklabels(["0", "0.25", "0.5", "0.75", "1"], fontsize=FS)
    ax.set_xlabel("P(proposed binder > inactive design)")
    ax.annotate("no information", xy=(0.5, len(B) - 0.45), xytext=(1.5, 0),
                textcoords="offset points", color=N["600"], fontsize=FS, va="center")
    _bare(ax)
    ax.spines["left"].set_visible(False)
    ax.tick_params(axis="y", length=0)
    fig.savefig(out / "h1metrics.svg", bbox_inches="tight", pad_inches=PAD_IN)
    fig.savefig(out / "h1metrics.png", bbox_inches="tight", pad_inches=PAD_IN, dpi=400)
    plt.close(fig)


# -------------------------------------------------------------- panel: h2legs
def panel_h2legs(out):
    C = pd.read_csv(IN / "h2_legs.csv")
    C = C.iloc[[0, 1, 2, 3, 4, 6, 5]].reset_index(drop=True)   # composite last
    fig, ax = plt.subplots(figsize=VAR_SIZE(86, 46))
    y = np.arange(len(C))[::-1]
    STRESS = [("p_af3", "o", "AF3 always"),
              ("p_best", "s", "Best refold"),
              ("p_nobif", "D", "Best refold,\nno BG+BIF")]
    ax.axvspan(0.05, 1.0, color=N["100"], lw=0, zorder=0)
    for yi, r in zip(y, C.itertuples()):
        c = SEQ_C if r.kind == "seq" else (STRUCT_C if r.kind == "struct" else TH["500"])
        ps = [getattr(r, k) for k, _, _ in STRESS]
        ax.plot(ps, [yi] * 3, "-", lw=0.5, color=c, alpha=0.55, zorder=2)
        for j, (p, (_, mk, _)) in enumerate(zip(ps, STRESS)):
            # AF3 drawn as an open ring: on the sequence legs it sits exactly under
            # the best-refold mark (no structure enters them) and must stay visible
            open_ = (j == 0)
            ax.plot([p], [yi], mk, ms=3.2 if open_ else 2.6,
                    mfc="none" if open_ else c, mec=c, mew=0.55 if open_ else 0,
                    zorder=3)
    ax.axvline(0.05, lw=0.5, color=N["500"], zorder=1)
    ax.set_xscale("log")
    ax.set_xlim(0.0055, 1.35)
    ax.set_xticks([0.01, 0.05, 0.1, 0.5, 1])
    ax.set_xticklabels(["0.01", "0.05", "0.1", "0.5", "1"], fontsize=FS)
    ax.set_yticks(y)
    ax.set_yticklabels(C.label, fontsize=FS)
    ax.set_ylim(-1.5, len(C) - 0.2)
    ax.set_xlabel("$P$, proposed binder vs inactive design")
    for lab, c in (("sequence only", SEQ_C), ("predicted structure", STRUCT_C)):
        ax.plot([], [], "o", ms=2.6, mfc=c, mec=c, mew=0, label=lab)
    for j, (_, mk, lab) in enumerate(STRESS):
        ax.plot([], [], mk, ms=3.2 if j == 0 else 2.6,
                mfc="none" if j == 0 else N["700"], mec=N["700"],
                mew=0.55 if j == 0 else 0, label=lab.replace("\n", " "))
    ax.legend(loc="upper left", bbox_to_anchor=(-0.02, -0.26), frameon=False,
              fontsize=FS, handletextpad=0.25, labelspacing=0.3, borderpad=0,
              ncol=3, columnspacing=0.8, handlelength=0.6)
    _bare(ax)
    ax.spines["left"].set_visible(False)
    ax.tick_params(axis="y", length=0)
    fig.savefig(out / "h2legs.svg", bbox_inches="tight", pad_inches=PAD_IN)
    fig.savefig(out / "h2legs.png", bbox_inches="tight", pad_inches=PAD_IN, dpi=400)
    plt.close(fig)


# ---------------------------------------------------------------- panel: arom
def panel_arom(out):
    d = json.load(open(IN / "arom.json"))
    fig, ax = plt.subplots(figsize=VAR_SIZE(86, 46))
    lib = np.array(d["lib"]) * 100
    nat = np.array(d["nat"]) * 100
    xs = np.linspace(0, 20, 400)

    def dens(v, bw):
        k = np.exp(-0.5 * ((xs[None, :] - v[:, None]) / bw) ** 2).sum(0)
        return k / k.max()

    # densities occupy the lower half, the 32 sit as dot rows above them
    dl, dn = dens(lib, 0.6) * 0.78, dens(nat, 0.9) * 0.78
    ax.fill_between(xs, 0, dl, color=N["100"], lw=0, zorder=0)
    ax.plot(xs, dl, lw=0.6, color=N["400"], zorder=1)
    ax.plot(xs, dn, lw=0.8, color=TH["500"], ls="--", zorder=2)
    ax.annotate("ordered library", xy=(2.6, 0.52), xytext=(0, 0), textcoords="offset points",
                color=N["600"], fontsize=FS, va="bottom", ha="center")
    ax.annotate("59 natural RBPs", xy=(13.4, 0.40), xytext=(0, 0),
                textcoords="offset points", color=TH["500"], fontsize=FS,
                va="bottom", ha="left")
    # Csy4 as a tick into the density field, not a floating glyph
    cx = d["csy4"] * 100
    ax.plot([cx, cx], [0, 0.90], "-", lw=0.7, color=TH["700"], zorder=3)
    ax.annotate("Csy4", xy=(cx, 0.90), xytext=(0, 1.5), textcoords="offset points",
                color=TH["700"], fontsize=FS, va="bottom", ha="center")

    rows = [("Flow binders", np.array(d["binders"]) * 100, 1.74, MAIN, True),
            ("Inactive designs", np.array(d["non"]) * 100, 1.26, N["600"], False)]
    rng = np.random.default_rng(3)
    for lab, v, yy, c, filled in rows:
        ax.plot(v, yy + rng.uniform(-0.085, 0.085, len(v)), "o", ms=2.3,
                mfc=c if filled else "none", mec=c, mew=0.55, alpha=0.95, zorder=4)
        ax.plot([v.mean(), v.mean()], [yy - 0.21, yy + 0.21], "-", lw=1.0,
                color=c, zorder=5)
        # label parked right of the data, where both rows are empty
        ax.annotate(f"{lab}\n$n$ = {len(v)}, mean {v.mean():.1f}%",
                    xy=(16.2, yy), xytext=(0, 0), textcoords="offset points",
                    color=c, fontsize=FS, va="center", ha="left", linespacing=1.3)
    ax.set_xlim(-0.45, 20)
    ax.set_ylim(0, 2.1)
    ax.set_xticks([0, 5, 10, 15, 20])
    ax.set_yticks([])
    ax.set_xlabel("Aromatic residues (% of sequence)")
    _bare(ax)
    ax.spines["left"].set_visible(False)
    fig.savefig(out / "arom.svg", bbox_inches="tight", pad_inches=PAD_IN)
    fig.savefig(out / "arom.png", bbox_inches="tight", pad_inches=PAD_IN, dpi=400)
    plt.close(fig)


def main(out):
    out.mkdir(parents=True, exist_ok=True)
    panel_scan(out)
    panel_h1metrics(out)
    panel_h2legs(out)
    panel_arom(out)
    print("wrote", out)


if __name__ == "__main__":
    if _VARIANT:
        main(Path(sys.argv[4]))
    else:
        main(Path(sys.argv[1]) if len(sys.argv) > 1
             else REPO / "output" / "Fig5fgh_h1null_h2hold")
