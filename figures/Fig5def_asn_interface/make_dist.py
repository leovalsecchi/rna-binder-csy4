"""Asparagine distributions: 12k library, flow inactive designs, flow binders, natural RBPs.

Two bare panels, same four groups in the same order and colours:
  dist_total      asparagine as a % of the whole sequence
  dist_interface  asparagine as a % of the interface residues

Everything is a PERCENTAGE, not a count, because the natural proteins run 61-1301 aa
against the designs' fixed 80, and their interfaces are 41 residues against ~21, so raw
counts are not comparable across the four groups. For the designs a count is recoverable:
sequence % x 0.8 = residues out of 80.

Large groups (12k sequences, 600 interfaces) are drawn as a density; the small ones
(15, 17, 59) as their actual points, jittered, because with n that small a density would
invent shape that is not in the data. Every row carries a mean tick.

    python code/Fig5def_asn_interface/make_dist.py [outdir]
    python code/Fig5def_asn_interface/make_dist.py W_MM H_MM FONT_PT OUTDIR
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
from common import palette, style  # noqa: E402
from common.style import mm  # noqa: E402

IN = REPO / "inputs" / "asn"
PAD_IN = 0.04
_V = ((float(sys.argv[1]), float(sys.argv[2]), float(sys.argv[3]))
      if len(sys.argv) >= 5 else None)
style.apply()
FS = plt.rcParams["font.size"]
if _V:
    FS = _V[2]
    plt.rcParams["font.size"] = FS
SIZE = (lambda w, h: (mm(_V[0]), mm(_V[1])) if _V else (mm(w), mm(h)))

N = palette.NEUTRAL
TH = palette.THEME
BIND, NON, LIB, NAT = palette.MAIN, N["600"], N["400"], TH["400"]


def _rows(ax, rows, xmax, xlabel):
    rng = np.random.default_rng(4)
    xs = np.linspace(0, xmax, 400)
    for i, (lab, v, c, big) in enumerate(rows):
        y = len(rows) - 1 - i
        v = np.asarray(v, float)
        v = v[np.isfinite(v)]
        if big:                                   # density, scaled into the row
            bw = max(xmax / 60, 0.25)
            d = np.exp(-0.5 * ((xs[None, :] - v[:, None]) / bw) ** 2).sum(0)
            d = d / d.max() * 0.62
            ax.fill_between(xs, y - 0.02, y - 0.02 + d, color=c, lw=0, alpha=0.5, zorder=2)
            ax.plot(xs, y - 0.02 + d, lw=0.6, color=c, zorder=3)
        else:                                     # the actual points
            ax.plot(np.clip(v, 0, xmax), y + 0.16 + rng.uniform(-0.1, 0.1, len(v)), "o",
                    ms=2.1, mfc=c, mec=c, mew=0, alpha=0.85, zorder=3)
        ax.plot([v.mean(), v.mean()], [y - 0.1, y + 0.52], "-", lw=1.0, color=c, zorder=4)
        ax.annotate(f"{lab}   $n$ = {len(v):,}   mean {v.mean():.2f}%   "
                    f"{100 * (v == 0).mean():.0f}% with none",
                    xy=(xmax * 0.995, y + 0.60), xytext=(0, 0), textcoords="offset points",
                    color=c, fontsize=FS, ha="right", va="bottom")
    ax.set_xlim(-xmax * 0.012, xmax)
    ax.set_ylim(-0.25, len(rows) - 0.1)
    ax.set_yticks([])
    ax.set_xlabel(xlabel)
    for s in ("top", "right", "left"):
        ax.spines[s].set_visible(False)
    ax.spines["bottom"].set_linewidth(0.5)
    ax.tick_params(width=0.5, length=1.8, pad=1.2, labelsize=FS)
    ax.tick_params(axis="y", length=0)


def panel_total(out):
    m = pd.read_csv(IN / "lib12k_aaN.csv.gz")
    T = pd.read_csv(IN / "asn32.csv")
    NA = pd.read_csv(IN / "asn_nat59.csv")
    fig, ax = plt.subplots(figsize=SIZE(86, 52))
    _rows(ax, [("Natural RBPs", NA.all_N * 100, NAT, False),
               ("Flow binders", T[T.binder == 1].sel_all_N * 100, BIND, False),
               ("12k library", m.aa_N * 100, LIB, True),
               ("Inactive designs", T[T.binder == 0].sel_all_N * 100, NON, False)],
          12.0, "Asparagine (% of sequence)")
    ax.annotate("designs are 80 aa, so 1.25% = one residue", xy=(0.5, -0.01),
                xycoords="axes fraction", xytext=(0, -22), textcoords="offset points",
                ha="center", va="top", color=N["500"], fontsize=FS)
    fig.savefig(out / "dist_total.svg", bbox_inches="tight", pad_inches=PAD_IN)
    fig.savefig(out / "dist_total.png", bbox_inches="tight", pad_inches=PAD_IN, dpi=400)
    plt.close(fig)


def panel_interface(out):
    L = pd.read_csv(IN / "asn_lib600.csv")
    T = pd.read_csv(IN / "asn32.csv")
    NA = pd.read_csv(IN / "asn_nat59.csv")
    fig, ax = plt.subplots(figsize=SIZE(86, 52))
    _rows(ax, [("Natural RBPs", NA.if_N * 100, NAT, False),
               ("Flow binders", T[T.binder == 1].sel_if_N * 100, BIND, False),
               ("12k library", L.if_N * 100, LIB, True),
               ("Inactive designs", T[T.binder == 0].sel_if_N * 100, NON, False)],
          22.0, "Asparagine (% of interface residues)")
    ax.annotate("interfaces: designs ~21 residues, naturals ~41", xy=(0.5, -0.01),
                xycoords="axes fraction", xytext=(0, -22), textcoords="offset points",
                ha="center", va="top", color=N["500"], fontsize=FS)
    fig.savefig(out / "dist_interface.svg", bbox_inches="tight", pad_inches=PAD_IN)
    fig.savefig(out / "dist_interface.png", bbox_inches="tight", pad_inches=PAD_IN, dpi=400)
    plt.close(fig)


def main(out):
    out.mkdir(parents=True, exist_ok=True)
    panel_total(out)
    panel_interface(out)
    print("wrote", out)


if __name__ == "__main__":
    main(Path(sys.argv[4]) if _V else
         (Path(sys.argv[1]) if len(sys.argv) > 1
          else REPO / "output" / "Fig5def_asn_interface"))
