"""Fig4g: the RNA epitope of the 32 screen hits against the whole 12k library.

Three copies of the Csy4 repeat hairpin, nucleotide circles filled by the mean
share of a design's protein-RNA contacts landing on that nucleotide:
library (n=11,968), hits (n=32), and hits minus library.

Topology is measured, not assumed. The 16-nt target is 4AL5 chain B residues
5-20 (ACUGCCGUAUAGGCAG, confirmed against the design CIFs, which number it
1-16). Base pairs are the five with an N1/N3 donor-acceptor distance under
3.3 A in 4AL5 itself:

    C6-G20  U7-A19  G8-C18  C9-G17  C10-G16   (chain B numbering)
    nt2-16  nt3-15  nt4-14  nt5-13  nt6-12    (target numbering)

so the fold is a 5-bp stem, a 5-nt GUAUA loop (nt7-11) and an unpaired A at
the 5' end. Nothing in the drawing is a cartoon guess.

Contact shares come from geom.py: heavy-atom pairs within 4.5 A, each
nucleotide's count divided by the design's total, so one design's 16 values sum
to 1 and the mean is a footprint comparable across designs regardless of how
large its interface is.

The library and hit panels share one colour scale; the difference panel has its
own, symmetric about zero. Both scales are drawn as frameless bars.

Data: inputs/hits32_geom_12k.csv.gz (12,000 rows, hit flag marks the 32).
Writes output/Fig4g_hit_epitope/hairpin_footprint.{svg,pdf,png}.
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
from matplotlib.patches import Circle
from common import paths, palette
from common.csy4 import SEQUENCE, PAIRS, LOOP, DANGLE

FIG = Path(__file__).resolve().parent.name
out = paths.fig_output(FIG)
_cfg = yaml.safe_load((paths.ROOT / "config" / "style.yaml").read_text())
plt.rcParams.update({
    "font.family": _cfg["font_family"], "font.size": _cfg["font_size_pt"],
    "pdf.fonttype": 42, "svg.fonttype": "none", "savefig.dpi": 600,
    "axes.grid": False, "text.color": palette.INK,
})
def mm(x): return x / 25.4
palette.register_colormaps()

INK = palette.INK
SEQ = matplotlib.colormaps["thesis_seq"]
DIV = matplotlib.colormaps["thesis_div"]

# ------------------------------------------------- topology: see common/csy4
R_NT = 0.34
# Loop arc: 5 nucleotides over a 160 deg span at r = 1.30 puts 0.89 between
# neighbours, against 1.0 up the stem - close enough to read as one chain, and
# wide enough that the circles do not touch. A tighter arc overlaps them.
LOOP_R, LOOP_SPAN = 1.30, (170, 10)


def layout():
    """nt index -> (x, y). Stem as two ladders, loop as an outward arc."""
    pos = {DANGLE: (-0.45, -1.0)}
    for i in range(2, 7):                       # 5' strand, bottom to top
        pos[i] = (-0.45, float(i - 2))
    for i in range(12, 17):                     # 3' strand, top to bottom
        pos[i] = (0.45, float(16 - i))
    ang = np.linspace(*np.deg2rad(LOOP_SPAN), len(LOOP))
    for k, i in enumerate(LOOP):                # loop bulges out past the stem
        pos[i] = (LOOP_R * np.cos(ang[k]), 4.0 + LOOP_R * np.sin(ang[k]))
    return pos


POS = layout()


def draw(ax, vals, cmap, vmin, vmax):
    """One hairpin. vals: nt index (1-based) -> number driving the fill."""
    norm = plt.Normalize(vmin, vmax)
    # sugar-phosphate backbone, through the circle centres in sequence order
    chain = np.array([POS[i] for i in range(1, 17)])
    ax.plot(chain[:, 0], chain[:, 1], "-", color=palette.NEUTRAL["400"],
            lw=0.9, zorder=1, solid_capstyle="round")
    # base pairs as rungs, inset so they meet the circles rather than cross them
    for a, b in PAIRS:
        (xa, ya), (xb, yb) = POS[a], POS[b]
        ax.plot([xa + R_NT, xb - R_NT], [ya, yb], "-",
                color=palette.NEUTRAL["500"], lw=0.9, zorder=1)
    for i in range(1, 17):
        x, y = POS[i]
        c = cmap(norm(vals[i]))
        ax.add_patch(Circle((x, y), R_NT, facecolor=c, edgecolor=INK,
                            lw=0.5, zorder=3))
        # letter in whichever ink survives on this fill
        lum = 0.299 * c[0] + 0.587 * c[1] + 0.114 * c[2]
        ax.text(x, y, SEQUENCE[i - 1], ha="center", va="center",
                fontsize=5.6, color="white" if lum < 0.55 else INK, zorder=4)
    ax.text(*(POS[DANGLE][0] - 0.52, POS[DANGLE][1]), "5′",
            ha="center", va="center", fontsize=6, color=palette.NEUTRAL["600"])
    ax.text(*(POS[16][0] + 0.52, POS[16][1]), "3′",
            ha="center", va="center", fontsize=6, color=palette.NEUTRAL["600"])
    ax.set_xlim(-1.74, 1.74)
    ax.set_ylim(-1.50, 5.80)
    ax.set_aspect("equal")
    ax.axis("off")


def main():
    df = pd.read_csv(paths.INPUTS / "hits32_geom_12k.csv.gz")
    nt = [f"af3_nt{i}" for i in range(1, 17)]
    # one design of the 12,000 has no AF3 protein-RNA contact at all
    df = df[df[nt].notna().all(axis=1)]
    n_lib, n_hit = int((df.hit == 0).sum()), int((df.hit == 1).sum())
    lib = df[df.hit == 0][nt].mean().values
    hit = df[df.hit == 1][nt].mean().values
    diff = hit - lib
    v = max(lib.max(), hit.max())
    d = np.abs(diff).max()

    fig = plt.figure(figsize=(mm(116), mm(74)))
    # Explicit figure coordinates: the three hairpins are narrow, and a
    # gridspec wide enough to separate the two colour scales leaves a band of
    # white between them.
    cols = [0.025, 0.350, 0.675]
    labels = [f"Library  n = {n_lib:,}", f"Candidates  n = {n_hit}", "Candidates − library"]
    for x0, lab, (vals, cmap, lo, hi) in zip(cols, labels, [
            (lib, SEQ, 0, v), (hit, SEQ, 0, v), (diff, DIV, -d, d)]):
        ax = fig.add_axes([x0, 0.21, 0.300, 0.78])
        draw(ax, {i + 1: vals[i] for i in range(16)}, cmap, lo, hi)
        # which hairpin this is; without it the panel cannot be read at all
        ax.text(0.5, -0.015, lab, transform=ax.transAxes, ha="center",
                va="top", fontsize=6.5, color=INK)

    # two frameless scales: one shared by the first two hairpins, one for the third
    cax = fig.add_axes([0.075, 0.095, 0.455, 0.028])
    cax.imshow(np.linspace(0, 1, 256)[None, :], aspect="auto", cmap=SEQ,
               extent=(0, v, 0, 1))
    cax.set_xticks([0, v / 2, v])
    cax.set_xticklabels(["0", f"{v/2:.2f}", f"{v:.2f}"], fontsize=6)
    cax.set_xlabel("Share of contacts", fontsize=6.5, labelpad=1.5)

    dax = fig.add_axes([0.700, 0.095, 0.250, 0.028])
    dax.imshow(np.linspace(0, 1, 256)[None, :], aspect="auto", cmap=DIV,
               extent=(-d, d, 0, 1))
    dax.set_xticks([-d, 0, d])
    dax.set_xticklabels([f"−{d:.02f}", "0", f"+{d:.2f}"], fontsize=6)
    for a in (cax, dax):
        a.set_yticks([])
        for s in a.spines.values():
            s.set_visible(False)
        a.tick_params(length=1.6, width=0.5, pad=1.2)

    for ext in ("svg", "pdf", "png"):
        fig.savefig(out / f"hairpin_footprint.{ext}", bbox_inches="tight")
    print(f"wrote {out}/hairpin_footprint.svg|pdf|png")
    top = np.argsort(lib)[::-1][:4] + 1
    print("library hottest nt:", [(int(i), SEQUENCE[i-1], round(lib[i-1], 3)) for i in top])
    print("max |hit - library|:", round(d, 4),
          "at nt", int(np.argmax(np.abs(diff)) + 1))


if __name__ == "__main__":
    main()
