"""Fig4h: the shape of the 32 screen hits against the 12k library.

  a  where the hits sit in the library's shape cloud - elongation against
     radius of gyration, library as a density, the 32 as dots in stack colour.
     Neither axis separates them: the hits' median within-stack percentile is
     51.7 for elongation and 51.4 for Rg, dead centre of the 33.6-66.5 null
     band in make_percentile.py. The raw medians differ a little (2.39 against
     2.66, 13.1 A against 13.6 A) only because both distributions are
     long-tailed to the right. The panel shows the hits inside the cloud; it
     does not claim they sit anywhere particular in it.
  b  the one shape difference worth a panel: strand content. The pooled
     comparison leads at full weight, with the per-stack rows under it.

Panel b leads with the pooled bars and then shows the stacks, because the
pooled number is the claim and the stacks are what qualifies it. Pooled, the
hits do look strand-rich:
25.0% of the library has any strand against 31.2% of the hits, and mean strand
content runs 0.054 against 0.085. But RFD3+LiMPNN almost never produces a
strand (1.6% of its library designs) and RFD3 is the stack most depleted among
the hits (4 of 32, against 26% of the library). Most of the pooled gap is that
mix shift, not strand. Ranking designs within their own stack and permuting
the hit labels 20,000 times inside each stack gives p = 0.31 (AF3) and
p = 0.22 (ESM) - nothing. Panel b shows the per-stack bars so the reader can
see the confound rather than be told about it.

Secondary structure is DSSP-style from geom.py on the AF3 complex; ESM agrees
throughout (strand 0.053 against 0.083 pooled). Helix fraction is not used as
an axis: the library is 80% helix at the median with strand at 0, so helix
saturates and cannot separate beta-rich from floppy - the same dead end
recorded in Fig3d's header.

Data: inputs/hits32_geom_12k.csv.gz.
Writes output/Fig4h_hit_geometry/shape_space.{svg,pdf,png}.
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
from matplotlib.colors import LogNorm
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
palette.register_colormaps()

INK = palette.INK
# View clipped to the 1st-99th percentile of the library on both axes. Every
# one of the 32 hits falls inside it (elongation 1.35-4.88, Rg 10.8-19.6 A),
# so nothing is hidden by the clip; the tail it drops is library-only.
XLIM, YLIM = (1.0, 9.0), (10.0, 24.0)


def panel_a(ax, df):
    lib, hit = df[df.hit == 0], df[df.hit == 1]
    hb = ax.hexbin(lib.af3_elong, lib.af3_rg, gridsize=52,
                   extent=(*XLIM, *YLIM), cmap="thesis_seq",
                   norm=LogNorm(vmin=1, vmax=None), mincnt=1, linewidths=0)
    for s in palette.STACK_ORDER:
        d = hit[hit["stack"] == palette.SHORT[s]]
        ax.plot(d.af3_elong, d.af3_rg, "o", ms=3.4,
                mfc=palette.STACKS[s], mec="white", mew=0.5, ls="none")
    ax.set_xlabel("Elongation")
    ax.set_ylabel("R$_g$ (Å)")
    ax.set_xlim(*XLIM)
    ax.set_ylim(*YLIM)
    ax.tick_params(length=2, width=0.6, pad=1.5)
    # bottom right: the only corner of this cloud that is empty
    cax = ax.inset_axes([0.60, 0.11, 0.34, 0.035])
    cb = ax.figure.colorbar(hb, cax=cax, orientation="horizontal")
    cb.outline.set_visible(False)
    cax.tick_params(length=1.6, width=0.5, pad=1.2, labelsize=5.5)
    cax.set_title("Designs", fontsize=6, pad=2, color=INK)


def panel_b(ax, df):
    rows, colors = [], []
    lib, hit = df[df.hit == 0], df[df.hit == 1]
    rows.append(("All", (lib.af3_sse_strand > 0).mean() * 100,
                 (hit.af3_sse_strand > 0).mean() * 100,
                 int((hit.af3_sse_strand > 0).sum()), len(hit)))
    colors.append(palette.NEUTRAL["500"])
    for s in palette.STACK_ORDER:
        short = palette.SHORT[s]
        l, h = lib[lib["stack"] == short], hit[hit["stack"] == short]
        rows.append((short, (l.af3_sse_strand > 0).mean() * 100,
                     (h.af3_sse_strand > 0).mean() * 100,
                     int((h.af3_sse_strand > 0).sum()), len(h)))
        colors.append(palette.STACKS[s])

    # The pooled row leads: it is the headline difference, drawn at full
    # weight. The five stack rows under it are the explanation, so they are
    # recessive - same encoding, less ink. barh puts row 0 at the bottom, so
    # the canonical stack order runs 5..1 downwards.
    y = np.array([6.4, 5.0, 4.0, 3.0, 2.0, 1.0])
    for k, (yi, (_, pl, ph, nk, n), c) in enumerate(zip(y, rows, colors)):
        lead = k == 0
        ax.barh(yi, pl, height=0.95 if lead else 0.60, color=c,
                alpha=0.55 if lead else 0.26, edgecolor="none")
        ax.plot([ph], [yi], "o", ms=5.5 if lead else 3.8, mfc=c,
                mec="white", mew=0.6, zorder=3)
        if lead:
            # both pooled numbers spelled out; this is the panel's headline
            ax.annotate(f"{pl:.1f}%", (pl, yi), textcoords="offset points",
                        xytext=(-4, 0), ha="right", va="center", fontsize=6.5,
                        color=INK)
            ax.annotate(f"{ph:.1f}%   {nk}/{n}", (ph, yi),
                        textcoords="offset points", xytext=(9, 0),
                        va="center", fontsize=6.5, color=INK)
        else:
            ax.annotate(f"{nk}/{n}", (ph, yi), textcoords="offset points",
                        xytext=(7, 0), va="center", fontsize=6, color=INK)
    ax.set_yticks(y)
    ax.set_yticklabels([r[0] for r in rows], fontsize=6.5)
    ax.get_yticklabels()[0].set_fontsize(7.5)
    for t, c in zip(ax.get_yticklabels(), colors):
        t.set_color(c if c != palette.NEUTRAL["500"] else INK)
    ax.set_xlabel("Designs with any β-strand (%)")
    ax.set_xlim(0, 62)
    ax.tick_params(length=2, width=0.6, pad=1.5)
    ax.spines["left"].set_visible(False)
    ax.tick_params(axis="y", length=0)
    # bar is the library, dot is the hits; said once, above the pooled pair
    # and clear of the two percentages printed on that row
    ax.set_ylim(0.35, 7.45)
    ax.annotate("Library", (rows[0][1], y[0]), textcoords="offset points",
                xytext=(-4, 17), ha="right", fontsize=6.5,
                color=palette.NEUTRAL["600"])
    ax.annotate("Hits", (rows[0][2], y[0]), textcoords="offset points",
                xytext=(11, 17), ha="left", fontsize=6.5, color=INK)


def main():
    df = pd.read_csv(paths.INPUTS / "hits32_geom_12k.csv.gz")
    df = df[df.af3_rg.notna()]
    fig, axes = plt.subplots(
        1, 2, figsize=(mm(_cfg["width_double_mm"]), mm(58)),
        gridspec_kw=dict(width_ratios=[1.15, 1], wspace=0.32))
    panel_a(axes[0], df)
    panel_b(axes[1], df)
    # not tight_layout: the inset colourbar in panel a is not compatible with it
    fig.subplots_adjust(left=0.07, right=0.985, top=0.95, bottom=0.14)
    for ext in ("svg", "pdf", "png"):
        fig.savefig(out / f"shape_space.{ext}", bbox_inches="tight")
    print(f"wrote {out}/shape_space.svg|pdf|png")

    lib, hit = df[df.hit == 0], df[df.hit == 1]
    for c in ("af3_elong", "af3_rg", "af3_sse_strand", "af3_sse_helix"):
        print(f"  {c:16s} library median {lib[c].median():6.2f} mean {lib[c].mean():6.3f}"
              f" | hits median {hit[c].median():6.2f} mean {hit[c].mean():6.3f}")
    inside = ((hit.af3_elong.between(*XLIM)) & (hit.af3_rg.between(*YLIM))).sum()
    print(f"hits inside the panel-a view: {inside}/{len(hit)}")


if __name__ == "__main__":
    main()
