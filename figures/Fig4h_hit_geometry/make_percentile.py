"""Fig4h, the negative result earned rather than asserted: across ten geometry
metrics, the 32 hits are where a random draw from the library would put them.

Two forms of the same comparison, both emitted:

  percentile_strip  one row per metric, one dot per hit at its percentile
                    within the library, with the band a random 32 would fall in.
  density_curves    one small panel per metric, the library's distribution with
                    the hits' drawn over it and the 32 values as a rug.

Percentiles are taken WITHIN each design's own stack. Pooling them would
manufacture a signal out of nothing: the stacks differ in shape, and the hits
are not drawn from the stacks in library proportion (RFD3 is 26% of the
library and 4 of 32 hits), so a pooled rank would read that mix shift as
geometry. This is the same within-stack ranking the original hit analysis used.

The band is the 2.5-97.5 percentile of the median of 32 uniform draws, from
40,000 simulations (33.6 - 66.5) - where the hits' median percentile has to
land for the metric to be saying nothing.

Eight of the ten land inside it. The two that do not are marginal and are
exactly the two weak trends the original 220-feature analysis already turned
up and discarded: interface atom pairs at 68.9 (band top 66.5, hits make
slightly larger interfaces) and helix fraction at 32.5 (band floor 33.6, hits
are slightly less helical). Ten metrics at a 95% band expect 0.5 outside by
chance, neither clears BH correction over the full feature set, and the helix
one is downstream of the stack mix the same way strand content is. They are
drawn, not hidden, and the caption should call them marginal rather than
findings.

Strand content is deliberately absent. 75% of the library has no strand at
all, so its ranks are one big tie and both a strip and a smooth curve would
draw that point mass as a feature. It gets its own binned panel in
make_shape.py instead.

scipy in ~/.venvs/binder_map is broken (dlopen on _propack), so the Gaussian
KDE here is six lines of numpy with a Silverman bandwidth rather than
scipy.stats.gaussian_kde.

Data: inputs/hits32_geom_12k.csv.gz.
Writes output/Fig4h_hit_geometry/{percentile_strip,density_curves}.{svg,pdf,png}.
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
LIB = palette.NEUTRAL["400"]
HIT = palette.THEME["700"]

METRICS = [
    ("af3_n_iface",          "Interface residues"),
    ("af3_iface_atom_pairs", "Interface atom pairs"),
    ("af3_rg",               "R$_g$ (Å)"),
    ("af3_elong",            "Elongation"),
    ("af3_sse_helix",        "Helix fraction"),
    ("af3_longest_helix",    "Longest helix (res)"),
    ("af3_n_nt_touched",     "Nucleotides touched"),
    ("af3_frac_base",        "Base contact fraction"),
    ("af3_nt_entropy",       "Contact entropy (bits)"),
    ("af3_if_KR",            "K+R at interface"),
]


def kde(x, grid):
    """Gaussian KDE, Silverman bandwidth. scipy is unusable in this venv."""
    x = np.asarray(x, float)
    spread = min(x.std(ddof=1),
                 (np.percentile(x, 75) - np.percentile(x, 25)) / 1.349)
    bw = 0.9 * max(spread, 1e-9) * len(x) ** (-0.2)
    d = (grid[:, None] - x[None, :]) / bw
    return np.exp(-0.5 * d * d).sum(1) / (len(x) * bw * np.sqrt(2 * np.pi))


def null_band(n_hit, nsim=40000, seed=0):
    """Where the median percentile of n_hit random designs falls, 95% of the time."""
    rng = np.random.default_rng(seed)
    return np.percentile(np.median(rng.random((nsim, n_hit)), axis=1),
                         [2.5, 97.5]) * 100


def main():
    df = pd.read_csv(paths.INPUTS / "hits32_geom_12k.csv.gz")
    is_hit = df.hit.values == 1

    # percentile within the design's own stack
    pct = {c: df.groupby("stack")[c].rank(pct=True).values * 100
           for c, _ in METRICS}
    lo, hi = null_band(int(is_hit.sum()))

    # ---------------------------------------------------------- strip form
    fig, ax = plt.subplots(figsize=(mm(_cfg["width_double_mm"]), mm(62)))
    rng = np.random.default_rng(1)
    ax.axvspan(lo, hi, color=palette.NEUTRAL["200"], zorder=0, lw=0)
    ax.axvline(50, color=palette.NEUTRAL["400"], lw=0.6, zorder=1)
    rows = []
    for k, (c, label) in enumerate(METRICS):
        y = len(METRICS) - k
        v = pct[c][is_hit]
        v = v[~np.isnan(v)]
        ax.plot(v, y + rng.uniform(-0.17, 0.17, len(v)), "o", ms=2.6,
                mfc=HIT, mec="none", alpha=0.55, zorder=3, ls="none")
        med = np.median(v)
        ax.plot([med, med], [y - 0.34, y + 0.34], "-", color=INK, lw=1.1,
                zorder=4, solid_capstyle="butt")
        rows.append((label, med, lo <= med <= hi))
    ax.set_yticks(range(1, len(METRICS) + 1))
    ax.set_yticklabels([m[1] for m in METRICS][::-1], fontsize=6.5)
    ax.set_xlabel("Percentile within own route")
    ax.set_xlim(-1, 101)
    ax.set_ylim(0.4, len(METRICS) + 0.75)
    ax.set_xticks([0, 25, 50, 75, 100])
    ax.tick_params(length=2, width=0.6, pad=1.5)
    ax.spines["left"].set_visible(False)
    ax.tick_params(axis="y", length=0)
    ax.annotate("median of a random 32", ((lo + hi) / 2, len(METRICS) + 0.62),
                ha="center", va="bottom", fontsize=6,
                color=palette.NEUTRAL["600"])
    fig.tight_layout()
    for ext in ("svg", "pdf", "png"):
        fig.savefig(out / f"percentile_strip.{ext}", bbox_inches="tight")

    # -------------------------------------------------------- density form
    fig, axes = plt.subplots(2, 5, figsize=(mm(_cfg["width_double_mm"]), mm(62)))
    for axk, (c, label) in zip(axes.ravel(), METRICS):
        lib = df.loc[~is_hit, c].dropna().values
        hit = df.loc[is_hit, c].dropna().values
        g0, g1 = np.percentile(lib, 0.5), np.percentile(lib, 99.5)
        # An integer metric with few levels combs under a KDE - the smoother
        # puts a spike on every integer. Bin those instead of smoothing them,
        # the same reasoning as the strand panel in make_shape.py.
        levels = np.unique(lib)
        discrete = np.allclose(levels, np.round(levels)) and len(levels) <= 30
        if discrete:
            edges = np.arange(g0 - 0.5, g1 + 1.5)
            dl, _ = np.histogram(lib, edges, density=True)
            dh, _ = np.histogram(hit, edges, density=True)
            grid = edges[:-1] + 0.5
            top = max(dl.max(), dh.max())
            axk.fill_between(grid, dl, step="mid", color=LIB, alpha=0.55, lw=0)
            axk.step(grid, dh, where="mid", color=HIT, lw=1.0)
        else:
            pad = 0.05 * (g1 - g0)
            grid = np.linspace(g0 - pad, g1 + pad, 220)
            dl, dh = kde(lib, grid), kde(hit, grid)
            top = max(dl.max(), dh.max())
            axk.fill_between(grid, dl, color=LIB, alpha=0.55, lw=0)
            axk.plot(grid, dh, color=HIT, lw=1.0)
        # the 32 values themselves, so the curve is never mistaken for data
        axk.plot(hit, np.full(len(hit), -0.055 * top), "|", color=HIT,
                 ms=3, mew=0.6, clip_on=False)
        axk.set_xlabel(label, fontsize=6.5, labelpad=1.5)
        axk.set_yticks([])
        axk.set_ylim(-0.10 * top, top * 1.08)
        axk.set_xlim(grid[0], grid[-1])
        axk.spines["left"].set_visible(False)
        axk.tick_params(length=2, width=0.6, pad=1.5, labelsize=6)
        axk.locator_params(axis="x", nbins=4)
    axes[0, 0].annotate("Library", (0.96, 0.88), xycoords="axes fraction",
                        ha="right", fontsize=6.5, color=palette.NEUTRAL["600"])
    axes[0, 0].annotate("Hits", (0.96, 0.68), xycoords="axes fraction",
                        ha="right", fontsize=6.5, color=HIT)
    fig.tight_layout(h_pad=1.6, w_pad=0.9)
    for ext in ("svg", "pdf", "png"):
        fig.savefig(out / f"density_curves.{ext}", bbox_inches="tight")

    print(f"wrote {out}/percentile_strip.* and density_curves.*")
    print(f"null band for the median of 32: {lo:.1f} - {hi:.1f}")
    for label, med, inside in rows:
        print(f"  {label:24s} hit median percentile {med:5.1f}"
              f"  {'inside' if inside else 'OUTSIDE'}")


if __name__ == "__main__":
    main()
