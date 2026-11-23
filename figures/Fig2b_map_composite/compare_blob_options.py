"""Decision aid: three ways of drawing the generator "territories" on the fold map.

Written 2026-11-17 for Leo, who asked what the options are for making the blobs
"more clean and standardized". It does NOT change the production panel; it draws
the three candidates side by side into working/ so the choice can be made by eye.

The problem being fixed
-----------------------
The shipped panel contours each generator at 5% of ITS OWN peak density, so the
areas are not comparable between generators and are not proportional to backbone
count. BoltzGen's territory is 129.1 units against RFD3's 40.6, but per 1,000
backbones that is 0.86 against 0.28 - the ranking survives, the ratio does not.
The folder README already concedes the shading is "only a reading aid".

The three panels
----------------
1. CURRENT      contour at 0.05 * that generator's own peak density.
2. OPTION 2     contour at the level enclosing HDR_MASS of that generator's own
                backbones (a highest-density region). Still per-generator, but
                now every outline means the same thing - "the smallest region
                holding 90% of this generator" - so areas are comparable and the
                caption is one sentence.
3. OPTION 5     no contours: one grid cell takes the colour of whichever
                generator holds most backbones in it, with opacity on log total
                count. Territory and overlap in a single honest layer, nothing
                to explain, but density within a generator is gone.

Run:
    ~/.venvs/binder_map/bin/python code/Fig2b_map_composite/compare_blob_options.py
Writes working/_blob_options/blob_options.{svg,png} and prints each generator's
territory area under options 1 and 2, with the per-1,000-backbone figures that
show whether the areas have become comparable.

scipy will not import in this venv (propack dlopen), so the Gaussian smoothing is
the same hand-rolled separable np.convolve the production script uses.
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import numpy as np                                                  # noqa: E402
import pandas as pd                                                 # noqa: E402
import matplotlib.pyplot as plt                                     # noqa: E402
from matplotlib.colors import to_rgb                                # noqa: E402
from common import palette, style                                   # noqa: E402

style.apply()

REPO = Path(__file__).resolve().parents[2]
SRC = REPO / "inputs" / "large" / "fold_map_scatter.parquet"
OUT = REPO / "working" / "_blob_options"

GB, SIGMA = 0.10, 3.0          # grid step in UMAP units, smoothing in grid cells
LEVEL_CURRENT = 0.05           # the shipped rule: 5% of the generator's own peak
HDR_MASS = 0.90                # option 2: enclose this share of its backbones
DOT, DOT_A = 0.6, 0.05         # reference scatter under panels 1 and 2

# One group per generator. Both BoltzGen stacks are one generator, as in the
# production panel; BG+BIF is 107 backbones and shares BG+LiMPNN's colour.
GROUP = {"RFD3+LiMPNN": "RFD3", "BG+LiMPNN": "BoltzGen", "BG+BIF": "BoltzGen",
         "BD1+LiMPNN": "BD1", "HD+LiMPNN": "HD"}
COL = {"RFD3": palette.stack_color("RFD3 + LigandMPNN"),
       "BoltzGen": palette.stack_color("BoltzGen + LigandMPNN"),
       "BD1": palette.stack_color("BoltzDesign + LigandMPNN"),
       "HD": palette.stack_color("HalluDesign + LigandMPNN")}
ORDER = ["RFD3", "BoltzGen", "BD1", "HD"]


def gaussian1d(a, sigma, axis):
    """Separable Gaussian blur; scipy.ndimage is unavailable in this venv."""
    r = int(4 * sigma + 0.5)
    k = np.exp(-0.5 * (np.arange(-r, r + 1) / sigma) ** 2)
    k /= k.sum()
    return np.apply_along_axis(lambda m: np.convolve(m, k, mode="same"), axis, a)


def smooth(h):
    return gaussian1d(gaussian1d(h, SIGMA, 0), SIGMA, 1)


def hdr_level(z, h, mass):
    """Density level enclosing `mass` of the backbones.

    Cells are ranked by smoothed density and real counts accumulated down that
    ranking, so the level is where the cumulative count first reaches the mass.
    That is the standard highest-density-region construction: the outline is the
    SMALLEST region holding `mass` of this generator.
    """
    order = np.argsort(z, axis=None)[::-1]
    cum = np.cumsum(h.ravel()[order])
    total = h.sum()
    if total <= 0:
        return z.max()
    hit = np.searchsorted(cum, mass * total)
    hit = min(hit, order.size - 1)
    return z.ravel()[order][hit]


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    d = pd.read_parquet(SRC, columns=["x", "y", "stack", "deg0"])
    d = d[~d["deg0"]]                       # degree-0 isolates are not plotted
    d["g"] = d["stack"].map(GROUP)
    print(f"{len(d):,} backbones plotted ({d['g'].value_counts().to_dict()})")

    ex = np.arange(d["x"].min() - 1, d["x"].max() + 1 + GB, GB)
    ey = np.arange(d["y"].min() - 1, d["y"].max() + 1 + GB, GB)
    xc, yc = (ex[:-1] + ex[1:]) / 2, (ey[:-1] + ey[1:]) / 2
    cell = GB * GB

    H, Z, N = {}, {}, {}
    for g in ORDER:
        m = d["g"] == g
        N[g] = int(m.sum())
        H[g], _, _ = np.histogram2d(d["x"][m], d["y"][m], bins=[ex, ey])
        Z[g] = smooth(H[g])

    fig, axes = plt.subplots(1, 3, figsize=(style.mm(178), style.mm(62)))
    for ax in axes:
        ax.set_xticks([]); ax.set_yticks([]); ax.set_aspect("equal")
        for s in ax.spines.values():
            s.set_visible(False)

    rows = []
    for ax, title, mode in ((axes[0], "Current: 5% of own peak", "peak"),
                            (axes[1], f"Option 2: encloses {HDR_MASS:.0%}", "hdr")):
        ax.scatter(d["x"], d["y"], s=DOT, c="#b3b8bc", alpha=DOT_A,
                   linewidths=0, rasterized=True)
        lv = {g: (LEVEL_CURRENT * Z[g].max() if mode == "peak"
                  else hdr_level(Z[g], H[g], HDR_MASS)) for g in ORDER}
        area = {g: float((Z[g] > lv[g]).sum()) * cell for g in ORDER}
        # widest first, so a big territory never hides a small one
        for g in sorted(ORDER, key=lambda q: -area[q]):
            ax.contourf(xc, yc, Z[g].T, levels=[lv[g], Z[g].max() * 1.01],
                        colors=[COL[g]], alpha=0.11, zorder=1)
            ax.contour(xc, yc, Z[g].T, levels=[lv[g]], colors=[COL[g]],
                       linewidths=0.5, alpha=0.6, zorder=2)
        ax.set_title(title, fontsize=6.0, color=palette.INK, pad=3)
        rows.append((title, area))

    # --- option 5: whichever generator holds the cell, opacity on log count ---
    ax = axes[2]
    stack = np.dstack([H[g] for g in ORDER])
    tot = stack.sum(axis=2)
    win = np.argmax(stack, axis=2)
    rgb = np.zeros(tot.shape + (4,))
    for i, g in enumerate(ORDER):
        rgb[win == i, :3] = to_rgb(COL[g])
    with np.errstate(divide="ignore"):
        a = np.log1p(tot)
    rgb[..., 3] = np.clip(a / np.percentile(a[tot > 0], 99), 0, 1)
    rgb[tot == 0, 3] = 0
    ax.imshow(np.transpose(rgb, (1, 0, 2)), origin="lower", interpolation="nearest",
              extent=[ex[0], ex[-1], ey[0], ey[-1]], aspect="equal", zorder=1)
    ax.set_title("Option 5: cell takes its majority generator", fontsize=6.0,
                 color=palette.INK, pad=3)

    # direct labels, no legend box (STYLE_NOTES). Under the first panel, where
    # the map has no ink - inside it they sat on BoltzGen's satellite contours.
    # Widths are MEASURED, not estimated at a fixed em-per-char: "BoltzGen
    # 149,933" and "BD1 113,175" differ enough in width that even spacing
    # collides them.
    # text_w returns the ink width in the same unit as the size, so points here;
    # the row is laid out in figure coordinates against the figure's width in
    # points. Even spacing collides "BoltzGen  149,933" with "BD1  113,175".
    KEY_PT, KEY_GAP_PT = 5.5, 14.0
    keys = [f"{g}  {N[g]:,}" for g in ORDER]
    fig_pt = fig.get_figwidth() * 72.0
    x = 0.01
    for g, k in zip(ORDER, keys):
        fig.text(x, 0.015, k, fontsize=KEY_PT, color=COL[g], va="bottom", ha="left")
        x += (style.text_w(k, KEY_PT) + KEY_GAP_PT) / fig_pt

    fig.tight_layout(pad=0.4)
    for ext in ("svg", "png"):
        # bbox_inches: the titles sit above the axes and are clipped without it
        fig.savefig(OUT / f"blob_options.{ext}", dpi=300, bbox_inches="tight")
    print("wrote", OUT / "blob_options.png")

    print("\nterritory area (UMAP units^2), and per 1,000 backbones:")
    print(f"{'generator':<10}{'n':>9}" + "".join(f"{t.split(':')[0]:>26}" for t, _ in rows))
    for g in ORDER:
        line = f"{g:<10}{N[g]:>9,}"
        for _, area in rows:
            line += f"{area[g]:>14.1f}{area[g] / N[g] * 1000:>12.2f}"
        print(line)
    print("\nIf the two right-hand numbers in a column track n, the areas are "
          "comparable; if they do not, that column's rule is normalising each "
          "generator against itself.")


if __name__ == "__main__":
    main()
