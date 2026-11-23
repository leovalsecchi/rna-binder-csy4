"""Variant of make_figure.py with the eight callout marks omitted.

The composite draws its own circles, so the baked dots and numbers are redundant
and sit underneath them. Everything else, including the extent and RECT, is
identical, so the two PNGs are registered.

ORIGINAL HEADER FOLLOWS
"""
"""Fig02c: the campaign map as an assembly component.

This writes the MAP ONLY: the thinned dot cloud, the eight numbered marks, and
the stack key. No structure panels, no leaders, no binder captions. The eight
renderings are separate files already, at
output/Fig2b_stack_binders/callout_*.png, and the zoom-ins are assembled by
hand in Illustrator.

Two variants, same frame and same marks, so they drop onto each other:
  map_only.*              dots alone
  map_only_territories.*  dots with each generator's density shaded and outlined

Neither carries the stack key: it is placed at assembly time with the
renderings. The only text is the one-line note on thinning at the foot, which
is a single text object if it needs to move.

Numbers on the map match the rendering filenames:
  1 callout_hd_a    2 callout_bd1_a   3 callout_rfd3_a  4 callout_hd_b
  5 callout_bg_a    6 callout_rfd3_b  7 callout_bd1_b   8 callout_bg_b

Three decisions, and why:

1. THE DOT CLOUD IS THINNED. 499,064 marks in a 110 mm square is a solid field
   in which an 8-point callout cannot be seen. The map is binned at 0.12 UMAP
   units and at most CAP backbones are kept per bin, so every sparse arm keeps
   all of its points and only the saturated masses lose any. The silhouette is
   unchanged because the parts that were thinned were already opaque. This is a
   drawing decision; no count anywhere comes from the thinned set.

2. THE DOTS GO OUT AS TRUE VECTOR. That is only possible because of the
   thinning: 41,072 circles is an SVG Illustrator will open, 499,064 is not.
   Every mark therefore stays selectable and recolourable.

3. TERRITORIES ARE A SEPARATE FILE. The smoothing that makes four masses read
   as territories also inflates BoltzGen's dozens of small islands, so its
   outlined area runs about five times the others' per backbone. See README.

Data: inputs/large/fold_map_scatter.parquet, inputs/map8_callouts.csv.
Writes output/Fig2b_map_composite/map_only{,_territories}.{svg,pdf,png}.
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
import matplotlib.patheffects as pe
from common import paths, palette

FIG = Path(__file__).resolve().parent.name
out = paths.fig_output(FIG)
_cfg = yaml.safe_load((paths.ROOT / "config" / "style.yaml").read_text())
plt.rcParams.update({
    "font.family": _cfg["font_family"], "pdf.fonttype": 42, "svg.fonttype": "none",
    "savefig.dpi": 600, "axes.grid": False,
})
def mm(x): return x / 25.4


def gaussian_filter(Z, sigma):
    """Separable gaussian blur; scipy in this venv fails to dlopen _propack."""
    r = int(4 * sigma)
    k = np.exp(-0.5 * (np.arange(-r, r + 1) / sigma) ** 2)
    k /= k.sum()
    Z = np.apply_along_axis(lambda v: np.convolve(v, k, "same"), 0, Z)
    return np.apply_along_axis(lambda v: np.convolve(v, k, "same"), 1, Z)


def hdr_level(z, h, mass):
    """Density level enclosing `mass` of this generator's backbones.

    2026-11-17 (Leo): replaces the old rule, which put each generator's contour
    at 5% of ITS OWN peak density. Under that rule the enclosed areas were not
    comparable between generators: a generator with a lower peak got a
    proportionally lower cut, so the area partly measured peak height rather
    than spread. The folder README conceded the shading was "only a reading aid".

    This is the standard highest-density region: rank cells by smoothed density,
    walk down that ranking accumulating REAL counts, and stop where the
    cumulative count first reaches `mass` of the total. The outline is then the
    smallest region holding `mass` of this generator - which means the same
    thing for every generator, so areas can be compared and the caption is one
    sentence.

    It does NOT equalise the areas, and should not: BoltzGen still needs ~3.6x
    more area per backbone than RFD3. That ratio is now a result - BoltzGen's
    backbones are more dispersed in fold space - rather than an artefact of
    where each peak happened to sit.

    `h` must be the count grid in the SAME orientation as `z`; the callers
    transpose the smoothed histogram, so they pass H.T.
    """
    order = np.argsort(z, axis=None)[::-1]
    total = h.sum()
    if total <= 0:
        return z.max()
    cum = np.cumsum(h.ravel()[order])
    hit = min(int(np.searchsorted(cum, mass * total)), order.size - 1)
    return float(z.ravel()[order][hit])


BIN, CAP = 0.12, 6                      # thinning grid and per-bin cap
GB, SIGMA, HDR_MASS = 0.10, 3.0, 0.90   # territory grid, smoothing, enclosed share
                                        # 0.90: the outline holds 90% of that
                                        # generator's backbones (see hdr_level)
FW = 170                                # mm, square-ish panel

MERGE = {"BG+BIF": "BoltzGen"}
RENAME = {"BG+LiMPNN": "BoltzGen"}
COL, ORDER = {}, []
for long in palette.STACK_ORDER:
    short = palette.SHORT[long]
    name = RENAME.get(short, MERGE.get(short, short))
    if name not in COL:
        COL[name] = palette.stack_color(long if short not in MERGE else
                                        "BoltzGen + LigandMPNN")
        ORDER.append(name)
INK, HALO = palette.INK, palette.SURFACE

d = pd.read_parquet(paths.INPUTS / "large" / "fold_map_scatter.parquet")
d["stack"] = d["stack"].astype(str)
NISO = int(d["deg0"].sum())
d = d[~d["deg0"]].copy()
d["group"] = d["stack"].replace({**RENAME, **MERGE})
N = d.groupby("group").size().to_dict()
X, Y, G = d["x"].values, d["y"].values, d["group"].values

ix = ((X - X.min()) / BIN).astype(np.int64)
iy = ((Y - Y.min()) / BIN).astype(np.int64)
key = ix * 1_000_000 + iy
order = np.argsort(key, kind="stable")
ks = key[order]
start = np.r_[0, np.flatnonzero(np.diff(ks)) + 1]
cnt = np.diff(np.r_[start, len(ks)])
rng = np.random.default_rng(0)
keep = np.concatenate([order[s:s + n] if n <= CAP
                       else rng.choice(order[s:s + n], CAP, replace=False)
                       for s, n in zip(start, cnt)])
rng.shuffle(keep)
print(f"{len(d):,} backbones thinned to {len(keep):,} marks "
      f"({100*len(keep)/len(d):.1f}%)")

X0, X1, Y0, Y1 = X.min(), X.max(), Y.min(), Y.max()
padx, pady = 0.03 * (X1 - X0), 0.03 * (Y1 - Y0)

ex = np.arange(X0 - padx, X1 + padx + GB, GB)
ey = np.arange(Y0 - pady, Y1 + pady + GB, GB)
xc, yc = (ex[:-1] + ex[1:]) / 2, (ey[:-1] + ey[1:]) / 2
DENS, AREA, LV = {}, {}, {}
for g in ORDER:
    m = G == g
    H, _, _ = np.histogram2d(X[m], Y[m], bins=[ex, ey])
    Z = gaussian_filter(H, SIGMA).T
    DENS[g] = Z
    LV[g] = hdr_level(Z, H.T, HDR_MASS)     # H.T: counts in Z's orientation
    AREA[g] = float((Z > LV[g]).sum()) * GB * GB

c = pd.read_csv(paths.INPUTS / "map8_callouts.csv").sort_values("n")
c["group"] = c["stack"].replace({**RENAME, **MERGE})
XY = d[["x", "y"]].values
for _, r in c.iterrows():
    dd = np.hypot(XY[:, 0] - r.x, XY[:, 1] - r.y).min()
    assert dd < 0.02, f"{r.panel} is {dd:.3f} from the nearest plotted backbone"
    print(f"  {int(r.n)} callout_{r.panel:7s} ({r.x:5.2f}, {r.y:5.2f})  "
          f"fold {int(r.cluster_r10)}")

ASP = (Y1 - Y0 + 2 * pady) / (X1 - X0 + 2 * padx)
FH = FW                                 # square panel
MAP_H = 0.93                            # no key, so the map fills the panel
MAP_W = MAP_H * FH / (ASP * FW)
RECT = [(1 - MAP_W) / 2, 0.045, MAP_W, MAP_H]


def build(territories: bool, stem: str):
    fig = plt.figure(figsize=(mm(FW), mm(FH)))
    ax = fig.add_axes(RECT)
    ax.set_xlim(X0 - padx, X1 + padx)
    ax.set_ylim(Y0 - pady, Y1 + pady)
    ax.set_aspect("equal")
    for sp in ax.spines.values():
        sp.set_visible(False)
    ax.set_xticks([]); ax.set_yticks([])

    if territories:                     # widest first, so compact ones stay visible
        for g in sorted(ORDER, key=lambda g: -AREA[g]):
            Z = DENS[g]
            lv = LV[g]
            ax.contourf(xc, yc, Z, levels=[lv, Z.max() * 1.01], colors=[COL[g]],
                        alpha=0.11, zorder=1)
            ax.contour(xc, yc, Z, levels=[lv], colors=[COL[g]], linewidths=0.5,
                       alpha=0.6, zorder=2)
    # true vector: 41,072 circles opens in Illustrator, 499,064 does not
    ax.scatter(X[keep], Y[keep], s=0.45, c=[COL[g] for g in G[keep]],
               alpha=0.50, linewidths=0, rasterized=False, zorder=3)

    # callout marks are NOT drawn in this variant; the composite circles them

    note = (f"{len(d):,} unique backbones, {NISO:,} unclustered isolates not "
            f"shown; drawn as {len(keep):,} thinned marks.")
    if territories:
        note += (" Shaded areas are smoothed density per generator, a reading "
                 "aid over the dots, not proportional to backbone count.")
    fig.text(0.5, 0.018, note, fontsize=5, color=palette.NEUTRAL["600"],
             ha="center", va="center")

    for ext in ("svg", "pdf", "png"):
        fig.savefig(out / f"{stem}.{ext}", transparent=False,
                    dpi=600 if ext == "png" else None)
    plt.close(fig)
    size = (out / f"{stem}.svg").stat().st_size / 1e6
    print(f"wrote {stem}.svg/pdf/png  ({size:.1f} MB svg)")


build(True, "map_only_territories_nomarks")
for g in ORDER:
    print(f"  territory {g:14s} {AREA[g]:7.1f} units^2, "
          f"{1000*AREA[g]/N[g]:5.2f} per 1,000 backbones "
          f"(encloses {HDR_MASS:.0%})")
print("->", out)
