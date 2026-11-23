"""Fig4i, the one positive result of the flow-FC analysis: the 32 ordered designs do not
resemble natural Csy4 in fold or sequence, but their interface CHEMISTRY sits part-way
along a direction the natural protein defines.

Two bare panels:

  ladder   five interface-chemistry features, each as a line across four groups —
           the 12k library, the 17 flow non-binders, the 15 flow binders, and Csy4 —
           plotted as a percentile within the library so the five share one axis.
  score    the Csy4-direction score: library distribution, the 32 as points split by
           flow outcome, and Csy4's own position.

Why percentile and not z: the five features have incomparable units (a fraction, a
residue count) and skewed library distributions, so a z would let mean helix length
dominate the panel. Percentile within the library puts them on one bounded axis and
is what "Csy4 sits at the 9.5th percentile on Ala" already means in the text.

The statistic is the MEAN WITHIN-LIBRARY PERCENTILE RANK of a group's members, not
the percentile of the group's mean. Three reasons: the library then sits at exactly
50 by construction, so every line starts from a common baseline and the panel is a
fan rather than five unexplained heights; ranks are immune to the skew and to the
ties that would otherwise wreck one of these features (af3_if_arom is 0 for 51.6% of
the library, so its median percentile is 0 and a median-based panel is unusable);
and it is the same statistic the permutation tests use, so the panel and the P value
measure the same thing. Error bars are 95% bootstrap intervals, 2,000 resamples of
the 15 / 17 designs. The library mark has no interval: it is the whole distribution,
not a sample of it.

Honest note on monotonicity: on this rank statistic the two aromatic features are
NOT strictly monotone — the non-binder step sits a touch BELOW the library (46.8 and
49.3) before the binders rise to 69.4 and 62.0. On group means it reads as a small
rise instead (0.029 -> 0.032). Either way that step is noise; what the panel shows
is the binders displaced from both the library and the non-binders toward Csy4. The
caption should say so rather than claim five monotone ladders.

Base-contact fraction is deliberately absent although it makes the same ladder
(0.365 / 0.384 / 0.449 / 0.561). It is the one feature of the set that does not
replicate on the second refolder — ESM gives 0.389 binders vs 0.403 non, i.e. the
reverse — so it is an AF3 artefact and does not belong on a panel that claims a
direction. Helix fraction and mean helix length do replicate on ESM; Ala and
aromatic content are sequence features and model-independent.

Csy4 = 4AL5 chain A with chain B nt 5-20, the exact design target (1.16 A C1' RMSD
onto the design RNA frame), run through the SAME geom.py pipeline as the 12k.
Source: ~/hits32_browse_flowFC/{csy4.py,csy4score.py}; the stack-stratified
permutation p for the score panel is 0.010 (20,000 permutations, labels shuffled
within stack).

Data: inputs/csy4_ladder_12k.csv.gz (12,000 designs + one Csy4 row).
Writes output/Fig4i_csy4_ladder/{ladder,score}.{svg,pdf,png}.
Nothing is written to final/ - assembly exports there. No panel letters drawn.
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


PAD_IN = 0.04

# --- size variant, for the Fig 5 assembler -------------------------------------
# Called as:  make_figure.py W_MM H_MM FONT_PT OUTDIR
# and then it redraws ONLY the two assembled panels (love_oriented, curves) at that
# geometry with its text already at FONT_PT, so the assembler can place them at
# scale 1.
#
# This CLI is not optional polish. Without it the assembler re-sets type in place:
# it rewrites every font-size from 5.5 to 5.5/scale (8.03 px at d's scale) so the
# label prints at the house size after the panel is shrunk. The x,y of each label
# are NOT recomputed, so a right-anchored y label grows 46% further LEFT - straight
# out of its own viewBox. rsvg will not draw outside the viewBox, so ink() cannot
# even see the overflow, measures a box that is too narrow, and the page clip cuts
# "Aromatic at interface" down to "matic at interface". Padding the canvas only
# buys a character or two. Redrawing is the fix, and it is what
# STYLE_NOTES / the Fig1e and Fig3e panels already do.
_VARIANT = ((float(sys.argv[1]), float(sys.argv[2]), float(sys.argv[3]))
            if len(sys.argv) >= 5 else None)

INK = palette.INK
LINE = palette.NEUTRAL["500"]
LIB = palette.NEUTRAL["400"]
NON = palette.NEUTRAL["700"]
BIND = palette.THEME["600"]
REF = palette.THEME["900"]
FS = _cfg["font_size_pt"]
if _VARIANT:                       # redraw at the assembler's geometry and type size
    FS = _VARIANT[2]
    plt.rcParams["font.size"] = FS
    out = Path(sys.argv[4])
    out.mkdir(parents=True, exist_ok=True)
VAR_SIZE = (lambda w_mm, h_mm: (mm(_VARIANT[0]), mm(_VARIANT[1])) if _VARIANT
            else (mm(w_mm), mm(h_mm)))

# label, column, and the y-nudge in percentile units that keeps the right-hand
# labels apart. Per-point, not one shared offset: a shared one crowds the pair at
# the top and strands the three at the bottom.


def pct_in(lib, v):
    """Percentile of value v within the library distribution."""
    lib = np.sort(np.asarray(lib, float))
    return 100.0 * np.searchsorted(lib, v) / len(lib)


def boot_rank(ranks, n=2000, seed=0):
    """95% bootstrap interval for a group's MEAN percentile rank."""
    rng = np.random.default_rng(seed)
    r = np.asarray(ranks, float)
    return np.percentile(rng.choice(r, (n, len(r)), replace=True).mean(1), [2.5, 97.5])


def load():
    df = pd.read_csv(paths.INPUTS / "csy4_ladder_12k.csv.gz")
    ref = df[df["stack"] == "natural"].iloc[0]
    lib = df[df["stack"] != "natural"].copy()
    assert len(lib) == 12000, len(lib)
    return lib, ref


def kde(x, grid):
    """Gaussian KDE, Silverman bandwidth. scipy in this venv is broken (_propack)."""
    x = np.asarray(x, float)
    spread = min(x.std(ddof=1), (np.percentile(x, 75) - np.percentile(x, 25)) / 1.349)
    bw = 0.9 * max(spread, 1e-9) * len(x) ** (-0.2)
    d = (grid[:, None] - x[None, :]) / bw
    return np.exp(-0.5 * d * d).sum(1) / (len(x) * bw * np.sqrt(2 * np.pi))


def panel_score(lib, ref):
    """The 32 ride ON the library's own density curve, each at its own score.

    A dot's height is therefore not a second variable: it is how common that score
    is in the library, so a design far up the left flank is ordinary and one out on
    the right tail is rare. The x range is clipped to the 0.2-99.8 percentile of the
    library, which drops only the empty right tail and stretches the part that
    carries all 32 designs. Binders are filled and non-binders open, so the two
    groups are not told apart by colour alone (STYLE_NOTES section 4).
    """
    fig, ax = plt.subplots(figsize=(mm(105), mm(46)))
    s = lib.csy4_score.dropna().to_numpy()
    tested = lib[lib.flow_FC.notna()]
    b = np.sort(tested.loc[tested.binder == 1, "csy4_score"].to_numpy())
    n = np.sort(tested.loc[tested.binder == 0, "csy4_score"].to_numpy())

    lo, hi = np.percentile(s, [0.2, 99.8])
    hi = max(hi, ref.csy4_score, b.max(), n.max()) + 0.08
    lo = min(lo, b.min(), n.min()) - 0.08
    g = np.linspace(lo, hi, 600)
    dens = kde(s, g)
    ax.fill_between(g, 0, dens, color=palette.NEUTRAL["200"], lw=0, zorder=1)
    ax.plot(g, dens, color=LIB, lw=0.8, zorder=2)
    top = dens.max()

    def on_curve(v):
        return np.interp(v, g, dens)

    # Csy4 first, so the dots sit over it
    ax.plot([ref.csy4_score, ref.csy4_score], [0, top * 1.05], color=REF, lw=0.9, zorder=3)
    ax.annotate("Csy4", (ref.csy4_score, top * 1.05), xytext=(0, 3),
                textcoords="offset points", fontsize=FS, color=REF, ha="center",
                va="bottom", annotation_clip=False)

    for vals, c, fill, lbl in ((n, NON, "white", "Inactive designs"),
                               (b, BIND, None, "Proposed binders")):
        y = on_curve(vals)
        ax.scatter(vals, y, s=11, facecolors=(c if fill is None else fill),
                   edgecolors=c, linewidths=0.7, zorder=5, clip_on=False)
        ax.plot([vals.mean()], [-top * 0.055], marker="^", ms=3, color=c, zorder=5,
                clip_on=False, lw=0)
        ax.annotate(f"{lbl} (n = {len(vals)})", (1.0, 0.92 if fill is None else 0.78),
                    xycoords="axes fraction", xytext=(4, 0), textcoords="offset points",
                    fontsize=FS, color=c, va="center", ha="left", annotation_clip=False)
    ax.annotate("$P$ = 0.010", (1.0, 0.64), xycoords="axes fraction", xytext=(4, 0),
                textcoords="offset points", fontsize=FS, color=INK, va="center",
                ha="left", annotation_clip=False)

    ax.set_xlabel("Csy4-direction score")
    ax.set_yticks([])
    ax.set_ylim(-top * 0.10, top * 1.30)
    ax.set_xlim(lo, hi)
    ax.tick_params(length=2, width=_cfg["line_width_pt"], direction="out")
    for sp in ("top", "right", "left"):
        ax.spines[sp].set_visible(False)
    ax.spines["bottom"].set_position(("data", -top * 0.10))
    fig.subplots_adjust(left=0.03, right=0.80, top=0.93, bottom=0.20)
    for ext in ("svg", "pdf", "png"):
        fig.savefig(out / f"score.{ext}", bbox_inches="tight")
    plt.close(fig)


def panel_curves(lib, ref):
    """Three overlapping densities on the score axis — library, non-binders, binders —
    each with its mean drawn as a bar up to its own curve.

    ONE BANDWIDTH FOR ALL THREE, and it is the Silverman width for n = 15 (the
    smallest group) using the library's spread. Per-group Silverman would give the
    n = 15 curve a bandwidth 4x the library's purely because it is small, so the
    binder curve would look broad for a reason that has nothing to do with binders.
    With a common width, a difference in shape is a difference in the data. The cost
    is that the library curve is slightly over-smoothed, which is harmless: its true
    shape is in the score panel.

    The rug under the axis is not decoration. A smooth curve drawn through 15 points
    implies a population; the rug says how many points there actually are.
    """
    fig, ax = plt.subplots(figsize=VAR_SIZE(105, 46))
    s = lib.csy4_score.dropna().to_numpy()
    tested = lib[lib.flow_FC.notna()]
    b = tested.loc[tested.binder == 1, "csy4_score"].to_numpy()
    n = tested.loc[tested.binder == 0, "csy4_score"].to_numpy()

    lo = min(s.min(), b.min(), n.min())
    hi = max(np.percentile(s, 99.8), ref.csy4_score, b.max(), n.max()) + 0.1
    lo = max(lo, np.percentile(s, 0.2)) - 0.1
    g = np.linspace(lo, hi, 600)

    spread = min(s.std(ddof=1), (np.percentile(s, 75) - np.percentile(s, 25)) / 1.349)
    BW = 0.9 * spread * min(len(b), len(n)) ** (-0.2)

    def dens(x):
        d = (g[:, None] - np.asarray(x, float)[None, :]) / BW
        return np.exp(-0.5 * d * d).sum(1) / (len(x) * BW * np.sqrt(2 * np.pi))

    # The library is background, not a third series: palest fill, hairline edge, and
    # it sits under everything. The comparison the panel is making is binders vs
    # non-binders; the library is the frame that comparison happens inside.
    series = [("Library (n = 12,000)", s, palette.NEUTRAL["300"], 0.6, True),
              ("Inactive designs (n = 17)", n, NON, 1.0, False),
              ("Proposed binders (n = 15)", b, BIND, 1.2, False)]
    peak = max(dens(v).max() for _, v, _, _, _ in series)
    for lbl, vals, c, lw, fill in series:
        y = dens(vals)
        if fill:
            ax.fill_between(g, 0, y, color=palette.NEUTRAL["100"], lw=0, zorder=0)
        ax.plot(g, y, color=c, lw=lw, zorder=3 if not fill else 1,
                solid_capstyle="round")
        # mean: a short tick riding the curve at that x, not a bar to the baseline
        mu = vals.mean()
        ym = np.interp(mu, g, y)
        h = peak * 0.18
        ax.plot([mu, mu], [ym - h, ym + h], color=c, lw=0.7, zorder=5 if not fill else 1,
                solid_capstyle="butt")
        # the designs themselves, riding their own group's curve
        if not fill:
            ax.scatter(vals, np.interp(vals, g, y), s=9, facecolors=c, edgecolors="white",
                       linewidths=0.4, zorder=6, clip_on=False)

    rug = peak * 0.04

    ax.plot([ref.csy4_score, ref.csy4_score], [0, peak * 1.06], color=REF, lw=0.9,
            zorder=2)
    ax.annotate("Csy4", (ref.csy4_score, peak * 1.06), xytext=(0, 3),
                textcoords="offset points", fontsize=FS, color=REF, ha="center",
                va="bottom", annotation_clip=False)
    # The library's curve is deliberately below 3:1 against white; its LABEL is not,
    # or the key would be unreadable (STYLE_NOTES section 9).
    for i, (lbl, vals, c, _, fill) in enumerate(series):
        ax.annotate(lbl, (1.0, 0.92 - 0.14 * i), xycoords="axes fraction", xytext=(4, 0),
                    textcoords="offset points", fontsize=FS,
                    color=palette.NEUTRAL["600"] if fill else c, va="center",
                    ha="left", annotation_clip=False)
    ax.annotate("$P$ = 0.010", (1.0, 0.92 - 0.14 * 3), xycoords="axes fraction",
                xytext=(4, 0), textcoords="offset points", fontsize=FS, color=INK,
                va="center", ha="left", annotation_clip=False)

    ax.set_xlabel("Csy4-direction score")
    ax.set_yticks([])
    ax.set_ylim(-rug * 1.3, peak * 1.30)
    ax.set_xlim(lo, hi)
    ax.tick_params(length=2, width=_cfg["line_width_pt"], direction="out")
    for sp in ("top", "right", "left"):
        ax.spines[sp].set_visible(False)
    ax.spines["bottom"].set_position(("data", -rug * 1.3))
    fig.subplots_adjust(left=0.03, right=0.78, top=0.93, bottom=0.20)
    for ext in ("svg", "pdf", "png"):
        fig.savefig(out / f"curves.{ext}", bbox_inches="tight", pad_inches=PAD_IN)
    plt.close(fig)
    print(f"curves: common bandwidth {BW:.3f}; means lib {s.mean():+.3f} "
          f"non {n.mean():+.3f} bind {b.mean():+.3f} Csy4 {ref.csy4_score:+.3f}")


def panel_love(lib, ref, oriented=False):
    """Love plot / standardised-difference dumbbell — the replacement for the ladder.

    `oriented=True` emits the one-direction version (`love_oriented`): every feature
    is multiplied by sign(Csy4 - library), so +x always means "toward Csy4" and all
    seven diamonds sit on the right. A row whose dots go negative is one where the
    designs move AWAY from the natural protein. The arrow on each row label says
    which way the RAW feature moves to get there, which is the only thing the sign
    flip throws away.

    NOT log2, which is what was asked for first and does not work on this data.
    A log2 fold change needs a ratio, and Csy4's value is exactly 0 for two of the
    seven features (no Asp/Glu and no Lys anywhere in its 30 interface residues), so
    those rows would be log2(0) = -inf and drop out - and they are the two most
    interesting rows. The ratios would also be unstable for the other near-zero
    composition features. A log does not fix direction anyway: log2 is still signed
    about 0, so the rows would point both ways exactly as before. Orientation is the
    thing that makes them agree, and SD units are what keep the near-zero features
    on the panel at all.

    One row per feature, x in LIBRARY SD UNITS with the library at zero. Two dots
    per row, non-binders and binders, joined so the shift reads as one move, each
    with a 95% bootstrap interval. Csy4's own position is a diamond on the same
    scale. Rows sorted by the binder effect.

    Why this instead of the four-column slopegraph: the slopegraph gave equal
    billing to a column that is 50 by construction (the library, zero variance) and
    one that is a single protein, implying a four-step trajectory among comparable
    groups. Here the library is the origin, Csy4 is one mark, and the only thing
    drawn as a sample-with-uncertainty is the thing that is one. Percentile rank
    also flattened effect size; SD units keep it.

    It also takes more rows for free, so the two features the ladder had to leave
    out are back: interface Asp/Glu and interface Lys. Both earn their place.
    Non-binders sit +0.74 SD on interface Asp/Glu while Csy4 is at -0.87 and the
    binders at -0.23, so the dumbbell straddles zero in the Csy4 direction. And on
    interface Lys BOTH groups sit above the library (+0.17, +0.20) while Csy4 is at
    -1.56: the one feature where the designs move AWAY from the natural protein,
    which is the Arg-not-Lys gap stated as a picture.
    """
    FEATS = [("Aromatic fraction", "arom"), ("Aromatic at interface", "af3_if_arom"),
             ("Lys at interface", "af3_if_K"), ("Asp/Glu at interface", "af3_if_neg"),
             ("Helix fraction", "af3_sse_helix"), ("Ala fraction", "aa_A"),
             ("Mean helix length", "af3_mean_helix_len")]
    tested = lib[lib.flow_FC.notna()]
    bi, no = tested[tested.binder == 1], tested[tested.binder == 0]

    rows = []
    for lbl, c in FEATS:
        mu, sd = lib[c].mean(), lib[c].std()
        dc = (ref[c] - mu) / sd
        sg = 1.0 if (not oriented or dc >= 0) else -1.0
        d = lambda v: sg * (np.asarray(v, float).mean() - mu) / sd
        ci = lambda v, seed: np.sort(sg * (np.percentile(
            np.random.default_rng(seed).choice(np.asarray(v, float),
                                               (2000, len(v)), replace=True).mean(1),
            [2.5, 97.5]) - mu) / sd)
        rows.append(dict(label=lbl + ("" if not oriented else (" ↑" if dc >= 0 else " ↓")),
                         non=d(no[c]), bind=d(bi[c]), csy4=sg * dc,
                         non_ci=ci(no[c], 1), bind_ci=ci(bi[c], 2)))
    R = pd.DataFrame(rows).sort_values("bind", ascending=False).reset_index(drop=True)

    fig, ax = plt.subplots(figsize=VAR_SIZE(105, 50))
    y = np.arange(len(R))[::-1]
    ax.axvline(0, color=palette.NEUTRAL["400"], lw=0.7, zorder=1)
    for yi, r in zip(y, R.itertuples()):
        ax.plot([r.non, r.bind], [yi, yi], color=palette.NEUTRAL["300"], lw=1.4,
                zorder=2, solid_capstyle="round")
        for lo_hi, c in ((r.non_ci, NON), (r.bind_ci, BIND)):
            ax.plot(lo_hi, [yi, yi], color=c, lw=0.6, zorder=3, solid_capstyle="butt")
        ax.scatter([r.non], [yi], s=16, facecolors="white", edgecolors=NON,
                   linewidths=0.9, zorder=5)
        ax.scatter([r.bind], [yi], s=16, facecolors=BIND, edgecolors=BIND,
                   linewidths=0.9, zorder=5)
        ax.scatter([r.csy4], [yi], s=20, marker="D", facecolors=REF, edgecolors=REF,
                   linewidths=0, zorder=4)

    ax.set_yticks(y)
    ax.set_yticklabels(R.label)
    ax.set_ylim(-0.7, len(R) - 0.3)
    ax.set_xlabel("Shift toward Csy4 (library SD)" if oriented
                  else "Difference from library (SD)")
    ax.tick_params(axis="y", length=0)
    ax.tick_params(axis="x", length=2, width=_cfg["line_width_pt"], direction="out")
    for sp in ("top", "right", "left"):
        ax.spines[sp].set_visible(False)
    key = [("Inactive designs (n = 17)", NON, "white"), ("Proposed binders (n = 15)", BIND, None),
           ("Csy4", REF, None)]
    for i, (lbl, c, fc) in enumerate(key):
        ax.annotate(lbl, (1.07, 0.95 - 0.12 * i), xycoords="axes fraction", xytext=(8, 0),
                    textcoords="offset points", fontsize=FS, color=c, va="center",
                    ha="left", annotation_clip=False)
        ax.scatter([1.07], [0.95 - 0.12 * i], transform=ax.transAxes, clip_on=False,
                   s=20 if i == 2 else 16, marker="D" if i == 2 else "o",
                   facecolors=(fc or c), edgecolors=c, linewidths=0.9, zorder=6)
    fig.subplots_adjust(left=0.22, right=0.76, top=0.97, bottom=0.17)
    for ext in ("svg", "pdf", "png"):
        fig.savefig(out / f"love_oriented.{ext}", bbox_inches="tight", pad_inches=PAD_IN)
    plt.close(fig)
    print(R.round(3).to_string(index=False))


def panel_features(lib, ref):
    """The five features of the Csy4-direction score, one small panel each, in RAW units.

    Step histograms, NOT kernel densities. Every one of these is a count over a
    denominator, so they live on a grid and carry point masses: aa_A and arom are
    counts over the 80-residue design (0.0125 steps, 44 and 22 distinct values), and
    the three interface shares are counts over the interface, with a spike at exactly
    zero holding 46.9% (Asp/Glu), 51.6% (aromatic) and 7.8% (Lys) of the library. A
    KDE smears those spikes into humps and invents density between grid points; the
    bar at zero IS the finding for two of these features. Same lesson as the Fig3d
    point mass.

    Binning: features on the 1/80 grid get one bar per level; the three ratio
    features get 26 equal bins from 0 to the library's 99.5th percentile, which is
    why their longest tail is not drawn (nothing of interest is out there - Csy4 is
    at 0 on two of them and the group means are all well inside).

    The 32 designs sit above the histogram as two dot rows rather than on the curve:
    a bar chart has no curve to ride, and at this panel width a strip is legible
    where a beeswarm would not be.
    """
    FEATS = [("Ala fraction", "aa_A"), ("Aromatic fraction", "arom"),
             ("Asp/Glu at interface", "af3_if_neg"), ("Lys at interface", "af3_if_K"),
             ("Aromatic at interface", "af3_if_arom")]
    tested = lib[lib.flow_FC.notna()]
    bi = tested[tested.binder == 1]
    no = tested[tested.binder == 0]

    fig, axes = plt.subplots(1, len(FEATS), figsize=(mm(178), mm(42)))
    for ax, (lbl, c) in zip(axes, FEATS):
        v = lib[c].dropna().to_numpy()
        u = np.unique(v)
        if len(u) <= 25:                      # exact grid: one bar per level
            step = np.min(np.diff(u))
            edges = np.append(u - step / 2, u[-1] + step / 2)
        else:
            edges = np.linspace(0, np.percentile(v, 99.5), 27)
        h, e = np.histogram(np.clip(v, edges[0], edges[-1]), bins=edges)
        h = h / h.max()
        ax.stairs(h, e, fill=True, color=palette.NEUTRAL["100"], zorder=1)
        ax.stairs(h, e, fill=False, color=palette.NEUTRAL["300"], lw=0.6, zorder=2)

        ax.plot([ref[c], ref[c]], [0, 1.08], color=REF, lw=0.9, zorder=3)
        for vals, y, col, fc in ((no[c], 1.26, NON, "white"), (bi[c], 1.46, BIND, None)):
            ax.scatter(vals, np.full(len(vals), y), s=7, facecolors=(fc or col),
                       edgecolors=col, linewidths=0.6, zorder=5, clip_on=False)
            ax.plot([vals.mean(), vals.mean()], [y - 0.10, y + 0.10], color=col,
                    lw=0.7, zorder=6, clip_on=False)
        ax.set_xlim(edges[0], edges[-1])
        ax.set_ylim(0, 1.62)
        ax.set_yticks([])
        ax.set_xlabel(lbl)
        ax.tick_params(axis="x", length=2, width=_cfg["line_width_pt"], direction="out")
        ax.locator_params(axis="x", nbins=4)
        for sp in ("top", "right", "left"):
            ax.spines[sp].set_visible(False)

    key = [("Library (n = 12,000)", palette.NEUTRAL["400"], palette.NEUTRAL["100"], "s"),
           ("Inactive designs (n = 17)", NON, "white", "o"),
           ("Proposed binders (n = 15)", BIND, None, "o"), ("Csy4", REF, None, "|")]
    x0 = 0.0
    for lbl, col, fc, mk in key:
        axes[0].scatter([x0], [-0.46], transform=axes[0].transAxes, clip_on=False,
                        s=18, marker=mk, facecolors=(fc or col), edgecolors=col,
                        linewidths=0.9, zorder=6)
        axes[0].annotate(lbl, (x0, -0.46), xycoords="axes fraction", xytext=(5, 0),
                         textcoords="offset points", fontsize=FS,
                         color=palette.NEUTRAL["600"] if mk == "s" else col,
                         va="center", ha="left", annotation_clip=False)
        x0 += (len(lbl) * 0.055 + 0.12)
    fig.subplots_adjust(left=0.02, right=0.99, top=0.97, bottom=0.34, wspace=0.22)
    for ext in ("svg", "pdf", "png"):
        fig.savefig(out / f"features.{ext}", bbox_inches="tight")
    plt.close(fig)
    print("features: " + " | ".join(
        f"{c} Csy4 {ref[c]:.3f} lib {lib[c].mean():.3f} non {no[c].mean():.3f} "
        f"bind {bi[c].mean():.3f}" for _, c in FEATS))


def main():
    lib, ref = load()
    if _VARIANT:                 # assembler asked for just the two placed panels
        panel_curves(lib, ref)
        panel_love(lib, ref, oriented=True)
        return
    panel_score(lib, ref)
    panel_curves(lib, ref)
    panel_love(lib, ref, oriented=True)
    panel_features(lib, ref)
    print("wrote", out)


if __name__ == "__main__":
    main()
