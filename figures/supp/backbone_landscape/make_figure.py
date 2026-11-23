"""backbone_landscape - the backbone fold landscape through selection.

ONE panel, THREE maps in a row, all on MAIN-TEXT FIG 2b's OWN EMBEDDING, so a
position means the same backbone here and in the main text, and a map differs
from its neighbour only by which backbones are still there. Each stage REMOVES
the backbones it filtered out rather than greying them.

2026-11-19 (Leo): "I want to use exactly the campaing umap, like it is in main
fig2b, then remove the dots that go filtered out and so on. I want to keep the
clustering we have in the main conserved. Also with the colors."

Done. How far, exactly:

  * THE EMBEDDING IS FIG 2b's. The previous version drew the union Foldseek
    layout and its README argued the swap was impossible. THAT ARGUMENT WAS WRONG
    ON ITS DECISIVE POINT AND IS WITHDRAWN. `inputs/large/fold_map_scatter.parquet`
    already carries a `stage` column, so all three drawn maps need NO JOIN AT
    ALL; the cut ipTM arm was placed on Euler by rebuilding that same table with
    one extra flag (`make_inputs.py`, `--fig2b`), and that flag is still read,
    because the arm is still measured.
  * IT IS PROVED, NOT CLAIMED. `fig2b_layers.parquet` is the Fig 2b table plus
    `in_counterfactual`, and `load()` ASSERTS that its x, y, stage, cluster_r10,
    stack and deg0 reproduce `fold_map_scatter.parquet` row for row. If the
    embedding ever drifts, the figure stops building.
  * THE CLUSTERING IS CONSERVED BECAUSE IT IS THE SAME CLUSTERING: the same
    551,516 nodes carrying the same Leiden `cluster_r10`. Nothing is re-clustered
    here, and this figure's own k=30 r8.0 / 917-fold clustering no longer draws
    anything - it is only quoted, for the coverage numbers, in `report()`.
  * THE COLOURS ARE FIG 2b's, including its BoltzGen merge: BoltzIF is an inverse
    folding step on BoltzGen backbones, so on this partition the two design routes
    are the same backbone (BG+BIF survives as 107/551,516 nodes, 0.02%) and the
    merged group keeps the BG+LiMPNN hex. Copied from Fig2b_map_callouts; no new
    hex (COLORS.md rules 1 and 2).
  * THE 52,452/551,516 (9.5%) ISOLATES ARE DROPPED, as Fig 2b drops them: degree 0
    in the thresholded fold graph, so UMAP's random init throws them onto the
    perimeter and their angular position carries no meaning.
  * THE MARKS ARE FIG 2b's TOO, as of 2026-11-19. `fig2b_marks()` READS `s` and
    `alpha` out of Fig2b_map_callouts' own scatter call instead of keeping a local
    copy, so the two cannot drift apart silently - the same prove-it-every-run
    principle as `load()`. Leo: "for the umap just take the umap from the main
    figure and use its code and display setting for low number binders." The old
    local DOT = 0.22 was a shrunken size invented here, and it left the 26k / 8k /
    6.6k maps sparse and faint. ONE marker size across every map - the thinning
    is the message and must not be flattered by growing the dots on the maps
    that have fewer of them.
  * FILTERED-OUT BACKBONES ARE STILL REMOVED, NOT GREYED. That is Leo's own
    instruction (above), and it is WHY MAP 1 HAS TO EXIST: with no grey field
    under maps 2-4, nothing else stands in for the pool they came from. The main
    text's own low-N idiom, Fig2b_map_composite's eight callouts, is a per-point
    treatment - a big ringed mark over the full coloured cloud - for EIGHT points;
    it does not transpose to 26,312. So the low-N fix here is Fig 2b's marker, not
    Fig 2b's callout. Nothing was dropped and nothing was greyed.
  * EVERY MAP STATES BOTH COUNTS, as of 2026-11-19. Leo, reading the backbone
    counts as design counts: "I'm not sure why there's 3 maps after the first one.
    The first subsample should be the 55,000 we filter to. And then the second
    subsample should be the one with 12,000 designs." It already was: map 2 IS the
    55,706 gated designs and map 3 IS the 12,000-design library, collapsed onto
    26,312 and 8,312 drawn backbones. Only the labels never said so. The design
    counts are derived in `designs()` from `funnel_real_gate.csv` and the coverage
    table, never typed in here, and the gate is cross-checked between the two.
  * MAP 4, THE ipTM COUNTERFACTUAL, IS CUT - 2026-11-20, Leo: "remove the iptm
    counterfactual and all the text below in grey." It was never a stage (the
    top 12,000 designs by ipTM alone; nothing passed through it) and Leo never
    asked for it. WHAT WENT WITH IT, recorded because it is a result and not a
    decoration: the ipTM arm contracts to 69.0% of the gated region while the
    ordered panel does not contract at all (52.7 against 53.2 units^2, 99.2%),
    and in the map that contraction was VISIBLE - the arm lost almost every
    peripheral island the other three keep. THE CLAIM APPEARS IN NO OTHER
    FIGURE OR CAPTION in this repo, so on 2026-11-20 it was put into THIS
    figure's caption as prose, marked STATED, NOT DRAWN - the idiom bif_effect
    uses for its chi-squared and repredictor_comparison for its ipTM-step
    rates. That was the COORDINATING AGENT'S call applied on Leo's behalf, not
    Leo's; the README records it that way and how to reverse it. The sentence
    reports the areas and stops: the VISIBLE half - the arm losing almost every
    peripheral island - is what a number cannot carry, and the caption must not
    be edited into implying the reader has seen it. The area is MEASURED on
    every run, by `outline(None, ...)` in `figure()` and printed by `report()`,
    so the caption's number and the data's cannot drift - but nothing draws it.
  * THE GREY NOTE STRIP UNDER THE ROW IS CUT with it, same instruction. Its
    load-bearing half - "is map 1 the whole campaign?" - is now on MAP 1's OWN
    LABEL, which states 499,064 of 551,516 backbones drawn beside the
    6,992,515 designs, and in the caption, which carries the isolate drop and
    "Nothing else is filtered out of it". STYLE_NOTES 0 was always against the
    strip: explanatory notes on the canvas belong to the caption.
  * THE FOLD-COVERAGE ARGUMENT IS NOT DRAWN. 917/917 (100%) folds for the ordered
    panel against 326/917 (35.6%) for the ipTM arm is measured on the k=30 r8.0
    map, a DIFFERENT partition from the one drawn here, so putting a k30 fold count
    beside a bbmap backbone count on the same map would mix partitions. The counts
    on the maps are the backbones you can see; the coverage numbers are re-derived
    in `report()` and belong to the caption, which names the partition once
    (STYLE_NOTES 0: anything the caption can carry is not drawn).

FIG 2b's CASCADE IS THE APPLIED GATE, not the unapplied `all_methods.py` one.
Settled 2026-11-19 against the files; see README. `stage` is the FURTHEST stage a
backbone reached, so the gated set is stage >= 1 (26,658), not stage == 1 (18,106).

No design-route colour means anything but design-route identity (COLORS.md rule
2); no hatch (rule 7); the one non-route mark is the outline, in MAIN, a single
series (rule 1).

Inputs are built on Euler by `make_inputs.py`, because the embeddings, the layer
flags and the cluster labels never leave the cluster.
"""
import json
import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from common import paths, style, palette

FIG = Path(__file__).resolve().parent.name
FIG2B = paths.INPUTS / "large" / "fold_map_scatter.parquet"   # the main-text table
LAYERS2B = paths.INPUTS / "large" / "fig2b_layers.parquet"    # + in_counterfactual
FIG2B_SRC = Path(__file__).resolve().parents[2] / "Fig2b_map_callouts" / "make_figure.py"
LAYERS = paths.INPUTS / "large" / "backbone_landscape_layers.parquet"
COVER = paths.INPUTS / "backbone_landscape_coverage.csv"
FUNNEL = paths.INPUTS / "funnel_real_gate.csv"   # the campaign + gate DESIGN counts
QUANT = paths.INPUTS / "backbone_landscape_quantities.json"
STEM = "landscape_stages"

FW = 178                # mm, double column: three maps in a row as ONE panel
N_MAPS = 3              # was 4 until 2026-11-20; see the ipTM note in the docstring
PAD_MM, GAP_MM = 1.5, 2.0
KEY_MM = 4.0            # strip under the row for the design-route key
KEY_PT = 5.5            # the key ships WITH the panel, as it does on
KEY_SW_MM, KEY_GAP_MM, KEY_SEP_MM = 1.3, 0.9, 4.0   # Fig2b_map_composite
LABEL_MM = 8.0          # strip over the row for the three-line per-map label
PAD = 0.03              # share of the data span left as margin, as in Fig 2b
GB, SIGMA, HDR_MASS = 0.10, 3.0, 0.90   # Fig2b_map_composite's own settings
N_FOLD = 917            # fold clusters in the gated pool, k=30 r8.0 (caption only)

# Fig 2b's own cascade, on Fig 2b's own partition, asserted on every run.
CASCADE = {"campaign": 551_516, "gated": 26_658, "ordered": 8_552}
N_ISOLATE = 52_452
# asserted on the union partition, which this figure no longer draws but whose
# numbers its README compares against
LAYER_SIZES = {"in_sample500k": 499_999, "in_passing60k": 31_763,
               "in_12k": 9_759, "in_counterfactual": 7_813}


def fig2b_marks():
    """Marker size and alpha out of MAIN-TEXT FIG 2b's OWN scatter call.

    Read from its source, not copied into a constant here. That file is a script,
    not a module, so importing it would rebuild the main-text map; and a copied
    number drifts silently the first time Fig 2b is retuned. Same principle as
    `load()`: prove the reuse on every run rather than claim it in a comment.
    """
    m = re.search(r"ax\.scatter\([^)]*?\bs=([0-9.]+).*?\balpha=([0-9.]+)",
                  FIG2B_SRC.read_text(), re.S)
    assert m, f"no scatter call carrying s= and alpha= in {FIG2B_SRC}"
    s, a = float(m.group(1)), float(m.group(2))
    assert 0 < s < 5 and 0 < a <= 1, (s, a)
    return s, a


DOT, ALPHA = fig2b_marks()      # pt^2 and alpha, identical in every map


# Fig 2b's colour rule, copied from Fig2b_map_callouts: fixed colour per design
# route, fixed order, the two BoltzGen routes collapsed in place because on this
# partition they are the same backbone. The merged group keeps the BG+LiMPNN hex;
# no new hex is invented (COLORS.md rules 1 and 2).
MERGE = {"BG+BIF": "BoltzGen"}
RENAME = {"BG+LiMPNN": "BoltzGen"}


def route_colors():
    col, order = {}, []
    for long in palette.STACK_ORDER:
        short = palette.SHORT[long]
        name = RENAME.get(short, MERGE.get(short, short))
        if name not in col:
            col[name] = palette.stack_color(
                long if short not in MERGE else "BoltzGen + LigandMPNN")
            order.append(name)
    return col, order


# --- borrowed from code/Fig2b_map_composite/make_figure.py, unchanged ---------
# Copied rather than imported: that file is a script, not a module - importing it
# would rebuild the main-text map. Copied rather than hoisted into common/, which
# would mean editing a main-text figure for a supplementary one.

def gaussian_filter(Z, sigma):
    """Separable gaussian blur; scipy in this venv fails to dlopen _propack, so
    every statistic in this file is hand-rolled on numpy."""
    r = int(4 * sigma)
    k = np.exp(-0.5 * (np.arange(-r, r + 1) / sigma) ** 2)
    k /= k.sum()
    Z = np.apply_along_axis(lambda v: np.convolve(v, k, "same"), 0, Z)
    return np.apply_along_axis(lambda v: np.convolve(v, k, "same"), 1, Z)


def hdr_level(z, h, mass):
    """Density level enclosing `mass` of this stage's backbones.

    The standard highest-density region: rank cells by smoothed density, walk down
    that ranking accumulating REAL counts, and stop where the cumulative count
    first reaches `mass` of the total. The outline is then the smallest region
    holding `mass` of this stage, which means the same thing for every stage, so
    the outlines are comparable and the caption is one sentence.

    `h` must be the count grid in the SAME orientation as `z`; the caller
    transposes the smoothed histogram, so it passes H.T.
    """
    order = np.argsort(z, axis=None)[::-1]
    total = h.sum()
    if total <= 0:
        return z.max()
    cum = np.cumsum(h.ravel()[order])
    hit = min(int(np.searchsorted(cum, mass * total)), order.size - 1)
    return float(z.ravel()[order][hit])
# -----------------------------------------------------------------------------


def load():
    """The Fig 2b table with `in_counterfactual` added, PROVED identical to the
    main-text one on every column the main text draws."""
    d = pd.read_parquet(LAYERS2B)
    ref = pd.read_parquet(FIG2B)
    assert len(d) == len(ref) == CASCADE["campaign"], (len(d), len(ref))
    shared = ["x", "y", "stage", "cluster_r10", "stack", "deg0"]
    key = ["x", "y", "cluster_r10", "stage"]
    a = d[shared].sort_values(key, kind="stable").reset_index(drop=True)
    b = ref[shared].sort_values(key, kind="stable").reset_index(drop=True)
    for c in shared:                    # category vs str: compare as strings
        assert (a[c].astype(str).values == b[c].astype(str).values).all(), \
            f"fig2b_layers.parquet does not reproduce the main-text map on {c}"
    assert {"campaign": len(d), "gated": int((d.stage >= 1).sum()),
            "ordered": int((d.stage == 2).sum())} == CASCADE
    assert int(d.deg0.sum()) == N_ISOLATE, int(d.deg0.sum())
    print(f"embedding verified against {FIG2B.name}: {len(d):,d} backbones, "
          f"all {len(shared)} drawn columns identical")
    return d


def main():
    style.apply()
    # 300 dpi turns a half-million-point cloud to mush at 42 mm; the dots are the
    # only rasterised thing here, the outline and the type stay vector.
    plt.rcParams["savefig.dpi"] = 600
    out = paths.fig_output(FIG)
    d = load()
    cov = pd.read_csv(COVER).set_index("arm")
    q = json.loads(QUANT.read_text())

    # The union-partition numbers this figure's README compares against, still
    # asserted even though they no longer draw anything.
    u = pd.read_parquet(LAYERS)
    assert len(u) == 526_467, len(u)
    assert dict(u[list(LAYER_SIZES)].sum()) == LAYER_SIZES
    assert (u.in_passing60k == u.in_pool_k30).all()
    assert (u.in_12k == u.in_ordered_k30).all()
    assert q["fold_pool"] == q["fold_ordered"] == N_FOLD

    areas, cf_area = figure(d, out, designs(cov))
    report(d, u, cov, q, areas, cf_area)


def stages(d):
    """The three maps, in reading order: drawn name, mask over the drawn rows, and
    the key naming that map's DESIGN count in `designs()`.

    The ipTM counterfactual was the fourth until 2026-11-20 and is CUT. It was
    never a stage - nothing passed through it - and Leo never asked for it; see
    the docstring for the one result that left with it.
    """
    return [("The campaign", d.stage >= 0, "campaign"),
            ("Five applied gates", d.stage >= 1, "gated"),
            ("Ordered panel", d.stage == 2, "ordered")]


def designs(cov):
    """The four DESIGN counts, read from the inputs, never typed in here.

    Four, not three: the ipTM arm no longer draws a map but report() still
    prints its contraction, so its design count is still needed.

    2026-11-19 (Leo), reading the backbone counts as design counts: "I'm not sure
    why there's 3 maps after the first one. The first subsample should be the
    55,000 we filter to. And then the second subsample should be the one with
    12,000 designs." It already is both: map 2 IS the 55,706 gated designs and
    map 3 IS the 12,000-design library; they collapse onto 26,312 and 8,312 drawn
    backbones, and only the labels never said so.

    The campaign total and the gate come off the campaign funnel, the two 12,000
    arms off the coverage table, and the gate is cross-checked between the two.
    """
    f = pd.read_csv(FUNNEL).set_index("step")["total"]
    n = {"campaign": int(f.iloc[0]), "gated": int(f.iloc[-1]),
         "ordered": int(cov.loc["ordered", "designs"]),
         "counterfactual": int(cov.loc["counterfactual", "designs"])}
    assert n["gated"] == int(cov.loc["pool", "designs"]), \
        (n["gated"], cov.loc["pool", "designs"])    # two files, one gate
    return n


def limits(d):
    """Fig 2b's own window: full extent of the drawn cloud, 3% pad, aspect equal.
    Not forced square - Fig 2b's window is not square, and registering with it is
    the whole point."""
    x0, x1, y0, y1 = d.x.min(), d.x.max(), d.y.min(), d.y.max()
    px, py = PAD * (x1 - x0), PAD * (y1 - y0)
    return (x0 - px, x1 + px), (y0 - py, y1 + py)


def outline(ax, hit, xlim, ylim):
    """The one non-route mark: the smallest smoothed-density region holding
    HDR_MASS of THIS stage's backbones. Recomputed per map, which is the point.

    `ax=None` MEASURES the region without drawing it. That is how the cut ipTM
    arm's area survives the loss of its map: the contraction is still computed
    on this same grid and printed by report(), rather than frozen into a number
    typed in a README - the same treatment metric_reliability gives the caption
    clause of its own cut panel.
    """
    ex = np.arange(xlim[0], xlim[1] + GB, GB)
    ey = np.arange(ylim[0], ylim[1] + GB, GB)
    H, _, _ = np.histogram2d(hit.x, hit.y, bins=[ex, ey])
    Z = gaussian_filter(H, SIGMA).T                 # H.T: counts in Z's orientation
    lv = hdr_level(Z, H.T, HDR_MASS)
    if ax is None:
        return float((Z > lv).sum()) * GB * GB
    # Fig2b_map_composite's own weight and alpha, and BENEATH the marks as it is
    # there. At 42 mm the 90% region of a strongly multi-modal cloud breaks into
    # ~25 components, so an opaque line over the dots reads as noise; receded, it
    # reads as the territory it is.
    ax.contour((ex[:-1] + ex[1:]) / 2, (ey[:-1] + ey[1:]) / 2, Z, levels=[lv],
               colors=[palette.MAIN], linewidths=0.5, alpha=0.6, zorder=1)
    return float((Z > lv).sum()) * GB * GB


def draw_key(fig, col, order, fh):
    """One small field per DESIGN ROUTE, under the row.

    The key ships WITH the panel for the same reason Fig2b_map_composite's does
    (2026-11-16, Leo): without it a standalone export shows three unexplained
    colours, and the colour is the one thing on these maps a reader cannot decode
    from the maps themselves. Labels are Fig 2b's own, so the two read as one
    system; no counts, because they differ per map and the maps carry their own.
    Widths are MEASURED, not estimated from a per-character guess.
    """
    from matplotlib.patches import Rectangle
    mmpt = 25.4 / 72.0
    tw = [style.text_w(g, KEY_PT) * mmpt for g in order]
    entry = [KEY_SW_MM + KEY_GAP_MM + w for w in tw]
    total = sum(entry) + KEY_SEP_MM * (len(entry) - 1)
    assert total < FW - 2 * PAD_MM, f"key is {total:.1f} mm, wider than the panel"
    kax = fig.add_axes([0, 0, 1, KEY_MM / fh])
    kax.set_xlim(0, FW)
    kax.set_ylim(0, 1)
    kax.axis("off")
    x = (FW - total) / 2
    sh = KEY_SW_MM / KEY_MM                     # square swatch, in axis units
    for g, w, e in zip(order, tw, entry):
        kax.add_patch(Rectangle((x, 0.5 - sh / 2), KEY_SW_MM, sh,
                                facecolor=col[g], lw=0, zorder=1))
        kax.text(x + KEY_SW_MM + KEY_GAP_MM, 0.50, g, fontsize=KEY_PT,
                 color=palette.INK, ha="left", va="center")
        x += e + KEY_SEP_MM


def figure(d, out, nd):
    col, order = route_colors()
    niso = int(d.deg0.sum())
    d = d[~d.deg0].copy()                           # Fig 2b drops the isolates
    d["group"] = d["stack"].astype(str).replace({**RENAME, **MERGE})
    ndrawn = len(d)
    ntot = niso + ndrawn
    xlim, ylim = limits(d)
    asp = (ylim[1] - ylim[0]) / (xlim[1] - xlim[0])

    w = (FW - 2 * PAD_MM - (N_MAPS - 1) * GAP_MM) / N_MAPS   # width of one map
    h = w * asp
    fh = h + KEY_MM + LABEL_MM
    fig = plt.figure(figsize=(style.mm(FW), style.mm(fh)))
    rng = np.random.default_rng(0)                  # fixed seed, as in Fig 2b
    areas = []
    for i, (lab, mask, key) in enumerate(stages(d)):
        ax = fig.add_axes([(PAD_MM + i * (w + GAP_MM)) / FW,
                           KEY_MM / fh, w / FW, h / fh])
        hit = d[mask]
        k = np.arange(len(hit))
        rng.shuffle(k)                              # no route drawn last on top
        ax.scatter(hit.x.values[k], hit.y.values[k], s=DOT,
                   c=[col[g] for g in hit.group.values[k]], alpha=ALPHA,
                   linewidths=0, rasterized=True, zorder=2)
        areas.append((lab, len(hit), outline(ax, hit, xlim, ylim)))
        print(f"  {lab:34s} {nd[key]:9,d} designs -> {len(hit):7,d} drawn backbones")

        # BOTH counts, every map: the design count the map represents and the
        # backbone count actually drawn. The share is of the drawn backbones,
        # which is what the eye is comparing.
        #
        # Map 1 states the backbone TOTAL as well, "499,064 of 551,516 drawn".
        # That is the surviving half of the grey note cut on 2026-11-20: the
        # question it answered - "these 499,000, these are the full campaign,
        # right?" - cost a round of confusion once, so the map that provokes it
        # carries the answer itself. The isolates and the "nothing else is
        # filtered out of it" are in the caption (README), which is where a
        # reader meets the rest.
        bb = (f"{len(hit):,d} of {ntot:,d} backbones drawn" if i == 0
              else f"{len(hit):,d} backbones ({100 * len(hit) / ndrawn:.1f}%)")
        txt = f"{lab}\n{nd[key]:,d} designs\n{bb}"
        for line in txt.split("\n"):        # measured, never estimated
            assert style.text_w(line, style.CFG["font_size_pt"]) * 25.4 / 72 < w, \
                f"label line wider than the {w:.1f} mm map: {line}"
        ax.annotate(txt, (0.5, 1.0), xycoords="axes fraction", xytext=(0, 2),
                    textcoords="offset points", va="bottom", ha="center",
                    linespacing=1.35, fontsize=style.CFG["font_size_pt"],
                    color=palette.INK)
        ax.set_xlim(*xlim)
        ax.set_ylim(*ylim)
        ax.set_aspect("equal")
        ax.set_xticks([])
        ax.set_yticks([])
        for sp in ax.spines.values():
            sp.set_visible(False)

    draw_key(fig, col, order, fh)
    # The grey note strip under the row was CUT on 2026-11-20 (Leo: "remove the
    # iptm counterfactual and all the text below in grey"). Its content did not
    # evaporate: the campaign-scope chain is on map 1's own label and in the
    # caption, and the HDR mass is in the caption. STYLE_NOTES 0 was always
    # against it - explanatory notes on the canvas belong to the caption.
    #
    # The ipTM arm's outline is still MEASURED, on the same grid as the three
    # drawn ones, so report() can print the contraction it used to show.
    cf_area = outline(None, d[d.in_counterfactual], xlim, ylim)
    # not style.save: bbox_inches="tight" would trim the authored 178 mm width,
    # and an assembled page lays this out from that width. SVG is emitted because
    # page assembly consumes SVG and matplotlib keeps labels as live text there
    # (svg.fonttype: none) where a converted PDF arrives as glyph outlines.
    for ext in ("pdf", "png", "svg"):
        fig.savefig(out / f"{STEM}.{ext}")
    plt.close(fig)
    print(f"wrote {STEM}.pdf/png/svg  {FW:.0f} x {fh:.1f} mm, "
          f"{N_MAPS} {w:.1f} x {h:.1f} mm maps")
    for lab, nn, a in areas:
        print(f"  {lab:34s} {nn:7,d}/{ndrawn:,d} "
              f"({100 * nn / ndrawn:5.1f}%)  outline {a:7.1f} UMAP units^2 "
              f"({HDR_MASS:.0%} of that stage)")
    return areas, cf_area


def report(d, u, cov, q, areas, cf_area):
    """Every number the caption carries, printed so RUNLOG can quote it."""
    drawn = d[~d.deg0]
    print("\n-- Fig 2b partition (drawn here), from its own `stage` column")
    print(f"  all {len(d):,d}  gated {int((d.stage >= 1).sum()):,d}  "
          f"ordered {int((d.stage == 2).sum()):,d}  "
          f"ipTM arm {int(d.in_counterfactual.sum()):,d}")
    print(f"  furthest-stage buckets {d.stage.value_counts().sort_index().to_dict()}"
          "  <- 18,106 is stage==1 ONLY; the gated set is stage>=1 = 26,658")
    print(f"  isolates {int(d.deg0.sum()):,d}/{len(d):,d} "
          f"({100 * int(d.deg0.sum()) / len(d):.1f}%); drawn {len(drawn):,d}")
    print(pd.crosstab(d["stack"], d["stage"]).to_string())
    cf, st = d.in_counterfactual, d.stage
    print(f"  ipTM arm vs the cascade: in the gated set {int((cf & (st >= 1)).sum()):,d}"
          f"/{int(cf.sum()):,d}, in the ordered panel {int((cf & (st == 2)).sum()):,d}"
          f"/{int(cf.sum()):,d}")
    print(f"  ordered-only backbones {int(((st == 2) & ~cf).sum()):,d}, "
          f"ipTM-only {int((cf & (st != 2)).sum()):,d}")

    print("\n-- union partition (not drawn; the numbers the README compares)")
    print("  layer sizes:", u.filter(like="in_").sum().to_dict())
    for a, b in (("in_sample500k", "in_passing60k"), ("in_sample500k", "in_12k"),
                 ("in_passing60k", "in_12k")):
        print(f"  {a} & {b}: {int((u[a] & u[b]).sum()):,d}")

    print("\n-- the diversity argument, on the k=30 r8.0 map (caption, not drawn)")
    print(cov.to_string())
    o, c = cov.loc["ordered"], cov.loc["counterfactual"]
    for lab, key in (("folds", "fold"), ("epitopes", "epi"),
                     ("(fold, epitope) pairs", "pair")):
        n = q[f"{key}_pool"]
        print(f"  {lab}: ordered {q[f'{key}_ordered']}/{n} "
              f"({100 * q[f'{key}_ordered'] / n:.1f}%) vs counterfactual "
              f"{q[f'{key}_cf']}/{n} ({100 * q[f'{key}_cf'] / n:.1f}%)")
    print(f"  mean ipTM ordered {o.mean_iptm:.4f} vs counterfactual "
          f"{c.mean_iptm:.4f}, cost {c.mean_iptm - o.mean_iptm:.4f}")
    print(f"  backbones: ordered only {q['nodes_ordered_only']:,d}, both "
          f"{q['nodes_both']:,d}, counterfactual only {q['nodes_cf_only']:,d}, "
          f"net {q['nodes_net']:,d}")
    print(f"  fold clusters reached by the ordered panel ONLY: "
          f"{q['fold_ordered_only']} (the counterfactual's folds are a strict "
          f"subset: counterfactual-only = {q['fold_cf_only']})")

    print("\n-- outline areas on THIS embedding. The union-partition figures "
          "39.5 / 47.7 / 27.8 were\n   measured on a different layout and do NOT "
          "carry over; these replace them.")
    for lab, _, a in areas:
        print(f"  {lab:34s} {a:7.1f} UMAP units^2")
    print(f"  {'ipTM counterfactual (NOT DRAWN)':34s} {cf_area:7.1f} UMAP units^2")
    print(f"  the ipTM arm contracts to {100 * cf_area / areas[1][2]:.1f}% "
          f"of the gated region; the ordered panel to "
          f"{100 * areas[2][2] / areas[1][2]:.1f}% - the one result the cut "
          f"fourth map carried. It is STATED, NOT DRAWN: these four areas are "
          f"the caption's, and the caption is the only place in the repo this "
          f"claim appears, so re-check it against this line if either moves.")

    # The three circulating numbers this figure was asked to re-derive.
    vals = set(q.values())
    for n, what in ((591, "fold clusters reached by the ordered panel only"),
                    (5_144, "backbones in the ordered panel and not the counterfactual"),
                    (1_946, "net backbone difference, ordered - counterfactual")):
        print(f"  {n}: {'REPRODUCES = ' + what if n in vals else 'DOES NOT REPRODUCE'}")
    for n in (528, 1_515):
        assert n not in vals
        print(f"  {n}: DOES NOT REPRODUCE - no quantity on this pool equals it")


if __name__ == "__main__":
    main()
