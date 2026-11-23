"""Fig03 metrics panels as violins over the full per-design distributions.

Replaces the bar version (means only), archived at
archive/fig03_campaign_benchmark/metrics_bars_2026-11-05.*.

The bars came from campaign_summary.csv, which holds one transcribed mean per
stack. Violins need the 6,992,515 per-design values, in
inputs/large/campaign_metrics_full.csv.gz - see inputs/MANIFEST.md for how it was
built, for the verified consensus definitions, and for the BG+BIF rescore it
applies.

COLOUR: strictly from config/palette.json via common.palette - no hand-picked
hexes (COLORS.md rule 1). Each stack carries its own fixed colour and the
canonical STACK_ORDER (rule 2), rather than the old "RFD3 in accent, everything
else grey" highlight, which conflicted with stack colours meaning stack identity.
x tick labels double as the direct labels rule 5 asks for, since gold, rose and
light purple converge in greyscale.

PANEL 1 IS A BAR, NOT A VIOLIN. "Designs" is a count: there is no distribution
behind it, and a violin would invent a spread that does not exist. Its bars are
labelled directly, so that panel carries no x axis at all.

LAYOUT IS HORIZONTAL. Five stacks x five panels meant 25 tick labels set at 45
degrees; laid on their side the stack names are written once, horizontally, in
the first panel, and STACK_ORDER reads top to bottom as rule 2 asks (so the y
axis is inverted - row 0 is the top). Metric names are panel titles rather than
rotated y labels.
"""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import collections
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager
from common import paths, style, palette

style.apply()
_have = {f.name for f in font_manager.fontManager.ttflist}
for _f in (style.CFG["font_family"], "Arial", "Helvetica", "Liberation Sans",
           "Nimbus Sans", "DejaVu Sans"):
    if _f in _have:
        matplotlib.rcParams["font.family"] = _f
        break

out = paths.fig_output("Fig2cd_benchmark")

INK = palette.INK
MUT = palette.NEUTRAL["600"]
EDGE = palette.NEUTRAL["700"]
FS = style.CFG["font_size_pt"]   # the one house label size. Ticks, values and
                                 # the key were FS / FS: at 5.5 that is
                                 # 4.5 / 3.5 pt, under the 5 pt floor. They keep
                                 # their muted colour; the size is one.
LW = style.CFG["line_width_pt"]

ORDER = [palette.SHORT[k] for k in palette.STACK_ORDER]   # canonical, never reordered
LABEL = {"RFD3+LiMPNN": "RFD3", "BD1+LiMPNN": "BD1", "HD+LiMPNN": "HD",
         "BG+BIF": "BG+BIF", "BG+LiMPNN": "BG+LiM"}
# (key, panel title, unit, (AF3 col, ESM col)). Order: count bar, GPU-h bar, then
# the split violins.
# Titles are SHORT by the house rule (STYLE_NOTES: a label is the quantity plus
# unit, about 17 characters at 5.5 pt; "Designs (M)" is that rule's own example).
# 2026-11-17: "Designs (millions)" and "Designs per GPU-h" collided with each
# other and with ipTM once the page re-set this panel's text larger - the strip
# went from 5 sub-panels to 6 in the same width, so the titles had to give.
PANELS = [("count", "Designs", "M", None),
          ("gpuh", "Designs/GPU-h", None, None),
          # ligand_RMSD, sc_DockQ, ipTM: the gate order of Figure 3's Venn
          # (ligand_RMSD, sc_DockQ, AA run, ipTM). Keep the two figures in step.
          ("rmsd", "ligand_RMSD", "\u00c5", ("rmsd_af3", "rmsd_esm")),
          ("dockq", "sc_DockQ", None, ("dockq_af3", "dockq_esm")),
          ("iptm", "ipTM", None, ("iptm_af3", "iptm_esm")),
          ("rg", "Rg", "\u00c5", ("rog_af3", "rog_esm"))]
COLS = [c for _, _, _, p in PANELS if p for c in p]
BIF = "BG+BIF"
BIF_TRUE = paths.INPUTS / "large" / "bif_true_per_backend.csv.gz"
BIF_FIX = {"dockq_af3": "dockq_af3_true", "dockq_esm": "dockq_esm_true",
           "rmsd_af3": "rmsd_af3_true", "rmsd_esm": "rmsd_esm_true"}


def bif_true_backbone(df):
    """BG+BIF stored DockQ/RMSD are vs a Boltz re-prediction (wrong). Swap in the
    true-backbone per-backend values; if the source file is absent, blank those
    four columns for BG+BIF (row left empty) rather than plot wrong numbers.
    ponytail: single swap point - change BIF_TRUE/BIF_FIX to change the source."""
    m = df["stack"] == BIF
    if not BIF_TRUE.exists():
        print(f"WARNING: {BIF_TRUE} missing - BG+BIF omitted from sc_DockQ and "
              "ligand_RMSD halves (stored values are the uncorrected re-prediction).")
        df.loc[m, list(BIF_FIX)] = np.nan
        return
    t = pd.read_csv(BIF_TRUE).drop_duplicates("design_id").set_index("design_id")
    ids = df.loc[m, "design_id"]
    for col, tcol in BIF_FIX.items():
        df.loc[m, col] = ids.map(t[tcol]).to_numpy()   # unmatched -> NaN, not stale
    print(f"BG+BIF true-backbone: {ids.isin(t.index).sum():,}/{m.sum():,} rows rescored")


def gpuh():
    """Designs per GPU-h by short stack key; fails loudly on any mismatch."""
    g = pd.read_csv(paths.INPUTS / "Fig3e_gpu_efficiency.csv", comment="#")
    g["short"] = g["stack"].map(palette.SHORT)
    if g["short"].isna().any() or set(g["short"]) != set(ORDER) or len(g) != len(ORDER):
        raise SystemExit(f"GPU-h stacks do not match strip stacks: {list(g['stack'])}")
    return (1.0 / (g["gpuh_gen"] + g["gpuh_refold"])).set_axis(g["short"])


def _half(ax, v, i, colour, top, alpha):
    kde = v if len(v) <= NSUB else rng.choice(v, NSUB, replace=False)
    body = ax.violinplot([kde], positions=[i], vert=False, widths=0.92,
                         showextrema=False)["bodies"][0]
    for p in body.get_paths():     # keep one side of the stack's centre line
        y = p.vertices[:, 1]
        p.vertices[:, 1] = np.minimum(y, i) if top else np.maximum(y, i)
        # read back from the drawn artist: top half's lowest edge / bottom half's
        # highest edge IS the centre line the two halves meet on.
        CENTRE[(ax, i)].append(p.vertices[:, 1].max() if top else p.vertices[:, 1].min())
    body.set_facecolor(colour); body.set_edgecolor("none"); body.set_alpha(alpha)
    q1, med, q3 = np.percentile(v, [25, 50, 75])
    yy = i - 0.17 if top else i + 0.17
    ax.plot([q1, q3], [yy, yy], color=INK, lw=1.0, solid_capstyle="butt", zorder=3)
    ax.plot([med], [yy], "o", ms=2.3, mfc=palette.SURFACE, mec=INK, mew=0.6, zorder=4)
    return np.percentile(v, 0.5), np.percentile(v, 99.5)


HALF_LABELS = []
BARC = {}                             # bar panel key -> row centres
CENTRE = collections.defaultdict(list)   # (violin ax, stack idx) -> [edge, edge]
NSUB = 60_000          # per stack, for the KDE only; stats use every row
rng = np.random.default_rng(0)


def _frame(ax, with_axis=True):
    for sp in ("top", "right", "left"):
        ax.spines[sp].set_visible(False)
    ax.spines["bottom"].set_visible(with_axis)
    ax.spines["bottom"].set_color(palette.NEUTRAL["300"])
    ax.tick_params(axis="y", length=0, labelsize=FS, colors=INK, pad=2)
    ax.tick_params(axis="x", width=LW, length=2, colors=MUT, labelsize=FS,
                   pad=1.5)


# Which generator is diffusion and which is hallucination. This is what panel d
# carried in its hatched legend; d was removed 2026-11-17 (Leo: "you can remove
# panel d now, its data is in the panel c now") so the distinction moves here.
# The two hallucination stacks are the LAST two of palette.STACK_ORDER, so the
# split is contiguous and one bracket per mode covers it. Asserted, not assumed:
# a reordering of STACK_ORDER would otherwise silently mislabel the brackets.
# Keys are palette.SHORT values (what ORDER holds), not the y-axis LABEL text.
MODE = {"RFD3+LiMPNN": "Diffusion", "BG+LiMPNN": "Diffusion",
        "BG+BIF": "Diffusion",
        "BD1+LiMPNN": "Hallucination", "HD+LiMPNN": "Hallucination"}


def _design_mode_brackets(ax, axw_mm):
    """Brackets left of the stack labels grouping diffusion vs hallucination."""
    import matplotlib.transforms as mtr
    missing = [s for s in ORDER if s not in MODE]
    assert not missing, f"no design mode for {missing}; MODE keys are palette.SHORT values"
    modes = [MODE[s] for s in ORDER]
    runs, start = [], 0
    for i in range(1, len(modes) + 1):
        if i == len(modes) or modes[i] != modes[start]:
            runs.append((modes[start], start, i - 1)); start = i
    assert len(runs) == len(set(m for m, _, _ in runs)), (
        f"design modes are not contiguous in STACK_ORDER: {modes}. "
        "Move the bracket to per-stack marks, or restore the grouping.")
    tr = mtr.blended_transform_factory(ax.transAxes, ax.transData)
    # Placed by MEASURED label width, not a fixed axes fraction: with six
    # sub-panels the first axes is narrow and -0.46 landed on the stack labels.
    PT = 25.4 / 72                  # mm per pt
    lab = (max(style.text_w(LABEL[s]) for s in ORDER) + 2 + 4) * PT   # label + tick pad + gap
    X = -lab / axw_mm               # axes fraction of the bracket line
    TICK = 1.2 / axw_mm
    for name, i0, i1 in runs:
        ax.plot([X, X], [i0 - 0.34, i1 + 0.34], transform=tr, color=MUT,
                lw=0.6, clip_on=False, solid_capstyle="butt")
        for yy in (i0 - 0.34, i1 + 0.34):
            ax.plot([X, X + TICK], [yy, yy], transform=tr, color=MUT, lw=0.6,
                    clip_on=False)
        ax.text(X - (FS * 0.7 * PT) / axw_mm, (i0 + i1) / 2, name, transform=tr, rotation=90,
                ha="center", va="center", fontsize=FS, color=MUT)


def check_alignment(fig, axes):
    """Numeric row alignment + half-label collision check; prints, then asserts."""
    pos = {tuple(np.round([a.get_position().y0, a.get_position().y1], 9)) for a in axes}
    ylim = {a.get_ylim() for a in axes}
    assert len(pos) == 1 and len(ylim) == 1, f"sub-panel rows differ: {pos} {ylim}"
    print("all 6 sub-panels share figure y-extent", next(iter(pos)), "and ylim", next(iter(ylim)))
    print(f"{'stack':12s}" + "".join(f"{k:>10s}" for k in BARC) + "   violin-centre per panel")
    for i, s in enumerate(ORDER):
        v = {a: CENTRE[(a, i)] for a in axes if (a, i) in CENTRE}
        print(f"{s:12s}" + "".join(f"{BARC[k][i]:10.3f}" for k in BARC) + "   "
              + " ".join(f"{e[0]:.3f}/{e[1]:.3f}" for e in v.values()))
        for k in BARC:
            assert abs(BARC[k][i] - i) < 1e-9, (k, s, BARC[k][i])
        for a, e in v.items():
            assert len(e) == 2 and all(abs(x - i) < 1e-9 for x in e), (s, e)
    # half key: must not overlap itself, any sub-panel title, or any axes text.
    fig.canvas.draw()
    rd = fig.canvas.get_renderer()
    others = [a.title for a in axes] + [t for a in axes for t in a.texts if t not in HALF_LABELS]
    others += [t for t in axes[0].get_yticklabels()]
    bb = [t.get_window_extent(rd) for t in HALF_LABELS]
    assert not bb[0].overlaps(bb[1]), "AF3 / ESMFold2 labels overlap"
    for t, b in zip(HALF_LABELS, bb):
        for o in others:
            if o.get_text():
                assert not b.overlaps(o.get_window_extent(rd)), (t.get_text(), "collides with", o)
    print("half key: no collisions; widths (mm):",
          {t.get_text(): round(b.width / fig.dpi * 25.4, 2) for t, b in zip(HALF_LABELS, bb)})


def main():
    src = paths.INPUTS / "large" / "campaign_metrics_full.csv.gz"
    df = pd.read_csv(src, usecols=["design_id", "stack"] + COLS)
    print(f"{len(df):,} designs")
    bif_true_backbone(df)
    print(df.groupby("stack")[COLS].mean().round(3).to_string())
    gp = gpuh()

    counts = df.groupby("stack").size()
    colors = [palette.stack_color(s) for s in ORDER]
    y = np.arange(len(ORDER))
    # DRAWN AT ITS PLACED SIZE. The page re-sets every panel's text to 5.5 pt
    # AFTER scaling it into its slot, so a panel drawn much wider than it is
    # placed has its text grow against a bounding box measured before the growth,
    # and the outer labels are clipped on the page though clean standalone.
    # pad_inches cannot fix that: the page crops back to ink.
    # 2026-11-17: this panel briefly shared its row with d and was drawn at
    # 112 mm for that reason. Leo then removed d - its GPU-hours data is the
    # "Designs/GPU-h" sub-panel here now - so c owns the full block again.
    PANEL_W_MM, PANEL_H_MM = 178.0, 40.0
    # 2026-11-17 (Leo): "scooch the violin plots over to the right a bit and
    # compress them a little bit, but have like a proper label to the left of
    # the first violin plot". A dedicated narrow COLUMN carries the key rather
    # than free text in a gutter: the gutter is only ~7.5 mm and "ESMFold2" at
    # 5.5 pt needs ~9.8 mm before its swatch, so free text would have run back
    # over the Designs/GPU-h value labels. Giving it a column also does the
    # scooching and the compressing by itself - the four violins go from 24.3
    # to 22.5 mm each (~7%) and all shift right, with no hand-placed offsets.
    KEY_RATIO = 0.45
    fig, axes_all = plt.subplots(
        1, len(PANELS) + 1, figsize=(style.mm(PANEL_W_MM), style.mm(PANEL_H_MM)),
        gridspec_kw={"width_ratios": [0.78, 1, KEY_RATIO] + [1] * (len(PANELS) - 2)})
    keyax = axes_all[2]
    keyax.set_axis_off()
    axes = [a for i, a in enumerate(axes_all) if i != 2]

    first_violin = next(a for a, P in zip(axes, PANELS) if P[3])
    for ax, (key, title, unit, cols) in zip(axes, PANELS):
        if cols is None:
            vals = [counts[s] / 1e6 if key == "count" else gp[s] for s in ORDER]
            # 2026-11-17 (Leo): "make these bars a bit slimmer". 0.92 left the five
            # bars nearly touching, which read as a stacked block; 0.70 keeps
            # the row pitch (and so the bracket geometry) untouched.
            bars = ax.barh(y, vals, height=0.70, color=colors, linewidth=0)
            BARC[key] = [b.get_y() + b.get_height() / 2 for b in bars]
            for i, v in enumerate(vals):       # direct labels instead of an axis
                ax.text(v + max(vals) * 0.04, i, f"{v:.2f}" if key == "count" else f"{v:.0f}", va="center",
                        fontsize=FS, color=MUT)
            ax.set_xlim(0, max(vals) * 1.26)
            ax.set_xticks([])
        else:
            lo, hi = [], []
            for i, s in enumerate(ORDER):
                for c, top, al in ((cols[0], True, 1.0), (cols[1], False, 0.5)):
                    v = df.loc[df["stack"] == s, c].to_numpy()
                    v = v[np.isfinite(v)]
                    if len(v) == 0:
                        continue
                    l, h = _half(ax, v, i, colors[i], top, al)
                    lo.append(l); hi.append(h)
            # 0.5-99.5 percentile window: tails would flatten the bodies.
            lo, hi = min(lo), max(hi)
            pad = 0.06 * (hi - lo)
            ax.set_xlim(max(0.0, lo - pad), hi + pad)
            if ax is first_violin:
                # Half key, now in its own column to the LEFT of this panel.
                # It shows the two TONES rather than naming "upper"/"lower":
                # every violin is its stack's colour at alpha 1.0 on top and
                # 0.5 beneath, so a neutral swatch at those same two alphas
                # reads across all five stacks without implying any one of
                # them. Stacked top-over-bottom so the key has the same
                # geometry as the thing it explains.
                keyax.text(0.0, 0.60, "Refolded by", transform=keyax.transAxes,
                           ha="left", va="bottom", fontsize=FS, color=MUT)
                # The two tones TOUCH at a shared edge, as the two halves of
                # every violin do. Drawn apart they read as two unrelated
                # swatches; drawn together the key is a miniature of the thing
                # it explains. KH is also the label pitch, so it cannot go
                # below ~2.3 mm without the two names colliding at 5.5 pt.
                KH = 0.065
                for txt, al, col in (("AF3", 1.0, INK), ("ESMFold2", 0.5, MUT)):
                    yk = 0.50 if al == 1.0 else 0.50 - KH
                    keyax.add_patch(plt.Rectangle(
                        (0.0, yk), 0.30, KH, transform=keyax.transAxes,
                        facecolor=INK, alpha=al, lw=0, clip_on=False))
                    HALF_LABELS.append(keyax.text(
                        0.38, yk + KH / 2, txt, transform=keyax.transAxes,
                        ha="left", va="center", fontsize=FS, color=col))

        ax.set_title(title + (f"  ({unit})" if unit else ""), fontsize=FS,
                     color=INK, loc="left", pad=4)
        ax.set_yticks(list(y))
        ax.set_yticklabels([LABEL[s] for s in ORDER] if ax is axes[0] else [])
        ax.set_ylim(len(ORDER) - 0.5, -0.5)     # canonical order, top to bottom
        _frame(ax, with_axis=cols is not None)

    # 0.32, was 0.18: six sub-panels need more gutter than five, and the stack
    # labels on the left of each violin need somewhere to go.
    fig.subplots_adjust(wspace=0.32)
    fig.canvas.draw()
    _design_mode_brackets(axes[0], axes[0].get_window_extent().width / fig.dpi * 25.4)
    check_alignment(fig, axes)
    for ext in ("svg", "pdf", "png"):
        # pad_inches, not a bare tight box: the page re-sets this panel's text
        # LARGER than it is drawn here, so a box that hugs the ink at this size
        # clips the outer labels ("BG+LiM" on the left, "ESMFold2" on the right)
        # once the panel is assembled.
        fig.savefig(out / f"metrics.{ext}", bbox_inches="tight", pad_inches=0.06)
    plt.close(fig)
    print("colours:", {s: palette.stack_color(s) for s in ORDER})
    print("->", out / "metrics.*")


if __name__ == "__main__":
    main()
