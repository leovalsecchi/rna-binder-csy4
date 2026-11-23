"""bif_effect - what the BoltzIF scoring defect did to the delivered panel.

BoltzGen + BoltzIF's stored dockq_* / sc_rmsd_* were scored against BoltzIF's own
Boltz RE-PREDICTION, not against the generated backbone. The *_true columns are the
correction. STORED = campaign history, what the gate actually saw; CORRECTED = what
those backbones really are. The two are never mixed inside one number here.

Three bare panels. (a) and (b) are the 2,397 rescored members of the 2,400-design
BoltzGen + BoltzIF slice of the delivered panel. (c) is the WHOLE ROUTE, all
1,435,617 designs, as a bare four-quadrant density - see panel_c().

⚠️ The decile curve and the quadrant density SWAPPED LETTERS on 2026-11-20 (Leo:
"move c up next to a and relabel it"), so the ordered-panel pair (a) and (b) now sit
side by side and the double-column route panel is (c). The book letters panels by
FILE-NAME order, so the letter lives in the emitted stem and nowhere else; nothing
either panel draws changed.

Encoding: one design route (BoltzGen + BoltzIF) and no route comparison, so there is
nothing for a route colour to distinguish - every panel is a single series in MAIN
(COLORS.md rule 1), and panel (c)'s density ramp is the house theme through
style.density_cmap(). Panel (c) carries no overlay at all (Leo, 2026-11-19: "remove
all the dots. Just make it a clean hex plot with the 4 quadrants"), so the ramp is
the only encoding in it. Panel (a) plots the POOLED AF3 + ESMFold2 convention
and panel (c) the gate's BEST-OF-TWO one; neither separates the two backends, so no
re-predictor colour (COLORS.md rule 7) is needed and there is no hatch anywhere in
this figure - the withdrawn texture encoding was never used here. Per-backend
values: README.

scipy is unusable in this venv (_propack dlopen failure), so the chi-square and its
p-value are hand-rolled below.
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from math import exp, factorial

import numpy as np
import pandas as pd
import pymupdf
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection

from common import paths, style, palette

FIG = Path(__file__).resolve().parent.name
LIB = paths.INPUTS / "csy4_binder_library_12k_v6_k30.csv"
TRUE = paths.INPUTS / "large" / "bif_true_per_backend.csv.gz"
# The route's STORED per-backend scores, all 1,480,140. The library CSV carries
# stored scores for the 2,400 ordered designs only, which is exactly why panel
# (c) could not hold a stored-FAIL point until this file existed.
ROUTE_STORED = paths.INPUTS / "large" / "bif_route_stored_per_backend.csv.gz"
STORED_COLS = ["rmsd_af3", "rmsd_esm", "dockq_af3", "dockq_esm"]
TRUE_COLS = ["rmsd_af3_true", "rmsd_esm_true", "dockq_af3_true", "dockq_esm_true"]
ROUTE = "BG+BIF"                 # the literal `stack` value in the CSV
RMSD_GATE, DOCKQ_GATE = 3.0, 0.7
ROUNDS = ["v2", "v4", "v5", "v6", "v7", "v8", "v9"]   # seven; there is no v3 here
N_DEC = 10


# --- data -------------------------------------------------------------------
def load():
    """The route's panel members, joined to their true-backbone scores."""
    lib = pd.read_csv(LIB)
    tru = pd.read_csv(TRUE)
    panel = lib[lib["stack"] == ROUTE].merge(tru, on="design_id", how="left",
                                             validate="one_to_one")
    tcols = ["rmsd_af3_true", "rmsd_esm_true", "dockq_af3_true", "dockq_esm_true"]
    scored = panel[tcols].notna().all(axis=1)
    # The three unscored rows are PRESENT with empty values, not absent, and they
    # are not one thing: see the README. Asserted so a changed source is noticed.
    assert len(panel) == 2400 and scored.sum() == 2397, (len(panel), int(scored.sum()))
    d = panel[scored].copy()
    d["stored_pass"] = gate(d.sc_rmsd_af3, d.sc_rmsd_esm, d.dockq_af3, d.dockq_esm)
    d["true_pass"] = gate(d.rmsd_af3_true, d.rmsd_esm_true,
                          d.dockq_af3_true, d.dockq_esm_true)
    d["flip"] = d.stored_pass & ~d.true_pass
    assert d.stored_pass.all(), "the panel is drawn from the gated set"
    return d


def load_route():
    """The WHOLE route, stored scores joined to corrected ones, ordered flagged.

    Panel (c) needs the stored-FAIL half of the plane and the ordered panel
    cannot hold one point there - it IS the gated set, so the gate wall is the
    edge of its data and two quadrants are empty by construction.

    ROUTE_STORED is one awk pass over the campaign master. That master's
    CONSENSUS sc_dockq/sc_rmsd are the true-backbone rescore, but its
    PER-BACKEND columns are the campaign's own stored values and were left
    alone on purpose ("the raw per-backend columns are kept as-is so the
    correction stays visible", export_campaign_metrics.py). Verified against
    the library's own stored columns on all 2,400 ordered designs: max |d| is
    one 4-dp rounding unit.
    """
    d = pd.read_csv(ROUTE_STORED).merge(pd.read_csv(TRUE), on="design_id",
                                        how="inner", validate="one_to_one")
    assert len(d) == 1_480_140, len(d)
    lib = pd.read_csv(LIB, usecols=["design_id", "stack"])
    ordered = set(lib.loc[lib["stack"] == ROUTE, "design_id"])
    assert len(ordered) == 2400, len(ordered)
    # Cross-checks against files built by OTHER scripts, taken here on the full
    # 1,480,140 rather than after the filter below, because that is the
    # population those files counted - comparing a narrowed one to them would
    # be comparing two different things and calling the gap an error.
    whole = dict(
        n=len(d),
        sc=int((d[["rmsd_af3", "rmsd_esm"]].min(axis=1) <= RMSD_GATE).sum()),
        gate=int(gate(d.rmsd_af3, d.rmsd_esm, d.dockq_af3, d.dockq_esm).sum()))
    assert whole["sc"] == 24_401, whole   # venn_regions_by_stack.csv, BG+BIF
    d = d[d[STORED_COLS].notna().all(axis=1)
          & d[TRUE_COLS].notna().all(axis=1)].copy()
    d["ordered"] = d.design_id.isin(ordered)
    # 1,480,140 route designs; 1,436,724 carry a true-backbone score and
    # 1,435,617 of those also carry a complete stored quadruple. All 2,397
    # rescored ordered designs survive, so panel (c)'s population CONTAINS
    # panel (a)'s exactly rather than merely overlapping it.
    assert len(d) == 1_435_617 and int(d.ordered.sum()) == 2397, (
        len(d), int(d.ordered.sum()))
    d.attrs["whole"] = whole
    return d


def gate(rmsd_a, rmsd_e, dockq_a, dockq_e):
    """The applied gate's two structural terms, either backend passes either term."""
    return (((rmsd_a <= RMSD_GATE) | (rmsd_e <= RMSD_GATE))
            & ((dockq_a >= DOCKQ_GATE) | (dockq_e >= DOCKQ_GATE)))


def pooled(d, af3, esm):
    """The AF3 + ESMFold2 POOLED convention: one value per design per backend."""
    return np.concatenate([d[af3].to_numpy(), d[esm].to_numpy()])


# --- (a) paired stored -> corrected, pooled over both backends --------------
def panel_a(d, out):
    fig, axes = plt.subplots(1, 2, figsize=(style.mm(85), style.mm(46)))
    # The median labels are nudged AWAY from the gate hairline, which both medians
    # start within a hair of: dy is (stored, corrected) in points.
    spec = [
        (axes[0], ("sc_rmsd_af3", "sc_rmsd_esm"), ("rmsd_af3_true", "rmsd_esm_true"),
         "ligand_RMSD (Å) ↓", RMSD_GATE, (0, 30), (-5, 5)),
        (axes[1], ("dockq_af3", "dockq_esm"), ("dockq_af3_true", "dockq_esm_true"),
         "sc_DockQ", DOCKQ_GATE, (0, 1), (5, -5)),
    ]
    for ax, s, t, label, thr, ylim, dy in spec:
        a, b = pooled(d, *s), pooled(d, *t)
        seg = np.stack([np.column_stack([np.zeros_like(a), a]),
                        np.column_stack([np.ones_like(b), b])], axis=1)
        # 4,794 paired segments: the pairing IS the point, so a line per pair,
        # rasterised so the PDF stays a few hundred kB rather than a few MB.
        ax.add_collection(LineCollection(seg, colors=palette.MAIN, linewidths=0.2,
                                         alpha=0.03, rasterized=True))
        ax.axhline(thr, color=palette.NEUTRAL["500"], lw=0.5, zorder=0)
        m = (np.median(a), np.median(b))
        ax.plot([0, 1], m, color=palette.INK, lw=1.0, marker="o", ms=2.5, zorder=3)
        fmt = "{:.2f}" if thr == DOCKQ_GATE else "{:.1f}"
        for x, v, ha, off in ((0, m[0], "right", -4), (1, m[1], "left", 4)):
            ax.annotate(fmt.format(v), (x, v), textcoords="offset points",
                        xytext=(off, dy[x]), ha=ha, va="center",
                        fontsize=style.CFG["font_size_pt"], zorder=4)
        ax.set_xlim(-0.35, 1.35)
        ax.set_ylim(*ylim)
        ax.set_xticks([0, 1])
        ax.set_xticklabels(["stored", "corrected"])
        ax.set_ylabel(label)
        ax.tick_params(axis="y", length=2, width=0.5)
        ax.tick_params(axis="x", length=0)
        ax.spines["bottom"].set_visible(False)
    fig.subplots_adjust(left=0.11, right=0.99, top=0.97, bottom=0.12, wspace=0.42)
    style.save(fig, out, "panel_a_paired_rescore")
    plt.close(fig)


# --- (c) stored vs corrected over the WHOLE ROUTE, per gate term -------------
# Leo, 2026-11-19: "B i want to see all quadrants. Just label the quadrants
# accordingly." On the ordered panel two of the four CANNOT hold a point. That
# panel is the gated set, so stored best-of-two ligand_RMSD topped out at exactly
# 3.0000 and stored sc_DockQ bottomed out at exactly 0.7000 - the gate wall was
# the edge of the data, and no labelling fixes an empty quadrant. Widening the
# population to the whole BoltzGen + BoltzIF route is what puts designs on the
# stored-FAIL side, so all four quadrants carry data.
#
# Leo, 2026-11-19: "remove all the dots. Just make it a clean hex plot with the
# 4 quadrants, and then the labels are enough." The INK overlay that marked the
# 2,397 ordered designs is gone, and with it every per-quadrant "ordered" count
# line - see the label block below. The counts are still derived and still
# printed by report_route(); they are a caption fact now, not a panel one.
#
# Each term on the gate's OWN convention - best-of-two across the two backends,
# min for ligand_RMSD and max for sc_DockQ - so one point is one design and the gate
# hairline is the line that design was actually judged against. The pooled
# convention of panel (a) would put stored values on the fail side of a line
# those designs passed, which reads as a rejection that never happened;
# best-of-two is the only convention a drawn gate can carry.
#
# aspect, limits, ticks: equal limits on both axes make y = x a true 45 deg
# reference and put the crossing of the two gate hairlines ON it. The ligand_RMSD
# ceiling is 70, not the ordered panel's 65: one route design sits at 66.1.
# sc_DockQ runs to 1.0 rather than the ordered panel's 0.95 - its natural ceiling,
# and 0.95 left the stored-PASS strip 26% of the axis where a 21.4 mm count
# line needs 30%. Measured; _fits() re-checks it on every run.
TERMS = [
    ("ligand_RMSD", "min", ("rmsd_af3", "rmsd_esm"),
     ("rmsd_af3_true", "rmsd_esm_true"), RMSD_GATE,
     "Stored ligand_RMSD (Å) ↓", "Corrected ligand_RMSD (Å) ↓",
     "log", (0.5, 70.0), [1, 3, 10, 30]),
    ("sc_DockQ", "max", ("dockq_af3", "dockq_esm"),
     ("dockq_af3_true", "dockq_esm_true"), DOCKQ_GATE,
     "Stored sc_DockQ", "Corrected sc_DockQ",
     "linear", (0.0, 1.0), [0, 0.25, 0.5, 0.75, 1.0]),
]

# Leo's four quadrants, as (stored passes this term, corrected passes this term,
# label, emphasised). "Falsely rejected" is his 2026-11-19 wording for the
# quadrant he first called "should have been ordered", and the rename keeps that
# label's logic: it is a counterfactual about designs that were NOT ordered, so
# it stays TRUE of every design in the quadrant however big the population
# grows. His other two originally ended "and ordered", which does NOT survive
# the widening: 7,840 of the 9,803 designs that the stored gate falsely accepted
# were never ordered, so the clause would be false of most of the quadrant. With
# the overlay gone there is nothing on the panel for an "ordered" count to refer
# to, so none of the four carries one - the caption does.
#
# The emphasis flag is on the quadrant that is the figure's result: designs the
# campaign discarded that the corrected scores accept. Weight, not colour.
QUADS = [
    (True, False, "Falsely accepted", False),
    (True, True, "Both accepted", False),
    (False, True, "Falsely rejected", True),
    (False, False, "Both rejected", False),
]

# (c) is the house DOUBLE-column width: eight labelled quadrants do not fit in
# 85 mm, MEASURED rather than guessed. The widest count line is 17.4 mm and
# sc_DockQ's stored-pass strip is 26% of its x axis - 13.7 mm at 85 mm against
# 20.6 mm here. _fits() re-measures every block on every run. Equal aspect
# makes the axes square, so the figure HEIGHT follows from the width instead
# of being chosen.
B_ADJ = dict(left=0.055, right=0.995, top=0.985, bottom=0.085, wspace=0.15)
B_AX_MM = (style.CFG["width_double_mm"] * (B_ADJ["right"] - B_ADJ["left"])
           / (2 + B_ADJ["wspace"]))
B_H_MM = B_AX_MM / (B_ADJ["top"] - B_ADJ["bottom"])
LBL_PAD = 0.012


def best_of_two(d, cols, how):
    """The gate's own quantity: min ligand_RMSD / max sc_DockQ over AF3 and ESMFold2."""
    return getattr(d[list(cols)], how)(axis=1).to_numpy()


def _frac(thr, lim, scale):
    """Where a gate threshold sits across the axes, 0 to 1."""
    lo, hi = lim
    if scale == "log":
        return (np.log(thr) - np.log(lo)) / (np.log(hi) - np.log(lo))
    return (thr - lo) / (hi - lo)


# The book's page geometry, from make_supp_book.py (A4 210 mm - 2 x 14 mm
# margin, 6 mm between panels in a row). A panel pair only shares a row while
# it fits in the text block, and the book reflows SILENTLY when it stops -
# which is how a late edit turns a side-by-side pair into two rows and gets
# blamed on something else. ROW_SLACK_MM is how much room we want to have.
BOOK_TEXT_MM, BOOK_GUTTER_MM, ROW_SLACK_MM = 182.0, 6.0, 1.5


def emitted_box(out, stem):
    """Page size of the EMITTED PDF in mm, asserting no text falls off it.

    ⚠️ The PDF, not the PNG, and not matplotlib's own idea of the text extent.
    On 2026-11-20 the x label of panel (c) was clipped 0.78 pt off the bottom of
    the MediaBox - losing the underscore of "ligand_RMSD" - and the 1.00 tick 2.81 pt
    off the right. Both rendered perfectly in the PNG, and matplotlib's renderer
    bbox put the label 0.5 pt INSIDE the canvas, because it measures the laid-out
    box and not the glyphs' own descender metrics. Only the artefact of record
    shows it, so that is what is measured here, on every run.
    """
    with pymupdf.open(out / f"{stem}.pdf") as doc:
        page = doc[0]
        R = page.rect
        for blk in page.get_text("rawdict")["blocks"]:
            for line in blk.get("lines", []):
                for sp in line["spans"]:
                    b = sp["bbox"]
                    over = max(b[3] - R.height, -b[1], b[2] - R.width, -b[0])
                    assert over <= 0, (
                        f"{stem}.pdf: {''.join(c['c'] for c in sp['chars'])!r} "
                        f"runs {over:.2f} pt off a {R.width:.1f} x {R.height:.1f} pt "
                        "page and will be clipped in print")
        return R.width * 25.4 / 72, R.height * 25.4 / 72


def _fits(text, width_mm, bold=False):
    """Nothing clipped (STYLE_NOTES checklist 3), measured rather than eyeballed.

    style.text_w returns ink width in the same unit as the size it is given, so
    points here; a quadrant's block has to clear its own side of the hairline.
    The 0.52 em/char stand-in is out by 10-33% per string, which is how a
    too-wide label reaches a PDF looking merely tight. Measured in the weight it
    is actually drawn in, bold being the wider.
    """
    w = max(style.text_w(ln, weight="bold" if bold else "normal")
            for ln in text.split("\n")) * 25.4 / 72
    assert w <= width_mm - 1.0, (
        f"{w:.1f} mm of label in a {width_mm:.1f} mm quadrant: {text!r}")
    return w


def panel_c(d, out):
    """Stored against corrected over the whole route, one point per design.

    All four quadrants carry data because the population is the ROUTE and not
    the gated set. Nothing is overlaid: the ramp and the four labelled quadrants
    are the whole panel. Every ordered design is a stored PASS by construction,
    so the two stored-FAIL quadrants hold none of them - "Falsely rejected" is
    the campaign's own discard pile, and that none of it was bought is a caption
    fact (report_route() still counts it).

    The counts are PER TERM and the two terms overlap, so no single one of them
    is the combined gate flip; report() prints the gate-level 2x2 beside them.
    """
    fig, axes = plt.subplots(1, 2, figsize=(style.mm(style.CFG["width_double_mm"]),
                                            style.mm(B_H_MM)))
    o = d.ordered.to_numpy()
    n_tot = len(d)
    counts = {}
    for ax, (name, how, s, t_, thr, xlab, ylab, scale, lim, ticks) in zip(axes, TERMS):
        a, b = best_of_two(d, s, how), best_of_two(d, t_, how)
        sp = (a <= thr) if how == "min" else (a >= thr)
        cp = (b <= thr) if how == "min" else (b >= thr)
        ax.set_xscale(scale), ax.set_yscale(scale)
        # 1.4M points is a density, not a scatter. House ramp and house count
        # scale (style.density_cmap / density_norm), shared with the main
        # text's 2-D densities so the two cannot drift apart. It also keeps the
        # sparse quadrants readable, which an alpha-stacked scatter does not
        # manage at this n - the quadrants that matter here hold ~1% each.
        ext = ([np.log10(lim[0]), np.log10(lim[1])] * 2 if scale == "log"
               else [lim[0], lim[1]] * 2)
        hb = ax.hexbin(a, b, gridsize=78, extent=ext, mincnt=1, linewidths=0,
                       xscale=scale, yscale=scale, cmap=style.density_cmap(),
                       rasterized=True, zorder=1)
        style.density_norm(hb)
        ax.plot(lim, lim, color=palette.NEUTRAL["400"], lw=0.5, ls=(0, (3, 2)),
                zorder=2)                                    # y = x
        ax.axvline(thr, color=palette.NEUTRAL["500"], lw=0.5, zorder=2)
        ax.axhline(thr, color=palette.NEUTRAL["500"], lw=0.5, zorder=2)
        f = _frac(thr, lim, scale)
        lo_w, hi_w = f * B_AX_MM, (1 - f) * B_AX_MM
        for s_pass, c_pass, lab, emph in QUADS:
            m = (sp if s_pass else ~sp) & (cp if c_pass else ~cp)
            k, k_o = int(m.sum()), int((m & o).sum())
            counts[(name, lab)] = (k, k_o)
            # Each block goes in its quadrant's OUTER corner, furthest from the
            # gate crossing, which is also where that quadrant is emptiest.
            low_x = s_pass if how == "min" else not s_pass
            low_y = c_pass if how == "min" else not c_pass
            x, ha = (LBL_PAD, "left") if low_x else (1 - LBL_PAD, "right")
            y, va = (LBL_PAD, "bottom") if low_y else (1 - LBL_PAD, "top")
            # Two lines, not three: the "N ordered" line went with the dots it
            # referred to. k_o is still counted for report_route().
            txt = f"{lab}\n{k:,}/{n_tot:,} ({100 * k / n_tot:.1f}%)"
            _fits(txt, lo_w if low_x else hi_w, emph)
            # A quadrant's own corner is not always empty - sc_DockQ's "both
            # rejected" corner is the densest part of the whole figure - so each
            # block sits on a soft SURFACE patch and stays readable over any bin
            # of the ramp. Legibility, not decoration: without it two of the
            # eight blocks are unreadable on near-black hexes.
            #
            # A patch and NOT a path-effects halo, which was tried first: a
            # stroked glyph has to be rendered as an outline, and that silently
            # turned these eight blocks into paths in the SVG. style.save()
            # writes SVG precisely so assembly gets live text it can restyle,
            # so the one element that must stay text cannot be the one we
            # outline. Verified: the SVG still carries every line as <text>.
            ax.annotate(txt, (x, y), xycoords="axes fraction", ha=ha, va=va,
                        fontsize=style.CFG["font_size_pt"], linespacing=1.35,
                        fontweight="bold" if emph else "normal", zorder=4,
                        bbox=dict(boxstyle="square,pad=0.22", linewidth=0,
                                  facecolor=palette.SURFACE, alpha=0.78))
        assert sum(k for (nm, _), (k, _) in counts.items() if nm == name) == n_tot

        ax.set_xlim(*lim), ax.set_ylim(*lim)
        ax.set_aspect("equal", adjustable="box")
        ax.set_xticks(ticks), ax.set_yticks(ticks)
        if scale == "log":
            ax.set_xticklabels([str(v) for v in ticks])
            ax.set_yticklabels([str(v) for v in ticks])
            ax.minorticks_off()
        else:
            # "{:g}", not the default "0.00 / 0.25 / ... / 1.00". The 1.00 tick
            # is CENTRED on the axes' right edge, 10.7 pt of type on a canvas
            # with 2.5 pt of margin there, so matplotlib drew it 2.81 pt over
            # the MediaBox and the PDF lost the last digit. "1" is 3.1 pt and
            # sits inside. It also makes this half agree with the log half,
            # which has always written its ticks bare (1 / 3 / 10 / 30), with
            # `repredictor_comparison` panel (a), and with STYLE_NOTES 7's
            # "round numbers elsewhere". No geometry moves, which is the point:
            # this panel's height is DERIVED from its width under equal aspect,
            # so any margin lever would have moved the documented box.
            ax.set_xticklabels([f"{v:g}" for v in ticks])
            ax.set_yticklabels([f"{v:g}" for v in ticks])
        # labelpad 2.5, not the 4.0 default. ⚠️ At 4.0 the x label's baseline sat
        # low enough that the UNDERSCORE of "Stored ligand_RMSD (Å) ↓" - which hangs
        # ~0.1 em below it - fell 0.78 pt off the bottom of the MediaBox and was
        # clipped out of the PDF, so the printed panel read "Stored ligand RMSD (Å) ↓".
        # The PNG was always correct, which is exactly why it survived a rename
        # sweep: measure the PDF, not the PNG. Set on BOTH halves so the two x
        # labels stay on one line; emitted_box() now asserts the whole thing.
        ax.set_xlabel(xlab, labelpad=2.5), ax.set_ylabel(ylab)
        ax.tick_params(axis="both", length=2, width=0.5)
    fig.subplots_adjust(**B_ADJ)
    # tight=False: the canvas is laid out in mm and the height was DERIVED from
    # the width, so a tight box would quietly resize what assembly places.
    # panel_b_stored_vs_corrected until 2026-11-20; the stem carries the letter.
    style.save(fig, out, "panel_c_stored_vs_corrected", tight=False)
    plt.close(fig)
    return counts


def by_round(d):
    t = (d.groupby("campaign")["flip"].agg(flip="sum", total="size")
         .reindex(ROUNDS))
    assert t.total.notna().all(), f"a round went missing: {t}"
    t["survive"] = t.total - t.flip
    return t.astype(int)


def chisq(obs):
    """Pearson chi-square of independence, and an exact p for even df.

    Hand-rolled: scipy in this venv raises on import (_propack dlopen). For even
    df the upper tail is closed-form, P(X >= x) = exp(-x/2) * sum_{i<df/2} (x/2)^i / i!
    """
    rt, ct = obs.sum(1, keepdims=True), obs.sum(0, keepdims=True)
    exp_ = rt @ ct / obs.sum()
    chi2 = float(((obs - exp_) ** 2 / exp_).sum())
    df = (obs.shape[0] - 1) * (obs.shape[1] - 1)
    assert df % 2 == 0, "the closed form below needs an even df"
    h = chi2 / 2
    p = exp(-h) * sum(h ** i / factorial(i) for i in range(df // 2))
    return chi2, df, p, float(exp_.min())


# --- (b) does high ipTM protect a design? -----------------------------------
def panel_b(d, out):
    """Flip rate by ipTM decile. ipTM never reads the designed backbone, so the
    curve is close to flat - which is the claim. `iptm_max` is the gate's own
    quantity; the per-backend conventions are in the README."""
    q = np.quantile(d.iptm_max, np.linspace(0, 1, N_DEC + 1))
    q[-1] = np.inf
    dec = np.searchsorted(q[1:-1], d.iptm_max, side="right")
    g = pd.DataFrame({"dec": dec, "flip": d.flip.to_numpy(),
                      "iptm": d.iptm_max.to_numpy()}).groupby("dec")
    t = g.agg(rate=("flip", "mean"), n=("flip", "size"), lo=("iptm", "min"))
    overall = d.flip.mean()

    fig, ax = plt.subplots(figsize=(style.mm(85), style.mm(38)))
    ax.axhline(100 * overall, color=palette.NEUTRAL["600"], lw=0.5, zorder=0)
    ax.plot(t.index + 1, 100 * t.rate, color=palette.MAIN, lw=0.8, marker="o", ms=3)
    # Reference line named at the left, where the curve is furthest above it; the
    # top decile's raw counts below its own point, where nothing else sits.
    ax.annotate(f"all designs {100 * overall:.1f}%", (1, 100 * overall),
                textcoords="offset points", xytext=(0, -7), ha="left", va="top",
                fontsize=style.CFG["font_size_pt"], color=palette.NEUTRAL["600"])
    top, bot = t.iloc[-1], t.iloc[0]
    for n_flip, n, xi, dy in ((round(bot.rate * bot.n), bot.n, 1, 5),
                              (round(top.rate * top.n), top.n, N_DEC, -9)):
        ax.annotate(f"{int(n_flip)}/{int(n)}", (xi, 100 * (n_flip / n)),
                    textcoords="offset points", xytext=(0, dy), ha="center",
                    fontsize=style.CFG["font_size_pt"])
    ax.set_xlim(0.4, N_DEC + 0.6)
    ax.set_xticks(range(1, N_DEC + 1))
    ax.set_xlabel("ipTM decile")
    ax.set_ylabel("Flip to fail (%)")
    ax.set_ylim(0, 100)
    ax.set_yticks([0, 25, 50, 75, 100])
    ax.tick_params(axis="both", length=2, width=0.5)
    fig.subplots_adjust(left=0.11, right=0.99, top=0.97, bottom=0.20)
    # panel_c_flip_by_iptm_decile until 2026-11-20; it moved up beside (a).
    style.save(fig, out, "panel_b_flip_by_iptm_decile")
    plt.close(fig)
    return t


# --- report ------------------------------------------------------------------
def report(d, chi, dec):
    chi2, df, p, min_exp = chi
    n = len(d)
    print(f"rescored {n}/2400; true_pass {d.true_pass.sum()}; "
          f"flip {d.flip.sum()}/{n} ({100 * d.flip.mean():.2f}%)")
    # The ordered panel's own per-term numbers and the overlap that keeps
    # them from being mistaken for the 81.9% combined flip. The panel itself
    # now draws the ROUTE; report_route() carries its quadrants.
    rc = d[["rmsd_af3_true", "rmsd_esm_true"]].min(axis=1) <= RMSD_GATE
    qc = d[["dockq_af3_true", "dockq_esm_true"]].max(axis=1) >= DOCKQ_GATE
    print(f"  ordered panel per term: fails corrected ligand_RMSD "
          f"{int((~rc).sum())}/{n} ({100 * (~rc).mean():.1f}%), sc_DockQ "
          f"{int((~qc).sum())}/{n} ({100 * (~qc).mean():.1f}%)")
    print(f"  ordered panel terms overlap: both fail {int((~rc & ~qc).sum())}, "
          f"ligand_RMSD only {int((~rc & qc).sum())}, sc_DockQ only {int((~qc & rc).sum())}, "
          f"union {int((~rc | ~qc).sum())} = the combined flip")
    print("  denominator variants, all 81.9%: "
          f"{d.flip.sum()}/{n} = {100 * d.flip.sum() / n:.2f}%, "
          f"{d.flip.sum() + 2}/{n + 2} = {100 * (d.flip.sum() + 2) / (n + 2):.2f}% "
          "(the two no_native_interface rows counted as flips), "
          f"{d.flip.sum() + 2}/2400 = {100 * (d.flip.sum() + 2) / 2400:.2f}%")
    for name, s, t in (("ligand_RMSD", ("sc_rmsd_af3", "sc_rmsd_esm"),
                        ("rmsd_af3_true", "rmsd_esm_true")),
                       ("sc_DockQ", ("dockq_af3", "dockq_esm"),
                        ("dockq_af3_true", "dockq_esm_true"))):
        a, b = pooled(d, *s), pooled(d, *t)
        lo = d[list(s)].min(axis=1) if name == "ligand_RMSD" else d[list(s)].max(axis=1)
        hi = d[list(t)].min(axis=1) if name == "ligand_RMSD" else d[list(t)].max(axis=1)
        print(f"  {name} median POOLED {np.median(a):.3f} -> {np.median(b):.3f} | "
              f"best-of-two {lo.median():.3f} -> {hi.median():.3f}")
    print(by_round(d)[["survive", "flip", "total"]].to_string())
    print(f"  chi2 {chi2:.4f} df {df} p {p:.4f} smallest expected {min_exp:.2f}")
    sp = pooled(d[~d.flip], "iptm_af3", "iptm_esm")
    fp = pooled(d[d.flip], "iptm_af3", "iptm_esm")
    print(f"  ipTM survive-flip mean gap: POOLED {sp.mean() - fp.mean():+.4f}", end="")
    for c in ("iptm_af3", "iptm_esm", "iptm_max"):
        gap = d.loc[~d.flip, c].mean() - d.loc[d.flip, c].mean()
        thr = np.quantile(d[c], 0.9)
        top = d[d[c] >= thr]
        print(f" | {c} {gap:+.4f}, top decile {100 * top.flip.mean():.2f}% "
              f"({top.flip.sum()}/{len(top)}, >= {thr:.4f})", end="")
    print()
    print("  flip rate by ipTM decile (%):",
          ", ".join(f"{100 * r:.1f}" for r in dec.rate))


def report_route(r, quad):
    """Panel (c): its population, its four quadrants, and the GATE-level 2x2.

    The per-term quadrants are what the panel draws. The gate-level 2x2 is what
    the caption needs, and the two are NOT the same event - the terms overlap.
    The two CHECK lines re-derive numbers that other files in this repo settled
    independently, so a changed input shows up here rather than in the caption.
    """
    n = len(r)
    print(f"\n(c) WHOLE ROUTE n={n:,} of 1,480,140 ({100 * n / 1_480_140:.1f}%); "
          f"ordered {int(r.ordered.sum()):,}")
    for spec in TERMS:
        print(f"  {spec[0]} best-of-two:")
        for _, _, lab, _ in QUADS:
            k, k_o = quad[(spec[0], lab)]
            print(f"      {lab:22s} {k:>9,}/{n:,} ({100 * k / n:6.3f}%)  "
                  f"ordered {k_o:,}")
    sg = gate(r.rmsd_af3, r.rmsd_esm, r.dockq_af3, r.dockq_esm)
    cg = gate(r.rmsd_af3_true, r.rmsd_esm_true, r.dockq_af3_true, r.dockq_esm_true)
    print("  GATE LEVEL, both structural terms (NOT what the panel draws):")
    for lab, m in (("Both rejected", ~sg & ~cg), ("Falsely accepted", sg & ~cg),
                   ("Both accepted", sg & cg), ("Falsely rejected", ~sg & cg)):
        print(f"      {lab:22s} {int(m.sum()):>9,}/{n:,} "
              f"({100 * m.sum() / n:6.3f}%)  ordered {int((m & r.ordered).sum()):,}")
    print(f"  ROUTE-WIDE flip rate {int((sg & ~cg).sum()):,}/{int(sg.sum()):,} "
          f"({100 * (sg & ~cg).sum() / sg.sum():.1f}%) - the ordered panel's "
          "81.9% is NOT this number and must not be quoted as it")
    print(f"  stored gate {int(sg.sum()):,} -> corrected gate {int(cg.sum()):,}")
    # Independent cross-checks, each against the population it was counted on.
    # load_route() already asserted the stored ligand_RMSD term on the FULL route
    # (24,401 = venn_regions_by_stack.csv); the corrected one lives here because
    # the corrected cascade counts only designs that HAVE a corrected score,
    # which is this population exactly.
    w = r.attrs["whole"]
    c_sc = int((r[["rmsd_af3_true", "rmsd_esm_true"]].min(axis=1) <= RMSD_GATE).sum())
    assert c_sc == 18_577, c_sc      # the README's corrected cascade, g1_sc_rmsd
    s_sc = int((r[["rmsd_af3", "rmsd_esm"]].min(axis=1) <= RMSD_GATE).sum())
    print(f"  CHECK stored ligand_RMSD pass {w['sc']:,} over all {w['n']:,} "
          f"= venn_regions_by_stack.csv; {s_sc:,} of them are drawn, the other "
          f"{w['sc'] - s_sc:,} carrying no complete pair of scores")
    print(f"  CHECK corrected ligand_RMSD pass {c_sc:,} = the corrected cascade's "
          f"g1_sc_rmsd; ratios {100 * w['gate'] / w['sc']:.1f}% stored and "
          f"{100 * int(cg.sum()) / c_sc:.1f}% corrected, as published")
    print(f"  CHECK stored both terms: {w['gate']:,} over the full route "
          f"against venn's 12,158, and {int(sg.sum()):,} drawn. It reconciles "
          f"exactly: {w['gate'] - 12_158:+,} from the campaign master's 4-dp "
          f"sc_DockQ rounding (designs landing on exactly 0.7000) and "
          f"{int(sg.sum()) - w['gate']:+,} with no true-backbone score to plot")


def check_emitted(out):
    """Every panel measured on its own emitted PDF, then the page layout checked.

    Two failures this catches, both of which reached a rendered page before:
      1. text clipped off a panel's MediaBox - emitted_box() asserts per panel;
      2. (a) and (b) no longer fitting one row of the book, which make_supp_book
         handles by SILENTLY reflowing (b) onto its own row, undoing Leo's
         "move c up next to a" without a word anywhere.
    The row assert is at the real limit, so the reflow cannot happen without the
    build stopping. It is NOT at BOOK_TEXT_MM - ROW_SLACK_MM, because the pair
    emits at 181.9 mm today and that would fail every run: (a) and (b) save with
    tight=True and their tight boxes run 2.6 and 3.3 mm past the authored 85 mm.
    Buying the slack means changing panel geometry nobody asked to change, so
    the thin margin is WARNED about on every run instead.
    """
    wa, ha = emitted_box(out, "panel_a_paired_rescore")
    wb, hb = emitted_box(out, "panel_b_flip_by_iptm_decile")
    wc, hc = emitted_box(out, "panel_c_stored_vs_corrected")
    row = wa + BOOK_GUTTER_MM + wb
    print(f"\nemitted: (a) {wa:.1f} x {ha:.1f}, (b) {wb:.1f} x {hb:.1f}, "
          f"(c) {wc:.1f} x {hc:.1f} mm; no text off any MediaBox")
    assert wc <= BOOK_TEXT_MM, (
        f"(c) is {wc:.1f} mm and the book's text block is {BOOK_TEXT_MM:.0f} mm: "
        "it will be scaled down to fit")
    assert row <= BOOK_TEXT_MM, (
        f"(a) + gutter + (b) = {row:.1f} mm against a {BOOK_TEXT_MM:.0f} mm text "
        "block: the book will drop (b) onto its own row and the side-by-side "
        "layout asked for on 2026-11-20 is gone. Narrow one of them, or emit it "
        "at its authored 85 mm with tight=False (which clips its annotations).")
    slack = BOOK_TEXT_MM - row
    if slack < ROW_SLACK_MM:
        print(f"  \u26a0\ufe0f  (a) + 6 mm + (b) = {row:.1f} mm leaves {slack:.1f} mm "
              f"in the {BOOK_TEXT_MM:.0f} mm text block, under the {ROW_SLACK_MM} mm "
              "we want. One longer tick label and the row breaks - see README.")
    else:
        print(f"  row (a)|(b) = {row:.1f} mm, {slack:.1f} mm of slack")


def main():
    style.apply()
    out = paths.fig_output(FIG)
    d = load()
    panel_a(d, out)
    dec = panel_b(d, out)
    # (c) is the whole route; (a) and (b) stay on the ordered panel.
    route = load_route()
    quad = panel_c(route, out)
    # The round-uniformity test LOST ITS PANEL when the route panel became the paired
    # scatter (Leo, 2026-11-19). It is still a result, so it is still derived
    # on every run and it is stated in the README's caption draft with its
    # statistic - it must not quietly vanish with the bars.
    chi = chisq(by_round(d)[["survive", "flip"]].to_numpy(float))
    check_emitted(out)
    report(d, chi, dec)
    report_route(route, quad)


if __name__ == "__main__":
    main()
