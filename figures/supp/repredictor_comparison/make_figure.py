"""Why Boltz-2 was dropped: two panels on the route-balanced 250,000-design sample.

Source: ~/rna_binder_plots/_data/refolder_agreement_designs.parquet (2026-11-03),
50,000 designs per design route, exactly balanced, all nine metric columns present
(3 re-predictors x ipTM, sc_DockQ, ligand_RMSD). Boltz-2 is in no other local table (the
filter store and campaign_metrics_full.csv.gz are AF3+ESM only), so the campaign-wide
version of this figure would cost another Euler pass. Because this sample is
route-BALANCED its rates are deliberately not campaign rates - report() prints the two
reweightings that bracket the campaign-wide value, for both gate quantities.

BG+BIF: the parquet's main dockq / sc_rmsd columns already hold the CORRECTED
true-backbone values (dockq_self_*_true, dockq_LRMSD_*_true from the 2026-10-21
rescore); the stored ones survive only as *_stored, unused here. So every panel is on
the corrected numbers, which is what Leo asked for ("the updated boltzgen+bif numbers
for all") and which matches Fig3a_backend_agreement's fix_bif(). The campaign's own
funnel was gated on the STORED values, so the chain report() prints is not the
campaign's funnel - see the README, which says so.

Panels
  (a) 3 x 3: rows ipTM / sc_DockQ / ligand_RMSD, columns the three re-predictor pairs. Density on
      a log count scale, y = x dotted, campaign thresholds dashed on both axes, and two
      numbers per cell: R2 (the share of variance one re-predictor's score explains in the
      other) and the share of designs on which the pair reaches the same verdict at that
      threshold. They disagree, and that disagreement is the figure's argument.
      Conventions, ramp and count scale are the main text's Fig3a_backend_agreement,
      whose three panels are this matrix's AF3-ESMFold2 column - see "reconciliation".
  (b) seed reproducibility per re-predictor: a box over the ICCs of 16 interface metrics
      (the `metric_reliability` figure's panel a, collapsed to one box per re-predictor).
      Different n and different designs from (a) - see README.

THE GATE PANEL IS GONE, 2026-11-20 (Leo: "Remove panel b. Try to get c next to a and
rename it"). It drew the gate as applied - g1 ligand_RMSD, g2 sc_DockQ, g3 seq run, g5 ipTM
- step by step over all 250,000 designs, with the count still standing after each, for
the campaign's either-backend rule and for each re-predictor alone. Its drawing code is
deleted; its ARITHMETIC is not. chain(), campaign_chain(), GATE and report()'s gate
block all stay, so every number the panel carried is still derived and printed on every
run and cannot drift out of the README - the same treatment bif_effect gave its
round-uniformity chi-square when that lost its panel. The last render is in
archive/repredictor_comparison_panel_b_gate_cascade_2026-11-20/, which also records
which of its two readings survive in other figures and which do not.

r vs R2. Fig3a_backend_agreement was switched to R2 on 2026-11-19 at Leo's request;
this panel follows it, so the book annotates one statistic and not two. Pearson r and
Spearman rho are still printed by report(). The superscript is a literal U+00B2 and
never mathtext: $R^2$ renders in DejaVu and would be the only non-Arial glyph in the
book (2026-11-16 font audit).

Reconciliation with Fig3a_backend_agreement. The middle column here draws the same
three measurements as that figure, and the two cannot be allowed to disagree:
  - ramp and count scale come from style.density_cmap / style.density_norm, which both
    scripts call, so neither can drift;
  - thresholds, axis ranges, tick sets, label text, the y = x rule and the statistic are
    copied from it;
  - the numbers differ only in their population, and by very little: campaign-wide
    (6.99 M designs, that figure) R2 = 0.063 / 0.245 / 0.323 against 0.065 / 0.248 /
    0.326 on this route-balanced 250,000 (ipTM / sc_DockQ / ligand_RMSD). report() prints both
    and the README tabulates them.

Re-predictor colours: palette.REFOLDERS, pipeline order palette.REFOLDER_ORDER, never
rank (COLORS.md rule 7). The hatch encoding this figure used until 2026-11-19 is
withdrawn - Leo asked twice for the dots and the stripes to go - and no mark here
carries a hatch.

scipy is broken in ~/.venvs/binder_map (_propack dlopen failure), so every statistic
here is hand-rolled: Pearson is numpy's corrcoef, Spearman a Pearson on pandas
average-tied ranks, and the box quartiles are numpy percentiles. No p-values: at
n = 250,000 every p underflows and carries no information.
"""
import itertools
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from common import paths, style, palette

FIG = Path(__file__).resolve().parent.name
SRC = Path.home() / "rna_binder_plots" / "_data" / "refolder_agreement_designs.parquet"
SEED_SRC = (Path.home() / "rna_binder_plots" / "seed_metric_reliability"
            / "seed_metric_reliability.csv")

N_SAMPLE = 250_000
PER_ROUTE = 50_000
K_SEEDS = 5
N_SEED_DESIGNS = 2_000    # designs per re-predictor in the seed experiment

# Pipeline order and colours come from the palette, never from a local list
# (COLORS.md rule 7: a re-predictor is an instrument, never a design route).
PRED = palette.REFOLDER_ORDER                      # AF3, Boltz-2, ESMFold2
PFX = {"AF3": "af3", "Boltz-2": "b2", "ESMFold2": "esm"}
PAIRS = list(itertools.combinations(PRED, 2))      # pipeline order, not by rank

# The campaign gate's three metric terms. Thresholds, labels, axis ranges and tick
# sets are Fig3a_backend_agreement's, unchanged. (key, label, lo, hi, thr, ge, log,
# ticks); ge=False means the metric passes BELOW its threshold.
M = [("iptm", "ipTM", 0.0, 1.0, 0.6, True, False, [0, .5, 1]),
     ("dockq", "sc_DockQ", 0.0, 1.0, 0.7, True, False, [0, .5, 1]),
     ("sc_rmsd", "ligand_RMSD (Å) ↓", 0.5, 100.0, 3.0, False, True, [1, 10, 100])]

# THE GATE AS APPLIED, in its own order, from the export's README section 4 ("The
# gate, exactly as applied") and validated there against the documented campaign
# chain 6,992,515 -> ... -> 55,706 on all five steps.
#
#   g1  ligand_RMSD  <= 3 A     EITHER backend (AF3 or ESMFold2)
#   g2  sc_DockQ >= 0.7     EITHER backend
#   g3  longest single-residue run <= 5     sequence property, NO backend
#   g4  footnb_010 <= 5 heavy atoms         BOTH backends  <- the one AND
#   g5  max(AF3, ESM) ipTM >= 0.6           EITHER backend
#
# \u26a0\ufe0f This was panel (b)'s x axis until 2026-11-20. The panel is gone; the chain is
# NOT - it is still derived and printed by report() on every run, which is the only
# place the quantifier asymmetry (either / either / no backend / either, with the
# missing g4 the one AND term) and the per-re-predictor counterfactuals now live.
# The tick labels kept a newline while they were two-line ticks; that is gone with
# the ticks, so they are plain one-line strings again.
# g4 is NOT a column in this export (footnb was joined only for
# cascade_by_stack.csv), which is why this chain ends at 2,894 and not at the
# campaign's own figure.
#
# (key, label, quantifier, threshold, ge, campaign flag). key "seq" has no
# per-predictor form: it reads the sequence, so all three lines take the same cut.
GATE = [("sc_rmsd", "ligand_RMSD \u2264 3 \u00c5", "either", 3.0, False, "g1"),
        ("dockq", "+ sc_DockQ \u2265 0.7", "either", 0.7, True, "g2"),
        ("seq", "+ seq run \u2264 5", "no backend", None, None, "g3"),
        ("iptm", "+ ipTM \u2265 0.6", "either", 0.6, True, "g5")]
# Thresholds must agree with M, which panel (a) draws, or the panel and the printed
# chain would describe different gates.
_MT = {k: (thr, ge) for k, _, _, _, thr, ge, _, _ in M}
for _k, _lab, _q, _thr, _ge, _f in GATE:
    assert _k == "seq" or _MT[_k] == (_thr, _ge), f"gate term {_k} disagrees with M"

# Designs generated and Boltz-2-eligible per route, from the export's README.txt
# sections 4-5. Used only for the reweightings printed by report().
GENERATED = {"RFD3+LiMPNN": 1_524_796, "BG+LiMPNN": 1_519_958, "BG+BIF": 1_480_140,
             "BD1+LiMPNN": 1_444_614, "HD+LiMPNN": 1_023_007}
ELIGIBLE = {"RFD3+LiMPNN": 1_524_724, "BG+LiMPNN": 1_476_307, "BG+BIF": 1_436_724,
            "BD1+LiMPNN": 1_444_520, "HD+LiMPNN": 846_924}
# Route order and short labels come from the palette, never from a local list.
ROUTES = list(palette.STACKS_SHORT)


def pearson(x, y):
    return float(np.corrcoef(x, y)[0, 1])


def spearman(x, y):
    return float(pd.Series(x).rank().corr(pd.Series(y).rank()))


def term(d, p, key, thr, ge):
    """One gate term for one re-predictor, as a boolean mask."""
    v = d[f"{PFX[p]}_{key}"]
    return v >= thr if ge else v <= thr


def chain(d, p):
    """Masks after each APPLIED gate step in turn, cumulative, for one re-predictor.

    This was panel (b) until 2026-11-20; it is now report() output only. The
    campaign never gated on a single re-predictor, so every one of these chains is
    a COUNTERFACTUAL - what that re-predictor alone would have cut - except the
    seq-run step, which reads the sequence and is therefore the same cut for all
    three. Boltz-2 never faced the gate at all.
    """
    m = np.ones(len(d), bool)
    out = []
    for key, _, _, thr, ge, flag in GATE:
        m = m & (d[flag].to_numpy() if key == "seq"
                 else term(d, p, key, thr, ge).to_numpy())
        out.append(m)
    return out


def campaign_chain(d):
    """The same steps for the gate AS APPLIED: each metric term on EITHER gated
    backend, which is what g1 / g2 / g5 are, plus the backend-free g3. g4 is not a
    column in this export, so this chain stops one term short of the campaign's
    own - see the README."""
    m = np.ones(len(d), bool)
    out = []
    for _, _, _, _, _, flag in GATE:
        m = m & d[flag].to_numpy()
        out.append(m)
    return out


def metric_gate(d, p):
    """The three METRIC terms at once for one re-predictor, without seq run. Kept
    because report()'s route reweightings and the README's per-term table are both
    on it."""
    m = np.ones(len(d), bool)
    for key, _, _, _, thr, ge, _, _ in M:
        m = m & term(d, p, key, thr, ge).to_numpy()
    return m


def fits(s, mm, where):
    assert style.text_w(s) / 72 * 25.4 < mm, f"{where}: {s!r} overruns {mm} mm"


def main():
    style.apply()
    out = paths.fig_output(FIG)

    cols = ["stack"] + [f"{q}_{k}" for q in PFX.values() for k, *_ in M]
    d = pd.read_parquet(SRC, columns=cols + ["g1", "g2", "g3", "g5", "passes_all"])
    assert len(d) == N_SAMPLE, f"expected {N_SAMPLE} rows, got {len(d)}"
    counts = d["stack"].value_counts()
    assert set(counts.index) == set(ROUTES) and (counts == PER_ROUTE).all(), \
        f"sample is not route-balanced: {counts.to_dict()}"
    assert d[cols[1:]].notna().all().all(), "sample should carry no NaN metric"
    # The gate flags shipped with the export must BE the three terms applied to the two
    # re-predictors the campaign gated on, or the printed campaign chain means nothing.
    # Recomputed here from the corrected columns: they agree on all 250,000 rows, which
    # also confirms the flags are on the corrected BG+BIF values and not the stored ones.
    for flag, got in (("g1", term(d, "AF3", "sc_rmsd", 3.0, False)
                             | term(d, "ESMFold2", "sc_rmsd", 3.0, False)),
                      ("g2", term(d, "AF3", "dockq", 0.7, True)
                             | term(d, "ESMFold2", "dockq", 0.7, True)),
                      ("g5", (np.maximum(d["af3_iptm"], d["esm_iptm"]) >= 0.6))):
        assert (d[flag] == got).all(), f"{flag} does not reproduce from the columns"
    # g3 (seq run <= 5) cannot be recomputed here - the export carries the flag but
    # not the sequence - so it is checked the one way the file allows: the export's
    # own passes_all must BE the four available terms, which is where the printed
    # campaign chain ends. g4 (foot) is in neither, hence 2,894 and not 56,179.
    assert (d["passes_all"] == (d.g1 & d.g2 & d.g3 & d.g5)).all(), \
        "passes_all is not g1 & g2 & g3 & g5"
    seed = pd.read_csv(SEED_SRC)

    panel_a(d, out)
    panel_b(seed, out)
    report(d, seed)


# --- (a) three metrics x three re-predictor pairs ------------------------------
def panel_a(d, out):
    """Rows are metrics, columns re-predictor pairs. Both axes of a cell carry the row's
    metric, so the metric is written ONCE per row in the left margin and each cell's
    own axis labels are just the two re-predictors - six short labels a row instead of six
    long ones (STYLE_NOTES sections 0 and 3). The density ramp and the count scale are
    style.density_cmap / style.density_norm, shared with Fig3a_backend_agreement."""
    # GAPX and L are set by the ligand_RMSD row, whose tick labels are "100": a gap that
    # fits "0.5" puts that row's y labels inside the cell to its left.
    W, H = 93.0, 96.5
    PW = PH = 20.0
    L, GAPX, B, GAPY = 13.0, 9.0, 8.5, 11.0            # mm
    cm = style.density_cmap()

    fig = plt.figure(figsize=(style.mm(W), style.mm(H)))
    for r, (k, lab, lo, hi, t, ge, log, ticks) in enumerate(M):
        y0 = B + (len(M) - 1 - r) * (PH + GAPY)        # first metric on top
        fits(lab, PH, f"row label {k}")
        fig.text(2.2 / W, (y0 + PH / 2) / H, lab, rotation=90, ha="center",
                 va="center", color=palette.INK)
        for c, (a, b) in enumerate(PAIRS):
            x0 = L + c * (PW + GAPX)
            ax = fig.add_axes([x0 / W, y0 / H, PW / W, PH / H])
            x = d[f"{PFX[a]}_{k}"].to_numpy(float)
            y = d[f"{PFX[b]}_{k}"].to_numpy(float)
            # Tails into the edge bins, exactly as Fig3a_backend_agreement clips them.
            xc, yc = np.clip(x, lo, hi), np.clip(y, lo, hi)
            ext = (np.log10(lo), np.log10(hi)) * 2 if log else (lo, hi) * 2
            # gridsize 50 across 0..1 is two ipTM quanta per hex, so the lattice
            # samples the data instead of beating against it; the log panel is not
            # quantised and takes 40. Both from Fig3a_backend_agreement.
            hb = ax.hexbin(xc, yc, gridsize=(40 if log else 50), extent=ext,
                           xscale="log" if log else "linear",
                           yscale="log" if log else "linear",
                           cmap=cm, mincnt=1, linewidths=0, rasterized=True)
            style.density_norm(hb)
            ax.plot([lo, hi], [lo, hi], color=palette.NEUTRAL["500"], lw=0.4,
                    ls=":", zorder=3)
            ax.axvline(t, color=palette.INK, lw=0.6, ls="--", zorder=4)
            ax.axhline(t, color=palette.INK, lw=0.6, ls="--", zorder=4)
            if log:
                ax.set_xscale("log")
                ax.set_yscale("log")
            ax.set_xlim(lo, hi)
            ax.set_ylim(lo, hi)
            ax.set_xticks(ticks)
            ax.set_yticks(ticks)
            ax.set_xticklabels([f"{v:g}" for v in ticks])
            ax.set_yticklabels([f"{v:g}" for v in ticks])
            ax.minorticks_off()
            fits(a, PW, "x label")
            fits(b, PH, "y label")
            ax.set_xlabel(a, labelpad=1.5)
            ax.set_ylabel(b, labelpad=1.5)
            ax.tick_params(axis="both", length=2, width=0.5)
            # Two numbers, both load-bearing and neither a title: how well the pair
            # RANKS designs together, and how often it reaches the same verdict at the
            # threshold drawn on the cell.
            pa, pb = (x >= t, y >= t) if ge else (x <= t, y <= t)
            same = f"same call {100 * float((pa == pb).mean()):.1f}%"
            fits(same, PW, "same-call annotation")
            ax.text(0.0, 1.145, f"R² = {pearson(x, y) ** 2:.2f}",
                    transform=ax.transAxes, ha="left", va="baseline")
            ax.text(0.0, 1.02, same, transform=ax.transAxes, ha="left",
                    va="baseline", color=palette.NEUTRAL["600"])
    style.save(fig, out, "panel_a_metric_pairs", tight=False)
    plt.close(fig)


# --- (b) seed reproducibility per re-predictor ---------------------------------
def panel_b(seed, out):
    """One box per re-predictor over the 16 interface metrics, replacing the 48 jittered
    dots this panel carried until 2026-11-19 ("show it as a candle plot, or at least
    dont show all the models with different dots"). Whiskers are the full range, there
    being 16 values and no reason to hide any as outliers; box = quartiles, rule =
    median. A box hides that the ORDER is unanimous, so the two counts that carry it
    are annotated: ESMFold2 is the highest of the three on 16 of 16 metrics, Boltz-2 the
    lowest on 13 of 16. AF3 is annotated with nothing because there is nothing to say -
    it is highest on none and lowest on 3."""
    # 52 mm, not 45: the two annotations are one line each, sitting under the boxes
    # clear of the lowest whisker cap (0.151), and at 45 mm they ran together.
    W, HH = 52.0, 36.0
    fig, ax = plt.subplots(figsize=(style.mm(W), style.mm(HH)))
    x = np.arange(len(PRED)).astype(float)
    wide = seed.pivot(index="metric", columns="refolder", values="icc")
    vals = []
    for p in PRED:
        v = seed.loc[seed["refolder"] == p, "icc"].to_numpy(float)
        assert len(v) == 16, f"{p}: expected 16 metrics, got {len(v)}"
        vals.append(v)
    bp = ax.boxplot(vals, positions=x, widths=0.52, whis=(0, 100), showfliers=False,
                    patch_artist=True,
                    medianprops=dict(color=palette.INK, lw=0.9, solid_capstyle="butt"),
                    whiskerprops=dict(color=palette.INK, lw=0.5),
                    capprops=dict(color=palette.INK, lw=0.5))
    for p, box in zip(PRED, bp["boxes"]):
        box.set(facecolor=palette.REFOLDERS[p], edgecolor=palette.INK, linewidth=0.5)
    slot = (0.99 - 0.19) * W / len(PRED)
    for i, p in enumerate(PRED):
        hi = int((wide.idxmax(axis=1) == p).sum())
        lo = int((wide.idxmin(axis=1) == p).sum())
        if hi == 0 and lo <= 3:          # the middle re-predictor: nothing to say
            continue
        s = f"highest {hi}/16" if hi else f"lowest {lo}/16"
        fits(s, slot, f"{p} unanimity annotation")
        # Below the box: the lowest ICC of any re-predictor is 0.151, so one line at
        # 0.03 clears every whisker cap.
        ax.text(x[i], 0.03, s, ha="center", va="bottom", color=palette.NEUTRAL["600"])
    ax.set_xticks(x)
    ax.set_xticklabels(PRED)
    ax.set_xlim(-0.6, len(PRED) - 0.4)
    ax.set_ylim(0, 1)
    ax.set_yticks([0, 0.5, 1])
    ax.set_ylabel(f"ICC ({K_SEEDS} seeds)")
    ax.tick_params(axis="both", length=2, width=0.5)
    fig.subplots_adjust(left=0.19, right=0.99, top=0.97, bottom=0.11)
    # Was panel_c_seed_icc until 2026-11-20. The book letters panels by FILE-NAME
    # order, so dropping the gate panel moves this one to (b) and the stem has to
    # say so; the descriptive half is unchanged because what it draws is unchanged.
    style.save(fig, out, "panel_b_seed_icc", tight=False)
    plt.close(fig)


# --- numbers the caption and the README carry -------------------------------
def report(d, seed):
    print(f"sample {len(d):,} designs, {PER_ROUTE:,} per route x {len(ROUTES)} routes")
    print("BG+BIF sc_DockQ/ligand_RMSD = CORRECTED true-backbone values (parquet main columns)")

    print("\n(a) 3 x 3, nine cells")
    for k, lab, lo, hi, t, ge, log, _ in M:
        for a, b in PAIRS:
            x = d[f"{PFX[a]}_{k}"].to_numpy(float)
            y = d[f"{PFX[b]}_{k}"].to_numpy(float)
            pa, pb = (x >= t, y >= t) if ge else (x <= t, y <= t)
            ag = pa == pb
            extra = ""
            if log:
                lx, ly = (np.log10(np.clip(v, 1e-3, None)) for v in (x, y))
                extra = f" [log10 r={pearson(lx, ly):.4f} R2={pearson(lx, ly) ** 2:.4f}]"
            print(f"  {k:8s} {a:8s} vs {b:9s} r={pearson(x, y):.4f} "
                  f"R2={pearson(x, y) ** 2:.4f} rho={spearman(x, y):.4f} | "
                  f"same call at {t}: {int(ag.sum()):,}/{len(d):,} ({100 * ag.mean():.2f}%)"
                  f" [both pass {int((pa & pb).sum()):,}, only {a} {int((pa & ~pb).sum()):,},"
                  f" only {b} {int((~pa & pb).sum()):,}, neither {int((~pa & ~pb).sum()):,}]"
                  f"{extra}")
    print("  Fig3a_backend_agreement, campaign-wide (6,992,515 designs), AF3 vs ESMFold2:")
    print("    r 0.251 / 0.495 / 0.568 -> R2 0.063 / 0.245 / 0.323 (ipTM / sc_DockQ / ligand_RMSD)")
    print("    the middle column above is the route-balanced version of those three")

    print("\nTHE GATE AS APPLIED (no panel since 2026-11-20) - g1 ligand_RMSD, g2 sc_DockQ,")
    print("    g3 seq run, g4 foot, g5 ipTM")
    print("    g4 (foot <= 5, BOTH backends) is not a column in this export, so every")
    print("    chain below stops one term short of the campaign's own.")
    for p in PRED:
        parts = []
        for k, lab, _, _, t, ge, _, _ in M:
            m = term(d, p, k, t, ge).to_numpy()
            parts.append(f"{k} {'>=' if ge else '<='}{t}: {int(m.sum()):,}/{N_SAMPLE:,} "
                         f"({100 * m.mean():.3f}%)")
        print(f"  {p:9s} each term ALONE: " + " | ".join(parts))
        print(f"  {'':9s} CHAIN: " + " -> ".join(
            f"{lab} {int(m.sum()):,}/{N_SAMPLE:,} ({100 * m.mean():.4f}%)"
            for (_, lab, _, _, _, _), m in zip(GATE, chain(d, p))))
        g = metric_gate(d, p)
        print(f"  {'':9s} three METRIC terms only, no seq run: {int(g.sum()):,}/"
              f"{N_SAMPLE:,} ({100 * g.mean():.4f}%) | median ipTM "
              f"{d[f'{PFX[p]}_iptm'].median():.4f}")
        # Route balance is not campaign balance: reweight by the two per-route
        # denominators the export's README gives, as a bracket on the campaign value.
        for name, wt in (("designs generated", GENERATED), ("eligible pool", ELIGIBLE)):
            per = pd.Series(g, index=d["stack"].to_numpy()).groupby(level=0).mean()
            w = sum(per[k2] * wt[k2] for k2 in wt) / sum(wt.values())
            print(f"  {'':9s}   reweighted by {name}: {100 * w:.4f}%")
        per = pd.Series(g, index=d["stack"].to_numpy()).groupby(level=0).mean()
        print(f"  {'':9s}   per route: "
              + "  ".join(f"{k2} {100 * per[k2]:.3f}" for k2 in ROUTES))
    for flag, what in (("g1", "ligand_RMSD <= 3, EITHER backend"),
                       ("g2", "sc_DockQ >= 0.7, EITHER backend"),
                       ("g3", "longest single-residue run <= 5, NO backend"),
                       ("g5", "max ipTM >= 0.6, EITHER backend")):
        print(f"  campaign term {flag} ({what}): {int(d[flag].sum()):,}/{N_SAMPLE:,}"
              f" ({100 * d[flag].mean():.4f}%)")
    print("  campaign CHAIN (every metric term on EITHER gated backend): " + " -> ".join(
        f"{lab} {int(m.sum()):,}/{N_SAMPLE:,} ({100 * m.mean():.4f}%)"
        for (_, lab, _, _, _, _), m in zip(GATE, campaign_chain(d))))
    print(f"  ends at passes_all = {int(d['passes_all'].sum()):,}/{N_SAMPLE:,}; the")
    print("  CAMPAIGN's own five-step chain is corpus-wide and is NOT this (README).")
    # The finding the gate panel used to draw: the ipTM term is applied LAST, and it
    # is the step at which the three re-predictors differ most.
    print("\n  THE FINDING - ipTM is where Boltz-2 differs, and only there:")
    for p in PRED:
        c = [int(m.sum()) for m in chain(d, p)]
        alone = int(term(d, p, "iptm", 0.6, True).sum())
        print(f"    {p:9s} ipTM alone {alone:,}/{N_SAMPLE:,} ({100 * alone / N_SAMPLE:.3f}%)"
              f" | at the ipTM STEP {c[-2]:,} -> {c[-1]:,}, keeps {100 * c[-1] / c[-2]:.1f}%")
    ia = {p: int(term(d, p, "iptm", 0.6, True).sum()) for p in PRED}
    fin = {p: int(chain(d, p)[-1].sum()) for p in PRED}
    print(f"    On ipTM alone Boltz-2 passes {ia['Boltz-2'] / ia['ESMFold2']:.2f}x ESMFold2;")
    print(f"    through the whole chain it ends on {fin['Boltz-2'] / fin['ESMFold2']:.2f}x.")
    print("    Its permissiveness is confined to the one metric that never reads the")
    print("    designed backbone - which is why the ipTM step barely cuts it.")
    print("  NOTE on wording: Boltz-2 does NOT rank designs least like the other two.")
    print("    By Spearman rho on ipTM the WEAKEST pair is AF3-ESMFold2 (0.247), the two")
    print("    the gate used; Boltz-2 correlates better with each (0.285, 0.347). The")
    print("    Boltz-2 claim is true of VERDICT agreement at 0.6 only - say which.")

    print(f"\n(b) seed ICC, {N_SEED_DESIGNS:,} designs per re-predictor, {K_SEEDS} seeds")
    w = seed.pivot(index="metric", columns="refolder", values="icc")
    for p in PRED:
        v = seed.loc[seed["refolder"] == p, "icc"]
        print(f"  {p:9s} min {v.min():.4f} q1 {v.quantile(.25):.4f} "
              f"median {v.median():.4f} q3 {v.quantile(.75):.4f} max {v.max():.4f} | "
              f"highest of the three on {int((w.idxmax(axis=1) == p).sum())}/16 metrics, "
              f"lowest on {int((w.idxmin(axis=1) == p).sum())}/16")
    print("  metrics where Boltz-2 is NOT the lowest:",
          list(w.index[w.idxmin(axis=1) != "Boltz-2"]))


if __name__ == "__main__":
    main()
