"""metric_reliability - how precisely can each interface metric be measured at all?

Two bare panels. The five refolding seeds are the instrument, not the subject: re-running
the same design under a different seed is how the measurement error is exposed, and the
quantity on the panels is the metric's own reliability. Renamed from seed_reliability
2026-11-19 at Leo's request ("more a precision / reliability test, the seed was just a
tool"); see README for why the slug is reliability and not precision.

Rebuilt 2026-11-18 from ~/rna_binder_plots/thesis_figures/fig_S2.py (2026-10-02), which
predates the style contract: it drew chart titles, panel letters and in-plot legend
boxes, and hand-picked re-predictor hexes.

Vocabulary, settled 2026-11-19: "re-predictor", never "refolder", in all prose, labels
and legends. The identifiers keep the old word on purpose and are not to be "fixed" -
`refolder` is a column in the source CSVs and `palette.REFOLDERS` is the palette's API.

Re-predictor encoding: palette v0.6 gives re-predictors their own colours
(palette.REFOLDERS, COLORS.md rule 7), so they neither borrow a design-route colour nor
need a texture. The earlier theme-step + hatch encoding was withdrawn 2026-11-19 at Leo's
request; see README. Order is palette.REFOLDER_ORDER, the pipeline's, never the result's.

scipy is unusable in this venv (_propack dlopen), so nothing here calls it: panel (b)
counts seeds with pandas, and the stored ICC and Spearman values are read, not recomputed
- the hand-rolled statistics are plain counts and min/max over the stored tables.
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from common import paths, style, palette

FIG = Path(__file__).resolve().parent.name
SRC = Path.home() / "rna_binder_plots" / "seed_metric_reliability"
N_DESIGNS = 2000          # designs per re-predictor, 5 seeds each
K_SEEDS = 5

# Fixed per re-predictor, everywhere. Order is the pipeline's, never the result's.
REFOLDERS = palette.REFOLDER_ORDER
# backend code in seed_metrics_merged_v2.csv -> re-predictor name
BACKEND = {"af3": "AF3", "b2": "Boltz-2", "esm": "ESMFold2"}
# AF3 carries six seeds (1-6); the stored ICCs use 2-6, so seed 1 is unused.
AF3_SEEDS = range(2, 7)
FAM_ORDER = ["SASA", "contacts", "BUNS", "H-bonds", "pi-stack"]

# Panel (a) is laid out in mm, not in subplot fractions, because the bottom
# margin has to hold a MEASURED rotated label and the gap between the sets of
# three has to be a real distance. Flipped to the full double-column width on
# 2026-11-20 (Leo: "Flip a to be horizontal full width of page. Give the sets
# of three a bit more space.").
FW_A = style.CFG["width_double_mm"]        # 178 mm
MARGIN_A = dict(l=10.0, r=2.0, t=4.0, b=19.5)   # mm; b holds the 45 deg labels
PLOT_H_A = 45.0            # mm of bar height for ICC 0..1
BAR_W = 0.22               # data units per bar: 3 x 0.22 leaves a 0.34 gap
ROT_A = 45                 # metric labels; 16 horizontal ones cannot fit
KEY_SW, KEY_GAP, KEY_SEP = 1.3, 0.9, 4.0   # mm, as in backbone_landscape's key

# Shorter on the panel, expanded in the caption (STYLE_NOTES section 3).
RELABEL = {
    "pi-stack count (no nt1)": "π-stack count",
    "pi-stack score (no nt1)": "π-stack score",
}


def bar_style(r):
    """No edge: every re-predictor colour is >= 4.7:1 on white, so COLORS.md's
    hairline-under-3:1 rule does not bite."""
    return dict(facecolor=palette.REFOLDERS[r])


def refolder_key(fig, y_mm=1.0):
    """Frameless key along the bottom, shared by the figure's panels (section 6).

    Entry widths are MEASURED with style.text_w, not laid out on a fixed
    per-entry pitch. The old fixed 0.30-of-figure pitch was tuned for an 85 mm
    canvas; on panel (a)'s 178 mm it would have scattered three short labels
    over 107 mm of empty paper.
    """
    fw, fh = (v * 25.4 for v in fig.get_size_inches())      # canvas in mm
    fs = style.CFG["font_size_pt"]
    tw = [style.text_w(r, fs) * 25.4 / 72 for r in REFOLDERS]
    total = sum(KEY_SW + KEY_GAP + w for w in tw) + KEY_SEP * (len(tw) - 1)
    assert total < fw, f"key is {total:.1f} mm, wider than the {fw:.0f} mm canvas"
    x = (fw - total) / 2
    for r, w in zip(REFOLDERS, tw):
        fig.patches.append(plt.Rectangle(
            (x / fw, y_mm / fh), KEY_SW / fw, KEY_SW / fh,
            transform=fig.transFigure, figure=fig, **bar_style(r)))
        fig.text((x + KEY_SW + KEY_GAP) / fw, y_mm / fh, r, fontsize=fs,
                 va="baseline", ha="left")
        x += KEY_SW + KEY_GAP + w + KEY_SEP


def main():
    style.apply()
    out = paths.fig_output(FIG)

    icc = pd.read_csv(SRC / "seed_metric_reliability.csv")
    dist = pistack_seed_counts()

    panel_a(icc, out)
    panel_b(dist, out)
    # Not a panel any more (the CV-vs-ipTM bars were cut 2026-11-19): one caption
    # clause carries all 45 values, so the range still has to be re-derived.
    report(icc, dist, pd.read_csv(SRC / "metric_std_vs_quality.csv"))


def pistack_seed_counts():
    """designs x (how many of the 5 seeds show a pi-stack), 0..5, per re-predictor.

    Derived by counting seeds in seed_metrics_merged_v2.csv - the stored
    pistack_presence_reliability.csv carries only the two endpoints. Plain
    pandas counting, no scipy (it cannot import in this venv).

    The assert IS the check: reproducing both stored endpoints from the per-seed
    rows is what proves the column and the seed subset are the ones the stored
    summary used.
    """
    d = pd.read_csv(SRC / "seed_metrics_merged_v2.csv")
    stored = pd.read_csv(SRC / "pistack_presence_reliability.csv").set_index("refolder")
    out = {}
    for bk, name in BACKEND.items():
        s = d[d["backend"] == bk]
        if bk == "af3":
            s = s[s["seed"].isin(AF3_SEEDS)]
        per = s.assign(hit=s["pis_count_nodangle"] > 0).groupby("base_design_id")["hit"]
        assert (per.size() == K_SEEDS).all(), f"{name}: not {K_SEEDS} seeds per design"
        k = per.sum()
        assert len(k) == N_DESIGNS, f"{name}: {len(k)} designs, expected {N_DESIGNS}"
        out[name] = k.value_counts().reindex(range(K_SEEDS + 1), fill_value=0)
        for col, want in (("ever", (k > 0).sum()), ("always", (k == K_SEEDS).sum())):
            got = round(stored.loc[name, col] / 100 * N_DESIGNS)
            assert got == want, f"{name} {col}: stored {got}, per-seed {want}"
    return pd.DataFrame(out)[REFOLDERS]


# --- (a) ICC per metric, grouped by family ---------------------------------
def panel_a(icc, out):
    """Vertical bars across the full double-column width, 2026-11-20.

    Was 16 horizontal rows on 85 mm. Flipped at Leo's request, which costs the
    y-axis its free label room and buys the sets of three real air: at 178 mm a
    metric gets a 10.4 mm slot, three 0.22-unit bars span 0.66 of it and the
    0.34 left over is a 3.5 mm gap between neighbouring sets - against 0.22
    (2.3 mm) on the old row pitch.

    The family separators and names moved with the flip. Hairlines were
    horizontal rules between row blocks with the family named at the right
    edge; they are now VERTICAL rules between column blocks with the family
    named ABOVE its block, horizontally. The right edge no longer exists as a
    place (it is one metric, not one family), and the top is the only clear
    strip left once the bottom carries 45 deg metric labels - so the names go
    there, unrotated, which is also the only orientation that does not add a
    second reading angle to the panel.
    """
    order = (icc.groupby(["family", "metric"], sort=False)["icc"].mean()
             .reset_index())
    order["family"] = pd.Categorical(order["family"], FAM_ORDER, ordered=True)
    order = order.sort_values(["family", "icc"], ascending=[True, False])
    metrics = list(order["metric"])
    fams = list(order["family"])
    labels = [RELABEL.get(m, m) for m in metrics]

    x = np.arange(len(metrics)).astype(float)
    m, fs = MARGIN_A, style.CFG["font_size_pt"]
    fh = m["t"] + PLOT_H_A + m["b"]
    pw = FW_A - m["l"] - m["r"]                     # plot width, mm
    xlim = (-0.6, len(metrics) - 0.4)
    pitch = pw / (xlim[1] - xlim[0])                # mm per metric slot

    # Everything that has to fit is MEASURED, never estimated (0.52 em/char is
    # out by 10-33% per string here). Three things can collide at 45 deg:
    mmpt = 25.4 / 72
    lh = fs * mmpt * 1.15                           # one line of type, mm
    diag = max(style.text_w(t, fs) * mmpt for t in labels) / 2 ** 0.5
    # (i) a label running off the bottom of the canvas,
    assert diag + lh + KEY_SW + 2 * KEY_GAP < m["b"], (
        f"longest metric label drops {diag + lh:.1f} mm below the axis; "
        f"the bottom margin is {m['b']:.1f} mm and the key needs part of it")
    # (ii) two neighbouring labels overlapping. At 45 deg they are parallel
    #      lines offset perpendicular by pitch/sqrt(2), so that gap - not the
    #      label width - is what has to clear one line of type.
    assert pitch / 2 ** 0.5 > lh, (
        f"metric labels collide: {pitch / 2 ** 0.5:.2f} mm perpendicular gap "
        f"against {lh:.2f} mm of type - widen the panel or rotate further")
    # (iii) a family name wider than the block it names (the two-metric blocks
    #       BUNS and pi-stacking are the tight ones).
    fam_label = {f: ("\u03c0-stacking" if f == "pi-stack" else f) for f in FAM_ORDER}
    for f in FAM_ORDER:
        n = fams.count(f)
        assert style.text_w(fam_label[f], fs) * mmpt < n * pitch, (
            f"family name {fam_label[f]!r} is wider than its {n}-metric block")

    fig = plt.figure(figsize=(style.mm(FW_A), style.mm(fh)))
    ax = fig.add_axes([m["l"] / FW_A, m["b"] / fh, pw / FW_A, PLOT_H_A / fh])
    for j, r in enumerate(REFOLDERS):
        v = (icc[icc["refolder"] == r].set_index("metric")
             .reindex(metrics)["icc"].values)
        ax.bar(x + (j - 1) * BAR_W, v, BAR_W * 0.9, **bar_style(r))

    # Family separators between the column blocks, and a name over each.
    for i in range(1, len(fams)):
        if fams[i] != fams[i - 1]:
            ax.axvline(i - 0.5, color=palette.GRID, lw=0.5, zorder=0)
    for f in FAM_ORDER:
        idx = [i for i, g in enumerate(fams) if g == f]
        ax.annotate(fam_label[f], (float(np.mean(idx)), 1.0),
                    xycoords=("data", "axes fraction"), xytext=(0, 2),
                    textcoords="offset points", fontsize=fs,
                    color=palette.NEUTRAL["600"], ha="center", va="bottom")

    ax.set_xticks(x)
    ax.set_xticklabels(labels, rotation=ROT_A, ha="right",
                       rotation_mode="anchor")
    ax.set_xlim(*xlim)
    ax.set_ylim(0, 1.0)
    ax.set_ylabel(f"ICC ({K_SEEDS} seeds)")
    ax.tick_params(axis="both", length=2, width=0.5)
    refolder_key(fig, y_mm=0.9)
    # not style.save: bbox_inches="tight" would trim the authored 178 mm width
    # the margins above were measured against.
    for ext in ("pdf", "png", "svg"):
        fig.savefig(out / f"panel_a_icc_by_metric.{ext}")
    plt.close(fig)
    print(f"panel a: {FW_A:.0f} x {fh:.1f} mm, {len(metrics)} metrics at "
          f"{pitch:.2f} mm pitch, bar {BAR_W * pitch:.2f} mm, gap between sets "
          f"{(1 - 3 * BAR_W) * pitch:.2f} mm, longest label drops "
          f"{diag + lh:.1f} of {m['b']:.1f} mm")


# --- (b) in how many of the 5 seeds does a pi-stack survive? ---------------
def panel_b(dist, out):
    """Distribution over 1..5 seeds. The 0 column (no stack in any seed) is not
    drawn: at 1,535-1,716 designs it is 5x the tallest bar here and would flatten
    the shape the panel is about. It is 2,000 - ever, and the caption gives it."""
    seeds = np.arange(1, K_SEEDS + 1)
    w = 0.27
    fig, ax = plt.subplots(figsize=(style.mm(85), style.mm(46)))
    for j, r in enumerate(REFOLDERS):
        ax.bar(seeds + (j - 1) * w, dist[r].reindex(seeds).values, w * 0.92,
               **bar_style(r))
        # Only the all-5 bars: 1 design in 2,000 is a bar you cannot see, and
        # that floor is the panel's point. Everything else goes in the caption.
        n = int(dist[r].loc[K_SEEDS])
        ax.annotate(f"{n:,}", (seeds[-1] + (j - 1) * w, n),
                    textcoords="offset points", xytext=(0, 2), ha="center",
                    fontsize=style.CFG["font_size_pt"])

    ax.set_xticks(seeds)
    ax.set_xlabel(f"Seeds with a \u03c0-stack (of {K_SEEDS})")
    ax.set_ylabel("Designs")
    ax.set_xlim(0.5, K_SEEDS + 0.5)
    ax.set_ylim(0, dist.loc[1].max() * 1.15)
    ax.tick_params(axis="both", length=2, width=0.5)
    # No key here: one frameless key serves the whole figure, on panel (a).
    fig.subplots_adjust(left=0.12, right=0.98, top=0.97, bottom=0.17)
    style.save(fig, out, "panel_b_pistack_reseed")
    plt.close(fig)


def report(icc, dist, cvq):
    """Numbers the caption has to carry, printed so RUNLOG can quote them."""
    fam = (icc.groupby("family")["icc"].agg(["min", "max"]).round(3)
           .reindex(FAM_ORDER))
    print(fam.to_string())
    best = (icc.loc[icc.groupby("metric")["icc"].idxmax(), "refolder"]
            .value_counts().to_dict())
    print("highest ICC per metric, by re-predictor:", best)
    def pct(n):
        v = n / N_DESIGNS * 100
        return f"{n:,}/{N_DESIGNS:,} ({v:.1f}%)" if v >= 1 else f"{n:,}/{N_DESIGNS:,} ({v:.2f}%)"

    for r in REFOLDERS:
        k = dist[r]
        ever, always = int(k[1:].sum()), int(k[K_SEEDS])
        print(f"{r}: never {pct(int(k[0]))}, ever {pct(ever)}, "
              f"all {K_SEEDS} {pct(always)}, "
              f"flips {(ever - always) / ever * 100:.1f}% of ever")
        print("   by seed count 1-5: "
              + ", ".join(f"{i}:{pct(int(k[i]))}" for i in range(1, K_SEEDS + 1))
              + f"  | exactly one seed = {k[1] / ever * 100:.1f}% of ever")

    # The cut panel (b)'s one surviving caption clause: all 45 values, not the 6
    # it drew. min/max by hand - no scipy in this venv, and none needed, the
    # Spearman rho itself is stored.
    rho = cvq["rho_cv_iptm"]
    assert (rho < 0).all(), "a non-negative rho_cv_iptm would break the clause"
    print(f"rho(CV, ipTM): {len(rho)} values = {cvq.metric.nunique()} metrics x "
          f"{cvq.refolder.nunique()} re-predictors, all negative, "
          f"{rho.min():.3f} to {rho.max():.3f}")


if __name__ == "__main__":
    main()
