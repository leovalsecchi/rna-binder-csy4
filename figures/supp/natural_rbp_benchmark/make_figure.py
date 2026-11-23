"""natural_rbp_benchmark - designed interfaces against natural RBP-RNA complexes across interface descriptors.

Three bare panels. The argument is a correction, not a list of deficits: panel (a) raises
the size confound, (b) applies it and shows that it destroys most of the apparent deficit,
and only the H-bond row survives - which is why the third panel then looks at the residues
that read bases. Do not flatten this into "designs make fewer H-bonds".

2026-11-19 (Leo): "I think I want to remove C." The H-bond distribution panel - raw per
complex above, per interface residue below - is gone. Its numbers are NOT: report() still
prints them, because the per-residue medians reproducing (b)'s 0.11x is the check that
this code path and the stored comparison table agree.

Letters: the book letters panels by filename order at assembly, so the base-reading panel
is drawn as (c) while its stem stays panel_d_*. The stem is deliberately not renamed -
an Illustrator assembly links panel files by name.

2026-11-19 (Leo, on the panels): "here you need to normalize for the interface residues.
but also we are analysing hbonds and very specific contacts for things we said before we
are not super sure about."

He is right on both counts, and the first was an internal contradiction: (a)'s whole
argument is that a raw extensive count is not comparable when the natural protein
interface is 2.18x larger, yet the H-bond and base-reading panels were raw per-complex
counts. Both were normalised PER PROTEIN INTERFACE RESIDUE, PER COMPLEX before
summarising - the same estimator (b) uses, never a ratio of medians - and the base-reading
panel says so on its axis.

On the second half: H-bonds and pi-stacking are the LEAST reproducible metrics measured
here (metric_reliability: ICC 0.163-0.514 and 0.151-0.519 across five refolding seeds,
against 0.634-0.928 for burial), so the H-bond result is carried by (b)'s single corrected
ratio and the caption's caveat, no longer by a panel of its own. The deficit is not
deleted - it is the one thing that survives the size correction - but it is not presented
as the firmest number in the figure.

Rebuilt 2026-11-18 from ~/rna_binder_plots/thesis_figures/fig_F8.py (2026-10-02), which
predates the style contract (4 chart titles, 4 panel letters, 3 in-plot legends) and
hand-picked every hex.

Colour here carries no route identity and needs none: natural is the reference
distribution, so it takes a NEUTRAL (COLORS.md rule 3), and the designs take MAIN
(rule 1). The one surviving metric is emphasised by WEIGHT and MARKER SIZE, never by a
fifth hue - STYLE_NOTES section 4. No re-predictor appears, so there is no hatch to remove
under palette v0.6 (COLORS.md rule 7); there never was any in this figure.

scipy is unusable in ~/.venvs/binder_map (its _propack extension fails to dlopen), so
every statistic here is hand-rolled on numpy/pandas. None of them is a test: the panels
report medians, means and counts only, which is the level of claim the seed-reliability
result supports.
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from common import paths, style, palette

FIG = Path(__file__).resolve().parent.name
SRC = Path.home() / "rbp_rna_db" / "subset59" / "nat_vs_design_findings"
SIZE_CSV = paths.INPUTS / "natural_rbp_benchmark_size_confound.csv"   # from make_inputs.py
SEL_CSV = paths.INPUTS / "natural_rbp_benchmark_selection.csv"        # from make_inputs.py

# Primary per-complex tables. The base-reading panel and report()'s H-bond check
# read these rather than the pre-aggregated nat_vs_design_findings artefacts
# (distributions.json, aa_base.csv), which are raw counts only and so cannot be
# normalised after the fact.
NAT_CSV = Path.home() / "rbp_rna_db" / "subset59" / "metrics_natural_59.csv"
DES_CSV = Path.home() / "Documents" / "metric_explorer" / "data" / "metrics_100k.csv"

N_NAT = 58          # non-redundant natural complexes (59 curated, 8JOZ dropped)
N_DES = 100_000

NAT_C = palette.NEUTRAL["500"]
DES_C = palette.MAIN
HERO = "n_h_bonds_interface"     # the one metric that survives normalisation
DEN = "n_contact_res_prot"       # the denominator panel (a) argues for

AAS = ["ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE",
       "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL"]

# Order: the survivor first, then the rest. Never sorted by value.
EXT_ORDER = ["n_h_bonds_interface", "n_contacts_base", "n_contacts_sugar",
             "n_heavy_contacts_4A", "n_contacts_backbone", "n_contacts_phosphate"]


def fs():
    return style.CFG["font_size_pt"]


def load_primaries():
    """The two per-complex tables, filtered exactly as normalize_metrics.py does."""
    nat = pd.read_csv(NAT_CSV)
    err = nat["dsasa_biotite_error"]     # 8JOZ drops: its biotite SASA failed
    nat = nat[err.isna() | (err.astype(str).str.strip() == "")]
    assert len(nat) == N_NAT, f"expected {N_NAT} natural complexes, got {len(nat)}"
    assert pd.to_numeric(nat[DEN], errors="coerce").min() > 0

    keep = {DEN, HERO} | {f"aa_base_{a}" for a in AAS}
    des = pd.read_csv(DES_CSV, low_memory=False, usecols=lambda c: c in keep)
    assert len(des) == N_DES, f"expected {N_DES:,} designs, got {len(des)}"
    # A design with no protein interface residue has no interface to normalise
    # by. Same `.where(den > 0)` drop as normalize_metrics.py, so the base-reading
    # panel rests on the same rows as (b).
    des = des[pd.to_numeric(des[DEN], errors="coerce") > 0]
    return nat, des


def per_res(df, col):
    """col / protein interface residues, PER COMPLEX. Never a ratio of medians."""
    return (pd.to_numeric(df[col], errors="coerce")
            / pd.to_numeric(df[DEN], errors="coerce"))


def main():
    style.apply()
    out = paths.fig_output(FIG)

    size = pd.read_csv(SIZE_CSV)
    sel = pd.read_csv(SEL_CSV)
    norm = pd.read_csv(SRC / "comparison_normalized.csv")
    nat, des = load_primaries()

    panel_a(size, out)
    panel_b(norm, out)
    aa = panel_d(nat, des, out)
    panel_e(sel, out)
    report(norm, nat, des, aa, sel)


# --- (a) the size confound: the objection, raised before the data is used ----
def panel_a(size, out):
    y = np.arange(len(size))[::-1].astype(float)

    fig, ax = plt.subplots(figsize=(style.mm(85), style.mm(32)))
    ax.barh(y, size["ratio_nat_over_des"], 0.46,
            facecolor=NAT_C, edgecolor=palette.INK, linewidth=0.4)
    ax.axvline(1.0, color=palette.INK, lw=0.7)

    for yy, r in zip(y, size.itertuples()):
        # Raw inputs beside the ratio (STYLE_NOTES section 5), and the n that
        # differs: dSASA is 53 of 58, the rest 58.
        n = "" if r.nat_n == N_NAT else f"  n = {r.nat_n}"
        ax.annotate(f"{r.ratio_nat_over_des:.2f}×   "
                    f"({r.nat_median:,.0f} vs {r.des_median:,.0f}){n}",
                    (r.ratio_nat_over_des, yy), textcoords="offset points",
                    xytext=(4, 0), va="center", fontsize=fs())

    ax.set_yticks(y)
    ax.set_yticklabels(size["label"])
    ax.set_xlim(0, 3.4)
    ax.set_ylim(-0.6, len(size) - 0.4)
    ax.set_xlabel("natural ÷ design (median)")
    ax.tick_params(axis="x", length=2, width=0.5)
    ax.tick_params(axis="y", length=0)
    fig.subplots_adjust(left=0.33, right=0.99, top=0.95, bottom=0.30)
    style.save(fig, out, "panel_a_size_confound")
    plt.close(fig)


# --- (b) what survives the correction ---------------------------------------
def panel_b(norm, out):
    ext = norm[norm["kind"] == "extensive"]
    # Clashes are AF3-model-vs-crystal confounded; normalisation does not fix
    # that, so the metric is EXCLUDED rather than normalised and shown.
    ext = ext[ext["metric"] != "n_interface_clashes"]
    ext = ext.set_index("metric").loc[EXT_ORDER].reset_index()
    ext["label"] = ext["label"].str.replace("(4A)", "(4 Å)", regex=False)
    y = np.arange(len(ext))[::-1].astype(float)

    fig, ax = plt.subplots(figsize=(style.mm(85), style.mm(42)))
    ax.axvline(1.0, color=palette.INK, lw=0.7, zorder=1)

    for yy, r in zip(y, ext.itertuples()):
        hero = r.metric == HERO
        ax.plot([r.raw_ratio, r.perres_ratio], [yy, yy],
                color=palette.NEUTRAL["300"], lw=0.8, zorder=2)
        ax.plot(r.raw_ratio, yy, "o", ms=3.0, mfc="white",
                mec=palette.NEUTRAL["500"], mew=0.7, zorder=3)
        ax.plot(r.perres_ratio, yy, "o", ms=5.0 if hero else 3.6,
                color=DES_C, zorder=4)
        ax.annotate(f"{r.perres_ratio:.2f}×", (r.perres_ratio, yy),
                    textcoords="offset points", xytext=(0, 4), ha="center",
                    fontsize=fs(), fontweight="bold" if hero else "normal")

    ax.set_yticks(y)
    ax.set_yticklabels(ext["label"])
    for lbl, m in zip(ax.get_yticklabels(), ext["metric"]):
        if m == HERO:
            lbl.set_fontweight("bold")
    ax.set_xscale("log")
    ax.set_xlim(0.055, 6.0)
    ax.set_ylim(-0.7, len(ext) - 0.3)
    ax.set_xticks([0.1, 0.25, 0.5, 1, 2, 4])
    ax.set_xticklabels(["0.1", "0.25", "0.5", "1", "2", "4"])
    ax.set_xlabel("design ÷ natural, per interface residue")
    ax.tick_params(axis="x", length=2, width=0.5)
    ax.tick_params(axis="y", length=0)
    fig.subplots_adjust(left=0.30, right=0.98, top=0.95, bottom=0.22)
    # Open = raw count, filled = corrected. Named once, under the row whose two
    # points sit furthest apart on the log axis, so the two words cannot touch.
    gap = (np.log10(ext["perres_ratio"]) - np.log10(ext["raw_ratio"])).abs()
    i = int(gap.idxmax())
    r0 = ext.iloc[i]
    for x, lab in ((r0.raw_ratio, "raw"), (r0.perres_ratio, "per residue")):
        ax.annotate(lab, (x, y[i]), textcoords="offset points",
                    xytext=(0, -9), ha="center", fontsize=fs(),
                    color=palette.NEUTRAL["600"])
    style.save(fig, out, "panel_b_normalisation")
    plt.close(fig)


# --- (d) which residues read bases ------------------------------------------
def panel_d(nat, des, out):
    # Normalised per complex before the mean, like every other panel. The mean
    # is kept rather than the median: aa_base_X is zero for most complexes, so
    # a median would be zero for all but arginine.
    # nat_raw/des_raw are carried but NOT drawn: they are the uncorrected means
    # the earlier version of this panel plotted, kept so report() can print both
    # factors side by side and the README's raw column is re-derived, not quoted.
    t = pd.DataFrame([{"res": a,
                       "nat_mean": per_res(nat, f"aa_base_{a}").mean(),
                       "des_mean": per_res(des, f"aa_base_{a}").mean(),
                       "nat_raw": pd.to_numeric(nat[f"aa_base_{a}"]).mean(),
                       "des_raw": pd.to_numeric(des[f"aa_base_{a}"]).mean()}
                      for a in AAS])
    top = (t.sort_values("nat_mean", ascending=False).head(12)
           .reset_index(drop=True))
    y = np.arange(len(top))[::-1].astype(float)
    h = 0.38

    fig, ax = plt.subplots(figsize=(style.mm(85), style.mm(58)))
    ax.barh(y + h / 2, top["nat_mean"], h, facecolor=NAT_C,
            edgecolor=palette.INK, linewidth=0.4)
    ax.barh(y - h / 2, top["des_mean"], h, facecolor=DES_C,
            edgecolor=palette.INK, linewidth=0.4)

    # Every row carries its factor now, not only the large ones: after the
    # correction several rows INVERT (designs above naturals), and a rule that
    # annotated only the deficits would hide exactly those.
    for yy, r in zip(y, top.itertuples()):
        ax.annotate(f"{r.nat_mean / max(r.des_mean, 1e-9):.1f}×",
                    (max(r.nat_mean, r.des_mean), yy),
                    textcoords="offset points", xytext=(3, 0), va="center",
                    fontsize=fs(), color=palette.NEUTRAL["600"])

    ax.set_yticks(y)
    ax.set_yticklabels(top["res"])
    ax.set_xlim(0, 1.55)
    ax.set_ylim(-0.7, len(top) - 0.3)
    ax.set_xlabel("Base contacts per interface residue (mean)")
    ax.tick_params(axis="x", length=2, width=0.5)
    ax.tick_params(axis="y", length=0)
    fig.subplots_adjust(left=0.17, right=0.98, top=0.98, bottom=0.14)
    # Direct labels on the top row rather than a legend box (section 6).
    r0 = top.iloc[0]
    ax.annotate(f"natural (n = {len(nat)})", (r0["nat_mean"] * 0.5, y[0] + h / 2),
                ha="center", va="center", fontsize=fs(), color="white")
    ax.annotate(f"designs (n = {len(des):,})", (r0["des_mean"] * 0.5, y[0] - h / 2),
                ha="center", va="center", fontsize=fs(), color="white")
    style.save(fig, out, "panel_d_base_reading_residues")
    plt.close(fig)
    return top


# --- (e, drawn as d) the benchmark set itself -------------------------------
def panel_e(sel, out):
    """How the 58 natural complexes were selected, and the second n the figure carries.

    Leo 2026-11-20: "Can you add a general overview of the dataset we used? Either the
    way we selected them or just their size distribution or something like that."

    The SELECTION, not the size distribution, because the size of this set is already
    the subject of panel (a) - it is the whole objection the figure is built on - while
    nothing on the page says what the set is or where 58 complexes came from. The
    content of the set goes in the caption, which is what a caption is for.

    This panel also puts the figure's TWO n's in one place: 58 is the set, and the 53
    on panel (a)'s ΔSASA bar is a per-metric subset of it, drawn here as a subordinate
    open bar under the final step rather than as a fifth step, because nothing was
    selected out at that point - a measurement is simply missing for five complexes.

    Every bar is 85 mm wide like the rest of the figure and takes the NEUTRAL the
    natural side wears in (a) and (c): no route appears here either (COLORS.md rule 3).
    """
    steps = sel.iloc[:4]
    n0 = int(steps.n.iloc[0])
    y = np.arange(len(steps))[::-1].astype(float) + 1.0
    y_sub = y[-1] - 0.62          # closer than a full row: a subset, not a step
    n_last, n_dsasa = int(steps.n.iloc[-1]), int(sel.n_dsasa.iloc[-1])

    # Every bar is annotated, so the axis has to hold bar + gap + its own label. The
    # label widths are MEASURED, not estimated at ~0.52 em/char, which is out by 10-33%
    # per string and would clip the longest one or waste a third of the axis.
    labels = [f"{n:,}" if i == 0 else f"{n:,}/{n0:,} ({100 * n / n0:.1f}%)"
              for i, n in enumerate(steps.n)]
    sub_label = f"{n_dsasa}/{n_last} ({100 * n_dsasa / n_last:.1f}%)"
    gap_pt = 3.0
    ax_w_pt = 85.0 / 25.4 * 72 * (0.98 - 0.36)      # right - left, below
    need = [(v, style.text_w(t, fs())) for v, t in
            zip(list(steps.n) + [n_dsasa], labels + [sub_label])]
    xmax = max(v / (1 - (gap_pt + w) / ax_w_pt) for v, w in need) * 1.01

    fig, ax = plt.subplots(figsize=(style.mm(85), style.mm(36)))
    for yy, r, lab in zip(y, steps.itertuples(), labels):
        last = r.Index == steps.index[-1]
        ax.barh(yy, r.n, 0.46, facecolor=NAT_C if last else palette.NEUTRAL["300"],
                edgecolor=palette.INK, linewidth=0.4)
        ax.annotate(lab, (r.n, yy), textcoords="offset points", xytext=(gap_pt, 0),
                    va="center", fontsize=fs(),
                    fontweight="bold" if last else "normal")
    ax.barh(y_sub, n_dsasa, 0.30, facecolor=palette.SURFACE, edgecolor=NAT_C,
            linewidth=0.6)
    ax.annotate(sub_label, (n_dsasa, y_sub), textcoords="offset points",
                xytext=(gap_pt, 0), va="center", fontsize=fs(),
                color=palette.NEUTRAL["600"])

    ax.set_yticks(list(y) + [y_sub])
    ax.set_yticklabels(list(steps.label) + ["of which ΔSASA measured"])
    for lbl in ax.get_yticklabels()[-1:]:
        lbl.set_color(palette.NEUTRAL["600"])
    ax.set_xlim(0, xmax)
    ax.set_ylim(y_sub - 0.45, y[0] + 0.55)
    ax.set_xticks([0, 100, 200, 300])
    ax.set_xlabel("Natural RBP–RNA complexes")
    ax.tick_params(axis="x", length=2, width=0.5)
    ax.tick_params(axis="y", length=0)
    fig.subplots_adjust(left=0.36, right=0.98, top=0.97, bottom=0.21)
    style.save(fig, out, "panel_e_dataset_overview")
    plt.close(fig)


def report(norm, nat, des, aa, sel):
    n0 = int(sel.n.iloc[0])
    print("selection of the benchmark set (drawn as (d)):")
    for r in sel.itertuples():
        print(f"  {r.label:30s} {r.n:4d}/{n0} ({100 * r.n / n0:5.1f}% of the search)"
              f"   {r.note}")
    n_last, n_dsasa = int(sel.n.iloc[-1]), int(sel.n_dsasa.iloc[-1])
    print(f"  ΔSASA measurable {n_dsasa}/{n_last} "
          f"({100 * n_dsasa / n_last:.1f}% of the set) - panel (a)'s ΔSASA bar\n")

    ext = norm[norm["kind"] == "extensive"].set_index("metric")
    print("(b) raw -> per-residue ratio (design / natural):")
    for m in EXT_ORDER:
        r = ext.loc[m]
        print(f"  {r['label']:22s} {r['raw_ratio']:.2f} -> {r['perres_ratio']:.2f}")

    # Not a panel since 2026-11-19 (Leo: "I think I want to remove C"), still the
    # check: the per-residue medians below must reproduce (b)'s 0.11x ratio, which
    # is what ties this code path to the stored comparison table.
    print(f"\nH-bond check (no panel) n = {len(nat)} natural, {len(des):,}/{N_DES:,} "
          f"({100 * len(des) / N_DES:.1f}%) designs with >= 1 interface residue")
    raw_n, raw_d = pd.to_numeric(nat[HERO]), pd.to_numeric(des[HERO])
    pr_n, pr_d = per_res(nat, HERO), per_res(des, HERO)
    print(f"  raw median        nat {raw_n.median():.0f}  des {raw_d.median():.0f}"
          f"   (design max {raw_d.max():.0f}, natural max {raw_n.max():.0f})")
    print(f"  per-residue median nat {pr_n.median():.3f}  des {pr_d.median():.3f}"
          f"   ratio {pr_d.median() / pr_n.median():.2f}")
    # Overlap, stated as a count: the 5th percentile of the natural
    # distribution is a far steadier threshold than its minimum at n = 58.
    p5 = pr_n.quantile(0.05)
    k = int((pr_d >= p5).sum())
    print(f"  designs reaching the natural 5th percentile ({p5:.3f}): "
          f"{k:,}/{len(pr_d):,} ({100 * k / len(pr_d):.1f}%)")
    z = int((raw_d == 0).sum())
    print(f"  designs with no interface H-bond at all: {z:,}/{len(raw_d):,} "
          f"({100 * z / len(raw_d):.1f}%)")

    print("\nbase contacts, natural ÷ design, raw -> per interface residue:")
    for r in aa.itertuples():
        print(f"  {r.res}  {r.nat_raw:6.2f}/{r.des_raw:5.2f} "
              f"{r.nat_raw / max(r.des_raw, 1e-9):5.1f}×  ->  "
              f"{r.nat_mean:.3f}/{r.des_mean:.3f} "
              f"{r.nat_mean / max(r.des_mean, 1e-9):5.1f}×")


if __name__ == "__main__":
    main()
