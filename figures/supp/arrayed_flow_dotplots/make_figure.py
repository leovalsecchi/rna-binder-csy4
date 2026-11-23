"""arrayed_flow_dotplots - the events behind `code/supp/arrayed_flow`.

That figure draws FlowJo's per-well summary of the LV2.33 / LV2.34.1 arrayed
re-test. This one draws the cells. Two bare panels:

  (a) the gating hierarchy, on one well, left to right as it is applied:
      Cells -> Single Cells -> the mTagBFP2 / mScarlet plot the rest of the
      figure is made of. Percentages are in-gate, of the parent.
  (b) six binder columns x their four targets, 24 density plots on one shared
      count scale and one shared pair of axes, each with the mScarlet+ gate
      drawn and its in-gate percentage printed.

GATE AND READER ARE NOT THIS FIGURE'S. Both are imported from
`LV_FlowCytometry/analysis/lv2_33_screen.py`, the lab's own analysis of these
files, through `make_inputs.py`. See that module's docstring for the hierarchy,
the dead FCS path in the original script, and the two checks that this figure's
numbers do not contradict `code/supp/arrayed_flow`.

Run: ~/.venvs/binder_map/bin/python3 code/supp/arrayed_flow_dotplots/make_figure.py
(`make_inputs.py` first; it caches the events and the per-well table.)

WHAT THE DESIGN FOLLOWS, and from where. Dot plots have conventions, so none of
this is invented here:
  * The gating hierarchy is drawn left to right with the in-gate percentage on
    every plot, and each axis names its analyte - the European Journal of
    Immunology's own figure policy for flow data (Eur J Immunol 2021;51:2703,
    which points at "Guidelines for the use of flow cytometry and cell sorting
    in immunological studies", 2nd ed. Eur J Immunol 2019, 3rd ed. 2021).
  * DENSITY shading, not opaque dots. At 60-80k events a dot plot saturates and
    hides exactly the subpopulation the figure is about; BD's and Cytobank's
    display guidance both say so, and recommend a density or contour display
    instead.
  * BIEXPONENTIAL axes, not plain log, so the events at and below zero are
    drawn rather than truncated - the MIFlowCyt-aligned facility guidance
    (Univ. Fribourg CAF) and the same Cytobank note. Implemented as asinh,
    which is what a logicle axis is; `COFACTOR` sets the width of the linear
    region.
  * ONE scale for every sample in the figure, both axes and the count ramp, so
    panels are comparable. Mixing scales between samples is the single thing
    all three sources warn about.
Not followed: their quadrant gates with four percentages (this assay has one
positivity gate on a transfected population, so one number is the honest one),
their "FL1"-style axis names (the fluorophore is what a reader needs; the
detector is in the caption), and any panel title - STYLE_NOTES section 1.
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib.patches import Rectangle

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from common import palette, paths, style  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parent))
import make_inputs as MI  # noqa: E402
S = MI.S                                   # lv2_33_screen: gate thresholds

FIG = Path(__file__).resolve().parent.name
WORK = paths.ROOT / "working" / FIG

# --- the biexponential axis ------------------------------------------------
# asinh(x / COFACTOR): linear inside +-COFACTOR, logarithmic outside, which is
# what a logicle/biexponential axis is. 200 puts the linear region over the
# detector's own noise - the untransfected control's 99.99th percentile in 561
# is 138-212 across its four wells - so the OFF population is spread out
# instead of piled on a log axis' left edge, and the mScarlet+ gate at 500 sits
# clear of it in the logarithmic part.
COFACTOR = 200.0
DECADES = [1e2, 1e3, 1e4, 1e5]
LABELLED = {1e3: "10$^3$", 1e5: "10$^5$"}   # 5.5 pt: label every other decade
GRID_N = 70                                 # density bins per axis


def tr(x):
    return np.arcsinh(np.asarray(x, float) / COFACTOR)


TICKS = [0.0] + [tr(v) for v in DECADES]
TICKLAB = ["0"] + [LABELLED.get(v, "") for v in DECADES]
# One pair of limits for all 25 plots of the two panels (the convention every
# display guide insists on). Set from the detector range, not per panel.
XLIM = (tr(-900), tr(2.62e5))
YLIM = (tr(-900), tr(2.62e5))

# --- sizes -----------------------------------------------------------------
# make_supp_book packs panels into rows no wider than the A4 text block
# (210 - 2 x 14 = 182 mm) with a 6 mm gutter, by the size each panel is EMITTED
# at. These two keep (a) and (b) on their own rows, stacked, inside one page.
BOOK_AVAIL, BOOK_GUTTER, BOOK_MARGIN, BOOK_TITLE = 182.0, 6.0, 14.0, 10.0
A4_H = 297.0
CELL_MM = 25.0                 # one square dot plot
A_W_MM, A_H_MM = 3 * CELL_MM + 16, CELL_MM + 11
B_W_MM = 4 * CELL_MM + 30      # 4 targets + row labels + colourbar
B_H_MM = 6 * CELL_MM + 12

# --- colours ---------------------------------------------------------------
# The ramp encodes COUNT, which is a data series, so it is the house density
# ramp and not a molecular colour (COLORS.md rule 8). The one molecular colour
# in the figure is the mScarlet gate, which is a gate on mScarlet.
MSC = palette.SYNBIO["reporter_alt"]
SUB = palette.NEUTRAL["500"]


# --------------------------------------------------------------------- data
def load():
    d = pd.read_csv(WORK / "well_stats.tsv", sep="\t").set_index("well")
    ev = np.load(WORK / "events.npz")
    onoff = pd.read_csv(paths.INPUTS / "lv2_34_1_flow_on_off.csv", index_col=0)
    idx = pd.read_csv(paths.INPUTS / "candidates32_index.tsv",
                      sep="\t").set_index("flow_well")
    rows = []
    for binder, kind, why in MI.choose(onoff):
        v = onoff.ONOFF_avg.get(binder, np.nan)
        if kind == "untransfected":
            name, note = "untransfected", "control"
        elif kind.startswith("Csy4"):
            # The plan's NAMES, as arrayed_flow prints them. 84 / 8.4 ng were
            # delivered; that belongs in the caption, not on a 25 mm panel.
            name, note = "Csy4-H29A", f"{kind.split()[1]} ng   {v:.0f}×"
        else:
            name = f"binder {int(idx.loc[binder, 'twist_id'])}"
            note = f"{kind.split()[0]}   {v:.2f}×"
        rows.append(dict(binder=binder, kind=kind, name=name, note=note,
                         onoff=v, why=why))
    return d, ev, pd.DataFrame(rows).set_index("binder")


def density(ax, x, y, cmap, norm=None):
    """One density plot. Uniform bins in the transformed axis, so the bins are
    even on the page and the shading is a count, not an artefact of the scale."""
    h, xe, ye = np.histogram2d(
        tr(x), tr(y), bins=GRID_N,
        range=[list(XLIM), list(YLIM)])
    m = ax.pcolormesh(xe, ye, np.ma.masked_less(h.T, 1), cmap=cmap, norm=norm,
                      rasterized=True, shading="flat")
    return m, h


def frame(ax, xlab=False, ylab=False):
    ax.set_xlim(*XLIM)
    ax.set_ylim(*YLIM)
    ax.set_xticks(TICKS)
    ax.set_yticks(TICKS)
    ax.set_xticklabels(TICKLAB if xlab else [])
    ax.set_yticklabels(TICKLAB if ylab else [])
    ax.tick_params(length=1.4, width=0.4, pad=1.2)
    for s in ax.spines.values():
        s.set_linewidth(0.4)
    ax.spines["top"].set_visible(True)
    ax.spines["right"].set_visible(True)


def pct(n, N):
    """Every rate as n/N (%), which is also what a reader needs to judge it."""
    return f"{n:,}/{N:,}\n{100 * n / N:.1f}%"


# ------------------------------------------------------------------ panel a
def panel_a(d, ev, cmap, norm, out):
    """The gate hierarchy, left to right, on the clearest ON well."""
    w = MI.GATE_WELL
    fsc, fsch, ssc = ev["gating/FSC-A"], ev["gating/FSC-H"], ev["gating/SSC-A"]
    bfp, rfp = ev[f"gating/{S.BFP}"], ev[f"gating/{S.RFP}"]
    cells = ((fsc >= 3e4) & (fsc <= 2.5e5) & (ssc >= 1e4) & (ssc <= 2.2e5))
    ratio = fsc / np.clip(fsch, 1, None)
    sing = cells & (ratio >= 1.0) & (ratio <= 1.9)

    fig, ax = plt.subplots(1, 3, figsize=(style.mm(A_W_MM), style.mm(A_H_MM)))

    a = ax[0]
    a.hexbin(fsc, ssc, gridsize=60, cmap=cmap, norm=norm, linewidths=0,
             extent=(0, 2.62e5, 0, 2.62e5), rasterized=True, mincnt=1)
    a.add_patch(Rectangle((3e4, 1e4), 2.2e5, 2.1e5, fill=False, lw=0.5,
                          ec=palette.INK))
    a.set_xlim(0, 2.62e5)
    a.set_ylim(0, 2.62e5)
    a.set_xlabel("FSC-A")
    a.set_ylabel("SSC-A")
    a.annotate("cells\n" + pct(int(cells.sum()), len(fsc)), (0.97, 0.95),
               xycoords="axes fraction", ha="right", va="top", linespacing=1.15)

    b = ax[1]
    sub = cells
    b.hexbin(fsc[sub], fsch[sub], gridsize=60, cmap=cmap, norm=norm,
             linewidths=0, extent=(0, 2.62e5, 0, 2.62e5), rasterized=True,
             mincnt=1)
    xs = np.array([0, 2.62e5])
    for k in (1.0, 1.9):
        b.plot(xs, xs / k, lw=0.5, color=palette.INK)
    b.set_xlim(0, 2.62e5)
    b.set_ylim(0, 2.62e5)
    b.set_xlabel("FSC-A")
    b.set_ylabel("FSC-H")
    b.annotate("single cells\n" + pct(int(sing.sum()), int(cells.sum())),
               (0.97, 0.95), xycoords="axes fraction", ha="right", va="top",
               linespacing=1.15)

    c = ax[2]
    density(c, bfp[sing], rfp[sing], cmap, norm)
    gates(c)
    frame(c, xlab=True, ylab=True)
    c.set_xlabel("mTagBFP2 (405)")
    c.set_ylabel("mScarlet (561)")
    tx = sing & (bfp > 1000)
    on = tx & (rfp > S.RFP_POS)
    c.annotate("transfected\n" + pct(int(tx.sum()), int(sing.sum())),
               (0.97, 0.03), xycoords="axes fraction", ha="right", va="bottom",
               color=SUB, linespacing=1.15)
    c.annotate("mScarlet$^+$\n" + pct(int(on.sum()), int(tx.sum())),
               (0.97, 0.97), xycoords="axes fraction", ha="right", va="top",
               color=MSC, linespacing=1.15)

    for a_ in ax[:2]:
        a_.tick_params(length=1.4, width=0.4, pad=1.2)
        a_.set_xticks([0, 1e5, 2e5])
        a_.set_yticks([0, 1e5, 2e5])
        a_.set_xticklabels(["0", "100k", "200k"])
        a_.set_yticklabels(["0", "100k", "200k"])
        for s in a_.spines.values():
            s.set_linewidth(0.4)
        a_.spines["top"].set_visible(True)
        a_.spines["right"].set_visible(True)

    fig.tight_layout(pad=0.25, w_pad=1.0)
    style.save(fig, out, "panel_a_gating")
    plt.close(fig)
    return dict(n_events=len(fsc), n_cells=int(cells.sum()),
                n_single=int(sing.sum()), n_tx=int(tx.sum()),
                n_on=int(on.sum()), well=w)


def gates(ax):
    """The two gates, each drawn where it cuts: transfected, then mScarlet+."""
    ax.axvline(tr(1000), lw=0.5, color=SUB, dashes=(2.5, 1.6))
    ax.plot([tr(1000), XLIM[1]], [tr(S.RFP_POS)] * 2, lw=0.6, color=MSC)


# ------------------------------------------------------------------ panel b
def panel_b(d, ev, meta, cmap, norm, out):
    """Six binder columns x four targets, one shared scale."""
    targets = list(MI.ROW_TARGET.items())          # A posA ... D Neg_noScaf
    TLAB = {"posA": "posA", "posB": "posB", "NegA": "NegA",
            "Neg_noScaf": "no array"}
    nr, nc = len(meta), len(targets)
    fig = plt.figure(figsize=(style.mm(B_W_MM), style.mm(B_H_MM)))
    gs = fig.add_gridspec(nr, nc, left=24 / B_W_MM, right=1 - 6.5 / B_W_MM,
                          top=1 - 5.0 / B_H_MM, bottom=6.5 / B_H_MM,
                          wspace=0.09, hspace=0.09)
    printed = []
    for i, (binder, m) in enumerate(meta.iterrows()):
        plate, col = int(binder[1]), int(binder.split("_C")[1])
        for j, (row, target) in enumerate(targets):
            well = f"Plate{plate}_{row}{col}"
            ax = fig.add_subplot(gs[i, j])
            density(ax, ev[f"{well}/bfp"], ev[f"{well}/rfp"], cmap, norm)
            gates(ax)
            frame(ax, xlab=(i == nr - 1), ylab=(j == 0))
            r = d.loc[well]
            n_on = round(r.pct_mscarlet / 100 * r.n_tx)
            ax.annotate(f"{r.pct_mscarlet:.1f}%", (0.95, 0.95),
                        xycoords="axes fraction", ha="right", va="top",
                        color=MSC)
            printed.append(dict(binder=binder, target=target, well=well,
                                pct=r.pct_mscarlet, n_on=n_on, n_tx=r.n_tx,
                                n_single=r.n_gated))
            if i == 0:
                ax.set_title(TLAB[target], pad=1.6,
                             fontsize=plt.rcParams["font.size"])
            if j == 0:
                # Two tiers: what the well is, then what arrayed_flow scores it.
                ax.set_ylabel(m["name"], rotation=0, ha="right", va="center",
                              labelpad=14 if m["note"] else 4)
                if m["note"]:
                    ax.annotate(m["note"], (0, 0.5), xycoords="axes fraction",
                                xytext=(-16, 0), textcoords="offset points",
                                ha="right", va="center", color=SUB)
    # One colourbar for all 24, because one scale is the point.
    cax = fig.add_axes([1 - 5.2 / B_W_MM, 6.5 / B_H_MM, 1.6 / B_W_MM,
                        14 / B_H_MM])
    cb = fig.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=cmap), cax=cax)
    cb.outline.set_linewidth(0.4)
    cb.ax.tick_params(length=1.4, width=0.4, pad=1.2)
    cb.set_label("events per bin", labelpad=1.5)
    fig.text(0.5 + 12 / B_W_MM, 1.2 / B_H_MM, "mTagBFP2 (405)", ha="center",
             va="bottom")
    fig.text(1.2 / B_W_MM, 0.5, "mScarlet (561)", rotation=90, ha="left",
             va="center")
    style.save(fig, out, "panel_b_wells")
    plt.close(fig)
    return pd.DataFrame(printed)


# -------------------------------------------------------------------- check
def count_scale(ev, meta):
    """ONE count scale for every plot, spanning where the occupied bins are.

    LogNorm(1, max) spends half the ramp on the handful of bins at the tail
    (`style.density_norm`'s finding, measured on Fig3a). Per-panel percentiles
    are wrong HERE, though: the comparison between wells is this figure's whole
    content, so the span is taken over the POOLED bins of all 24 and used
    unchanged everywhere.
    """
    allh = []
    for binder in meta.index:
        plate, col = int(binder[1]), int(binder.split("_C")[1])
        for row in MI.ROW_TARGET:
            w = f"Plate{plate}_{row}{col}"
            h, _, _ = np.histogram2d(tr(ev[f"{w}/bfp"]), tr(ev[f"{w}/rfp"]),
                                     bins=GRID_N, range=[list(XLIM), list(YLIM)])
            allh.append(h[h >= 1])
    c = np.concatenate(allh)
    lo, hi = max(1.0, np.percentile(c, 25)), np.percentile(c, 99.5)
    return LogNorm(lo, hi), c


def report(d, meta, printed, gating, c, norm):
    """Every number the README quotes, re-derived and printed on every run."""
    print(f"\ngating well {gating['well']}: {gating['n_events']:,} events -> "
          f"{gating['n_cells']:,} cells -> {gating['n_single']:,} single cells "
          f"-> {gating['n_tx']:,} transfected -> {gating['n_on']:,} mScarlet+")
    print(f"count scale: LogNorm({norm.vmin:.0f}, {norm.vmax:.0f}) over "
          f"{len(c):,} occupied bins pooled across the 24 wells "
          f"(max bin {c.max():.0f})")
    print("\nwells drawn (mScarlet+ % of transfected, this figure's gate):")
    p = printed.pivot(index="binder", columns="target", values="pct")
    p = p[["posA", "posB", "NegA", "Neg_noScaf"]].loc[meta.index]
    print(p.round(2).to_string())
    on = p[["posA", "posB"]].mean(1) / p[["NegA", "Neg_noScaf"]].mean(1)
    cmp = pd.DataFrame({"this figure": on, "arrayed_flow": meta.onoff})
    cmp["ratio"] = cmp["this figure"] / cmp["arrayed_flow"]
    print("\nON/OFF, avg(posA,posB)/avg(NegA,no scaffold), two readouts:")
    print(cmp.round(3).to_string())
    print(f"\ntransfected, of single cells: "
          f"{printed.n_tx.div(printed.n_single).mul(100).min():.1f}-"
          f"{printed.n_tx.div(printed.n_single).mul(100).max():.1f}% "
          f"over the 24 wells; single cells {printed.n_single.min():,}-"
          f"{printed.n_single.max():,}")
    # The untransfected control has to BE background in both channels, or the
    # row that anchors the figure is not a control.
    u = printed[printed.binder == "P6_C3"]
    assert u.pct.max() < 0.5, f"P6_C3 is not background: {u.pct.max():.3f}%"
    assert (u.n_tx / u.n_single).max() < 0.01, "P6_C3 is not untransfected"
    print(f"P6_C3 confirmed untransfected from the events: "
          f"{(u.n_tx / u.n_single * 100).max():.3f}% transfected, "
          f"{u.pct.max():.3f}% mScarlet+, "
          f"max n_tx {u.n_tx.max()}")
    return cmp


def emitted(out, stem):
    head = (out / f"{stem}.svg").read_text()[:600]
    import re
    g = [re.search(rf'{k}="([\d.]+)pt"', head) for k in ("width", "height")]
    assert all(g), f"{stem}.svg has no pt size"
    return tuple(float(m.group(1)) * 25.4 / 72.0 for m in g)


def main():
    style.apply()
    out = paths.fig_output(FIG)
    d, ev, meta = load()
    cmap = style.density_cmap()
    norm, c = count_scale(ev, meta)
    gating = panel_a(d, ev, cmap, norm, out)
    printed = panel_b(d, ev, meta, cmap, norm, out)
    cmp = report(d, meta, printed, gating, c, norm)

    sizes = {s: emitted(out, s) for s in ("panel_a_gating", "panel_b_wells")}
    for s, (w, h) in sizes.items():
        print(f"emitted {s}: {w:.1f} x {h:.1f} mm")
    (aw, ah), (bw, bh) = sizes["panel_a_gating"], sizes["panel_b_wells"]
    assert max(aw, bw) <= BOOK_AVAIL, "a panel is wider than the A4 text block"
    assert aw + BOOK_GUTTER + bw > BOOK_AVAIL, \
        "(a) and (b) would share a book row; (b) must sit under (a)"
    page = 2 * BOOK_MARGIN + BOOK_TITLE + ah + BOOK_GUTTER + bh
    assert page <= A4_H, f"the two panels need {page:.0f} mm on a 297 mm page"
    print(f"book page: {page:.0f} mm of {A4_H:.0f}")
    cmp.to_csv(out / "onoff_two_readouts.csv")


if __name__ == "__main__":
    main()
