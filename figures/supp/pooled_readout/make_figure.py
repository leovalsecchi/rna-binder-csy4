"""pooled_readout - how deeply the 12,000-member library was actually read, and
how hard the 32 depend on where the read floor was set.

NOT BLOCKED. The pooled ONT count tables have been on the NAS since 2026-10-19 at
LV_BinderDesign/LV2.1.1_RNAbinders_PoolCraft_V1/11_alignment/counts__*.csv - eight
tables (POS/NEG x context A/B x rep 1/2), 12,000 rows each. Everything here is
computed from those primary tables, not from the derived inputs/screen_points.csv
(which is reconstructed out of the 48 MB splicing_report.html and holds only the
members with >= 1 POS read).

The pipeline is reproduced rather than trusted: load() asserts that the published
gate (reads >= 75 in POS and in the matched NEG, splicing >= 0.07, fold change >= 10)
returns exactly 7 + 25 = 32 candidates off these tables. That assert IS the check on
every derived quantity below.

Sample encoding (COLORS.md rule 2): a sequencing sample is not a design route, so it
takes no route colour. Context carries the theme hue step, POS/NEG carries the line
style, and the tick label names both - no third channel and no legend box.

scipy is unavailable in ~/.venvs/binder_map (its _propack extension fails to dlopen),
so nothing here imports it; no statistics are needed on this figure in any case.
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from common import paths, style, palette

FIG = Path(__file__).resolve().parent.name
SRC = Path("/Volumes/wicky/users/lvalsecchi/Leo_Valsecchi/LV_BinderDesign"
           "/LV2.1.1_RNAbinders_PoolCraft_V1/11_alignment")

N_MEMBERS = 12000        # the delivered library. NOT the lab protocol's 24,682 example.
READ_FLOOR = 75          # the floor the ordered set was called at
MIN_SPLICING = 0.07
MIN_FC = 10.0
REPS = (1, 2)

# Panel widths, which ARE the page layout. make_supp_book.py packs panels greedily
# into rows no wider than its text block - A4_W 210 mm less 2 x MARGIN 14 = 182 mm, NOT
# the 178 mm that file's own comment claims - in filename order, and letters them
# itself. So whether two panels share a row is decided here and nowhere else. The
# EMITTED width is what it measures, and that is wider than the authored figsize,
# because style.save crops to the ink box and pads it (+3.3 mm on (a), +5.9 mm on
# (b) at 85 mm authored). Hence (a) and (b) authored under 85 mm, so that emitted +
# 6 mm gutter + emitted still fits the row. (c) goes the other way: its generous
# fractional side margins are mostly CROPPED, so 188 mm authored emits 176 mm and
# fills the row beneath the pair rather than stopping two thirds across.
# Leo 2026-11-19: "pull B next to A and C can be below them."
# Emitted now: (a) 84.4 + 6 + (b) 85.2 = 175.6 mm, (c) 175.9 mm. Both rows clear the
# block with room to spare, so a MARGIN tweak upstream cannot reflow the page.
W_A, W_B, W_C = 80.5, 78.0, 188.0

# Pipeline order, never result order. Context -> theme step, role -> line style.
SAMPLES = ["POS_A", "NEG_A", "POS_B", "NEG_B"]
CTX_STEP = {"A": "900", "B": "500"}
ROLE_DASH = {"POS": (), "NEG": (2.2, 1.2)}


def scolor(s):
    return palette.THEME[CTX_STEP[s.split("_")[1]]]


def sdash(s):
    return ROLE_DASH[s.split("_")[0]]


# --- data ------------------------------------------------------------------
def pooled(sample):
    """Replicates 1 + 2 summed per library member, as the screen's own report does."""
    parts = [pd.read_csv(SRC / f"counts__{sample}_rep{r}.csv",
                         usecols=["transcript", "intact_count", "removed_count", "total"])
             .set_index("transcript") for r in REPS]
    m = parts[0] + parts[1]
    assert len(m) == N_MEMBERS, f"{sample}: {len(m)} members, expected {N_MEMBERS}"
    assert (m.total == m.intact_count + m.removed_count).all(), f"{sample}: total != intact+removed"
    m["splicing"] = np.where(m.total > 0, m.removed_count / m.total, np.nan)
    return m


def load():
    per_sample = {s: pooled(s) for s in SAMPLES}

    # POS against its matched NEG, one row per (context, member). The fold change is
    # the report's own rule: when no NEG read was spliced, NEG splicing = 1/NEG reads.
    frames = []
    for ctx in ("A", "B"):
        j = per_sample[f"POS_{ctx}"].join(per_sample[f"NEG_{ctx}"], rsuffix="_neg")
        neg_spl = np.where(j.removed_count_neg > 0, j.splicing_neg,
                           np.where(j.total_neg > 0, 1.0 / j.total_neg, np.nan))
        j["fold_change"] = j.splicing / neg_spl
        j["context"] = ctx
        frames.append(j)
    pairs = pd.concat(frames)

    n = gated(pairs, READ_FLOOR).context.value_counts()
    assert (n.get("A"), n.get("B")) == (7, 25), (
        f"gate at floor {READ_FLOOR} gives A={n.get('A')} B={n.get('B')}, "
        "expected 7 + 25 = 32 (gated_Q1_75Reads_OrderedHits.xlsx)")
    return per_sample, pairs


def gated(pairs, floor):
    return pairs[(pairs.total >= floor) & (pairs.total_neg >= floor)
                 & (pairs.splicing >= MIN_SPLICING) & (pairs.fold_change >= MIN_FC)]


# --- (a) reads per library member ------------------------------------------
def sample_key(fig, y=0.012):
    """Frameless key along the bottom (STYLE_NOTES section 6).

    Direct labels were tried and do not fit: the four curves converge on 1 read, and
    POS_B crosses NEG_B in the band where they are otherwise furthest apart, so a
    label on each line lands on the wrong curve or on its neighbour's label.
    """
    for i, s in enumerate(SAMPLES):
        x = 0.14 + i * 0.21
        fig.add_artist(plt.Line2D([x, x + 0.035], [y + 0.008] * 2, lw=0.8,
                                  color=scolor(s), dashes=sdash(s),
                                  transform=fig.transFigure))
        fig.text(x + 0.043, y + 0.002, s, fontsize=style.CFG["font_size_pt"],
                 va="baseline", ha="left")


def panel_b(per_sample, out):
    fig, ax = plt.subplots(figsize=(style.mm(W_B), style.mm(52)))
    for s in SAMPLES:
        v = np.sort(per_sample[s].total.values)[::-1]
        v = v[v > 0]                      # a log axis cannot carry the unseen members
        ax.plot(np.arange(1, len(v) + 1), v, lw=0.8, color=scolor(s),
                dashes=sdash(s), solid_capstyle="butt")

    ax.axhline(READ_FLOOR, color=palette.NEUTRAL["600"], lw=0.5, dashes=(1.5, 1.5))
    # Named in the plot rather than given a 75 tick: at 5.5 pt a 75 tick label
    # touches the 100 decade label 0.12 of a decade above it.
    ax.annotate(f"{READ_FLOOR}-read floor", (N_MEMBERS * 0.015, READ_FLOOR),
                textcoords="offset points", xytext=(0, 2), ha="left",
                fontsize=style.CFG["font_size_pt"], color=palette.NEUTRAL["600"])
    ax.set_yscale("log")
    ax.set_xlim(0, N_MEMBERS * 1.02)
    ax.set_ylim(1, 10 ** 5)
    ax.set_yticks([1, 10, 100, 1000, 10000, 100000])
    ax.set_yticklabels(["1", "10", "100", "1,000", "10,000", "100,000"])
    ax.set_yticks([], minor=True)
    ax.set_xticks([0, 3000, 6000, 9000, 12000])
    ax.set_xticklabels(["0", "3,000", "6,000", "9,000", "12,000"])
    ax.set_xlabel("Library members, ranked by reads")
    ax.set_ylabel("Reads")
    ax.tick_params(axis="both", which="major", length=2, width=0.5)
    fig.subplots_adjust(left=0.13, right=0.98, top=0.97, bottom=0.24)
    sample_key(fig)
    style.save(fig, out, "panel_b_reads_per_member")
    plt.close(fig)


# --- (b) library recovery --------------------------------------------------
def recovery(per_sample):
    rows = []
    for s in SAMPLES:
        t = per_sample[s].total
        rows.append(dict(sample=s, seen=int((t > 0).sum()),
                         above=int((t >= READ_FLOOR).sum()), reads=int(t.sum())))
    return pd.DataFrame(rows).set_index("sample")


def fill_key(fig, y=0.012):
    """Filled = seen at all, open = above the floor. Frameless, under the panel."""
    pairs = ((palette.NEUTRAL["700"], palette.NEUTRAL["700"], "any read"),
             (palette.SURFACE, palette.NEUTRAL["700"], f"≥ {READ_FLOOR} reads"))
    for i, (fc, ec, lab) in enumerate(pairs):
        x = 0.22 + i * 0.30
        fig.patches.append(plt.Rectangle((x, y), 0.022, 0.020, figure=fig,
                                         transform=fig.transFigure, facecolor=fc,
                                         edgecolor=ec, linewidth=0.6))
        fig.text(x + 0.030, y + 0.004, lab, fontsize=style.CFG["font_size_pt"],
                 va="baseline", ha="left")


def panel_a(per_sample, out):
    r = recovery(per_sample)
    x = np.arange(len(SAMPLES), dtype=float)
    w = 0.34

    fig, ax = plt.subplots(figsize=(style.mm(W_A), style.mm(48)))
    for i, s in enumerate(SAMPLES):
        c = scolor(s)
        ax.bar(x[i] - w / 2, r.seen[s], w, facecolor=c,
               edgecolor=palette.INK, lw=0.4)
        ax.bar(x[i] + w / 2, r.above[s], w, facecolor=palette.SURFACE,
               edgecolor=c, lw=0.6)
        # The count only: the denominator is the same 12,000 for every bar and is
        # drawn once as the reference line, so repeating it eight times is the
        # label creep STYLE_NOTES section 0 is about. Rates go in the caption.
        for dx, v in ((-w / 2, r.seen[s]), (w / 2, r.above[s])):
            ax.annotate(f"{v:,}", (x[i] + dx, v), textcoords="offset points",
                        xytext=(0, 2), ha="center",
                        fontsize=style.CFG["font_size_pt"])

    ax.axhline(N_MEMBERS, color=palette.NEUTRAL["400"], lw=0.5)
    ax.annotate(f"{N_MEMBERS:,} ordered", (x[0] - w, N_MEMBERS),
                textcoords="offset points", xytext=(0, 3), ha="left",
                fontsize=style.CFG["font_size_pt"], color=palette.NEUTRAL["600"])
    ax.set_xticks(x)
    ax.set_xticklabels(SAMPLES)
    ax.set_ylabel("Library members")
    ax.set_ylim(0, N_MEMBERS * 1.22)
    ax.set_yticks([0, 3000, 6000, 9000, 12000])
    ax.set_yticklabels(["0", "3,000", "6,000", "9,000", "12,000"])
    ax.tick_params(axis="both", length=2, width=0.5)
    fig.subplots_adjust(left=0.13, right=0.98, top=0.97, bottom=0.26)
    fill_key(fig)
    style.save(fig, out, "panel_a_library_recovery")
    plt.close(fig)


# --- (c) candidate count against the read floor ----------------------------
FLOORS = np.concatenate([np.arange(1, 100, 1), np.arange(100, 505, 5)])


def sweep(pairs):
    rows = []
    for f in FLOORS:
        g = gated(pairs, f)
        n = g.context.value_counts()
        rows.append(dict(floor=int(f), A=int(n.get("A", 0)), B=int(n.get("B", 0)),
                         total=len(g)))
    return pd.DataFrame(rows).set_index("floor")


def panel_c(pairs, out):
    s = sweep(pairs)
    fig, ax = plt.subplots(figsize=(style.mm(W_C), style.mm(48)))

    # All three series end at 0 and B tracks the total once A is exhausted, so a
    # label at the end of a line lands on the other two. Frameless key instead,
    # the same idiom as panels (a) and (b).
    series = [("total", palette.MAIN, (), "Both contexts"),
              ("A", palette.THEME["900"], (2.2, 1.2), "Context A"),
              ("B", palette.THEME["500"], (1.2, 1.2), "Context B")]
    # Leo 2026-11-20: "colour the volume under the curve the colour of the hits it
    # represents". Each curve fills its own area in its own colour.
    #
    # OVERLAP, decided rather than left to chance: total = A + B exactly, so the
    # total area CONTAINS both context areas and the two context areas overlap each
    # other below floor ~120, where both still hold candidates. A stacked fill would
    # partition the total exactly (A's share, then B's) but its upper boundary would
    # not be the drawn B line, so the band labelled B would not be the curve labelled
    # B. Drawn separately translucent instead, in a FIXED back-to-front order -
    # total, then B, then A - so the paint order never depends on the data: the three
    # fills are nested backdrops under their own lines, and the region where both
    # contexts contribute is simply the darker one. Alpha is low enough that two
    # layers (0.19) stay under the 75-read rule and its annotated point, which are
    # lifted above every fill by zorder.
    for col, c, dash, _ in series:
        ax.fill_between(s.index, s[col], color=c, alpha=0.10, lw=0, zorder=0)
    for col, c, dash, _ in series:
        ax.plot(s.index, s[col], lw=0.9, color=c, dashes=dash, zorder=3)
    for i, (_, c, dash, lab) in enumerate(series):
        x = 0.16 + i * 0.27
        fig.add_artist(plt.Line2D([x, x + 0.035], [0.020] * 2, lw=0.9, color=c,
                                  dashes=dash, transform=fig.transFigure))
        fig.text(x + 0.043, 0.014, lab, fontsize=style.CFG["font_size_pt"],
                 va="baseline", ha="left")

    n_at_floor = int(s.total[READ_FLOOR])
    ax.plot([READ_FLOOR], [n_at_floor], "o", ms=3, color=palette.MAIN, zorder=6)
    ax.annotate(f"{n_at_floor} ordered\nat {READ_FLOOR} reads",
                (READ_FLOOR, n_at_floor), textcoords="offset points", xytext=(5, 6),
                ha="left", fontsize=style.CFG["font_size_pt"], color=palette.INK,
                zorder=6)
    # Above the fills, not merely after them: a fill drawn later would otherwise
    # wash the floor rule out, which is the one mark the panel is read against.
    ax.axvline(READ_FLOOR, color=palette.NEUTRAL["600"], lw=0.5, dashes=(1.5, 1.5),
               zorder=5)

    ax.set_xlim(0, FLOORS[-1] * 1.02)
    ax.set_ylim(0, s.total.max() * 1.12)
    # No 75 tick: it sits 25 units from 100 and the two labels touch at 5.5 pt.
    # The dashed line and its annotation carry the floor instead.
    ax.set_xticks([0, 100, 200, 300, 400, 500])
    ax.set_xlabel("Read floor, POS and matched NEG (reads)")
    ax.set_ylabel("Candidates")
    ax.tick_params(axis="both", length=2, width=0.5)
    fig.subplots_adjust(left=0.11, right=0.98, top=0.97, bottom=0.30)
    style.save(fig, out, "panel_c_candidates_vs_read_floor")
    plt.close(fig)
    return s


# --- the numbers the caption and README carry ------------------------------
def report(per_sample, pairs, s):
    r = recovery(per_sample)
    print("recovery, out of", f"{N_MEMBERS:,}", "ordered members")
    for smp in SAMPLES:
        print(f"  {smp:6s} any read {r.seen[smp]:6,}/{N_MEMBERS:,} "
              f"({100*r.seen[smp]/N_MEMBERS:5.1f}%)   >= {READ_FLOOR} "
              f"{r.above[smp]:6,}/{N_MEMBERS:,} ({100*r.above[smp]/N_MEMBERS:5.1f}%)"
              f"   {r.reads[smp]:>12,} reads assigned")
    print(f"  assigned reads, all four samples: {r.reads.sum():,}")
    print(f"  100x coverage target at {N_MEMBERS:,} members = "
          f"{100*N_MEMBERS:,} reads per sample; "
          f"met by {(r.reads >= 100*N_MEMBERS).sum()}/4 samples")

    print("candidates against read floor")
    for f in (1, 25, 50, 75, 100, 150, 200, 300, 500):
        print(f"  floor {f:3d}: A {s.A[f]:2d}  B {s.B[f]:2d}  total {s.total[f]:2d}")

    # Monotone, and nested: raising the floor only ever removes a candidate.
    base = set(gated(pairs, READ_FLOOR).index)
    assert s.total.is_monotonic_decreasing, "candidate count is not monotone in the floor"
    for f in (100, 150, 200, 300):
        new = set(gated(pairs, f).index) - base
        assert not new, f"floor {f} ADDS {len(new)} candidates absent at {READ_FLOOR}"
        print(f"  of the {len(base)} ordered, {len(base & set(gated(pairs, f).index))} "
              f"still pass at floor {f}")
    print("  no floor above 75 admits a candidate the 75-read set does not contain")


def main():
    style.apply()
    out = paths.fig_output(FIG)
    per_sample, pairs = load()
    # Leo 2026-11-19: recovery is now (a) and reads-per-member (b); the
    # assembler letters panels by filename order, so the stems carry it.
    panel_a(per_sample, out)
    panel_b(per_sample, out)
    s = panel_c(pairs, out)
    report(per_sample, pairs, s)


if __name__ == "__main__":
    main()
