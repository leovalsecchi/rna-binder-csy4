"""Fig03a: horizontal filter funnel, filters in the order RMSD, DockQ, AA run,
ipTM, foot clash. Flat top; ribbons show each stack's share of the survivors
at every stage. Two variants: depth on a log count scale, and linear.

Data: inputs/funnel_real_gate_v2.csv, from funnel_order_v2.py (same gate as
fig07, reordered; reproduces 55,706).
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
plt.rcParams.update({"font.family": _cfg["font_family"], "pdf.fonttype": 42,
                     "svg.fonttype": "none", "savefig.dpi": 600})
def mm(x): return x / 25.4

INK, MUTED = palette.INK, palette.NEUTRAL["500"]
STACK_COL = {palette.SHORT[k]: palette.STACKS[k] for k in palette.STACK_ORDER}
ORDER = list(STACK_COL)
FLOOR = 2e4          # count at zero depth; depth = log10(n / FLOOR)

d = pd.read_csv(paths.INPUTS / "funnel_real_gate_v2.csv")
N = d["total"].to_numpy(float)
assert int(N[-1]) == 55706, "gate no longer reproduces the passing set"
xk = np.arange(len(N), dtype=float)


def _pchip(xk, yk, xs):
    """Monotone cubic Hermite (Fritsch-Carlson); scipy is broken in the venv."""
    h = np.diff(xk); delta = np.diff(yk) / h
    m = np.zeros_like(yk)
    with np.errstate(divide="ignore", invalid="ignore"):
        m[1:-1] = np.where(delta[:-1] * delta[1:] > 0,
                           2.0 / (1.0 / delta[:-1] + 1.0 / delta[1:]), 0.0)
    m[0], m[-1] = delta[0], delta[-1]
    i = np.clip(np.searchsorted(xk, xs) - 1, 0, len(h) - 1)
    t = (xs - xk[i]) / h[i]; t2, t3 = t * t, t * t * t
    return ((2*t3 - 3*t2 + 1) * yk[i] + (t3 - 2*t2 + t) * h[i] * m[i]
            + (-2*t3 + 3*t2) * yk[i+1] + (t3 - t2) * h[i] * m[i+1])


FR = d[ORDER].to_numpy(float); FR /= FR.sum(1, keepdims=True)   # stack shares per stage
xs = np.linspace(0, xk[-1], 800)
_i = np.clip(np.searchsorted(xk, xs) - 1, 0, len(xk) - 2)
CUMK = np.hstack([np.zeros((len(N), 1)), np.cumsum(FR, 1)])           # band edges at stages


def draw(scale, stem):
    if scale == "log":
        dep = np.log10(N / FLOOR)
        ys = -_pchip(xk, dep, xs)
        ticks = [(1e5, "$10^5$"), (1e6, "$10^6$")]
        tick_y = lambda v: -np.log10(v / FLOOR)
        ylab, LH = "designs (log)", 0.10
    else:
        dep = N / 1e6                       # depth in millions, linear
        ys = -np.interp(xs, xk, dep)        # straight segments: PCHIP overshoots on this drop
        ticks = [(v * 1e6, f"{v}M") for v in (2, 4, 6)]
        tick_y = lambda v: -v / 1e6
        ylab, LH = "designs", 0.12
    # Band edges between stages follow the outline's own progress s in [0, 1],
    # so every edge is monotone wherever the outline is (no swirls). Multiplying
    # the outline by separately eased shares was not.
    yk = -dep
    s_ = ((ys - yk[_i]) / (yk[_i + 1] - yk[_i]))[:, None]
    BK = yk[:, None] * CUMK
    B = BK[_i] + (BK[_i + 1] - BK[_i]) * s_
    fig, ax = plt.subplots(figsize=(mm(180), mm(72)))

    for k, s_ in enumerate(ORDER):          # stack ribbons, share of local depth
        ax.fill_between(xs, B[:, k], B[:, k + 1],
                        color=STACK_COL[s_], lw=0, zorder=1)
        if k:
            ax.plot(xs, B[:, k], color="white", lw=0.4, zorder=2)
    ax.plot(xs, ys, color=INK, lw=1.0, zorder=3)
    ax.plot([0, xk[-1]], [0, 0], color=INK, lw=1.0, zorder=3)
    ax.plot([0, 0], [0, ys[0]], color=INK, lw=1.0, zorder=3)

    for i in range(1, len(N)):              # filter boundaries and names
        ax.plot([xk[i]] * 2, [0, -dep[i]], color="white", lw=0.9, zorder=2)
        ax.text((xk[i - 1] + xk[i]) / 2, LH, d["step"].iloc[i].replace("Longest AA", "AA")
                .replace(" (max)", ""), fontsize=7, color=INK, ha="center", va="bottom")
    for i in range(len(N)):                 # survivors just below the curve
        ax.text(xk[i] + 0.05, -dep[i] - LH, f"{int(N[i]):,}", fontsize=7, color=INK,
                fontweight="semibold", ha="left", va="top")

    for v, lab in ticks:                    # scale bar, left
        y = tick_y(v)
        ax.plot([-0.12, -0.06], [y, y], color=MUTED, lw=0.6)
        ax.text(-0.15, y, lab, fontsize=6.5, color=MUTED, ha="right", va="center")
    ax.plot([-0.09, -0.09], [0, -dep[0]], color=MUTED, lw=0.6)
    ax.text(-0.45, -dep[0] / 2, ylab, fontsize=6.5, color=MUTED, rotation=90,
            ha="center", va="center")

    # stack legend with passing-set shares, at the spout
    H = dep[0] * 0.6
    yy = np.linspace(-0.1 * dep[0], -0.1 * dep[0] - H, len(ORDER))
    for k, s_ in enumerate(ORDER):
        ax.add_patch(plt.Rectangle((xk[-1] + 0.45, yy[k] - 0.06 * dep[0]), 0.12,
                                   0.12 * dep[0], facecolor=STACK_COL[s_], lw=0))
        ax.text(xk[-1] + 0.63, yy[k], f"{s_}  {100*FR[-1, k]:.1f}%", fontsize=6.5,
                color=INK, ha="left", va="center")
    ax.text(xk[-1] + 0.45, LH, "Passing set mix", fontsize=7, color=INK,
            ha="left", va="bottom", fontweight="semibold")

    ax.set_xlim(-0.6, xk[-1] + 1.6)
    ax.set_ylim(-dep[0] * 1.12, dep[0] * 0.14)
    ax.axis("off")
    for ext in ("svg", "pdf", "png"):
        fig.savefig(out / f"{stem}.{ext}", bbox_inches="tight")
    plt.close(fig)
    print("wrote", stem)


def draw_lanes(stem, floor=300):
    """One log-depth lane per stack, same scale for all.

    Stacking shares inside a log outline makes the shares look skewed (a band's
    thickness is share x log(total), so RFD3 seems to lose almost nothing).
    Here each stack's own count gets its own log depth, so attrition per stack
    is honest and directly comparable between lanes.
    """
    C = d[ORDER].to_numpy(float)
    D = np.log10(C / floor)
    gap = 0.45
    fig, ax = plt.subplots(figsize=(mm(180), mm(110)))
    top = 0.0
    for k, s_ in enumerate(ORDER):
        ys = top - _pchip(xk, D[:, k], xs)
        ax.fill_between(xs, top, ys, color=STACK_COL[s_], lw=0, zorder=1)
        ax.plot(xs, ys, color=INK, lw=0.7, zorder=3)
        for i in range(1, len(N)):
            ax.plot([xk[i]] * 2, [top, top - D[i, k]], color="white", lw=0.8, zorder=2)
        ax.text(-0.08, top - D[0, k] / 2, s_, fontsize=7, color=STACK_COL[s_],
                fontweight="semibold", ha="right", va="center")
        ax.text(0.05, top - D[0, k] - 0.06, f"{int(C[0, k]):,}", fontsize=6.5,
                color=INK, ha="left", va="top")
        ax.text(xk[-1] + 0.08, top - D[-1, k] / 2,
                f"{int(C[-1, k]):,}  ({100*C[-1, k]/C[0, k]:.2f}%)", fontsize=6.5,
                color=INK, ha="left", va="center")
        top -= D[0].max() + gap
    for i in range(1, len(N)):
        ax.text((xk[i - 1] + xk[i]) / 2, 0.12, d["step"].iloc[i].replace("Longest AA", "AA")
                .replace(" (max)", ""), fontsize=7, color=INK, ha="center", va="bottom")
    # one shared decade bar, top right of the first lane
    ax.plot([xk[-1] + 1.15] * 2, [-0.2, -1.2], color=MUTED, lw=0.8)
    ax.text(xk[-1] + 1.2, -0.7, "10x", fontsize=6.5, color=MUTED, ha="left", va="center")
    ax.set_xlim(-1.1, xk[-1] + 1.5)
    ax.set_ylim(top + gap - 0.4, 0.6)
    ax.axis("off")
    for ext in ("svg", "pdf", "png"):
        fig.savefig(out / f"{stem}.{ext}", bbox_inches="tight")
    plt.close(fig)
    print("wrote", stem)


draw("log", "filter_funnel_horizontal")
draw_lanes("filter_funnel_horizontal_lanes")


def draw_survival(stem):
    """Per-stack survival, % of that stack's own designs, on a log axis.

    Log depth from a floor makes every stack look like it keeps most of its
    depth; on a log axis with decade gridlines, loss reads as decades dropped.
    """
    C = d[ORDER].to_numpy(float)
    S = 100 * C / C[0]
    fig, ax = plt.subplots(figsize=(mm(150), mm(80)))
    for k, s_ in enumerate(ORDER):
        ax.plot(xk, S[:, k], color=STACK_COL[s_], lw=1.6, marker="o", ms=3.2,
                zorder=3)
    ends = list(S[-1])                       # direct labels, nudged apart
    order = np.argsort(ends)[::-1]
    ly = {}
    prev = None
    for k in order:
        y = ends[k] if prev is None else min(ends[k], prev / 1.45)
        ly[k] = prev = y
    for k, s_ in enumerate(ORDER):
        ax.text(xk[-1] + 0.12, ly[k], f"{s_}  {S[-1, k]:.2f}%  ({int(C[-1, k]):,})",
                fontsize=6.5, color=STACK_COL[s_], ha="left", va="center")
    ax.set_yscale("log")
    ax.set_ylim(0.05, 130)
    ax.set_yticks([100, 10, 1, 0.1])
    ax.set_yticklabels(["100%", "10%", "1%", "0.1%"])
    ax.grid(axis="y", color=palette.GRID, lw=0.6, zorder=0)
    ax.set_xticks(xk)
    ax.set_xticklabels(["All"] + [t.replace("Longest AA", "AA").replace(" (max)", "")
                                  for t in d["step"].iloc[1:]], fontsize=7)
    ax.set_xlim(-0.2, xk[-1] + 0.1)
    ax.set_ylabel("Designs surviving, % of route", fontsize=7)
    ax.tick_params(labelsize=7, length=2.5)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    fig.tight_layout()
    for ext in ("svg", "pdf", "png"):
        fig.savefig(out / f"{stem}.{ext}", bbox_inches="tight")
    plt.close(fig)
    print("wrote", stem)


draw_survival("stack_survival_log")


def draw_combo(stem, ratio=2.0, height_mm=46, box_mm=0):
    """Top: total funnel on a log axis. Bottom: a thin, separate strip with the
    stack mix at every stage (shares, 100%); shares labelled at start and end.

    Both share x; each knot is the pool AFTER that filter, so the filter names
    sit once, under the strip.
    """
    dep = np.log10(N / FLOOR)
    # Three rows: an empty grey box at the very top to paste the ChimeraX
    # renders into, then the funnel, then the stack strip. The box sits ABOVE
    # the filter labels, which are pinned to the top of the funnel axes, so the
    # gap between row 0 and row 1 has to clear them.
    funnel_mm = height_mm * ratio / (ratio + 1.0)
    strip_mm = height_mm / (ratio + 1.0)
    if box_mm:
        fig, (abox, af, ab) = plt.subplots(
            3, 1, figsize=(mm(160), mm(height_mm + box_mm)), sharex=True,
            gridspec_kw={"height_ratios": [box_mm, funnel_mm, strip_mm],
                         "hspace": 0.30})
        # abox only reserves space; the grey is drawn later as figure-level
        # rectangles so the columns run continuously down past the strip
        # instead of stopping at each axes boundary.
        abox.set_xticks([]); abox.set_yticks([])
        abox.patch.set_alpha(0.0)
        for side in ("top", "right", "bottom", "left"):
            abox.spines[side].set_visible(False)
    else:
        fig, (af, ab) = plt.subplots(
            2, 1, figsize=(mm(160), mm(height_mm)), sharex=True,
            gridspec_kw={"height_ratios": [funnel_mm, strip_mm], "hspace": 0.16})
    # funnel: log axis, light fill, hairline outline, faint stage guides
    # y increases upward: more designs = higher, so the funnel descends left to
    # right and 10^6 sits above 10^5, the way a count axis normally reads.
    ys = _pchip(xk, dep, xs)
    af.fill_between(xs, 0, ys, color=palette.THEME["200"], lw=0, zorder=1)
    af.plot(xs, ys, color=palette.THEME["700"], lw=0.9, zorder=3)
    af.plot([0, xk[-1]], [0, 0], color=palette.THEME["700"], lw=0.9, zorder=3)
    for i in range(len(N)):
        if i == 0:                      # tuck under the curve so it clears the
            x, y = xk[0] + 0.16, dep[0] - 0.08      # first filter's title
            ha, va = "left", "top"
        else:
            x, y = xk[i], dep[i] + 0.10
            ha, va = "center", "bottom"
        af.text(x, y, f"{int(N[i]):,}", fontsize=6.5, color=INK, ha=ha, va=va)
    maj = [1e5, 1e6]
    mnr = [m * 10**e for e in (4, 5, 6) for m in range(2, 10)
           if FLOOR < m * 10**e <= N[0]]
    af.set_yticks([np.log10(v / FLOOR) for v in maj])
    # Plain comma-grouped numbers, not mathtext: with svg.fonttype="none"
    # mathtext emits its own DejaVu Sans font reference, which any viewer or
    # Illustrator without that font substitutes badly. Arial has no superscript
    # digits either, so the exponent form is out. These match the count labels.
    af.set_yticklabels([f"{int(v):,}" for v in maj])
    af.set_yticks([np.log10(v / FLOOR) for v in mnr], minor=True)
    af.set_ylim(0.0, dep[0] + 0.42)
    af.set_ylabel("Designs remaining", fontsize=7)
    af.tick_params(labelsize=7, length=2.5, bottom=False)
    af.tick_params(which="minor", length=1.5)
    for side in ("top", "right", "bottom"):
        af.spines[side].set_visible(False)
    if not box_mm:                      # the grey columns already separate stages
        for x in xk[1:]:
            af.axvline(x, color=palette.GRID, lw=0.6, zorder=0)

    # stack strip: shares at each stage, linear between stages (monotone edges)
    cum = np.vstack([np.interp(xs, xk, CUMK[:, k]) for k in range(len(ORDER) + 1)]).T
    for k, s_ in enumerate(ORDER):
        ab.fill_between(xs, 1 - cum[:, k], 1 - cum[:, k + 1], color=STACK_COL[s_],
                        lw=0, zorder=1)
    if not box_mm:
        for x in xk[1:-1]:
            ab.axvline(x, color="white", lw=0.6, zorder=2)
    # Labels are spread evenly over the strip height rather than tracking each
    # band's centre: the bands differ by more than an order of magnitude
    # (RFD3 50.1% vs HD 1.4%), so centre-tracking crowds the thin ones together.
    # Order still matches the stacking order top to bottom.
    n_lab = len(ORDER)
    even = [1 - (k + 0.5) / n_lab for k in range(n_lab)]
    for side, j, ha, dx in (("start", 0, "right", -0.06), ("end", -1, "left", 0.06)):
        cen = even
        for k, s_ in enumerate(ORDER):
            lab = f"{100*FR[j, k]:.1f}%"
            if side == "end":
                lab = f"{s_}  {lab}"
            ab.text(xk[j] + dx, cen[k], lab, fontsize=6, color=STACK_COL[s_],
                    ha=ha, va="center", clip_on=False)
    ab.set_ylim(0, 1)
    ab.set_yticks([])
    for side in ("top", "right", "left", "bottom"):
        ab.spines[side].set_visible(False)
    ab.tick_params(axis="x", bottom=False, labelbottom=False)
    # Filter names sit at the MIDPOINT of the interval they act on: the filter
    # between knot i-1 and knot i is what takes the pool from one to the other.
    # "All designs" is deliberately NOT drawn here: at a midpoint layout it
    # overlaps the first gate, and the starting pool is already labelled by the
    # 6,992,515 count at xk[0].
    for i in range(1, len(N)):
        lab = d["step"].iloc[i].replace("Longest AA", "AA").replace(" (max)", "")
        af.text((xk[i - 1] + xk[i]) / 2, 1.02, lab, fontsize=7, color=INK,
                ha="center", va="bottom", clip_on=False,
                transform=af.get_xaxis_transform())
    ab.set_xlim(-0.05, xk[-1] + 0.05)

    if box_mm:
        # One grey column per filter, running from the top of the render box
        # down to just past the stack strip. Drawn in figure coordinates and at
        # zorder 0 so it passes behind the axes and across the gaps between them.
        from matplotlib.patches import Rectangle
        inv = fig.transFigure.inverted()
        fx = lambda xd: inv.transform(af.transData.transform((xd, 0)))[0]
        # the funnel and strip axes default to an opaque white patch, which
        # would hide the columns behind them
        af.patch.set_alpha(0.0)
        ab.patch.set_alpha(0.0)
        y_top = abox.get_position().y1
        y_bot = ab.get_position().y0                  # flush with the strip's bottom
        GAP = 0.010                                    # data units, white seam
        for i in range(1, len(N)):
            x0, x1 = fx(xk[i - 1] + GAP), fx(xk[i] - GAP)
            fig.add_artist(Rectangle((x0, y_bot), x1 - x0, y_top - y_bot,
                                     transform=fig.transFigure,
                                     facecolor=palette.NEUTRAL["100"],
                                     edgecolor="none", zorder=0))
        # The stack strip is filled, so the seams between columns have to be
        # repainted on top of it (figure artist above the axes, zorder 3).
        for i in range(1, len(N) - 1):
            x0, x1 = fx(xk[i] - GAP), fx(xk[i] + GAP)
            fig.add_artist(Rectangle((x0, ab.get_position().y0), x1 - x0,
                                     ab.get_position().height,
                                     transform=fig.transFigure,
                                     facecolor="white", edgecolor="none", zorder=3))
    for ext in ("svg", "pdf", "png"):
        fig.savefig(out / f"{stem}.{ext}", bbox_inches="tight")
    plt.close(fig)
    print("wrote", stem)


draw_combo("Fig3a_funnel", 2.0, 46)
draw_combo("Fig3a_funnel_3to1", 3.0, 60)
draw_combo("Fig3a_funnel_3to1_box", 3.0, 60, box_mm=30)
draw("linear", "filter_funnel_horizontal_linear")
