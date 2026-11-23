"""Fig3b: four-set schematic Venn of the filter criteria (replaces the funnel).

Sets: ligand_RMSD, sc_DockQ, AA run, ipTM. Foot clash is NOT a set; it is a
separate strip below (80,983 all-four-pass -> 25,277 removed -> 55,706).
Areas are SCHEMATIC (regions span 1,646 to 2,901,047). The four-ellipse
layout is the standard one (two pairs of rotated ellipses, as in
matplotlib-venn), which has all 16 regions; it is stretched wider than tall
(affine, topology unchanged) to keep the panel short.
Data: inputs/venn_regions.csv. Renders: selection/Fig3a_box_*.png.
"""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.image as mimg
from matplotlib.patches import Polygon, Rectangle
from common import paths, palette, style

FIG = Path(__file__).resolve().parent.name
out = paths.fig_output(FIG)
style.apply()
INK, MUTED = palette.INK, palette.NEUTRAL["600"]

# 2026-11-19 (Leo): "add colors to the venn diagramm: (different from the ones
# we have for design method and not only blue). Also show me a version where
# the tiles are darker the more binders they include."
#
#   criteria  one hue per CRITERION, fills and outlines and the criterion's own
#             name; overlaps blend, which is what a Venn's overlaps are for
#   count     every one of the 15 regions shaded by how many designs are in it,
#             on a log ramp, outlines left neutral so only the shading speaks
#
# The four hues are petrol, teal, green and gold, all from config/palette.json
# and deliberately NOT the stack colours (rust, two purples, gold, pink): the
# criteria are not design methods and must not read as them. Gold is the one
# that sits anywhere near a stack colour (BD1's #bd8703); it is kept because
# the panel carries no stack mark at all, so nothing here can be mistaken for
# one, and four distinct hues is the point.
MODE = (sys.argv[1] if len(sys.argv) > 1 else "criteria")
assert MODE in ("criteria", "count"), MODE
STEM = "Fig3b_filter_venn" + ("" if MODE == "criteria" else "_count")
CRIT_C = [palette.THEME["600"], palette.SYNBIO["enzyme"],
          palette.SYNBIO["reporter"], palette.SYNBIO["enzyme_alt"]]
FS = style.CFG["font_size_pt"]
PT = 25.4 / 72                                  # mm per point
W = style.CFG["width_double_mm"]
H = 60.0

# ---- data + checks ---------------------------------------------------------
d = pd.read_csv(paths.INPUTS / "venn_regions.csv")
KEYS = ["sc_rmsd", "dockq", "aa_run", "iptm"]
N = {tuple(int(r[k]) for k in KEYS): int(r["n"]) for _, r in d.iterrows()}
assert len(N) == 16 and sum(N.values()) == 6_992_515, sum(N.values())
ALL4 = N[(1, 1, 1, 1)]
assert ALL4 == 80_983
FOOT_REMOVED, FINAL_N = 25_277, 55_706
assert ALL4 - FOOT_REMOVED == FINAL_N
for k in sorted(N): print(k, f"{N[k]:,}")
print("sum", f"{sum(N.values()):,}", " all-four", f"{ALL4:,}")

# ---- ellipses (matplotlib-venn 4-set layout, unit coords, then stretched) ---
SX, SY = 112.0, 80.0                            # mm per layout unit, x and y
CEN = [(.350, .400), (.450, .500), (.544, .500), (.644, .400)]
ANG = [140, 140, 40, 40]
EW, EH = .72, .45
NAMES = ["ligand_RMSD", "sc_DockQ", "AA run", "ipTM"]
THRESH = ["\u2264 3 \u00c5", "\u2265 0.7", "\u2264 5", "\u2265 0.6"]   # <= 3 A, >= 0.7, <= 5, >= 0.6 (verified vs filter JSON)
t = np.linspace(0, 2 * np.pi, 720)


def unit_edge(i):
    a = np.deg2rad(ANG[i]); ex, ey = EW / 2 * np.cos(t), EH / 2 * np.sin(t)
    return np.c_[CEN[i][0] + ex * np.cos(a) - ey * np.sin(a),
                 CEN[i][1] + ex * np.sin(a) + ey * np.cos(a)]


U = [unit_edge(i) for i in range(4)]
allu = np.vstack(U)
VX0, VY0 = 74.0 - SX * (allu[:, 0].min() + allu[:, 0].max()) / 2, 3.0 - SY * allu[:, 1].min()
to_mm = lambda p: np.c_[VX0 + SX * p[:, 0], VY0 + SY * p[:, 1]]
edge_l = [to_mm(u) for u in U]                  # outlines in mm


def inside(i, x, y):
    ux, uy = (x - VX0) / SX - CEN[i][0], (y - VY0) / SY - CEN[i][1]
    a = np.deg2rad(ANG[i])
    u = ux * np.cos(a) + uy * np.sin(a); v = -ux * np.sin(a) + uy * np.cos(a)
    return (u / (EW / 2)) ** 2 + (v / (EH / 2)) ** 2


# region label anchors: point in region farthest from every ellipse edge
gx, gy = np.meshgrid(np.linspace(0, W, 900), np.linspace(0, H, 400))
mem = np.stack([inside(i, gx, gy) < 1 for i in range(4)], -1)
edge = np.vstack(edge_l)
anchor, room = {}, {}
for key in N:
    if key == (0, 0, 0, 0):
        continue
    m = np.all(mem == np.array(key, bool), -1)
    assert m.any(), f"region {key} empty: layout lost a region"
    px, py = gx[m], gy[m]
    # clearance of the actual text box (w x h mm) from every outline: >1 = clear.
    # A circle-distance criterion wastes the long axis of thin lens regions.
    bw = style.text_w(f"{N[key]:,}", FS) * PT / 2 + 0.3
    bh = FS * PT * 0.72 / 2 + 0.3
    best, bd = None, -1
    for j in range(0, len(px), 3):
        sc_ = np.maximum(np.abs(edge[:, 0] - px[j]) / bw, np.abs(edge[:, 1] - py[j]) / bh).min()
        if sc_ > bd:
            best, bd = (px[j], py[j]), sc_
    anchor[key], room[key] = best, bd

# ---- canvas ----------------------------------------------------------------
fig = plt.figure(figsize=(style.mm(W), style.mm(H)))
ax = fig.add_axes([0, 0, 1, 1])
ax.set_xlim(0, W); ax.set_ylim(0, H); ax.axis("off")

if MODE == "criteria":
    for e, c in zip(edge_l, CRIT_C):
        ax.add_patch(Polygon(e, facecolor=c, alpha=0.22, lw=0, zorder=1))
    for e, c in zip(edge_l, CRIT_C):
        ax.add_patch(Polygon(e, fill=False, edgecolor=c, lw=0.9, zorder=3))
else:
    # Shade the 15 regions by their count. Drawn as a RASTER rather than as
    # boolean polygons: the membership test is already a grid, so a finer grid
    # of the same test gives every region's area for free and no polygon
    # library is needed. 600 dpi on the save keeps it crisp in print.
    from matplotlib.colors import LogNorm as _LN
    RES = 2400
    rx, ry = np.meshgrid(np.linspace(0, W, RES), np.linspace(0, H, int(RES * H / W)))
    rmem = np.stack([inside(i, rx, ry) < 1 for i in range(4)], -1)
    code = (rmem * np.array([8, 4, 2, 1])).sum(-1)                  # 0..15
    look = np.full(16, np.nan)
    for key, v in N.items():
        if key != (0, 0, 0, 0):
            look[int(key[0]) * 8 + int(key[1]) * 4 + int(key[2]) * 2 + int(key[3])] = v
    img = np.ma.masked_invalid(look[code])
    cm_ = plt.get_cmap("YlGnBu")
    cm_.set_bad("none")
    ax.imshow(img, origin="lower", extent=(0, W, 0, H), cmap=cm_,
              norm=_LN(max(1, int(np.nanmin(look))), np.nanmax(look)),
              interpolation="nearest", zorder=1, rasterized=True)
    for e in edge_l:
        ax.add_patch(Polygon(e, fill=False, edgecolor=palette.NEUTRAL["700"],
                             lw=0.7, zorder=3))

    def _tile_ink(key):
        """Ink or white, by the tile's own luminance: the two ~2.9M regions
        fill near-navy and swallowed a dark label."""
        rgb = cm_(_LN(max(1, int(np.nanmin(look))), np.nanmax(look))(N[key]))[:3]
        f = [(c / 12.92 if c <= .04045 else ((c + .055) / 1.055) ** 2.4) for c in rgb]
        return "white" if .2126 * f[0] + .7152 * f[1] + .0722 * f[2] < .42 else INK
for key, (x, y) in anchor.items():
    big = key == (1, 1, 1, 1)
    s = f"{N[key]:,}"
    print(f"  {key} {s:>10}  clearance x{room[key]:.2f}"
          + ("   <-- TIGHT" if room[key] < 1.15 else ""))
    ax.text(x, y, s, fontsize=FS + (1 if big else 0),
            color=(_tile_ink(key) if MODE == "count" else INK), ha="center",
            va="center", fontweight="bold" if big else "normal", zorder=5)
# none-of-four: outside every set, beside the lower-left of the Venn
ax.text(2, 7.0, f"{N[(0, 0, 0, 0)]:,}", fontsize=FS, color=INK, ha="left", va="center")
ax.text(2, 3.8, "none of the four", fontsize=FS, color=MUTED, ha="left", va="center")


# ---- renders with brackets -------------------------------------------------
def load(stem):
    im = mimg.imread(paths.ROOT / "selection" / f"Fig3a_box_{stem}.png")
    ys, xs = np.where(im[..., 3] > 0.05)
    pad = 20
    c = im[max(ys.min() - pad, 0):ys.max() + pad, max(xs.min() - pad, 0):xs.max() + pad]
    print(f"  render {stem}: source {im.shape[1]}x{im.shape[0]} px, cropped {c.shape[1]}x{c.shape[0]} px")
    return c                                     # full resolution: no decimation


def put(stem, xc, yc, w):
    im = load(stem)
    h = w * im.shape[0] / im.shape[1]
    if h > w:                       # fit a w x w box, keep aspect
        w, h = w * w / h, w
    ax.imshow(im, extent=(xc - w / 2, xc + w / 2, yc - h / 2, yc + h / 2),
              zorder=2, interpolation="lanczos", interpolation_stage="rgba")
    ax.set_xlim(0, W); ax.set_ylim(0, H); ax.set_aspect("equal")


def tangent(i, x, y):
    """Unit tangent (mm space) of ellipse i at its outline point (x, y), from the
    parametrisation P(t) = C + R(a) (EW/2 cos t, EH/2 sin t), stretched by (SX, SY)."""
    a = np.deg2rad(ANG[i])
    ux, uy = (x - VX0) / SX - CEN[i][0], (y - VY0) / SY - CEN[i][1]
    t0 = np.arctan2((-ux * np.sin(a) + uy * np.cos(a)) / (EH / 2),
                    (ux * np.cos(a) + uy * np.sin(a)) / (EW / 2))
    eu, ev = -EW / 2 * np.sin(t0), EH / 2 * np.cos(t0)      # dP/dt in the ellipse frame
    d = np.array([SX * (eu * np.cos(a) - ev * np.sin(a)), SY * (eu * np.sin(a) + ev * np.cos(a))])
    return d / np.hypot(*d)


def attach(i, y_want, side, x_far):
    """Rim point of ellipse i on `side` (-1 left, +1 right) at the valid y nearest
    y_want, valid = the horizontal ray out to x_far crosses no other ellipse."""
    e = edge_l[i]
    for dy in np.arange(0, 14, 0.25):
        for y in (y_want - dy, y_want + dy):
            xs_ = []
            for a, b in zip(e, np.roll(e, -1, 0)):
                if (a[1] - y) * (b[1] - y) < 0:
                    xs_.append(a[0] + (y - a[1]) * (b[0] - a[0]) / (b[1] - a[1]))
            if not xs_:
                continue
            x = min(xs_) if side < 0 else max(xs_)
            ray = np.linspace(x, x_far, 60)
            if all((inside(k, ray[2:], np.full(58, y)) >= 1).all() for k in range(4) if k != i):
                return x, y
    raise AssertionError(f"no clear horizontal leader for set {i}")


RW = 15.0
STANDOFF, BARH = 1.4, 1.6      # mm clear of the outline (perpendicular); bar half-length
XL, XR = 9.0, 138.0
GAP = 0.8
SPEC = [(0, XL, 24.0, "scrmsd", -1), (1, XL, 52.0, "dockq", -1),
        (2, XR, 52.0, "aarun", 1), (3, XR, 24.0, "iptm", 1)]
for i, xc, y_want, stem, side in SPEC:
    x_far = xc + (RW / 2 + GAP) * (1 if side < 0 else -1)   # render edge facing the Venn
    x, y = attach(i, y_want, side, x_far)
    put(stem, xc, y, RW)
    # bracket: bar parallel to the rim, standing STANDOFF mm clear of the outline measured
    # PERPENDICULAR to it (nearest-point distance), with the horizontal leader running
    # from the bar to the render. A horizontal offset still lands on a sloped rim.
    allp = np.vstack(edge_l)
    xb = x
    while np.hypot(*(e_i := edge_l[i] - np.array([xb, y])).T).min() < STANDOFF:
        xb += side * 0.05
    tg = tangent(i, *edge_l[i][np.hypot(*e_i.T).argmin()])
    bar = np.array([[xb - BARH * tg[0], y - BARH * tg[1]], [xb + BARH * tg[0], y + BARH * tg[1]]])
    pts = bar[0] + (bar[1] - bar[0]) * np.linspace(0, 1, 40)[:, None]
    clr = min(np.hypot(*(allp - q).T).min() for q in pts)
    assert clr >= STANDOFF - 0.1, f"{NAMES[i]}: bar within {clr:.2f} mm of an outline"
    ax.plot([xb, x_far], [y, y], color=INK, lw=0.6, zorder=4)          # horizontal leader
    ax.plot(bar[:, 0], bar[:, 1], color=INK, lw=0.6, solid_capstyle="butt", zorder=4)
    print(f"  {NAMES[i]}: bar clearance {clr:.2f} mm, bar x {xb:.2f} (rim x {x:.2f})")
    ax.text(xc, y - RW / 2 - 0.8, NAMES[i], fontsize=FS + 1,
            color=(CRIT_C[i] if MODE == "criteria" else INK), ha="center",
            va="top", fontweight="semibold", zorder=5)
    # The set NAME is centred on xc and is wider than its threshold line, so it -
    # not THRESH - is what reaches the canvas edge. `ligand_RMSD` (2026-11-21,
    # was `sc_RMSD`) made it the leftmost ink in the whole figure, 1.8 mm clear.
    _nw = style.text_w(NAMES[i], FS + 1) * PT / 2
    assert 0 <= xc - _nw and xc + _nw <= W, f"{NAMES[i]}: set name runs off the canvas"
    # gate threshold, same size as the name, wrapped under it (gate terms: filter JSON in
    # passing_v6_final.csv; ligand_RMSD/DockQ pass on EITHER backend, ipTM on max(AF3, ESM) - caption)
    ax.text(xc, y - RW / 2 - 0.8 - (FS + 1) * PT * 1.25, THRESH[i], fontsize=FS + 1, color=INK,
            ha="center", va="top", zorder=5)
    assert style.text_w(THRESH[i], FS + 1) * PT / 2 <= xc if side < 0 else xc + style.text_w(THRESH[i], FS + 1) * PT / 2 <= W
    print(f"  {NAMES[i]}: leader y {y:.1f} mm, from x {x:.1f} to {x_far:.1f}")

# ---- foot clash coda (NOT a set): one horizontal bar, 80,983 -> 55,706 + 25,277 ----
# Compact column right of the Venn; proximity does the work, no connector.
CX0, CW = 149.0, 28.0
sc = CW / ALL4
w_pass = FINAL_N * sc
yb0, bh_ = 38.0, 4.0
ax.text(CX0, yb0 + bh_ + 4.6, f"{ALL4:,}", fontsize=FS + 1, color=INK, ha="left", va="bottom", fontweight="bold")
ax.text(CX0, yb0 + bh_ + 1.2, "pass all four criteria", fontsize=FS, color=INK, ha="left", va="bottom")
# BOTH sections carry the same stroke. A stroke is centred on the path, so it
# grows the patch by lw/2 on every side: an lw=0 filled section beside an
# lw=0.7 hatched one is drawn 0.7 pt shorter and the two stop lining up
# (Leo, 2026-11-18: "sections of the bar are not the same height").
BAR_LW = 0.7
ax.add_patch(Rectangle((CX0, yb0), w_pass, bh_, facecolor=palette.THEME["700"],
                       edgecolor=palette.THEME["700"], lw=BAR_LW, zorder=2))
ax.add_patch(Rectangle((CX0 + w_pass, yb0), CW - w_pass, bh_, facecolor="none",
                       edgecolor=palette.THEME["700"], lw=BAR_LW, hatch="////", zorder=2))
ax.text(CX0, yb0 - 1.0, f"{FINAL_N:,}", fontsize=FS + 1, color=INK, ha="left", va="top", fontweight="bold")
t_pass, t_rem = "pass foot clash", "removed"
ax.text(CX0, yb0 - 1.0 - (FS + 1) * PT * 1.25, t_pass, fontsize=FS, color=INK, ha="left", va="top")
ax.text(CX0 + CW, yb0 - 1.0, f"{FOOT_REMOVED:,}", fontsize=FS, color=MUTED, ha="right", va="top")
ax.text(CX0 + CW, yb0 - 1.0 - (FS + 1) * PT * 1.25, t_rem, fontsize=FS, color=MUTED, ha="right", va="top")
assert style.text_w(t_pass, FS) * PT + style.text_w(t_rem, FS) * PT + 2 < CW, "coda labels collide"
assert CX0 + CW <= W
# the criterion itself: render + name, below the counts
put("footclash", CX0 + 7.5, 17.0, RW)
ax.text(CX0 + 16.5, 17.0 + 1.5, "Foot clash", fontsize=FS + 1, color=INK, ha="left", va="center", fontweight="semibold")
ax.text(CX0 + 16.5, 17.0 - 1.7, "\u2264 5", fontsize=FS + 1, color=INK, ha="left", va="center")
assert CX0 + 16.5 + style.text_w("Foot clash", FS + 1) * PT <= W, "foot clash label clipped"

# explicit dpi on EVERY format: embedded rasters are resampled to the save dpi
for ext in ("svg", "pdf", "png"):
    fig.savefig(out / f"{STEM}.{ext}", dpi=600)
print("size mm", W, H, "mode", MODE, "->", STEM)
