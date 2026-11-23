"""Fig 3c: Alpha/Beta/Loop ternary plots, 4 generators x (before, after) filtering."""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import numpy as np, pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import to_rgb
from common import paths, style, palette

FIG = Path(__file__).resolve().parent.name
style.apply()
out = paths.fig_output(FIG)

# ---- data -----------------------------------------------------------------
d = pd.read_parquet(paths.INPUTS / "sse_stage.parquet")
d["gen"] = d["stack"].map({"RFD3_LiMPNN": "RFD3", "BG_LiMPNN": "BoltzGen", "BG_BIF": "BoltzGen",
                           "BD1_LiMPNN": "BD1", "HD_LiMPNN": "HD"})
d["loop"] = 1 - d.helix - d.strand
GENS = ["RFD3", "BoltzGen", "BD1", "HD"]
LONG = {"RFD3": "RFD3 + LigandMPNN", "BoltzGen": "BoltzGen + LigandMPNN",
        "BD1": "BoltzDesign + LigandMPNN", "HD": "HalluDesign + LigandMPNN"}
assert [g for g in GENS] == [g for s in palette.STACK_ORDER for g in GENS if LONG[g] == s]  # order = STACK_ORDER
# stage is EXCLUSIVE (0 gate-fail, 1 gate, 2 library): after filtering = stage >= 1
EXPECT = {"RFD3": (144015, 14312), "BoltzGen": (149933, 9895), "BD1": (113175, 1661), "HD": (91941, 444)}
sets = {}
for g in GENS:
    s = d[d.gen == g]
    sets[g, 0], sets[g, 1] = s, s[s.stage >= 1]
    got = (len(sets[g, 0]), len(sets[g, 1]))
    print(f"{g:9s} before {got[0]:>8,}  after {got[1]:>7,}")
    assert got == EXPECT[g], (g, got, EXPECT[g])

# ---- geometry (mm; one axes spans the whole canvas) ------------------------
# 2026-11-19 (Leo): "remove top triangle lane". The ALL-BACKBONES row is gone;
# only the passed-gate row is drawn. Its unfiltered counterpart is still
# computed - it is the dashed 100% outline on every triangle, which is now the
# only place the before-gate distribution appears. H 62 -> 34 for the one row.
W, H = style.CFG["width_double_mm"], 34
SIDE, TH = 26.0, 26.0 * 3**.5 / 2
X0, PITCH = 17.0, 41.5
BASE = {0: 33.3, 1: 6.5}                      # baseline y of row 0 (before) / row 1 (after)

def xy(a, l, x0, y0):                         # Beta (0,0), Loop (1,0), Alpha apex
    return x0 + SIDE * (l + a / 2), y0 + TH * a

def density(a, l, sigma, nb=140):
    x, y = l + a / 2, a * 3**.5 / 2
    Z, xe, ye = np.histogram2d(x, y, bins=nb, range=[[0, 1], [0, 3**.5 / 2]])
    r = int(4 * sigma); k = np.exp(-.5 * (np.arange(-r, r + 1) / sigma) ** 2); k /= k.sum()
    Z = np.apply_along_axis(lambda v: np.convolve(v, k, "same"), 0, Z)
    Z = np.apply_along_axis(lambda v: np.convolve(v, k, "same"), 1, Z)
    return Z, (xe[:-1] + xe[1:]) / 2, (ye[:-1] + ye[1:]) / 2

def hdr(Z, mass):                             # density level enclosing `mass` of the points
    z = np.sort(Z.ravel())[::-1]; c = np.cumsum(z) / z.sum()
    return z[np.searchsorted(c, mass)]

def at_pts(Z, s):                             # smoothed density at each data point's bin
    nb_ = Z.shape[0]
    ix = np.clip((s.loop.values + s.helix.values / 2) * nb_, 0, nb_ - 1).astype(int)
    iy = np.clip(s.helix.values * nb_, 0, nb_ - 1).astype(int)
    return Z[ix, iy]

def lvl100(Z, s):                             # "100%" level = least-dense data point, so every point is inside
    return at_pts(Z, s).min()

TONES = (.15, .4, .8)                         # outer 100% / 90% / inner 50%
fig = plt.figure(figsize=(style.mm(W), style.mm(H)))
ax = fig.add_axes([0, 0, 1, 1]); ax.set_xlim(0, W); ax.set_ylim(0, H); ax.axis("off")
ink, fs = palette.INK, style.CFG["font_size_pt"]
SIGMA = {"RFD3": 2.5, "BoltzGen": 2.5, "BD1": 3.5, "HD": 6.0}   # bins of 1/140; by the generator's smaller set
CHECK = {}
DASH_MASS = .99                                 # unfiltered-set region drawn dashed in the after row
lbl = {0: "All backbones", 1: "Passed gate"}

for j, g in enumerate(GENS):
    col = palette.stack_color(LONG[g])
    x0 = X0 + j * PITCH
    for r in (1,):                            # row 0 (all backbones) removed
        y0 = BASE[r]
        s = sets[g, r]
        tri = np.array([xy(0, 0, x0, y0), xy(0, 1, x0, y0), xy(1, 0, x0, y0), xy(0, 0, x0, y0)])
        # contour fill needs clipping to the triangle (blur bleeds past the edges)
        patch = plt.Polygon(tri[:3], closed=True, fc="none", ec=ink, lw=.5, zorder=4)
        ax.add_patch(patch)
        # SAME sigma for both rows of a generator (so the pair is comparable); wider for generators with small n
        sigma = SIGMA[g]
        Z, gx, gy = density(s.helix.values, s.loop.values, sigma)
        GX, GY = np.meshgrid(gx, gy, indexing="ij")
        l100 = lvl100(Z, s)
        cs = ax.contourf(x0 + SIDE * GX, y0 + SIDE * GY, Z, levels=[l100, hdr(Z, .9), hdr(Z, .5), 1e9],
                         colors=[tuple(1 - (1 - np.array(to_rgb(col))) * f) for f in TONES], zorder=2)
        cs.set_clip_path(patch)
        if r == 1:                            # support of the unfiltered set: 99% region, same sigma, clipped
            sb = sets[g, 0]
            Zb, gx, gy = density(sb.helix.values, sb.loop.values, sigma)
            lvb = lvl100(Zb, sb)
            cd = ax.contour(x0 + SIDE * GX, y0 + SIDE * GY, Zb, levels=[lvb],
                            colors=[palette.NEUTRAL["500"]], linewidths=.4, linestyles="dashed", zorder=3)
            cd.set_clip_path(patch)
            a100 = Z >= l100
            CHECK[g] = ((at_pts(Zb, s) < lvb).mean() * 100, (a100 & (Zb < lvb)).sum() / a100.sum() * 100,
                        a100.sum() / (Zb >= lvb).sum())
        # vertex names, every triangle
        ax.text(x0 - 1, y0, "Beta", ha="right", va="center", fontsize=fs, color=ink)
        ax.text(x0 + SIDE + 1, y0, "Loop", ha="left", va="center", fontsize=fs, color=ink)
        ax.text(x0 + SIDE / 2, y0 + TH + 1.2, "Alpha", ha="center", va="bottom", fontsize=fs, color=ink)
        # direct label: generator over its n. Both used to be split across the
        # two rows (name on the top one, n on each); with one row they stack here.
        top = y0 + TH
        ax.text(x0 - 4, top - .3, g, ha="left", va="top", fontsize=fs, color=col, fontweight="bold")
        ax.text(x0 - 4, top - 3.6, f"n = {len(s):,}", ha="left", va="top", fontsize=fs,
                color=ink, linespacing=1)
        # ticks + numbers: ONE triangle only (RFD3, after row)
        if j == 0 and r == 1:
            for t in (0.25, .5, .75):
                s_ = f"{t:.2f}".rstrip("0") if t != .5 else "0.5"
                px, py = xy(0, 1 - t, x0, y0)      # bottom edge: Beta fraction t (Loop = 1 - t)
                ax.plot([px, px], [py, py - .8], c=ink, lw=.4); ax.text(px, py - 1.1, s_, ha="center", va="top", fontsize=fs - .5, color=ink)
                px, py = xy(t, 0, x0, y0)           # left edge: Alpha fraction t
                ax.plot([px, px - .7], [py, py + .4], c=ink, lw=.4)
                px2, py2 = xy(1 - t, t, x0, y0)     # right edge: Loop fraction t
                ax.plot([px2, px2 + .7], [py2, py2 + .4], c=ink, lw=.4)
                ax.text(px - 1.0, py, s_, ha="right", va="center", fontsize=fs - .5, color=ink)
                ax.text(px2 + 1.0, py2, s_, ha="left", va="center", fontsize=fs - .5, color=ink)

# ---- key: frameless, one line below the grid ----
mmw = lambda t: style.text_w(t, fs) / 72 * 25.4            # text width in mm
KY, x, sw = 1.3, 3.0, 1.6
def put(t):
    global x
    ax.text(x, KY, t, ha="left", va="center", fontsize=fs, color=ink); x += mmw(t) + 1.0
def swatch(c):
    global x
    ax.add_patch(plt.Rectangle((x, KY - sw / 2), sw, sw, fc=c, ec="none", zorder=5)); x += sw + .5
for c, t in ((palette.NEUTRAL["700"], "50%"), (palette.NEUTRAL["400"], "90%"), (palette.NEUTRAL["200"], "100%")):
    swatch(c); put(t)
x += .8
ax.plot([x, x + 3.0], [KY] * 2, c=palette.NEUTRAL["500"], lw=.4, ls="dashed"); x += 3.6
put("unfiltered")
print("key ends at mm:", round(x, 1))

for r in (1,):
    ax.text(3, BASE[r] + TH / 2, lbl[r], rotation=90, ha="center", va="center", fontsize=fs, color=ink)

plt.rcParams["savefig.pad_inches"] = 0   # emitted size = ink extent
style.save(fig, out, "panel_c")
from PIL import Image
im = Image.open(out / "panel_c.png")
print("generator   passed outside dashed 100%   after-100% area outside   after-100% / before-100% area")
for g, (p, a_, r_) in CHECK.items():
    print(f"  {g:9s}  {p:.2f}%   {a_:.2f}%   {r_:.2f}")
mm_ = lambda p: p / 300 * 25.4
print("emitted mm: %.1f x %.1f" % (mm_(im.size[0]), mm_(im.size[1])))
assert mm_(im.size[1]) <= 34
fig.savefig(out / "panel_c.svg", bbox_inches="tight")   # assembly reads SVG; same bbox/pad as pdf/png
