"""Fig3cd_pools: panel c = 50 largest backbone clusters, pool vs ipTM vs ordered
backbones; panel d = stack mix of the three pools as pies with arrows.
Emits output/Fig3cd_pools/panel_c.{svg,pdf,png} and panel_d.{...}.
Panel c source: inputs/fold_cluster_arms_k30.csv only (see README.md)."""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, Rectangle, Wedge
import matplotlib.patheffects as pe
from common import paths, palette, style

FIG = Path(__file__).resolve().parent.name
out = paths.fig_output(FIG)
style.apply()
HALO = [pe.withStroke(linewidth=1.6, foreground='white')]
FS = style.CFG["font_size_pt"]
MM = 25.4 / 72  # mm per pt
def pt(mm): return mm / MM

def emit(fig, stem):
    for ext in ("svg", "pdf", "png"):
        fig.savefig(out / f"{stem}.{ext}", dpi=300)   # no bbox tight: exact size
    print(stem, "size mm:", [round(v * 25.4, 1) for v in fig.get_size_inches()])

# ------------------------------------------------------------------ panel c
A = pd.read_csv(paths.INPUTS / "fold_cluster_arms_k30.csv")
assert len(A) == 917 and A.pool.sum() == 31763 and A.iptm.sum() == 7813 and A.library.sum() == 9759
assert (A.iptm <= A.pool).all() and (A.library <= A.pool).all()
print(f"pool {A.pool.sum():,}  ipTM {A.iptm.sum():,}  library {A.library.sum():,}  (k30 headline)")
N = 50
top = A.sort_values("pool", ascending=False, kind="stable").head(N)
ids = list(top.cluster); pc = top.pool.values; ic = top.iptm.values; lc = top.library.values
for i, a_, b_, c_ in zip(ids, pc, ic, lc):
    print(f"cluster {i:>3}  pool {a_:>4}  ipTM {b_:>4}  ordered {c_:>4}")
print("drawn pool", pc.sum(), "ipTM", ic.sum(), "ordered", lc.sum())
assert (ic <= pc).all() and (lc <= pc).all(), "an inset exceeds its pool column"

W, H = style.CFG["width_double_mm"], 44
fig = plt.figure(figsize=(W / 25.4, H / 25.4))
L, B, R, T = 11, 5.5, 1, 1
PW = W - L - R
ax = fig.add_axes([L / W, B / H, PW / W, (H - B - T) / H])
x = np.arange(N)
slot = PW / (N + 0.2)   # mm per column slot
YMAX = 800
CM = matplotlib.colormaps["thesis_div"]
C_IP, C_LIB = CM(0.12), CM(0.88)   # cool = ipTM, warm = ordered, as in Fig3d
bw, iw = 0.86, 0.36
assert iw * slot >= 0.8, f"inset {iw*slot:.2f} mm too thin"
ax.bar(x, pc, width=bw, color=palette.NEUTRAL["200"], edgecolor=palette.NEUTRAL["500"], linewidth=0.4)
ax.bar(x - iw / 2, ic, width=iw, color=C_IP, linewidth=0)
ax.bar(x + iw / 2, lc, width=iw, color=C_LIB, linewidth=0)
# every inset inside its pool column, in x as well as height
assert (iw <= bw / 2 * 2 - 1e-9) and ((ic <= pc).all()) and ((lc <= pc).all())
ax.set_xticks(x, [str(i) for i in ids])
# consecutive ids must not collide at the house size
gap = slot - max(style.text_w(str(i)) for i in ids) * MM
print(f"x-label gap {gap:.2f} mm (slot {slot:.2f} mm, inset {iw*slot:.2f} mm)")
assert gap > 0.6, "x labels collide"
ax.set_ylabel("Backbones", labelpad=1.5)
ax.set_xlim(-0.6, N - 0.4)
ax.set_ylim(0, YMAX)
ax.set_yticks([0, 200, 400, 600])
print("y range 0 -", YMAX)
ax.tick_params(length=2, width=0.5, pad=1.5)
ax.tick_params(axis="x", pad=1, labelsize=FS)
# ---- coverage, asserted
multi = A.pool > 1
cov = {k: (int((A[k] > 0).sum()), int(((A[k] > 0) & multi).sum()), int(((A[k] > 0) & ~multi).sum())) for k in ("iptm", "library")}
assert cov["iptm"] == (326, 206, 120) and cov["library"] == (917, 256, 661), cov
assert multi.sum() == 256 and (~multi).sum() == 661
assert ((A.iptm > 0) & (A.library == 0)).sum() == 0   # nothing is ipTM-only
miss = (A.iptm == 0)
assert miss.sum() == 591 and (miss & ~multi).sum() == 541
print("coverage (folds, multi, single): ipTM", cov["iptm"], " ordered", cov["library"],
      "| ipTM misses", miss.sum(), "of which singletons", (miss & ~multi).sum())

# ---- VIEW A: 917-fold matrix, row-major in pool-size order, drawn in mm on an overlay axes
Asort = A.sort_values("pool", ascending=False, kind="stable")
reached_ip = (Asort.iptm.values > 0)
NC, NR = 46, 20
gx0, gy0, cw, ch = 124.0, 8.0, 1.1, 0.7     # mm; grid 50.6 x 14 mm
oa = fig.add_axes([0, 0, 1, 1]); oa.set_xlim(0, W); oa.set_ylim(H, 0); oa.axis("off"); oa.patch.set_alpha(0)
C_BOTH = palette.NEUTRAL["400"]
rows_txt = [["."] * NC for _ in range(NR)]
for n, ri in enumerate(reached_ip):
    r_, c_ = divmod(n, NC)
    oa.add_patch(Rectangle((gx0 + c_ * cw, gy0 + r_ * ch), cw, ch, fc=C_BOTH if ri else C_LIB,
                           ec="white", lw=0.25))
    rows_txt[r_][c_] = "B" if ri else "o"
print("cell grid (row-major, largest fold first; B = reached by both, o = ordered only):")
for row in rows_txt: print("  " + "".join(row))
assert sum(r_.count("B") for r_ in rows_txt) == 326 and sum(r_.count("o") for r_ in rows_txt) == 591
oa.add_patch(plt.Polygon([(gx0, gy0), (gx0 + NC * cw, gy0), (gx0 + NC * cw, gy0 + ch), (gx0 + 4 * cw, gy0 + ch),
                          (gx0 + 4 * cw, gy0 + 2 * ch), (gx0, gy0 + 2 * ch)], closed=True, fill=False,
                         ec=palette.INK, lw=0.5))
T_A = "All 917 backbone clusters, largest first (outlined: 50 drawn)"
oa.text(gx0, gy0 - 1.2, T_A, fontsize=FS, ha="left", va="bottom", color=palette.INK)
assert gx0 + style.text_w(T_A) * MM <= W
for k_, (lab, col) in enumerate([("Reached by both selections", C_BOTH), ("Reached by ordered only", C_LIB)]):
    xx, yy = 90.0, gy0 + 7.0 + k_ * 3.4
    oa.add_patch(Rectangle((xx, yy - 0.8), 1.6, 1.6, fc=col, ec="none"))
    oa.text(xx + 2.4, yy, lab, fontsize=FS, ha="left", va="center", color=palette.INK)
    assert xx + 2.4 + style.text_w(lab) * MM <= gx0 - 0.8, "key hits grid"
x0, y0 = 10.0, 770
for k, (lab, col, ec) in enumerate([("Passed the gate", palette.NEUTRAL["200"], palette.NEUTRAL["500"]),
                                    ("ipTM-ranked selection", C_IP, C_IP),
                                    ("Diversity-optimised (ordered)", C_LIB, C_LIB)]):
    y = y0 - k * 80
    ax.add_patch(Rectangle((x0, y - 25), 0.6, 50, fc=col, ec=ec, lw=0.4))
    ax.text(x0 + 1.1, y, lab, va="center", ha="left", fontsize=FS, color=palette.INK)
    assert L + (x0 + 1.1 + 0.6) / (N + 0.2) * PW + style.text_w(lab) * MM <= W, "label clipped"
emit(fig, "panel_c")
plt.close(fig)

# ------------------------------------------------------------------ panel d
M = pd.read_csv(paths.INPUTS / "stackmix_k30.csv").set_index("stack")
M = M.loc[[palette.SHORT[k] for k in palette.STACK_ORDER]]   # fixed stack order
assert M["pool"].sum() == 55706 and M["ipTM"].sum() == 12000 and M["library"].sum() == 12000
print("pie totals: pool", M.pool.sum(), "ipTM", M.ipTM.sum(), "library", M.library.sum())

W, H = 90, 40     # emitted at 90 mm wide so the page places it at ~scale 1 (no font rescale)
fig = plt.figure(figsize=(W / 25.4, H / 25.4))
a = fig.add_axes([0, 0, 1, 1]); a.set_xlim(0, W); a.set_ylim(H, 0); a.axis("off")
a.set_aspect("equal")
cols = [palette.stack_color(s) for s in M.index]
LH = 2.3   # mm between label lines
# White type on every stack colour: all five are dark enough to carry it. If a
# colour is ever lightened this fails loudly rather than going quietly illegible.
def _lum(hex_):
    def ch(u):
        u = int(hex_.lstrip("#")[u:u + 2], 16) / 255
        return u / 12.92 if u <= 0.04045 else ((u + 0.055) / 1.055) ** 2.4
    return 0.2126 * ch(0) + 0.7152 * ch(2) + 0.0722 * ch(4)
assert all(_lum(c) < 0.40 for c in cols), "a stack colour is too light for white type"

PCT_R = 0.62      # label radius as a fraction of the pie radius
PCT_PAD = 0.6     # mm of air the text needs inside the slice edges
def pie(cx, cy, r, col):
    ang = -90.0   # 12 o'clock; y axis is inverted, so increasing angle runs clockwise on screen
    for v, c in zip(M[col].values, cols):
        d = 360.0 * v / M[col].sum()
        a.add_patch(Wedge((cx, cy), r, ang, ang + d, fc=c, ec="white", lw=0.4))
        # the share, set inside the slice - but only where the slice can hold it.
        # The room across a slice at the label radius is its chord there, so the
        # test is measured, not a percentage cutoff: the small stacks drop out
        # on their own rather than by a hand-picked threshold.
        txt = f"{100 * v / M[col].sum():.0f}%"
        rl = PCT_R * r
        if 2 * rl * np.sin(np.radians(d / 2)) >= style.text_w(txt) * MM + PCT_PAD:
            mid = np.radians(ang + d / 2)
            a.text(cx + rl * np.cos(mid), cy + rl * np.sin(mid), txt, ha="center",
                   va="center", fontsize=FS, color="white", zorder=3)
        ang += d
L_POOL = ["Passed the gate", f"{M.pool.sum():,}"]
L_IP = ["ipTM-ranked", f"{M.ipTM.sum():,}"]
L_LIB = ["Diversity-optimised", "(ordered)", f"{M.library.sum():,}"]   # wraps to two name lines
RP, RS = 10.0, 8.3
MARG, GAPL = 0.5, 1.5
wl = max(style.text_w(t_) * MM for t_ in L_LIB + L_IP)
lx = W - MARG - wl                       # common left edge of the label column
sx = lx - GAPL - RS                      # selection pie centre
px = MARG + RP
ty1 = MARG + RS; ty2 = ty1 + 2 * RS + 1.0
py = (ty1 + ty2) / 2
KY = H - 1.4
# horizontal budget, all from measured widths
print(f"x budget: pool {px-RP:.1f}-{px+RP:.1f}, arrows {px+RP+1:.1f}-{sx-RS-1.5:.1f} ({sx-RS-1.5-px-RP-1:.1f} mm), "
      f"selection pies {sx-RS:.1f}-{sx+RS:.1f}, labels {lx:.1f}-{W-MARG:.1f} (width {wl:.1f})")
assert sx - RS - 1.5 - (px + RP + 1.0) >= 15, "arrows too short"
assert ty2 + RS + 0.6 < KY - 0.8, "lower pie hits the key row"
pie(px, py, RP, "pool")
for j, t_ in enumerate(L_POOL):          # beneath the pool pie: the arrows leave from its right side
    a.text(px, py + RP + 1.0 + j * LH, t_, ha="center", va="top", fontsize=FS, color=palette.INK)
assert py + RP + 1.0 + LH + 2.0 < KY - 0.8, "pool label hits the key row"
for ty, col, lines in ((ty1, "ipTM", L_IP), (ty2, "library", L_LIB)):
    pie(sx, ty, RS, col)
    y0_ = ty - LH * (len(lines) - 1) / 2
    for j, t_ in enumerate(lines):
        a.text(lx, y0_ + j * LH, t_, ha="left", va="center", fontsize=FS, color=palette.INK)
    a.add_patch(FancyArrowPatch((px + RP + 1.0, py + (ty - py) * 0.3), (sx - RS - 1.5, ty),
                                arrowstyle="-|>", mutation_scale=pt(2.2), lw=0.6, color=palette.INK,
                                shrinkA=0, shrinkB=0))
SW, GAP_SW, GAP_E = 1.6, 1.0, 3.0
ent = [style.text_w(s_) * MM + SW + GAP_SW for s_ in M.index]
tot = sum(ent) + GAP_E * (len(ent) - 1)
print(f"key row {tot:.1f} mm of {W} mm (one line)")
assert tot <= W - 2, "key does not fit on one line; wrap to two rows"
kx = (W - tot) / 2
for s_, c, e in zip(M.index, cols, ent):
    a.add_patch(Rectangle((kx, KY - 0.8), SW, SW, fc=c, ec="none"))
    a.text(kx + SW + GAP_SW, KY, s_, va="center", ha="left", fontsize=FS, color=palette.INK)
    kx += e + GAP_E
emit(fig, "panel_d")

# ------------------------------------------------------------------ VIEW B: picks per fold vs fold size
W, H = 85, 45
fig = plt.figure(figsize=(W / 25.4, H / 25.4))
L, B, R, T = 11, 7.5, 1, 1
ax = fig.add_axes([L / W, B / H, (W - L - R) / W, (H - B - T) / H])
rng = np.random.default_rng(0)
ZLANE = 0.4          # y = 0 picks live in a lane drawn at this height on the log axis
SEP = 0.72           # break between the zero lane and the log axis proper
def yplot(v):
    v = np.asarray(v, float)
    # zeros: jittered inside the lane (jitter is vertical only, purely for visibility)
    z = ZLANE * 10 ** rng.uniform(-0.1, 0.1, size=len(v))
    return np.where(v > 0, v, z)
for col_, c_, mk, lab, ms in (("iptm", C_IP, "o", "ipTM-ranked", 2.2), ("library", C_LIB, "^", "Diversity-optimised (ordered)", 2.2)):
    ax.scatter(A.pool.values, yplot(A[col_].values), s=ms ** 2, marker=mk, color=c_, alpha=0.4, linewidths=0, zorder=3)
nz = {k: int((A[k] == 0).sum()) for k in ("iptm", "library")}
print("zero picks:", nz, "| ipTM zeros among singletons", int(((A.iptm == 0) & ~multi).sum()))
assert nz["iptm"] == 591 and nz["library"] == 0
ax.set_xscale("log"); ax.set_yscale("log")
ax.set_xlim(0.8, 800); ax.set_ylim(0.25, 400)
ax.axhline(SEP, color=palette.NEUTRAL["500"], lw=0.4, ls=(0, (2, 2)), zorder=1)
ax.set_xticks([1, 10, 100, 659], ["1", "10", "100", "659"])
ax.set_yticks([ZLANE, 1, 10, 100], ["0", "1", "10", "100"])
ax.minorticks_off()
ax.set_xlabel("Cluster size (pool backbones)", labelpad=1.5)
ax.set_ylabel("Selected backbones", labelpad=1.5)
ax.tick_params(length=2, width=0.5, pad=1.5)
ax.plot([1, 659], [1, 659], color=palette.NEUTRAL["500"], lw=0.4, zorder=2)   # y = x: whole fold taken
# direct labels
for k_, (lab, c_, mk) in enumerate((("ipTM-ranked selection", C_IP, "o"), ("Diversity-optimised (ordered)", C_LIB, "^"))):
    ya = 0.93 - k_ * 0.1
    ax.plot([0.03], [ya], mk, transform=ax.transAxes, color=c_, ms=2.4, mec="none", clip_on=False)
    ax.text(0.06, ya, lab, transform=ax.transAxes, fontsize=FS, va="center", ha="left", color=palette.INK)
# label set ALONG the diagonal, just above it, rotated to its on-screen angle
fig.canvas.draw()
p0, p1 = ax.transData.transform((3, 3)), ax.transData.transform((30, 30))
ang = np.degrees(np.arctan2(p1[1] - p0[1], p1[0] - p0[0]))
px = ax.transData.transform((11, 11))
off = 1.3 / 25.4 * fig.dpi                      # 1.3 mm perpendicular, up-left of the line
px = px + off * np.array([-np.sin(np.radians(ang)), np.cos(np.radians(ang))])
lx, ly = ax.transData.inverted().transform(px)
ax.text(lx, ly, "whole cluster taken", rotation=ang, rotation_mode="anchor", ha="center", va="bottom",
        fontsize=FS, color=palette.NEUTRAL["600"])
emit(fig, "panel_picks")
