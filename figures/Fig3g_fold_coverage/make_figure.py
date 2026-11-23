"""Fig3g_fold_coverage: all 917 backbone clusters as an EQUAL-cell matrix (largest first) coloured by ipTM coverage,
plus a thin strip splitting the 31,763 backbones by the same colours.
Source: inputs/fold_cluster_arms_k30.csv only. Emits output/Fig3g_fold_coverage/panel.{svg,pdf,png}."""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from common import paths, palette, style

FIG = Path(__file__).resolve().parent.name
out = paths.fig_output(FIG)
style.apply()
FS = style.CFG["font_size_pt"]
MM = 25.4 / 72

# ---- data + asserts
A = pd.read_csv(paths.INPUTS / "fold_cluster_arms_k30.csv")
assert len(A) == 917 and A.pool.sum() == 31763 and A.iptm.sum() == 7813 and A.library.sum() == 9759
multi = A.pool > 1
assert ((A.iptm > 0).sum(), ((A.iptm > 0) & multi).sum(), ((A.iptm > 0) & ~multi).sum()) == (326, 206, 120)
assert ((A.library > 0).sum(), ((A.library > 0) & multi).sum(), ((A.library > 0) & ~multi).sum()) == (917, 256, 661)
miss = A.iptm == 0
assert miss.sum() == 591 and (miss & ~multi).sum() == 541 and A.pool[miss].sum() == 657
assert ((A.iptm > 0) & (A.library == 0)).sum() == 0
S = A.sort_values("pool", ascending=False, kind="stable").reset_index(drop=True)
TOP = 50

W, H = style.CFG["width_double_mm"], 45
NC, P, GAP = 71, 2.49, 0.25            # columns, cell pitch (mm), white gap between cells
NR = -(-917 // NC); assert NC * P <= W
gx, gy = (W - NC * P) / 2, 4.2
fig = plt.figure(figsize=(W / 25.4, H / 25.4))
ax = fig.add_axes([0, 0, 1, 1]); ax.set_xlim(0, W); ax.set_ylim(H, 0); ax.axis("off")
C_BOTH, C_ORD = palette.NEUTRAL["400"], palette.DIVERGING[-1]   # warm end = ordered, as in Fig3cd panel c
n_both = n_ord = 0
for n, ip in enumerate(S.iptm.values > 0):
    r, c = divmod(n, NC)
    ax.add_patch(Rectangle((gx + c * P + GAP / 2, gy + r * P + GAP / 2), P - GAP, P - GAP, fc=C_BOTH if ip else C_ORD, ec="none"))
    n_both += ip; n_ord += not ip
assert (n_both, n_ord) == (326, 591)
ink = palette.INK
def T(x, y, s, **k):
    k.setdefault('color', ink); ax.text(x, y, s, fontsize=FS, **k); return style.text_w(s) * MM
def key(x, y, col, s, outline=False):
    ax.add_patch(Rectangle((x, y - 0.8), 1.6, 1.6, fc="none" if outline else col, ec=ink if outline else "none", lw=0.6))
    return x + 2.4 + T(x + 2.4, y, s, ha="left", va="center")
ky = 2.0
xe = key(gx, ky, C_BOTH, "326 backbone clusters reached by both selections")
xe = key(xe + 6, ky, C_ORD, "591 backbone clusters reached by ordered only")
# reading-order cue, right end of the key row: largest -> smallest, along rows
lg, sm, rw = "Largest", "smallest", ""
AL = 12.0
xc = W - gx - (style.text_w(lg) + style.text_w(sm)) * MM - AL - 2.0
assert xc > xe + 6, "cue collides with key"
xc += T(xc, ky, lg, ha="left", va="center") + 1.0
ax.annotate("", xy=(xc + AL, ky), xytext=(xc, ky), arrowprops=dict(arrowstyle="-|>", lw=0.5, color=ink, mutation_scale=4, shrinkA=0, shrinkB=0))
xc += AL + 1.0
T(xc, ky, sm, ha="left", va="center")
assert xe <= W, "key clipped"
# mass strip: same two colours, backbones instead of folds
assert A.pool[A.iptm == 0].max() == 5
# ---- caption line (the mass fact) and strip B (clusters)
tot = int(S.pool.sum()); miss_bb = int(S.pool[S.iptm == 0].sum()); assert (tot, miss_bb) == (31763, 657)
assert A.pool[A.iptm == 0].max() == 5
assert round(100 * (917 / 326 - 1)) == 181 and round(100 * miss_bb / tot) == 2
SW, bh = NC * P, 4.2
WHITE = palette.SURFACE
def strip(y, n_grey, n_all):
    wg = SW * n_grey / n_all
    ax.add_patch(Rectangle((gx, y), wg, bh, fc=C_BOTH, ec="none"))
    ax.add_patch(Rectangle((gx + wg, y), SW - wg, bh, fc=C_ORD, ec="none"))
    return wg
cy = gy + NR * P - 1.3   # (no caption line: the 657-backbone fact lives in the figure caption)
yb = cy + 2.0 + 1.3
wg = strip(yb, 326, 917)
tb1 = "326 of 917 backbone clusters reached by both selections"
tb2 = "591 backbone clusters reached by ordered only: +181% more clusters"
assert style.text_w(tb1) * MM < wg - 2 and style.text_w(tb2) * MM < SW - wg - 2
T(gx + 1.0, yb + bh / 2, tb1, ha="left", va="center")
T(gx + wg + 1.0, yb + bh / 2, tb2, ha="left", va="center", color=WHITE)
assert yb + bh <= H and cy + 2.0 <= yb, "caption collides with strip B"
assert yb + bh + 0.4 <= H
for ext in ("svg", "pdf", "png"):
    fig.savefig(out / f"panel.{ext}", dpi=300)
print("size mm:", [round(v * 25.4, 1) for v in fig.get_size_inches()])
print("cell pitch mm", P, "drawn", P-GAP, "| rows", NR, "| strip B bottom", yb+bh)
