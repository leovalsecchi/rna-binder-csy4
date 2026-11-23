"""Fig3a (new): AF3 vs ESMFold2 agreement for ipTM, DockQ and ligand_RMSD.

Per-backend columns only (_af3/_esm); the unsuffixed iptm/sc_dockq/sc_rmsd are
min/min/max consensus and would collapse the scatter. 2-D histogram, log counts.
Gate (funnel_order_v2.py): ipTM >= 0.6, DockQ >= 0.7, ligand_RMSD <= 3 A.
Dashed lines = thresholds on both axes; one % per quadrant (share of drawn points, sums to 100).
The gate itself is "either backend passes" (ipTM: max) = union of 3 quadrants: caption, not drawn.
BG+BIF DockQ/RMSD come from inputs/large/bif_true_per_backend.csv.gz (fix_bif).
Emits exactly 178 x 37 mm as PDF, PNG and SVG (live text).
Dev speed-up: FIG3A_CACHE=<file.pkl> caches the merged frame (not needed for a normal run).
"""
import os
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from common import palette, paths, style

FIG = Path(__file__).resolve().parent.name
style.apply()
out = paths.fig_output(FIG)
W_MM, H_MM = 178, 37

TRUE_BIF = paths.INPUTS / "large" / "bif_true_per_backend.csv.gz"


def fix_bif(d):
    """ONE place for the BG+BIF DockQ/RMSD source. Returns (d, drop): drop marks BG+BIF
    rows with no true-backbone value, excluded from the DockQ/RMSD sub-panels. The 43,416
    are genuinely unrescoreable on Euler: 42,447 no_native_interface, 899 unpaired
    (no_limpnn_uuid), 70 ref_error. Swap this function to change the source."""
    bif = d["stack"].eq("BG+BIF")
    if not TRUE_BIF.exists():
        print(f"WARNING: {TRUE_BIF.name} missing -> excluding ALL {bif.sum():,} BG+BIF rows "
              "from DockQ/RMSD (stored values are wrong for that stack)")
        return d, bif
    t = pd.read_csv(TRUE_BIF).drop_duplicates("design_id").set_index("design_id")
    d = d.copy()
    for k in ("dockq", "rmsd"):
        for b in ("af3", "esm"):
            d.loc[bif, f"{k}_{b}"] = d.loc[bif, "design_id"].map(t[f"{k}_{b}_true"]).to_numpy()
    drop = bif & d[["dockq_af3", "dockq_esm", "rmsd_af3", "rmsd_esm"]].isna().any(axis=1)
    print(f"BG+BIF: {bif.sum():,} rows; true per-backend values drawn for {int((bif&~drop).sum()):,}; "
          f"{int(drop.sum()):,} unrescoreable, excluded from DockQ/RMSD")
    return d, drop


cache = os.environ.get("FIG3A_CACHE")
if cache and os.path.exists(cache):
    d, drop = pd.read_pickle(cache)
else:
    cols = ["design_id", "stack", "iptm_af3", "iptm_esm", "dockq_af3", "dockq_esm",
            "rmsd_af3", "rmsd_esm", "bif_rescored"]
    d = pd.read_csv(paths.INPUTS / "large" / "campaign_metrics_full.csv.gz", usecols=cols)
    d, drop = fix_bif(d)
    d = d.drop(columns="design_id")
    if cache: pd.to_pickle((d, drop), cache)
print(f"rows: {len(d):,}")

# (key, label, lo, hi, thr, passes_if_ge, log, ticks, excl_bif)
M = [("iptm", "ipTM", 0, 1, 0.6, True, False, [0, .5, 1], False),
     ("dockq", "DockQ", 0, 1, 0.7, True, False, [0, .5, 1], True),
     ("rmsd", "RMSD (Å)", 0.5, 100, 3.0, False, True, [1, 10, 100], True)]
# The ramp and the count scale now live in code/common/style.py, shared with
# code/supp/repredictor_comparison, which draws the same AF3-ESMFold2 cells
# among its nine. Two copies would eventually disagree there. The four-round
# history of the colours and the reasoning behind the percentile scale are
# recorded on those two functions.
cm = style.density_cmap("fig3a_blue")
MUTED = palette.NEUTRAL["500"]


fig = plt.figure(figsize=(style.mm(W_MM), style.mm(H_MM)))
PW, PH, X0, GAP, Y0 = 43.0, 23.0, 10.0, 15.0, 8.0     # mm
for i, (k, lab, lo, hi, t, ge, log, ticks, ex) in enumerate(M):
    x0 = X0 + i * (PW + GAP)
    ax = fig.add_axes([x0 / W_MM, Y0 / H_MM, PW / W_MM, PH / H_MM])
    s = d[~drop] if ex else d
    x, y = s[f"{k}_af3"].to_numpy(float), s[f"{k}_esm"].to_numpy(float)
    ok = np.isfinite(x) & np.isfinite(y); x, y = x[ok], y[ok]
    n = len(x)
    pa, pe = (x >= t, y >= t) if ge else (x <= t, y <= t)
    q = dict(both=int((pa & pe).sum()), af3=int((pa & ~pe).sum()),
             esm=int((~pa & pe).sum()), none=int((~pa & ~pe).sum()))
    dis = q["af3"] + q["esm"]
    r = np.corrcoef(x, y)[0, 1]
    rho = pd.Series(x).rank().corr(pd.Series(y).rank())
    extra = (f"; Pearson on log10 = {np.corrcoef(np.log10(np.clip(x,1e-3,None)), np.log10(np.clip(y,1e-3,None)))[0,1]:.3f}"
             if log else "")
    print(f"{k}: points drawn n={n:,} (excluded {len(d)-n:,}) Pearson={r:.3f} Spearman={rho:.3f}{extra}"
          f" | opposite sides of {t}: {dis:,} = {100*dis/n:.2f}%"
          f" | both {q['both']:,} ({100*q['both']/n:.2f}%), AF3 only {q['af3']:,} ({100*q['af3']/n:.2f}%),"
          f" ESM only {q['esm']:,} ({100*q['esm']/n:.2f}%), neither {q['none']:,} ({100*q['none']/n:.2f}%)"
          f" | gate pass (any) {n-q['none']:,}")
    # 2026-11-19 (Leo): "turn the Fig3a into hexes". Hexagonal bins instead of
    # the square pcolormesh grid: a hex lattice has no axis-aligned rows, so the
    # square grid's faint vertical/horizontal banding over quantised scores goes
    # away and a diagonal cloud reads as a cloud rather than as a staircase.
    #
    # The grid is chosen against the 0.01 quantisation of ipTM: 50 columns
    # across 0..1 is exactly TWO quanta per hex, so the lattice samples the data
    # instead of beating against it. It is also coarse enough that the hexes
    # READ as hexes - at 43 mm wide a 100-column grid prints 0.43 mm cells and
    # the texture disappears into a smooth wash, which is the square grid's look
    # under another name. The log panel is not quantised and takes 40.
    # Counts still colour on a log scale; empty hexes are not drawn (mincnt=1).
    xc, yc = np.clip(x, lo, hi), np.clip(y, lo, hi)       # tail -> edge bins
    ext = ((np.log10(lo), np.log10(hi)) * 2 if log else (lo, hi) * 2)
    hb = ax.hexbin(xc, yc, gridsize=(40 if log else 50), extent=ext,
                   xscale="log" if log else "linear",
                   yscale="log" if log else "linear",
                   cmap=cm, mincnt=1, linewidths=0, rasterized=True)

    # Count scale: see style.density_norm - spans p25 to p99.5 of the occupied
    # hexes, because LogNorm(1, max) spent half the ramp where ~1-4% of the
    # points are and put every panel on a different number of decades.
    _lo, _hi, _pos = style.density_norm(hb)
    _c = hb.get_array()
    _c = np.asarray(_c[~np.ma.getmaskarray(_c)] if np.ma.isMaskedArray(_c) else _c, float)
    print(f"  {k} norm {_lo:,.0f}-{_hi:,.0f}  ramp position: "
          f"median {_pos(np.median(_c)):.2f}  p90 {_pos(np.percentile(_c, 90)):.2f}")
    ax.plot([lo, hi], [lo, hi], color=MUTED, lw=0.4, ls=":", zorder=3)
    ax.axvline(t, color=palette.INK, lw=0.6, ls="--", zorder=4)
    ax.axhline(t, color=palette.INK, lw=0.6, ls="--", zorder=4)
    if log: ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlim(lo, hi); ax.set_ylim(lo, hi)
    ax.set_xticks(ticks); ax.set_yticks(ticks)
    ax.set_xticklabels([f"{v:g}" for v in ticks]); ax.set_yticklabels([f"{v:g}" for v in ticks])
    ax.minorticks_off()
    arrow = " ↓" if not ge else ""
    ax.set_xlabel(f"AF3 {lab}{arrow}", labelpad=1.5); ax.set_ylabel(f"ESMFold2 {lab}{arrow}", labelpad=1.5)
    # 2026-11-19 (Leo): R^2, not r - the share of variance one backend's score
    # explains in the other, which is the honest number for "how much do they
    # agree". r = 0.25 / 0.50 / 0.57 becomes 0.06 / 0.25 / 0.33. Written with a
    # literal U+00B2, NOT mathtext: $R^2$ renders in DejaVu and would be the
    # only non-Arial glyph in the book (see the 2026-11-16 font audit). The
    # stdout line above still prints Pearson and Spearman unchanged.
    ax.text(0.0, 1.04, f"R\u00b2 = {r ** 2:.2f}", transform=ax.transAxes,
            ha="left", va="bottom", color=palette.INK)
    # one % per quadrant, in its outer corner; largest-remainder rounding so the four sum to 100.0
    raw = {kk: 1000 * v / n for kk, v in q.items()}
    pct = {kk: int(v) for kk, v in raw.items()}
    for kk in sorted(raw, key=lambda kk: raw[kk] - pct[kk], reverse=True)[:1000 - sum(pct.values())]:
        pct[kk] += 1
    print(f"  {k} quadrant %: " + ", ".join(f"{kk} {pct[kk]/10:.1f}" for kk in q) + f" (sum {sum(pct.values())/10:.1f})")
    hx = {True: 0.96 if ge else 0.04, False: 0.04 if ge else 0.96}   # side where AF3 passes / fails
    qw = {True: (1 - (t - lo) / (hi - lo)) if ge else (t - lo) / (hi - lo)}
    if log:
        f = np.log(t / lo) / np.log(hi / lo); qw = {True: 1 - f if ge else f}
    qw[False] = 1 - qw[True]
    for ap in (True, False):   # label must fit inside its quadrant (mm)
        assert style.text_w("99.9%") / 72 * 25.4 < qw[ap] * PW - 1.5, (k, ap)
    for (ap, ep), key in {(True, True): "both", (True, False): "af3", (False, True): "esm", (False, False): "none"}.items():
        cx, cy = hx[ap], (0.93 if ep == ge else 0.07)
        ax.text(cx, cy, f"{pct[key]/10:.1f}%", transform=ax.transAxes,
                ha="right" if cx > 0.5 else "left", va="top" if cy > 0.5 else "bottom",
                color=palette.INK, bbox=dict(fc="white", ec="none", alpha=0.6, pad=0.3))
fig.savefig(out / "Fig3a_backend_agreement.pdf")
fig.savefig(out / "Fig3a_backend_agreement.png")
fig.savefig(out / "Fig3a_backend_agreement.svg")
print(f"canvas {W_MM} x {H_MM} mm -> {out}")
