"""Fig5c: interface chemistry, the two metrics that separate binders.

Assembles the two ChimeraX zooms with the numbers behind them.

  a  the repulsive case: GLU4 carboxylate on A1's phosphate, FC 0.84
  b  the clean case: ARG12 salt-bridges C6's phosphate, FC 1.56
  c  the acid-phosphate flag, binders against non-binders against the library
  d  the five features that define the Csy4 direction, Csy4 against the library
  e  flow fold change for all 32, flagged designs marked

Flow fold change is LV2.34.1 `ONOFF_avg` = mean(PosA, PosB) / mean(FullNeg,
NegA), which is what `flow_FC` in the candidate index is; verified identical to
4 decimal places. This is NOT the LV2.33 arrayed screen.

Binder = FC >= 1.2, a cut fixed before analysis, giving 15 binders / 17 not.
Everything here is from HANDOVER_METRICS.md; see the README for what does not
survive and for the robustness caveats, which are the honest headline.
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import yaml
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from PIL import Image

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "code"))
from common import paths, palette                                 # noqa: E402

FIG = Path(__file__).resolve().parent.name
out = paths.fig_output(FIG)
_cfg = yaml.safe_load((paths.ROOT / "config" / "style.yaml").read_text())
plt.rcParams.update({
    "font.family": _cfg["font_family"], "pdf.fonttype": 42, "svg.fonttype": "none",
    "savefig.dpi": 600, "axes.grid": False,
})
def mm(x): return x / 25.4

AN = Path.home() / "Downloads" / "design29_742db4d3" / "an32"
CUT = 1.2
INK, GREY = palette.INK, palette.NEUTRAL["500"]
OK, BAD = palette.THEME["700"], palette.DIVERGING[-1]    # NEED, see README

m = pd.read_pickle(AN / "m32.pkl")
g = pd.read_pickle(AN / "geom_wide.pkl")
for c in ("af3_if_K", "af3_if_R", "af3_if_arom", "af3_if_neg"):
    m[c] = g[c].values
fl = pd.read_csv(paths.INPUTS / "lv2_34_1_flow_on_off.csv")
d = m.merge(fl[["design_id", "ONOFF_avg"]], on="design_id").rename(
    columns={"ONOFF_avg": "FC"})
d["binder"] = d.FC >= CUT
d["flag"] = d.red_acid_phos_per_ifaceres_af3 > 0
LIB_FLAG = float((m.red_acid_phos_per_ifaceres_af3 > 0).mean())
nb, nn = int(d.binder.sum()), int((~d.binder).sum())
fb, fn = int(d[d.binder].flag.sum()), int(d[~d.binder].flag.sum())
print(f"{len(d)} designs, {nb} binders / {nn} not at FC >= {CUT}")
print(f"acid flag: binders {fb}/{nb}, non-binders {fn}/{nn}, "
      f"library {100*LIB_FLAG:.1f}%")

FEAT = [("aa_A", "alanine fraction", 0.069, -1),
        ("arom", "aromatic fraction", 0.085, +1),
        ("af3_if_arom", "aromatic at interface", 0.100, +1),
        ("af3_if_neg", "D/E at interface", 0.000, -1),
        ("af3_if_K", "Lys at interface", 0.000, -1)]

fig = plt.figure(figsize=(mm(178), mm(108)))

# ---------------------------------------------------- a, b: the two zooms
for k, (name, title, sub) in enumerate([
        ("acid_zoom", "a   Carboxylate on phosphate",
         "BG+BIF d503930b   FC 0.84   flagged"),
        ("clean_zoom", "b   Arginine on phosphate",
         "BG+LiMPNN e8cb617b   FC 1.56   clean")]):
    ax = fig.add_axes([0.035 + k * 0.335, 0.43, 0.30, 0.50])
    p = out / f"{name}.png"
    if p.exists():
        im = np.asarray(Image.open(p).convert("RGBA"))
        ys, xs = np.nonzero(im[..., 3] > 8)
        ax.imshow(im[ys.min():ys.max(), xs.min():xs.max()],
                  interpolation="antialiased")
    else:
        ax.text(0.5, 0.5, "render missing", ha="center", va="center",
                transform=ax.transAxes, fontsize=6, color=GREY)
    ax.set_xticks([]); ax.set_yticks([]); ax.set_frame_on(False)
    fig.text(0.035 + k * 0.335, 0.945, title, fontsize=7, color=INK,
             ha="left", va="bottom", weight="bold")
    fig.text(0.035 + k * 0.335, 0.405, sub, fontsize=5.6,
             color=palette.NEUTRAL["600"], ha="left", va="bottom")

# ------------------------------------------------------- c: the kill flag
cx = fig.add_axes([0.715, 0.62, 0.255, 0.30])
vals = [100 * fb / nb, 100 * fn / nn, 100 * LIB_FLAG]
labs = [f"proposed binders\nn={nb}", f"inactive designs\nn={nn}", "ordered library"]
cx.bar(range(3), vals, color=[OK, BAD, GREY], width=0.62, linewidth=0)
for i, (v, n) in enumerate(zip(vals, [f"{fb}/{nb}", f"{fn}/{nn}", ""])):
    cx.text(i, v + 1.6, (n if n else f"{v:.1f}%"), ha="center", va="bottom",
            fontsize=5.6, color=INK)
cx.set_xticks(range(3)); cx.set_xticklabels(labs, fontsize=5.4,
                                            color=palette.NEUTRAL["600"])
cx.set_ylim(0, 52); cx.set_yticks([0, 25, 50])
cx.set_yticklabels(["0", "25", "50%"], fontsize=5.6,
                   color=palette.NEUTRAL["600"])
for sp in ("top", "right"):
    cx.spines[sp].set_visible(False)
for sp in ("left", "bottom"):
    cx.spines[sp].set_color(palette.NEUTRAL["400"]); cx.spines[sp].set_linewidth(0.6)
cx.tick_params(length=2, width=0.6, pad=1.5, color=palette.NEUTRAL["400"])
cx.set_ylabel("carry a repulsive\nacid-phosphate contact", fontsize=5.6,
              color=palette.NEUTRAL["600"], labelpad=2, linespacing=1.5)
cx.text(0, 1.06, "c   The kill flag", transform=cx.transAxes, fontsize=7,
        color=INK, ha="left", va="bottom", weight="bold")
cx.text(1.0, 1.06, "Fisher p = 0.0077", transform=cx.transAxes, fontsize=5.4,
        color=palette.NEUTRAL["600"], ha="right", va="bottom")

# -------------------------------------------- d: the Csy4 direction features
dx = fig.add_axes([0.165, 0.135, 0.30, 0.215])
ys = np.arange(len(FEAT))[::-1]
for y, (col, lab, csy4, sign) in zip(ys, FEAT):
    lo, hi = m[col].quantile([0.05, 0.95])
    dx.plot([lo, hi], [y, y], lw=2.6, color=palette.NEUTRAL["300"],
            solid_capstyle="round", zorder=1)
    dx.plot([m[col].mean()], [y], "o", ms=3.2, color=GREY, zorder=2)
    dx.plot([csy4], [y], "D", ms=4.0, color=OK, zorder=3)
dx.set_yticks(ys)
dx.set_yticklabels([f"{lab}  {'less' if s < 0 else 'more'}"
                    for _, lab, _, s in FEAT], fontsize=5.6, color=INK)
dx.set_xlim(-0.02, 0.33)
dx.set_xticks([0, 0.1, 0.2, 0.3])
dx.set_xticklabels(["0", "0.1", "0.2", "0.3"], fontsize=5.6,
                   color=palette.NEUTRAL["600"])
for sp in ("top", "right", "left"):
    dx.spines[sp].set_visible(False)
dx.spines["bottom"].set_color(palette.NEUTRAL["400"])
dx.spines["bottom"].set_linewidth(0.6)
dx.tick_params(length=2, width=0.6, pad=1.5, color=palette.NEUTRAL["400"])
dx.set_xlabel("fraction", fontsize=5.6, color=palette.NEUTRAL["600"], labelpad=1.5)
fig.text(0.035, 0.385, "d   What natural Csy4 does differently", fontsize=7,
         color=INK, ha="left", va="bottom", weight="bold")
fig.text(0.035, 0.360, "diamond Csy4, dot library mean, bar 5-95% of the library",
         fontsize=5.2, color=palette.NEUTRAL["600"], ha="left", va="bottom")

# ------------------------------------------------------ e: the 32 measured
ex = fig.add_axes([0.575, 0.135, 0.395, 0.215])
o = d.sort_values("FC").reset_index(drop=True)
xs = np.arange(len(o))
ex.bar(xs, o.FC, width=0.68, linewidth=0,
       color=[BAD if f else (OK if b else GREY)
              for f, b in zip(o.flag, o.binder)])
ex.axhline(CUT, lw=0.6, color=palette.NEUTRAL["600"], ls=(0, (3, 2)))
ex.text(0.3, CUT + 0.03, f"binder cut {CUT}", fontsize=5.2,
        color=palette.NEUTRAL["600"], va="bottom")
ex.set_ylim(0, 1.72); ex.set_yticks([0, 0.5, 1.0, 1.5])
ex.set_yticklabels(["0", "0.5", "1.0", "1.5"], fontsize=5.6,
                   color=palette.NEUTRAL["600"])
ex.set_xticks([])
for sp in ("top", "right"):
    ex.spines[sp].set_visible(False)
for sp in ("left", "bottom"):
    ex.spines[sp].set_color(palette.NEUTRAL["400"]); ex.spines[sp].set_linewidth(0.6)
ex.tick_params(length=2, width=0.6, pad=1.5, color=palette.NEUTRAL["400"])
ex.set_ylabel("flow ON/OFF fold change", fontsize=5.6,
              color=palette.NEUTRAL["600"], labelpad=2)
ex.set_xlabel("the 32 ordered designs, ranked", fontsize=5.6,
              color=palette.NEUTRAL["600"], labelpad=1.5)
fig.text(0.575, 0.385, "e   Every design measured", fontsize=7, color=INK,
         ha="left", va="bottom", weight="bold")
for i, (c, t) in enumerate([(BAD, "flagged"), (GREY, "below cut"),
                            (OK, "binder, clean")]):
    ex.text(0.33 + i * 0.23, 1.02, "■", transform=ex.transAxes,
            fontsize=6, color=c, ha="left", va="bottom")
    ex.text(0.355 + i * 0.23, 1.02, t, transform=ex.transAxes, fontsize=5.2,
            color=palette.NEUTRAL["600"], ha="left", va="bottom")

for ext in ("svg", "pdf", "png"):
    fig.savefig(out / f"interface_chemistry.{ext}", transparent=False,
                dpi=600 if ext == "png" else None)
plt.close(fig)
print("->", out / "interface_chemistry.png")
