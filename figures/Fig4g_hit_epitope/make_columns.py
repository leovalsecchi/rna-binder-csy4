"""Fig4g, column form: the same epitope comparison as make_hairpin.py, read
along the sequence instead of around the fold.

Per nucleotide, the mean share of a design's protein-RNA contacts, split into
the part reaching the base and the part reaching the sugar-phosphate backbone.
Library bar beside hit bar; the hit bar carries a bootstrap 95% interval,
because n = 32 and the whole point of the panel is whether the two differ by
more than sampling noise.

The split is base against backbone, not phosphate / sugar / base. geom.py
stores per-nucleotide totals and the base part of each (`nt{i}` and
`nt{i}_base`); the phosphate and sugar shares exist only summed over the whole
hairpin, so a three-way split per nucleotide is not in the data. Those two
whole-hairpin numbers are printed by this script for the caption instead of
being drawn.

Stem / loop / stem brackets under the axis come from the measured pairing in
make_hairpin.py, not from a guess at the fold.

Data: inputs/hits32_geom_12k.csv.gz.
Writes output/Fig4g_hit_epitope/epitope_columns.{svg,pdf,png}.
Nothing is written to final/ - assembly exports there.
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
from common.csy4 import SEQUENCE, LOOP, BRACKETS

FIG = Path(__file__).resolve().parent.name
out = paths.fig_output(FIG)
_cfg = yaml.safe_load((paths.ROOT / "config" / "style.yaml").read_text())
plt.rcParams.update({
    "font.family": _cfg["font_family"], "font.size": _cfg["font_size_pt"],
    "pdf.fonttype": 42, "svg.fonttype": "none", "savefig.dpi": 600,
    "axes.grid": False, "text.color": palette.INK,
    "axes.labelcolor": palette.INK, "axes.edgecolor": palette.INK,
    "xtick.color": palette.INK, "ytick.color": palette.INK,
    "axes.linewidth": _cfg["line_width_pt"],
    "axes.spines.top": False, "axes.spines.right": False,
})
def mm(x): return x / 25.4

INK = palette.INK
LIB_BASE, LIB_BB = palette.NEUTRAL["600"], palette.NEUTRAL["300"]
HIT_BASE, HIT_BB = palette.THEME["700"], palette.THEME["300"]
NBOOT = 4000


def main():
    df = pd.read_csv(paths.INPUTS / "hits32_geom_12k.csv.gz")
    # One design of the 12,000 has no AF3 protein-RNA contact at all, so its
    # shares are 0/0. It is dropped rather than carried as NaN, which silently
    # poisons every mean taken with numpy.
    cols = [f"af3_nt{i}" for i in range(1, 17)]
    keep = df[cols].notna().all(axis=1)
    if (~keep).any():
        print(f"dropped {int((~keep).sum())} design(s) with no AF3 RNA contact")
    df = df[keep].reset_index(drop=True)
    tot = df[cols].values
    base = df[[f"af3_nt{i}_base" for i in range(1, 17)]].values
    bb = tot - base                                    # sugar + phosphate
    is_hit = df.hit.values == 1

    lib_b, lib_k = base[~is_hit].mean(0), bb[~is_hit].mean(0)
    hit_b, hit_k = base[is_hit].mean(0), bb[is_hit].mean(0)

    rng = np.random.default_rng(0)
    idx = rng.integers(0, is_hit.sum(), size=(NBOOT, is_hit.sum()))
    boot = tot[is_hit][idx].mean(axis=1)               # NBOOT x 16
    lo, hi = np.percentile(boot, [2.5, 97.5], axis=0)
    hit_tot = hit_b + hit_k

    x = np.arange(16)
    w = 0.38
    fig, ax = plt.subplots(figsize=(mm(_cfg["width_double_mm"]), mm(52)))

    ax.bar(x - w / 2, lib_k, w, color=LIB_BB, edgecolor="none")
    ax.bar(x - w / 2, lib_b, w, bottom=lib_k, color=LIB_BASE, edgecolor="none")
    ax.bar(x + w / 2, hit_k, w, color=HIT_BB, edgecolor="none")
    ax.bar(x + w / 2, hit_b, w, bottom=hit_k, color=HIT_BASE, edgecolor="none")
    ax.errorbar(x + w / 2, hit_tot, yerr=[hit_tot - lo, hi - hit_tot],
                fmt="none", ecolor=INK, elinewidth=0.6, capsize=1.3,
                capthick=0.6)

    ax.set_xticks(x)
    ax.set_xticklabels([f"{SEQUENCE[i]}{i+1}" for i in range(16)], fontsize=6)
    ax.set_ylabel("Share of contacts")
    ax.set_xlim(-0.75, 15.75)
    ax.tick_params(length=2, width=0.6, pad=1.5)

    # direct labels instead of a legend box (STYLE_NOTES 6), parked over the
    # tallest pair where there is headroom
    top = int(np.argmax(lib_k + lib_b))
    ax.annotate("Library", (top - w / 2, lib_k[top] + lib_b[top]),
                textcoords="offset points", xytext=(-13, 7), fontsize=6.5,
                color=LIB_BASE, ha="center")
    ax.annotate("Hits", (top + w / 2, hi[top]), textcoords="offset points",
                xytext=(13, 5), fontsize=6.5, color=HIT_BASE, ha="center")
    # base / backbone key, frameless and in the panel's own empty corner
    # (STYLE_NOTES 6): the stacking order cannot be direct-labelled on a bar
    # only a few mm wide.
    ymax = max(hi.max(), (lib_b + lib_k).max())
    for k, (lab, col) in enumerate((("base", HIT_BASE), ("backbone", HIT_BB))):
        yk = ymax * (0.92 - 0.11 * k)
        ax.add_patch(plt.Rectangle((12.0, yk), 0.32, ymax * 0.045,
                                   facecolor=col, edgecolor="none"))
        ax.text(12.5, yk + ymax * 0.022, lab, fontsize=6, va="center",
                ha="left", color=INK)

    # stem / loop / stem, from the measured pairing; the three runs are
    # contiguous, so they need a visible break or they read as one rule
    y = -0.030 * ymax / 0.253
    for a, b, lab in BRACKETS:
        ax.plot([a, b], [y, y], color=palette.NEUTRAL["500"], lw=0.8,
                clip_on=False, solid_capstyle="butt")
        ax.text((a + b) / 2, y * 1.9, lab, ha="center", va="top", fontsize=6,
                color=palette.NEUTRAL["600"], clip_on=False)
    ax.set_ylim(0, None)

    fig.tight_layout()
    for ext in ("svg", "pdf", "png"):
        fig.savefig(out / f"epitope_columns.{ext}", bbox_inches="tight")
    print(f"wrote {out}/epitope_columns.svg|pdf|png")

    loop_nt = [i - 1 for i in LOOP]
    print("n: library %d  hits %d" % ((~is_hit).sum(), is_hit.sum()))
    print("loop share of all contacts:  library %.3f   hits %.3f"
          % (tot[~is_hit][:, loop_nt].sum(1).mean(),
             tot[is_hit][:, loop_nt].sum(1).mean()))
    for name, col in (("phosphate", "frac_phos"), ("sugar", "frac_sug"),
                      ("base", "frac_base")):
        print(f"whole-hairpin {name:10s} library {df[~is_hit][f'af3_{col}'].mean():.3f}"
              f"   hits {df[is_hit][f'af3_{col}'].mean():.3f}")
    outside = [i + 1 for i in range(16)
               if not (lo[i] <= (lib_b + lib_k)[i] <= hi[i])]
    print("nt whose library mean falls outside the hits' 95% interval:",
          outside or "none")


if __name__ == "__main__":
    main()
