"""Fig4d: screen hit rate per stack = ordered hits / designs in the 12k library.

One horizontal bar per stack, longest (best rate) at the top. Sorted by value
on Leo's request (2026-11-12), which overrides the fixed STACK_ORDER rule
(STYLE_NOTES 9) for this panel only; colour still follows the stack. Bars are
annotated with the raw counts (STYLE_NOTES 5), the axis carries the rate.

Inputs: inputs/csy4_binder_library_12k_v6_k30.csv (the 12,000-design library)
and inputs/candidates32_mapped.csv (the 32 ordered hits).
"""
import sys, csv
from collections import Counter
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import matplotlib.pyplot as plt
from common import paths, style, palette

FIG = Path(__file__).resolve().parent.name
style.apply()
out = paths.fig_output(FIG)


def stack_of(name):
    for short in palette.STACKS_SHORT:
        if f"_{short.replace('+', '_')}_" in name:
            return short
    raise ValueError(name)


def count(fname):
    with open(paths.INPUTS / fname) as f:
        return Counter(stack_of(r["design_id"]) for r in csv.DictReader(f))


lib = count("csy4_binder_library_12k_v6_k30.csv")
hit = count("candidates32_mapped.csv")
assert sum(lib.values()) == 12000 and sum(hit.values()) == 32, (lib, hit)

rows = sorted(lib, key=lambda s: hit[s] / lib[s])        # ascending: barh row 0 is the bottom
rate = [100 * hit[s] / lib[s] for s in rows]

def draw(size_mm=(55, 65), fs=6.5, outdir=out):
    """The panel. `size_mm`/`fs` exist so the Fig 4 assembly can redraw it at
    the size it actually lands on the page with its text already at the page
    size - scaled into its slot instead, the text re-set to T grows against
    the bars and the longest stack name is clipped off the left edge
    (2026-11-13, same failure as Fig2's gpu panel)."""
    fig, ax = plt.subplots(figsize=(style.mm(size_mm[0]), style.mm(size_mm[1])))
    ax.barh(rows, rate, height=0.68,
            color=[palette.stack_color(s) for s in rows], linewidth=0)
    for y, s in enumerate(rows):
        ax.text(rate[y] + 0.006, y, f"{hit[s]}/{lib[s]:,}", va="center",
                ha="left", fontsize=fs)
    ax.set_xlim(0, 0.45)
    ax.set_xticks([0, 0.1, 0.2, 0.3, 0.4])
    ax.set_xlabel("Hit rate (%)", fontsize=fs)
    ax.tick_params(axis="y", length=0, labelsize=fs)
    ax.tick_params(axis="x", labelsize=fs)
    ax.spines["left"].set_visible(False)
    outdir.mkdir(parents=True, exist_ok=True)
    fig.savefig(outdir / "Fig4d_stack_hit_rate.svg", bbox_inches="tight")
    if outdir == out:
        style.save(fig, out, "Fig4d_stack_hit_rate")
    plt.close(fig)


if __name__ == "__main__" and len(sys.argv) == 5:
    # W_MM H_MM FONT_PT OUTDIR, for page_assembly.variant()
    draw((float(sys.argv[1]), float(sys.argv[2])), float(sys.argv[3]),
         Path(sys.argv[4]))
elif __name__ == "__main__":
    draw()
    print({s: (hit[s], lib[s], round(100 * hit[s] / lib[s], 3)) for s in rows[::-1]})
