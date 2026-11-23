#!/usr/bin/env python3
"""Fig 3c: which stack produced each of the 32 ordered candidates.

Source of truth is the ordered-hit set itself,
working/structures/hits32/index.tsv (32 rows, one per ordered design),
copied from LV2.1.1 12_analysis/hits32_browse. Slice order and colour are the
fixed thesis stack order, so this panel reads against every other stack panel.
"""
import sys
from collections import Counter
from pathlib import Path

import matplotlib.pyplot as plt

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from common import palette, style                               # noqa: E402

HITS = Path(__file__).resolve().parents[2] / "working/structures/hits32"
OUT = Path(__file__).resolve().parents[2] / "output/Fig4ef_ordered32"


def counts():
    rows = HITS.joinpath("index.tsv").read_text().splitlines()[1:]
    names = [r.split("\t")[0] for r in rows if r.strip()]
    c = Counter(palette.SHORT[k] if k in palette.SHORT else k
                for k in (stack_of(n) for n in names))
    return c, len(names)


def stack_of(name):
    for short in palette.STACKS_SHORT:
        key = short.replace("+", "_")
        if f"_{key}_" in name:
            return short
    raise ValueError(name)


def main():
    style.apply()
    c, n = counts()
    assert n == 32, f"expected 32 ordered designs, found {n}"
    labels = [s for s in (palette.SHORT[k] for k in palette.STACK_ORDER) if c[s]]
    vals = [c[s] for s in labels]
    cols = [palette.stack_color(s) for s in labels]

    fig, ax = plt.subplots(figsize=(style.mm(78), style.mm(62)))
    wedges, _ = ax.pie(vals, colors=cols, startangle=90, counterclock=False,
                       wedgeprops=dict(edgecolor="white", linewidth=1.2))
    ax.set_aspect("equal")

    # direct labels on leader lines, no legend: the reader never has to
    # match a colour swatch to a name.
    for w, lab, v in zip(wedges, labels, vals):
        a = (w.theta1 + w.theta2) / 2
        import math
        x, y = math.cos(math.radians(a)), math.sin(math.radians(a))
        side = 1 if x >= 0 else -1
        ax.annotate(f"{lab}\n{v} of {n}  ({100*v/n:.0f}%)",
                    xy=(0.82 * x, 0.82 * y), xytext=(1.28 * side, 1.18 * y),
                    ha="left" if side > 0 else "right", va="center",
                    fontsize=6.5, color=palette.INK,
                    arrowprops=dict(arrowstyle="-", color=palette.NEUTRAL["500"],
                                    lw=0.6, shrinkA=0, shrinkB=2,
                                    connectionstyle="arc3,rad=0"))
    ax.set_xlim(-2.0, 2.0)
    OUT.mkdir(parents=True, exist_ok=True)
    style.save(fig, OUT, "Fig4e_stack_pie")
    print({s: c[s] for s in labels}, "total", n)


if __name__ == "__main__":
    main()
