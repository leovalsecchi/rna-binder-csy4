"""Standalone scaffold glyphs as SVGs, no titles, for panel A assembly.

Four files: each scaffold as a positive (Csy4 hairpin tips in the site colour)
and as a negative (grey tips, same geometry). F30 variation uses the short stem.
"""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import matplotlib.pyplot as plt
from common import style, paths, palette, scaffold_glyphs

style.apply()
TEAL, PLUM = palette.THEME["700"], palette.THEME["400"]
POS, NEG = palette.DIVERGING[-1], palette.NEUTRAL["500"]   # hairpin tip colours

out = paths.fig_output("Fig1f_scaffold_splicing") / "elements"
out.mkdir(parents=True, exist_ok=True)

def panel(stem, W, H, draw):
    fig = plt.figure(figsize=(style.mm(W), style.mm(H)))
    ax = fig.add_axes([0, 0, 1, 1]); ax.set_xlim(0, W); ax.set_ylim(H, 0); ax.axis("off")
    draw(ax)
    fig.savefig(out / f"{stem}.svg", bbox_inches="tight", transparent=True)
    plt.close(fig)
    print("wrote", out / f"{stem}.svg")

for tag, tip in (("positive", POS), ("negative", NEG)):
    panel(f"f30_variation_{tag}", 26.0, 16.0,
          # Short 5'/3' overhangs: the strand spans the four junctions
          # (13 +/- 7.5) plus 1.1 either side, not the full panel width.
          lambda ax, t=tip: scaffold_glyphs.array(ax, 13.0, 18.5, TEAL, t, 4.4, 21.6))
    panel(f"circular_scaffold_{tag}", 22.0, 22.0,
          # The ring meets the transcript at one point, so its strand is short.
          lambda ax, t=tip: scaffold_glyphs.ring(ax, 11.0, 22.0, PLUM, t, 7.3, 14.7))
