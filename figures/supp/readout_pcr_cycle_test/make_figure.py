"""Splicing-readout PCR run at 24 and at 29 cycles, on the same ten samples.

Two bare panels, one gel each. Same template set, same primer pair
(PR10148 + PR10137, ~850 bp expected), same master-mix recipe; the only
difference between the two runs is the cycle number, which is why 29 cycles
was adopted for the readout.

Was the left half of panel (c) of the retired `library_build` figure; split
out 2026-11-19 on Leo's request.

Lanes are NOT annotated, deliberately. `LV1.14.xlsx` and `LV1.15.xlsx` each
list ten samples (LV1.13.1 to .10) in lanes 1 to 10, but each capture holds
two gels, one above the other, which neither sheet describes. Assigning lanes
would be a guess.
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

import matplotlib.image as mpimg
import matplotlib.pyplot as plt

from common import paths, style, palette

FIG = Path(__file__).resolve().parent.name
GELS = Path("/Volumes/wicky/users/lvalsecchi/Leo_Valsecchi/LV_Gels")

INK = palette.INK
FS = 5.5
W, H = 85.0, 92.0           # authored panel size, mm; two of these per page

PANELS = [
    ("24 cycles", "gel_hairpinTests_24cycles_inverted.png", "cycles24"),
    ("29 cycles", "gel_hairpinTests_29cycles_day2_inverted.png", "cycles29"),
]
CROP = (320, 60, 1260, 1040)


def show(ax, path):
    """Crop and flip so the wells sit at the top and migration runs down.

    That flip is the only change made to the pixels; the raw captures have the
    wells at the bottom.
    """
    img = mpimg.imread(path)
    x0, y0, x1, y1 = CROP
    ax.imshow(img[y0:y1, x0:x1][::-1], cmap="gray", aspect="equal",
              extent=(x0, x1, y1, y0), interpolation="none")
    ax.set_xlim(x0, x1)
    ax.set_ylim(y1, y0)
    ax.axis("off")


def main():
    style.apply()
    out = paths.fig_output(FIG)
    for lab, fname, stem in PANELS:
        fig = plt.figure(figsize=(style.mm(W), style.mm(H)))
        ax = fig.add_axes([0.02, 0.02, 0.96, 0.90])
        show(ax, GELS / fname)
        # The condition, not a title: it is the one thing that differs between
        # the two panels, so it belongs on the canvas (STYLE_NOTES section 1).
        fig.text(0.02, 0.945, lab, ha="left", va="bottom", fontsize=FS,
                 color=INK)
        style.save(fig, out, stem, tight=False)
        plt.close(fig)
    print("  LV1.14 24 cycles vs LV1.15 29 cycles, ten samples each,")
    print("  PR10148 + PR10137, ~850 bp expected, KAPA HiFi, 50 uL")
    print("wrote", out)


if __name__ == "__main__":
    main()
