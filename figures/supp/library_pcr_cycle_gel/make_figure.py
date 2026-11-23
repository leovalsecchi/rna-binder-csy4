"""Oligo-pool PCR cycle series, and the two bands that were cut out of it.

One bare panel: the LV2.2 gel, three cycle counts pulled from one master mix,
with the band that was gel-extracted and taken forward and the band kept as a
backup each named at the right. Was panel (b) of the retired `library_build`
figure; split out 2026-11-19 on Leo's request, with the two band labels
rewritten so the extracted band and the backup band are told apart.

The identities are the journal's (LV2.4) and the cycle counts the worksheet's
(`LV_PCR/LV2.2_Oligo_RNA_Binders/LV2.2.xlsx`): one master mix in three tubes,
pulled at 16, 20 and 24 cycles. README has the quotes.
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
MUTED = palette.NEUTRAL["600"]
FS = 5.5

W, H = 115.0, 66.0          # authored panel size, mm
GEL = dict(path=GELS / "gel_LV2.2_oligopool_2026-09-30_inverted.png",
           crop=(200, 170, 1060, 1000), well_y=945, band_y=260,
           lanes={"ladder": 288.0, "16": 492.0, "20": 706.0, "24": 928.0})

# What each band was. Journal LV2.4 (Zymoclean): the 20-cycle band was the one
# carried into the BsmBI oligo-pool assembly, the 24-cycle band was cut as well
# and stored. The 16-cycle lane was never cut, which is said on its lane label.
# The taken-forward label sits ABOVE the backup one: both leaders then leave
# their band going up-right, in lane order, and cannot cross.
BANDS = [
    ("20", "Excised, taken forward", "20 cycles \u00b7 31.3 ng/\u00b5L", 0.52),
    ("24", "Excised, kept as backup", "24 cycles \u00b7 87.0 ng/\u00b5L", 0.26),
]
TEXT_X = 0.60               # left edge of the label column, figure fraction


def show(ax, spec):
    """Crop and flip so the wells sit at the top and migration runs down.

    That flip is the only change made to the pixels; the raw capture has the
    wells at the bottom.
    """
    img = mpimg.imread(spec["path"])
    x0, y0, x1, y1 = spec["crop"]
    ax.imshow(img[y0:y1, x0:x1][::-1], cmap="gray", aspect="equal",
              extent=(x0, x1, y1, y0), interpolation="none")
    ax.set_xlim(x0, x1)
    ax.set_ylim(y1, y0)
    ax.axis("off")
    return x0 + x1, y0 + y1


def main():
    style.apply()
    out = paths.fig_output(FIG)
    fig = plt.figure(figsize=(style.mm(W), style.mm(H)))
    ax = fig.add_axes([0.01, 0.03, 0.55, 0.82])
    _, sy = show(ax, GEL)
    wy = sy - GEL["well_y"]
    band_y = sy - GEL["band_y"]

    # Two lines per lane, so the lane names stay on one baseline and the
    # 16-cycle lane can say why it carries no band label. The blank second
    # line is a no-break space, which keeps the block heights equal.
    for lab, lx in GEL["lanes"].items():
        txt = lab if lab == "ladder" else f"{lab} cycles"
        second = "not excised" if lab == "16" else "\u00a0"
        ax.annotate(f"{txt}\n{second}", (lx, wy), xytext=(0, 7),
                    textcoords="offset points", ha="center", va="bottom",
                    fontsize=FS, color=INK, linespacing=1.4,
                    annotation_clip=False,
                    arrowprops=dict(arrowstyle="-", lw=0.4, color=MUTED))

    for lane, head, sub, fy in BANDS:
        lx = GEL["lanes"][lane]
        # The leader stops just LEFT of the label block, between its two
        # lines: run it to the text anchor instead and it strikes through
        # whichever line it arrives at.
        ax.annotate("", (lx, band_y), xytext=(TEXT_X - 0.018, fy),
                    textcoords=fig.transFigure, annotation_clip=False,
                    arrowprops=dict(arrowstyle="-", lw=0.4, color=MUTED,
                                    shrinkA=0, shrinkB=0))
        fig.text(TEXT_X, fy + 0.012, head, ha="left", va="bottom",
                 fontsize=FS, color=INK)
        fig.text(TEXT_X, fy - 0.012, sub, ha="left", va="top",
                 fontsize=FS, color=MUTED)

    style.save(fig, out, "pcr_cycle_gel", tight=False)
    plt.close(fig)
    print(f"  one master mix, three tubes: 16 / 20 / 24 cycles, 302 bp expected")
    print(f"  taken forward: 20 cycles, 31.3 ng/uL;  backup: 24 cycles, 87.0 ng/uL")
    print("wrote", out)


if __name__ == "__main__":
    main()
