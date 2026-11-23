"""The two readout gels that were run but whose captures are not archived.

Two bare panels, one per missing capture. Each draws a dashed empty frame in
place of the image and states what the written record still gives for that
gel, so the gap is visible rather than silently cropped out of the set.

Was the right half of panel (c) of the retired `library_build` figure; split
out 2026-11-19 on Leo's request. The frames used to say only "image not on
the NAS"; they now carry the worksheet and journal facts that survive without
the picture.

`LV_Gels/` holds four captures in total (the primer test, the LV2.2 oligo
pool, and the two hairpin-test gels). Neither of these two is among them.
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

from common import paths, style, palette

FIG = Path(__file__).resolve().parent.name

INK = palette.INK
MUTED = palette.NEUTRAL["600"]
LINEG = palette.NEUTRAL["400"]
FS = 5.5
W, H = 85.0, 46.0

# Each entry: file stem, the gel, and the record that survives it.
# LV1.27 from LV_PCR/LV1.27.xlsx; LV2.18 from journal LV2.18 (13.10) and the
# Nanodrop readings recorded with it.
RECORDS = [
    ("csy4_retest", "LV1.27, Csy4 retest readout PCR", [
        "Capture not on the NAS.",
        "",
        "LV1.27.xlsx: five reactions (LV1.27.1 to .5),",
        "KAPA HiFi 50 µL, PR10148 + PR10137,",
        "~850 bp expected, 24 cycles.",
        "Lanes 1 to 4 and 10.",
    ]),
    ("double_band", "LV2.18, the unexplained double band", [
        "Capture not on the NAS.",
        "",
        "Journal LV2.18: 1.25% gel, 250 mL TAE,",
        "25 µL SYBR Safe. Eight samples —",
        "pLV009 to pLV012 at 24 and at 29 cycles.",
        "Both bands extracted for 8/8 samples (100%).",
        "24-cycle yields 87.3, 61.9, 84.0, 67.5 ng/µL,",
        "all noted “purity not optimal”; 29-cycle set",
        "not measured.",
    ]),
]


def main():
    style.apply()
    out = paths.fig_output(FIG)
    for stem, head, lines in RECORDS:
        fig = plt.figure(figsize=(style.mm(W), style.mm(H)))
        ax = fig.add_axes([0.02, 0.02, 0.96, 0.90])
        ax.set_xlim(0, 1)
        ax.set_ylim(0, 1)
        ax.axis("off")
        ax.add_patch(Rectangle((0, 0), 1, 1, facecolor="none", edgecolor=LINEG,
                               lw=0.5, linestyle=(0, (3, 3))))
        ax.text(0.5, 0.5, "\n".join(lines), ha="center", va="center",
                fontsize=FS, color=MUTED, linespacing=1.55)
        # Which gel the empty frame stands for: the panel's only ink in ink.
        fig.text(0.02, 0.945, head, ha="left", va="bottom", fontsize=FS,
                 color=INK)
        style.save(fig, out, stem, tight=False)
        plt.close(fig)
    print("  two captures missing: LV1.27 (5 reactions) and LV2.18 (8 samples)")
    print("  LV_Gels/ holds 4 captures, neither of these among them")
    print("wrote", out)


if __name__ == "__main__":
    main()
