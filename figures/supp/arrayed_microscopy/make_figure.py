"""Arrayed microscopy of the 32 ordered candidates, LV2.34 (48 h, LV2.33).

Three bare panels from the crops `make_inputs.py` caches:

  a  all three channels merged, across all 140 wells - 4 targets x
     (32 candidates + 3 controls)
  b  the same fields in mScarlet alone - the readout, which the merge buries
     under mNeonGreen because nearly every transfected cell carries it
  c  the same six wells again at full resolution, each channel and the merge

(a) and (b) are the same field of the same well. (b) exists because the merge
cannot carry the result on its own: mNeonGreen marks every transfected cell
and saturates, so a three-colour tile reads green whatever mScarlet is doing.
The transfection control is the blue channel of (a) and the third column of (c).

Every field is the 480 x 480 px (311.7 um) window of its well holding the most
mTagBFP edge, with the middle of the well excluded outright - see
make_inputs.py, which also records why a centre crop and a whole-well tile were
both tried and both discarded. Placement never consults mScarlet or
mNeonGreen: both are reporters, so scoring on either would let the rule reach
for a field because of the readout.

Channel colours come from the palette, never from here: mScarlet is
`synbio.reporter_alt` and mNeonGreen is `synbio.reporter`, both matched to the
real emission so these panels and the construct maps agree (COLORS.md rule 8).
The palette has no mTagBFP entry, so the blue channel takes theme step 400,
the lightest theme blue that still reads as blue against black.

One display range per channel per set, as `inputs/micrographs/README.md`
requires: scaling each tile to its own maximum would put the result in the
contrast rather than in the cells. Ranges are printed on each run and recorded
in this figure's README.

    python3 code/supp/arrayed_microscopy/make_figure.py
"""
import csv
import sys
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from common import paths, style, palette                  # noqa: E402

FIG = Path(__file__).resolve().parent.name
IN = paths.INPUTS / "micrographs" / "arrayed"

INK = palette.INK
FS = 5.5
GAP_MM = 0.25             # hairline between tiles

UM_PX = 0.649403          # inputs/micrographs/arrayed/acquisition.txt
FIELD_PX = 480
FIELD_UM = FIELD_PX * UM_PX

MSCARLET, MNEON, BFP = 0, 2, 1            # component order in the ND2 is R,B,G
CHANNELS = [
    ("mScarlet", MSCARLET, palette.SYNBIO["reporter_alt"]),
    ("mNeonGreen", MNEON, palette.SYNBIO["reporter"]),
    ("mTagBFP", BFP, palette.THEME["400"]),
]

ROWS = ["A", "B", "C", "D"]
ROW_LABEL = {"A": "posA", "B": "posB", "C": "NegA", "D": "no array"}
CONTROL_COLS = [("P6_C3", "untransf."), ("P6_C4", "Csy4 1"), ("P6_C5", "Csy4 2")]
CLOSEUPS = [("P6_C4", "A"), ("P6_C5", "A"), ("P6_C3", "A"),
            ("P4_C5", "A"), ("P1_C4", "A"), ("P4_C5", "D")]


def variant_columns():
    """P<plate>_C<col> for V1..V32, in variant order."""
    return [(f"P{(v - 1) // 6 + 1}_C{(v - 1) % 6 + 1}", str(v))
            for v in range(1, 33)]


def load():
    tiles = dict(np.load(IN / "tiles.npz"))
    close = dict(np.load(IN / "closeups.npz"))
    wells = {(r["binder"], r["row"]): r
             for r in csv.DictReader((IN / "wells.csv").open())}
    return tiles, close, wells


def display_range(arrays, lo_pct=25, hi_pct=99.9):
    """One (lo, hi) for a whole set, from its pooled pixels.

    The low percentile is where that set stops being background: 25th for the
    3x block-averaged tiles, whose averaging lifts the floor, and lower for the
    full-resolution close-ups, where background really is the bulk of the
    frame. The counts are listed in this figure's README.
    """
    pool = np.concatenate([a.ravel() for a in arrays])
    return float(np.percentile(pool, lo_pct)), float(np.percentile(pool, hi_pct))


def grey(a, lo, hi):
    return np.clip((a.astype(np.float32) - lo) / max(hi - lo, 1.0), 0, 1)


def tint(g, hexcol):
    rgb = np.array([int(hexcol[i:i + 2], 16) / 255 for i in (1, 3, 5)])
    return g[..., None] * rgb


HUE = {ci: col for _, ci, col in CHANNELS}


def single(a, ci, ranges):
    """One channel in its own emission colour, not grey.

    Grey is the usual convention for a single channel, but it is not what was
    asked for here and it throws away the one cue that ties a panel to its
    channel without reading the label.
    """
    return tint(grey(a[:, :, ci], *ranges[ci]), HUE[ci])


def merge(a, ranges):
    """Additive three-colour merge, each channel on its own display range."""
    out = np.zeros(a.shape[:2] + (3,), np.float32)
    for _, ci, col in CHANNELS:
        out += tint(grey(a[:, :, ci], *ranges[ci]), col)
    return np.clip(out, 0, 1)


def scalebar(ax, field_um, npx, bar_um=100, text=None):
    """Scale bar, bottom right, in image pixel coordinates."""
    px = bar_um / field_um * npx
    pad = 0.06 * npx
    h = max(0.022 * npx, 1.0)
    ax.add_patch(Rectangle((npx - pad - px, npx - pad - h), px, h,
                           facecolor="white", edgecolor="none", zorder=5))
    if text:
        ax.text(npx - pad - px / 2, npx - pad - h - 0.03 * npx, text,
                color="white", fontsize=FS, ha="center", va="bottom", zorder=5)


def bare(ax):
    ax.set_xticks([]), ax.set_yticks([])
    for s in ax.spines.values():
        s.set_visible(False)


def save(fig, outdir, stem):
    """600 dpi, not the house 300, and SVG for page assembly.

    A vector save resamples every imshow raster to the figure's dpi, so a tile
    drawn 4.6 mm wide from 160 px would be thrown away at 300 dpi (54 px). See
    the matplotlib vector-image-dpi note.
    """
    for ext in ("pdf", "png", "svg"):
        fig.savefig(outdir / f"{stem}.{ext}", bbox_inches="tight", dpi=600)


# ---------------------------------------------------------------- panels a, b
def array_panel(tiles, imgfn, outdir, stem):
    cols = variant_columns() + [(None, None)] + list(CONTROL_COLS)
    ncol, nrow = len(cols), len(ROWS)
    tile_mm, gut_mm, foot_mm, head_mm = 4.55, 13.0, 12.0, 1.0
    w = gut_mm + ncol * tile_mm
    h = head_mm + nrow * tile_mm + foot_mm
    fig = plt.figure(figsize=(style.mm(w), style.mm(h)))
    for i, row in enumerate(ROWS):
        for j, (binder, label) in enumerate(cols):
            if binder is None or f"{binder}_{row}" not in tiles:
                continue
            # A hairline gap on every side. Butted tiles antialias unevenly -
            # every second boundary picks up a white line and the array reads
            # as pairs of wells rather than a grid.
            ax = fig.add_axes([(gut_mm + j * tile_mm) / w,
                               (foot_mm + (nrow - 1 - i) * tile_mm) / h,
                               (tile_mm - GAP_MM) / w, (tile_mm - GAP_MM) / h])
            img = imgfn(tiles[f"{binder}_{row}"])
            ax.imshow(img, cmap=None if img.ndim == 3 else "gray",
                      vmin=None if img.ndim == 3 else 0,
                      vmax=None if img.ndim == 3 else 1,
                      interpolation="nearest")
            bare(ax)
            if i == 0 and j == 0:
                scalebar(ax, FIELD_UM, img.shape[0])
            if i == nrow - 1:
                ax.set_xlabel(label, fontsize=FS, color=INK, labelpad=1.5,
                              rotation=90 if len(label) > 2 else 0,
                              ha="center", va="top")
            if j == 0:
                ax.set_ylabel(ROW_LABEL[row], fontsize=FS, color=INK,
                              rotation=0, ha="right", va="center", labelpad=4)
    outdir.mkdir(parents=True, exist_ok=True)
    save(fig, outdir, stem)
    plt.close(fig)


# ------------------------------------------------------------------ panel c
def closeup_panel(close, wells, ranges, outdir, stem):
    short = dict(CONTROL_COLS)          # the long control names do not fit
    labels = []
    for binder, row in CLOSEUPS:
        rec = wells[(binder, row)]
        name = short.get(binder, rec["variant"])
        fc = rec["onoff_avg"]
        tag = name if not fc else f"{name}  FC {float(fc):.2f}"
        labels.append(f"{tag}\n{ROW_LABEL[row]}")

    cols = [(n, ci, col) for n, ci, col in CHANNELS] + [("merge", None, None)]
    nrow = len(CLOSEUPS)
    tile_mm, gut_mm, hdr_mm = 26.0, 24.0, 4.0
    w = gut_mm + len(cols) * tile_mm
    h = hdr_mm + nrow * tile_mm
    fig = plt.figure(figsize=(style.mm(w), style.mm(h)))
    for i, (binder, row) in enumerate(CLOSEUPS):
        a = close[f"{binder}_{row}"]
        for j, (cname, ci, _) in enumerate(cols):
            ax = fig.add_axes([(gut_mm + j * tile_mm) / w,
                               (nrow - 1 - i) * tile_mm / h,
                               (tile_mm - GAP_MM) / w, (tile_mm - GAP_MM) / h])
            img = merge(a, ranges) if ci is None else single(a, ci, ranges)
            ax.imshow(img, interpolation="nearest")
            bare(ax)
            if i == 0:
                ax.set_title(cname, fontsize=FS, color=INK, pad=2)
            if i == 0 and j == 0:
                scalebar(ax, FIELD_UM, a.shape[0], 100, "100 µm")
            if j == 0:
                ax.set_ylabel(labels[i], fontsize=FS, color=INK, rotation=0,
                              ha="right", va="center", labelpad=4,
                              linespacing=1.4)
    outdir.mkdir(parents=True, exist_ok=True)
    save(fig, outdir, stem)
    plt.close(fig)


def main():
    style.apply()
    tiles, close, wells = load()
    outdir = paths.fig_output(FIG)
    got = next(iter(tiles.values())).shape
    print(f"{len(tiles)} tiles {got}, {len(close)} close-ups "
          f"{next(iter(close.values())).shape}, field {FIELD_UM:.1f} um")

    ranges, close_ranges = {}, {}
    for name, ci, _ in CHANNELS:
        ranges[ci] = display_range([a[:, :, ci] for a in tiles.values()], 25)
        close_ranges[ci] = display_range([a[:, :, ci] for a in close.values()], 1)
        print(f"display range {name:<11} tiles {ranges[ci][0]:7.0f} -"
              f" {ranges[ci][1]:7.0f}   close-ups {close_ranges[ci][0]:7.0f} -"
              f" {close_ranges[ci][1]:7.0f}")

    array_panel(tiles, lambda a: merge(a, ranges), outdir, "panel_a_merge_array")
    array_panel(tiles, lambda a: single(a, MSCARLET, ranges),
                outdir, "panel_b_mscarlet_array")
    closeup_panel(close, wells, close_ranges, outdir, "panel_c_controls")
    print(f"-> {outdir}")


if __name__ == "__main__":
    main()
