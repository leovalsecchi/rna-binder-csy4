"""Cut the LV2.34 arrayed micrographs down to what the figure draws.

Reads the six ND2 plates on the NAS (~45 GB, read-only) and writes ~25 MB into
`inputs/micrographs/arrayed/`, so `make_figure.py` never touches the NAS.

Per well, and no field anybody picked:

* a **whole-well** map - the full 7344 x 7344 stitched scan, sampled every 16th
  pixel and block-averaged to 153 px. Nothing cropped.
* the 480 x 480 px (311.7 um) field in that well holding the **most cells**,
  at least EXCLUDE_MM from the well centre and under ARTEFACT_CAP, scored as
  mTagBFP edge over the whole-well map (see `densest_field`). Scoring on the
  transfection marker and never on mScarlet is the point: the rule cannot
  reach for a field because of the readout. This is the figure's tile, and the
  six close-up wells keep theirs at full resolution.

Two earlier versions of this were wrong and are worth not repeating. **A centre
crop** falls into the cell-free hole, a millimetre or more across, that sits in
the middle of many of these wells - plate 4 C6 is the clearest - and scored
0.0000 bright-mTagBFP area in wells the flow run puts at 43% BFP-positive. **A
whole-well tile** avoids that but cannot resolve a cell: at 153 px one pixel is
31 um, so a block mean averages the sparse mScarlet cells away and a block max
turns every tile into speckle. A density-chosen field is the one that both
misses the hole and still shows cells.

    python3 code/supp/arrayed_microscopy/make_inputs.py
"""
import csv
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
sys.path.insert(0, str(Path(__file__).resolve().parent))

from common import paths                                  # noqa: E402
import nd2lite                                            # noqa: E402

ND2_DIR = Path("/Volumes/wicky/expdata/Microscopy")
OUT = paths.INPUTS / "micrographs" / "arrayed"

CROP = 480          # field, px -> 311.7 um at 0.649403 um/px
TILE = 160          # montage tile: CROP block-averaged 3x
STEP = 16           # whole-well subsampling: 7344 -> 459 px
BIN = 3             # 459 -> 153 px whole-well map
BFP = 1             # component index of 390_mTagBFP
UM_PER_PX = 0.649403                      # asserted against the file below

# No field may be centred within this of the well centre. Measured, not
# guessed: over the 136 transfected wells, the probability that a point is
# cell-free (mTagBFP below a quarter of its own well's 90th percentile) is
# 0.10 at the dead centre, climbs to a peak of 0.60 at r = 0.75 mm, and falls
# back to its outer-field floor of ~0.29 only beyond r = 1.5 mm. 1.2 mm is
# where the risk has come back down to within a few points of that floor,
# while still leaving the whole 1.2-2.2 mm annulus of a 4.77 mm scan to
# choose from. The clearest single example is plate 4 well C6.
EXCLUDE_MM = 1.2

# A field is rejected if its mean mTagBFP exceeds this multiple of the well's
# own median. Measured, not guessed: with the exclusion disc but no cap, the
# four untransfected wells take the four highest field/well median ratios of
# all 140 - 1.78, 2.57, 2.80, 2.88 - against a median of 1.22 across the set
# and a maximum of 1.55 among the 136 transfected wells. They have no cells,
# so the only thing an edge scorer can find in them is a flare, and they came
# out looking brightly transfected. At 1.4 every untransfected field lands at
# or below its own well median (0.91, 0.94, 1.00, 1.15) and the six
# transfected wells the cap does constrain simply take their next-densest
# field. No well is left without a candidate.
ARTEFACT_CAP = 1.4

# Below this fraction of cell-scale texture a well is taken to hold no cells.
# Measured: the four untransfected wells are the four lowest of all 140, at
# 0.0252, 0.0307, 0.0309 and 0.0315, and the lowest transfected well is 0.0730
# - more than double. 0.05 sits in the middle of that gap.
CELL_TEXTURE_MIN = 0.05

# Cell-plate rows are the four targets (LV2.33 as-performed, "Cell plate maps").
# The flow table's column names are in brackets.
TARGETS = {
    "A": ("posA", "pLV13"),         # PosA
    "B": ("posB", "pLV14"),         # PosB
    "C": ("NegA", "pLV16"),         # NegA
    "D": ("no array", "pLV17"),  # FullNeg
}

# Plate/column -> what is in the well. Columns 1-6 of plates 1-5 are V1..V30 and
# plate 6 columns 1-2 are V31/V32, exactly as planned. Plate 6 columns 3-5 are
# NOT what the plan drew - see README, "Where the files contradict the plan".
CONTROLS = {
    (6, 3): "untransfected",
    (6, 4): "Csy4 (plate 6 col 4)",
    (6, 5): "Csy4 (plate 6 col 5)",
}

# Close-ups: the two ends of the fold-change range plus every control, all on
# the posA target, and the top design again on the no-scaffold target.
CLOSEUPS = [
    ("P4_C5", "A"), ("P1_C4", "A"),
    ("P6_C4", "A"), ("P6_C5", "A"), ("P6_C3", "A"),
    ("P4_C5", "D"),
]


def _box(a, r):
    """Mean over a (2r+1) square, edge-padded, via a summed-area table."""
    n = a.shape[0]
    c = np.pad(a, ((r, r), (r, r)), mode="edge").cumsum(0).cumsum(1)
    c = np.pad(c, ((1, 0), (1, 0)))
    k = 2 * r + 1
    return (c[k:k + n, k:k + n] - c[0:n, k:k + n]
            - c[k:k + n, 0:n] + c[0:n, 0:n]) / (k * k)


def cell_texture(a):
    """Fraction of map pixels standing >20% above their own neighbourhood.

    One map pixel is ~31 um, about one cell, so cells are point-like here and
    a flare or a bubble is not: this separates "has cells" from "has a bright
    artefact", which brightness alone cannot.
    """
    return float((a > 1.2 * _box(a, 5)).mean())


def densest_field(bfp_coarse, crop, step):
    """Top-left of the `crop`-px window holding the most mTagBFP *edge*.

    Scored on the gradient magnitude, not on brightness, and summed with an
    area table so it is the same deterministic rule in every well. mScarlet is
    never consulted, so the rule cannot reach for a field because of the
    readout.

    Brightness was the first scorer and it was wrong: several of these wells
    carry a saturated artefact a few hundred microns across - a bubble or a
    fluorescent speck - and a solid bright disc beats any field of cells on
    total brightness. It won in plate 2 B6 and filled that tile with a white
    blob. A disc is flat, so it has edge only at its rim, while a field of
    cells has edge at every cell boundary. Checked by eye against the
    brightness scorer on six wells: edge lands in cell texture in all six and
    brightness lands on the artefact in four.
    """
    a = bfp_coarse.astype(np.float64)
    grad = (np.abs(np.diff(a, axis=1, prepend=a[:, :1]))
            + np.abs(np.diff(a, axis=0, prepend=a[:1, :])))
    win = max(crop // step, 1)
    s = np.pad(grad.cumsum(0).cumsum(1), ((1, 0), (1, 0)))
    tot = (s[win:, win:] - s[:-win, win:] - s[win:, :-win] + s[:-win, :-win])

    # Keep off the middle of the well. The cell-free region there is not a
    # thing the scorer should have to dodge by luck, so it is excluded
    # outright: any window whose centre falls inside EXCLUDE_MM of the well
    # centre is removed before the argmax. See EXCLUDE_MM for the radius.
    n = bfp_coarse.shape[0]
    cy = cx = (n - 1) / 2.0

    # A well with no cells has no densest field, and an edge scorer turned
    # loose on one finds the brightest flare near the wall - which is how the
    # untransfected control first came out looking brightly transfected. Test
    # for cells first, and if there are none take a fixed position instead:
    # nothing in such a well is more representative than anything else, so the
    # honest choice is a stated constant, not an argmax over artefacts.
    if cell_texture(a) < CELL_TEXTURE_MIN:
        off = int(round(EXCLUDE_MM * 1000.0 / UM_PER_PX / step)) + win // 2
        return int((cy - win / 2) * step), int((cx + off - win / 2) * step)

    yy, xx = np.mgrid[0:tot.shape[0], 0:tot.shape[1]]
    rad = np.hypot(yy + win / 2 - cy, xx + win / 2 - cx) * step * UM_PER_PX / 1000.0

    # A window whose mean mTagBFP is far above the well's own median is a
    # broad bright artefact - a flare, a bubble, a speck - not a dense field
    # of cells, and in a well with NO cells it is the only thing an edge
    # scorer can find. See ARTEFACT_CAP.
    s2 = np.pad(a.cumsum(0).cumsum(1), ((1, 0), (1, 0)))
    mean = (s2[win:, win:] - s2[:-win, win:]
            - s2[win:, :-win] + s2[:-win, :-win]) / (win * win)
    # Also keep off any unscanned corner of the stitch: a window holding
    # zeros puts a black bar in the panel and drags the whole set's display
    # range down to 0, which washes out every other tile.
    z = np.pad((a <= 0).astype(np.float64).cumsum(0).cumsum(1), ((1, 0), (1, 0)))
    zeros = (z[win:, win:] - z[:-win, win:] - z[win:, :-win] + z[:-win, :-win])
    ok = ((rad >= EXCLUDE_MM) & (mean <= ARTEFACT_CAP * np.median(a))
          & (zeros == 0))
    tot = np.where(ok, tot, -np.inf)
    if not np.isfinite(tot).any():
        raise ValueError("exclusion radius and artefact cap leave no field")
    iy, ix = np.unravel_index(int(np.argmax(tot)), tot.shape)
    return iy * step, ix * step


def plate_files():
    out = {}
    for p in range(1, 7):
        hits = sorted(ND2_DIR.glob(f"LV2_34_{p}_Plate{p}_*.nd2"))
        if len(hits) != 1:
            raise SystemExit(f"plate {p}: expected 1 ND2, found {len(hits)}")
        out[p] = hits[0]
    return out


def fold_change():
    rows = csv.DictReader((paths.INPUTS / "lv2_34_1_flow_on_off.csv").open())
    # "stack" is the upstream CSV's column name; it is a generator plus an
    # inverse folder, which is a design ROUTE. Written out under that name.
    return {r["binder"]: (r["ONOFF_avg"], r["design_id"], r["stack"]) for r in rows}


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    # run_all.sh runs every make_*.py, and this one reads 4 GB off the NAS at
    # ~11 MB/s. Skip when the cache is already newer than the plates.
    cached = [OUT / f for f in ("tiles.npz", "wholewell.npz", "closeups.npz",
                                "wells.csv", "acquisition.txt")]
    if "--force" not in sys.argv and all(p.exists() for p in cached):
        newest = max(f.stat().st_mtime for f in plate_files().values())
        if min(p.stat().st_mtime for p in cached) > newest:
            print(f"cache is current, skipping ({OUT}); --force to rebuild")
            return
    # The whole-well maps are an input to field selection, not an output of
    # it, and they cost one strided read per well - the slowest thing this
    # script does, and unusable when the NAS is busy. Reuse them whenever they
    # are newer than the plates; --rebuild-maps forces them.
    newest = max(f.stat().st_mtime for f in plate_files().values())
    wwp = OUT / "wholewell.npz"
    reuse = {}
    if "--rebuild-maps" not in sys.argv and wwp.exists() \
            and wwp.stat().st_mtime > newest:
        reuse = dict(np.load(wwp))
        print(f"reusing {len(reuse)} whole-well maps from {wwp.name}")

    fc = fold_change()
    want = {k: None for k in CLOSEUPS}
    tiles, whole, meta, fields = {}, {}, [], {}

    for plate, path in sorted(plate_files().items()):
        n = nd2lite.ND2(path)
        assert abs(n.um_per_px - UM_PER_PX) < 1e-6, (
            f"{path.name} is {n.um_per_px} um/px, not the {UM_PER_PX} the "
            "exclusion radius and every stated field width assume")
        for frame, well in enumerate(n.frame_wells()):
            row, col = well[0], int(well[1:])
            binder = f"P{plate}_C{col}"
            variant = CONTROLS.get((plate, col))
            if variant is None:
                variant = f"V{(plate - 1) * 6 + col}"
            target, plasmid = TARGETS[row]
            key = f"{binder}_{row}"
            if key in reuse:
                whole[key] = reuse[key]
            else:
                coarse = n.coarse(frame, STEP)
                k = coarse.shape[0] // BIN * BIN
                whole[key] = (coarse[:k, :k]
                              .reshape(k // BIN, BIN, k // BIN, BIN, n.comp)
                              .mean(axis=(1, 3)).astype(np.uint16))
            # Always from the STORED map, never from the finer `coarse`, so a
            # rebuild that reuses the cache places fields identically to one
            # that does not.
            map_step = n.height // whole[key].shape[0]
            y0, x0 = densest_field(whole[key][:, :, BFP], CROP, map_step)
            field = n.crop(frame, CROP, y0, x0)
            b = CROP // TILE
            tiles[key] = (field.reshape(TILE, b, TILE, b, n.comp)
                               .mean(axis=(1, 3)).astype(np.uint16))
            fields[key] = (y0, x0)
            if (binder, row) in want:
                want[(binder, row)] = field
            onoff, design_id, route = fc.get(binder, ("", "", ""))
            meta.append(dict(plate=plate, well=well, row=row, col=col,
                             binder=binder, variant=variant, target=target,
                             target_plasmid=plasmid, onoff_avg=onoff,
                             design_id=design_id, design_route=route,
                             field_y=y0, field_x=x0))
            print(f"  plate {plate} {well:>3}  {variant:<22} {target}", flush=True)

    missing = [k for k, v in want.items() if v is None]
    if missing:
        raise SystemExit(f"close-up wells not found: {missing}")

    np.savez_compressed(OUT / "tiles.npz", **tiles)
    np.savez_compressed(OUT / "wholewell.npz", **whole)
    np.savez_compressed(OUT / "closeups.npz",
                        **{f"{b}_{r}": want[(b, r)] for b, r in CLOSEUPS})
    with (OUT / "wells.csv").open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(meta[0]))
        w.writeheader()
        w.writerows(meta)

    n = nd2lite.ND2(plate_files()[1])
    (OUT / "acquisition.txt").write_text(
        f"source      {ND2_DIR}/LV2_34_<n>_Plate<n>_HEKwt_32RNAcandidateBinders_48h.nd2\n"
        f"acquired    2026-11-06, 48 h post transfection (LV2.33, 2026-11-04)\n"
        f"objective   Nikon Plan Apo lambdaD 10x OFN25, NA 0.45, Ti2 widefield\n"
        f"camera      Photometrics Kinetix, 1x1 binning, 16-bit HDR\n"
        f"calibration {n.um_per_px:.6f} um/px\n"
        f"scan        {n.width} x {n.height} px stitched per well "
        f"({n.width * n.um_per_px / 1000:.2f} x {n.height * n.um_per_px / 1000:.2f} mm)\n"
        f"channels    {', '.join(n.channels())}\n"
        f"whole-well  every {STEP}th px then {BIN}x binned -> "
        f"{whole[next(iter(whole))].shape[0]} px\n"
        f"field       {CROP} x {CROP} px "
        f"({CROP * n.um_per_px:.1f} x {CROP * n.um_per_px:.1f} um), densest "
        f"mTagBFP-edge field per well, centre at least {EXCLUDE_MM} mm from "
        f"the well centre, mean mTagBFP under {ARTEFACT_CAP}x the well "
        f"median; tile is that field at {TILE} px\n"
        f"            per-well top-left px are in wells.csv (field_y, field_x)\n"
        f"wells       {len(meta)}\n")
    print(f"\n{len(tiles)} tiles, {len(CLOSEUPS)} close-ups -> {OUT}")


if __name__ == "__main__":
    main()
