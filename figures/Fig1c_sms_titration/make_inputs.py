"""Rebuild this figure's two inputs from Matteo's own files on the share.

`make_panel.py` never reads his tree - it reads `inputs/`, like every other
figure script (the rule in `common/cartoon.py`: his working copies are revised
continuously). This script is the bridge, run by hand when the inputs need
regenerating, and it is also the record of where they came from.

    python code/Fig1c_sms_titration/make_inputs.py          # verify only
    python code/Fig1c_sms_titration/make_inputs.py --write  # rewrite inputs

Writes, both under `inputs/`:

  sms_titration_flow.csv        16 rows, 8 constructs x 2 doses
  micrographs/sms_titration/    5 PNG tiles, wells B3 B4 B5 C1 C4

Experiment E.MLS11.104.1 / MLS11.104.1, "Csy4mRNAsensingMultipleHairpins".

Flow: his `Norm mScarlet` and `FoldChange` columns, copied unchanged. The
workbook holds only well ids, so the construct names come from the x axis of
his bar chart, in order. His FoldChange column is what confirms the mapping:
it reads 1 at `F30-2xNegCtrlLoop` and 2.098 / 74.97 / 135.65 at the hairpin,
x2 and x8 - exactly the 2x, 75x and 135x brackets drawn on that chart.

Micrographs: cut from his eight-panel montage, not from the raw 16-bit wells
beside it. The montage already carries ONE display range across all eight
panels, which is the property the figure depends on; scaling each panel to its
own maximum would put the titration in the contrast rather than in the cells.
Panels xy09-xy16 are the 10 ng wells B3-C4 in the same order as the flow table.
"""
import argparse
import itertools
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from PIL import Image

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from common import paths                                      # noqa: E402

MLS = Path("/Volumes/wicky/users/mlampis/Matteo_Lampis")
XLSX = (MLS / "Matteo_Lampis_FlowCytometry" / "E.MLS11.104.1"
        / "E_MLS11_104_1_Analysis.xlsx")
MONTAGE = (MLS / "Matteo_Lampis_Microscopy" / "MLS11.104.1"
           / "MLS11_104_1_Plate3_HEK293H_Csy4MultipleHairpins_Part3" / "Montages"
           / "mScarletZoomIn_xy09_10_11_12_13_14_15_16.png")

CSV_OUT = paths.INPUTS / "sms_titration_flow.csv"
TILE_OUT = paths.INPUTS / "micrographs" / "sms_titration"

# his chart's x axis, same order in both dose groups
CONSTRUCT = ["Csy4-H29A-Only", "F30-2xNegCtrlLoop", "F30-2xCsy4Hairpin",
             "(F30-2xCsy4Hairpin)x2", "(F30-2xCsy4Hairpin)x3",
             "(F30-2xCsy4Hairpin)x4", "(F30-2xCsy4Hairpin)x5",
             "(F30-2xCsy4Hairpin)x8"]
UNITS = [0, 1, 1, 2, 3, 4, 5, 8]
GROUPS = [("100 ng", ["A1", "A2", "A3", "A4", "A5", "A6", "B1", "B2"]),
          ("10 ng",  ["B3", "B4", "B5", "B6", "C1", "C2", "C3", "C4"])]
# the montage panel order, and which five the figure uses
TILES = ["B3_only", "B4_negctrl", "B5_1x", None, "C1_3x", None, None, "C4_8x"]


def flow_table():
    d = pd.read_excel(XLSX)
    d["well"] = d.iloc[:, 0].astype(str).str.extract(r"Plate1_([A-H]\d+)\.fcs")
    d = d.dropna(subset=["well"]).set_index("well")
    rows = []
    for dose, wells in GROUPS:
        for w, c, u in zip(wells, CONSTRUCT, UNITS):
            r = d.loc[w]
            fc = r["FoldChange"]
            rows.append(dict(
                dose=dose, well=w, construct=c, f30_units=u,
                csy4_sites=0 if u == 0 or "NegCtrl" in c else u * 2,
                norm_mscarlet=round(float(r["Norm mScarlet"]), 6),
                fold_change=round(float(fc), 4) if isinstance(fc, (int, float))
                else None))
    t = pd.DataFrame(rows)
    # his own column anchors the construct order; fail loudly if it ever moves
    ten = t[t.dose == "10 ng"].set_index("construct").fold_change
    assert ten["F30-2xNegCtrlLoop"] == 1, "FoldChange is not 1 at the negative"
    assert round(ten["(F30-2xCsy4Hairpin)x8"]) == 136, "x8 is not his 135x bar"
    return t


def tiles():
    im = Image.open(MONTAGE).convert("RGB")
    col = np.asarray(im).astype(float).mean(axis=(0, 2))
    runs = [(r[0][0], r[-1][0]) for k, r in
            ((k, list(g)) for k, g in itertools.groupby(enumerate(col < 150),
                                                        key=lambda t: t[1])) if k]
    wide = np.median([e - s for s, e in runs])
    runs = [r for r in runs if r[1] - r[0] > wide * 0.5]   # drop edge slivers
    assert len(runs) == len(TILES), f"{len(runs)} panels, expected {len(TILES)}"
    return {name: im.crop((s, 0, e + 1, im.height))
            for (s, e), name in zip(runs, TILES) if name}


def main(write):
    t = flow_table()
    cut = tiles()
    if write:
        TILE_OUT.mkdir(parents=True, exist_ok=True)
        t.to_csv(CSV_OUT, index=False)
        for name, tile in cut.items():
            tile.save(TILE_OUT / f"{name}.png")
        print(f"wrote {CSV_OUT.name} ({len(t)} rows) and {len(cut)} tiles")
        return
    old = pd.read_csv(CSV_OUT)
    same = old.equals(t[old.columns])
    print(f"{CSV_OUT.name}: {'matches' if same else 'DIFFERS FROM'} the share")
    for name, tile in cut.items():
        p = TILE_OUT / f"{name}.png"
        ok = p.exists() and np.array_equal(np.asarray(Image.open(p).convert("RGB")),
                                           np.asarray(tile))
        print(f"  {name}.png: {'matches' if ok else 'DIFFERS'}  {tile.size}")
    print("re-run with --write to regenerate")


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--write", action="store_true",
                    help="rewrite the inputs instead of only checking them")
    main(ap.parse_args().write)
