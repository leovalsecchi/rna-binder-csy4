#!/usr/bin/env python3
"""ChimeraX renders for panel (a): the 4AL5 complex, and the hairpin taken out of it.

Both images are saved from ONE camera, so the second reads as the first with the
protein removed rather than as a different picture. The camera is the repo's
canonical hairpin frame (`Fig2b_map_callouts/render_shared_frame.canonical_orientation`:
helix axis vertical, open end down, flipped-out U10 toward the viewer), so this panel
shows the same molecule in the same pose as Fig 2a and Fig 2b.

Source: the deposited 4AL5 mmCIF. `prepare()` writes a cleaned copy - chain A protein
and chain B RNA, waters, glycerol and the deposited hydrogens dropped - into
working/structures/target_definition/, which is also what makes the repo's chain and
frame helpers work on it (they identify the RNA chain as the one whose residue names
are all in {A,U,G,C}, and chain B carries GOL/HOH in the deposition).

    python code/supp/target_definition/render_target.py              # 4 azimuths
    python code/supp/target_definition/render_target.py --azimuth 30 # one, for the panel

ChimeraX cannot render offscreen on macOS, so this steals focus while it runs and
--nogui silently saves nothing (code/common/CHIMERAX_RENDERING.md).
"""
import argparse
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "Fig2b_map_callouts"))

import numpy as np
import biotite.structure as struc
import biotite.structure.io.pdbx as pdbx

from common import palette, paths                                    # noqa: E402
from render_shared_frame import CHIMERAX, canonical_orientation      # noqa: E402

FIG = Path(__file__).resolve().parent.name
SRC = Path.home() / "4al5.cif"
WORK = paths.ROOT / "working" / "structures" / FIG
CLEAN = WORK / "4al5_protein_rna.cif"

RNA_RES = {"A", "U", "G", "C"}
CHAIN = (4, 21)          # every residue of the campaign target file, drawn whole
CORE = (5, 20)           # the 16 of them that carry a base; residues 4 and 21 are a
                         # lone O3' and a bare 5'-phosphate (verified in make_figure.py)

# The hairpin takes THE MAIN FIGURES' COLOURS, not `synbio.rna` (changed
# 2026-11-19 on Leo's comment that the colour read "weird" here). Fig 2a renders
# this same molecule in this same canonical frame as two-tone petrol -
# `Fig2a_pipeline/render_hairpin.py` BACKBONE/BASES, theme 700 and 300 - so the
# supplement and the main text now show one hairpin, not two. `synbio.rna` is
# theme 400 and was applied to cartoon AND slab alike, so the whole hairpin came
# out one flat pale tone with the base slabs invisible against the backbone.
#
# This is a deliberate departure from the synbio molecular-identity rule
# (COLORS.md numbers it 8; the file has no literal rule 6, and this script's
# earlier comment cited the number the rule carried before the list was
# renumbered). That rule sets `rna` = theme 400 "the target is context", but the
# target is the SUBJECT of this figure, and no main-figure structure render uses
# `synbio.rna`: Fig 2a uses theme 700/300, Fig 2b theme 600, Fig 5c neutral 400.
# The warm diverging end `#ae4e1a` is the Csy4 site colour on the SCHEMATIC
# glyphs only (rule 4); it stays out of structure renders, where COLORS.md rule 4
# and CHIMERAX_RENDERING.md both say it has no categorical use.
#
# Protein stays neutral (context, rule 3) and the silhouette stays neutral
# because the protein shares it. No hex is written by hand: every value is read
# from config/palette.json through common.palette.
C_RNA = palette.THEME["700"]      # backbone/cartoon, = Fig2a BACKBONE
C_BASE = palette.THEME["300"]     # base slabs,       = Fig2a BASES
C_PROTEIN = palette.NEUTRAL["300"]
C_EDGE = palette.SYNBIO["outline"]


def prepare():
    """Chain A protein + chain B RNA, no solvent, no glycerol, no hydrogens."""
    WORK.mkdir(parents=True, exist_ok=True)
    f = pdbx.CIFFile.read(SRC)
    a = pdbx.get_structure(f, model=1)
    keep = (((a.chain_id == "A") & struc.filter_amino_acids(a))
            | ((a.chain_id == "B") & np.isin(a.res_name, list(RNA_RES))))
    a = a[keep & (a.element != "H")]
    out = pdbx.CIFFile()
    pdbx.set_structure(out, a)
    out.write(CLEAN)
    rna = a[a.chain_id == "B"]
    ids, names = struc.get_residues(rna)
    print(f"{CLEAN.name}: {len(a)} heavy atoms, protein {int((a.chain_id=='A').sum())}, "
          f"RNA {len(rna)} over residues {ids[0]}-{ids[-1]} "
          f"({''.join(map(str, names))})")
    return a


def build(out, width, azimuths):
    centre, v, ang = canonical_orientation(CLEAN)
    rna = f"#1/B:{CORE[0]}-{CORE[1]}"      # the only residues with a base to slab
    L = [f"open {CLEAN}", f"open {paths.ROOT}/config/thesis_colors.cxc",
         "set bgColor white", "hide atoms", "hide cartoons"]
    for ax, d in zip("xyz", -centre):
        L.append(f"move {ax} {d:.4f} models #1")
    L.append(f"turn {v[0]:.6f},{v[1]:.6f},{v[2]:.6f} {ang:.4f} center 0,0,0 models #1")
    # `nucleotides ... tube/slab` only acts on DISPLAYED atoms, and `fill` resets
    # slab colour, so show first, style, then colour (CHIMERAX_RENDERING.md).
    # The whole chain is cartooned, so nothing of the campaign file is trimmed
    # away; residues 4 and 21 carry no base and no traceable backbone, so what
    # draws is the 16 complete nucleotides. They are NOT shown as loose atoms:
    # a bare 5'-phosphate drawn as four balls reads as a bound small molecule.
    #
    # ⚠️ Colouring the slabs is RESIDUE-then-revert, never an atom subset. A slab
    # takes its WHOLE RESIDUE's colour, so `color ... target a` paints a slab
    # NOTHING - that atom-subset idiom works on a molecular SURFACE, not here,
    # and it silently left every slab at the backbone colour until 2026-11-20.
    # So: paint the residue pale (ribbon AND slab go C_BASE), then put the ribbon
    # back with `target c`, which reverts ONLY the cartoon - in tube/slab the
    # base IS the slab and the sugar+phosphate ARE the ribbon.
    L += ["cartoon #1/B", f"show {rna} atoms", f"nucleotides {rna} tube/slab",
          f"color #1/B {C_RNA}",
          f"color {rna} {C_BASE}", f"color {rna} {C_RNA} target c",
          "cartoon #1/A", f"color #1/A {C_PROTEIN}",
          f"graphics silhouettes true width 2 color {C_EDGE}",
          "lighting soft", "camera ortho", "view", "zoom 0.92"]
    prev = 0
    for az in azimuths:
        if az - prev:
            L.append(f"turn y {az - prev}")
        prev = az
        # Same camera for both saves: only the protein's visibility changes.
        L += [f"save {out}/complex_y{az:03d}.png width {width} supersample 3 "
              "transparentBackground true",
              "hide #1/A cartoons",
              f"save {out}/hairpin_y{az:03d}.png width {width} supersample 3 "
              "transparentBackground true",
              "cartoon #1/A"]
    L.append("exit")
    return "\n".join(L)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--width", type=int, default=1600)
    ap.add_argument("--azimuth", default="0,30,60,90")
    a = ap.parse_args()
    out = paths.fig_output(FIG) / "renders"
    out.mkdir(parents=True, exist_ok=True)
    prepare()
    script = WORK / "render_target.cxc"
    script.write_text(build(out, a.width,
                            [int(x) for x in a.azimuth.split(",")]) + "\n")
    r = subprocess.run([CHIMERAX, "--exit", "--cmd", f"open {script}"],
                       capture_output=True, text=True)
    for ln in (r.stdout + r.stderr).strip().splitlines()[-12:]:
        print("  cx|", ln)
    pngs = sorted(p for p in out.glob("*_y*.png") if not p.name.startswith("._"))
    print(f"wrote {len(pngs)} png into {out}")


if __name__ == "__main__":
    main()
