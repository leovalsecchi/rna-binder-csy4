#!/usr/bin/env python3
"""Render the Figure 2b interaction examples in ChimeraX.

One .cxc per design plus three turntable views, following
code/common/CHIMERAX_RENDERING.md. Colours come from config/palette.json.

Chain roles are detected by CONTENT, never by chain id: RFD3 design files carry
the RNA in chain A and the binder in chain B, the opposite of every other
stack, so a hardcoded "#1/A = binder" silently colours the RNA as the binder.

Interface H-bonds are drawn with ChimeraX `hbonds` so the panel shows the
interaction rather than just the pose.

Usage:  python render_examples.py <cif_dir> <out_dir>
"""
import subprocess
import sys
from pathlib import Path

import numpy as np
import biotite.structure.io.pdbx as pdbx

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from common import palette                                     # noqa: E402

CHIMERAX = "/Applications/ChimeraX-1.11.1.app/Contents/MacOS/ChimeraX"
RNA_RES = {"A", "U", "G", "C"}
SHORT = {"RFD3_LiMPNN": "RFD3+LiMPNN", "BG_LiMPNN": "BG+LiMPNN",
         "BG_BIF": "BG+BIF", "BD1_LiMPNN": "BD1+LiMPNN",
         "HD_LiMPNN": "HD+LiMPNN"}


def roles(cif):
    """(protein_chain, rna_chain) decided by residue content."""
    a = pdbx.get_structure(pdbx.CIFFile.read(str(cif)), model=1)
    a = a[a.element != "H"]
    rna, pro = [], []
    for c in np.unique(a.chain_id):
        names = set(map(str, np.unique(a[a.chain_id == c].res_name)))
        (rna if names <= RNA_RES else pro).append(str(c))
    if len(rna) != 1 or len(pro) != 1:
        raise ValueError(f"{cif.name}: rna={rna} prot={pro}")
    return pro[0], rna[0]


def stack_of(design_id):
    for k, v in SHORT.items():
        if f"_{k}_" in design_id:
            return v
    raise ValueError(f"no stack in {design_id}")


def cxc_for(cif, out_dir):
    cif = cif.resolve()
    out_dir = out_dir.resolve()
    pro, rna = roles(cif)
    stem = cif.stem
    col = palette.stack_color(stack_of(stem))
    ink = palette.NEUTRAL["700"]
    rna_surf = palette.NEUTRAL["200"]
    lines = [
        f"open {cif}",
        "set bgColor white",
        "hide atoms",
        "cartoon",
        f"color #1/{pro} {col}",
        # The RNA bases must be visible or stacking and groove contacts are
        # invisible. Verified: `nucleotides fill` and `nucleotides ladder` do
        # NOTHING here after `hide atoms` (both produce a byte-identical image
        # to plain cartoon). Showing the RNA atoms as sticks is what works.
        f"color #1/{rna} {palette.THEME['600']}",
        f"show #1/{rna} atoms",
        f"style #1/{rna} stick",
        f"size #1/{rna} stickRadius 0.20",
        # binder side chains that actually touch the RNA
        f"show #1/{pro} & #1/{rna} :< 4.5 atoms",
        f"style #1/{pro} & #1/{rna} :< 4.5 stick",
        f"size #1/{pro} stickRadius 0.18",
        f"hbonds #1/{pro} restrict #1/{rna} reveal true color {ink} "
        "radius 0.06 dashes 0",
        "lighting soft",
        f"graphics silhouettes true width 1.2 color {ink}",
        "view",
        "zoom 0.85",
    ]
    for ang in (0, 120, 240):
        if ang:
            lines.append(f"turn y 120")
        lines.append(
            f"save {out_dir}/{stem}_y{ang:03d}.png width 1800 "
            "supersample 3 transparentBackground true")
    lines.append("close session")
    return "\n".join(lines)


def main():
    cif_dir, out_dir = Path(sys.argv[1]).resolve(), Path(sys.argv[2]).resolve()
    out_dir.mkdir(parents=True, exist_ok=True)
    cifs = sorted(p for p in cif_dir.glob("*.cif") if not p.name.startswith("._"))
    if not cifs:
        sys.exit(f"no cifs in {cif_dir}")
    blocks, skipped = [], []
    for c in cifs:
        try:
            blocks.append(cxc_for(c, out_dir))
        except Exception as e:
            skipped.append(f"{c.name}: {e}")
    script = out_dir / "render_all.cxc"
    script.write_text("\n".join(blocks) + "\nexit\n")
    print(f"{len(blocks)} designs -> {script}")
    for s in skipped:
        print("SKIPPED", s)
    # --nogui parses and runs but saves nothing (no OpenGL context).
    r = subprocess.run([CHIMERAX, "--exit", "--cmd", f"open {script}"],
                       capture_output=True, text=True)
    tail = (r.stdout + r.stderr).strip().splitlines()
    for ln in tail[-25:]:
        print("  cx|", ln)
    pngs = sorted(out_dir.glob("*.png"))
    print(f"wrote {len(pngs)} png (expect {3*len(blocks)})")


if __name__ == "__main__":
    main()
