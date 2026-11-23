#!/usr/bin/env python3
"""Render the Figure 2b interaction examples under the house convention in
code/common/CHIMERAX_RENDERING.md: design vs AF3 refold overlay, binder in its
stack colour, refold in neutral 400 at 35% transparency, RNA as a light
surface, soft lighting, neutral-700 silhouettes, transparent background.

Everything the doc's Gotchas warn about is handled explicitly:
  - chain roles detected by CONTENT (RFD3 design files are flipped)
  - RNA residue RANGES read per file, never assumed: one design here numbers
    its RNA 2-17 while its own refold uses 1-16, so a naive 1-16 align
    mispairs without erroring
  - no ';' inline comments in the .cxc (a ';' line aborts the rest of the file)
  - sidechains come from the refold only; the design is backbone-only
  - `set bgColor white` kept even with transparentBackground
  - `view` called ONCE, then zoom, then `turn y 45` between saves

Superposition is on the RNA (the DockQ convention in the doc's table), because
these panels are about where the binder sits on the hairpin, not fold error.

Usage: python render_convention.py <out_dir> [--width 1800]
"""
import argparse
import subprocess
import sys
from pathlib import Path

import numpy as np
import biotite.structure as struc
import biotite.structure.io.pdbx as pdbx

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from common import palette                                      # noqa: E402

CHIMERAX = "/Applications/ChimeraX-1.11.1.app/Contents/MacOS/ChimeraX"
BASE = Path(__file__).resolve().parents[2] / "working/structures/fig02b_examples"
RNA_RES = {"A", "U", "G", "C"}
CANON = "ACUGCCGUAUAGGCAG"
SHORT = {"RFD3_LiMPNN": "RFD3+LiMPNN", "BG_LiMPNN": "BG+LiMPNN",
         "BG_BIF": "BG+BIF", "BD1_LiMPNN": "BD1+LiMPNN", "HD_LiMPNN": "HD+LiMPNN"}

PICKS = [
    ("a_backbone_minorgroove", "csy4free_RFD3_LiMPNN_len80_8b13398e_design50_seq0"),
    ("b_head_loop",            "csy4free_BG_LiMPNN_len80_1caaae08_design23_seq0"),
    ("c_majorgroove_stem",     "csy4free_RFD3_LiMPNN_len80_383063df_design6_seq0"),
    ("d_pi_stacking",          "csy4free_RFD3_LiMPNN_len80_998adbc2_design46_seq3"),
]


def probe(path):
    """(protein_chain, rna_chain, first_canonical_resid, last) read from file."""
    a = pdbx.get_structure(pdbx.CIFFile.read(str(path)), model=1)
    a = a[a.element != "H"]
    rna_c = [str(c) for c in np.unique(a.chain_id)
             if set(map(str, np.unique(a[a.chain_id == c].res_name))) <= RNA_RES]
    pro_c = [str(c) for c in np.unique(a.chain_id) if str(c) not in rna_c]
    if len(rna_c) != 1 or len(pro_c) != 1:
        raise ValueError(f"{path.name}: rna={rna_c} prot={pro_c}")
    rna = a[a.chain_id == rna_c[0]]
    ids, names = struc.get_residues(rna)
    off = "".join(map(str, names)).find(CANON)
    if off < 0:
        raise ValueError(f"{path.name}: canonical frame not found")
    return pro_c[0], rna_c[0], int(ids[off]), int(ids[off + len(CANON) - 1])


def stack_of(did):
    for k, v in SHORT.items():
        if f"_{k}_" in did:
            return v
    raise ValueError(did)


def block(tag, did, out_dir, width, frames):
    dpath = BASE / "cifs" / f"{did}.cif"
    rpath = BASE / "af3" / f"{did}_AF3_seed1_sample0.cif"
    dpro, drna, d0, d1 = probe(dpath)
    rpro, rrna, r0, r1 = probe(rpath)
    col = palette.stack_color(stack_of(did))
    L = [
        f"open {dpath}",
        f"open {rpath}",
        "set bgColor white",
        "hide atoms",
        "cartoon",
        f"align #2/{rrna}:{r0}-{r1}@C1' toAtoms #1/{drna}:{d0}-{d1}@C1' move structures",
        f"hide #2/{rrna} cartoons",
        f"color #1/{dpro} {col}",
        f"color #2/{rpro} {palette.NEUTRAL['400']}",
        "transparency #2 35 target c",
        f"surface #1/{drna}:{d0}-{d1}",
        f"color #1/{drna} #d7dcdf target s",
        f"hide #1/{drna} cartoons",
        "lighting soft",
        f"graphics silhouettes true width 1.5 color {palette.NEUTRAL['700']}",
        "view",
        "zoom 0.8",
    ]
    for i in range(frames):
        if i:
            L.append("turn y 45")
        L.append(f"save {out_dir}/{tag}_y{45*i:03d}.png width {width} "
                 "supersample 3 transparentBackground true")
    L.append("close session")
    return "\n".join(L)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("out_dir")
    ap.add_argument("--width", type=int, default=1800)
    ap.add_argument("--frames", type=int, default=8)
    a = ap.parse_args()
    out = Path(a.out_dir).resolve()
    out.mkdir(parents=True, exist_ok=True)
    script = out / "render_convention.cxc"
    script.write_text("\n".join(block(t, d, out, a.width, a.frames)
                                for t, d in PICKS) + "\nexit\n")
    print(f"{len(PICKS)} panels -> {script}")
    r = subprocess.run([CHIMERAX, "--exit", "--cmd", f"open {script}"],
                       capture_output=True, text=True)
    for ln in (r.stdout + r.stderr).strip().splitlines()[-20:]:
        print("  cx|", ln)
    pngs = [p for p in out.glob("*.png") if not p.name.startswith("._")]
    print(f"wrote {len(pngs)} png (expect {a.frames*len(PICKS)})")


if __name__ == "__main__":
    main()
