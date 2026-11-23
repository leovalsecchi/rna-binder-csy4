#!/usr/bin/env python3
"""Figure 2b interaction panels in ONE SHARED RNA FRAME.

Every complex is superposed onto a single reference hairpin and drawn with the
same camera, so the RNA is in the identical position and orientation in every
panel and only the binder moves. Without this you cannot see that one binder
reads the major groove while another sits on the loop: each panel would have
its own arbitrary RNA pose.

How the shared frame is built: every model, design AND refold, is aligned
independently onto the SAME reference RNA (panel a's design) by C1' atoms. Two
independent fits to one reference land in one frame, which avoids needing to
move a design and its refold as a group. Only the reference RNA is drawn; all
other RNA copies are hidden, so the surface in every panel is literally the
same atoms.

House convention from code/common/CHIMERAX_RENDERING.md is kept: binder in its
stack colour, refold in neutral 400 at 35% transparency, RNA a #d7dcdf surface,
soft lighting, neutral-700 silhouettes width 1.5, white bg, transparent PNG,
supersample 3. `view` is called once and never re-run between frames.

RNA residue ranges are read per file (one design numbers its RNA 2-17 while its
refold uses 1-16) and chain roles are detected by content (RFD3 design files
have the chains flipped).

Usage: python render_shared_frame.py <out_dir> [--width 1800] [--frames 8]
"""
import argparse
import subprocess
import sys
from pathlib import Path

import numpy as np
import biotite.structure as struc
import biotite.structure.io.pdbx as pdbx
import biotite.structure.io.pdb as pdbio

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from common import palette                                      # noqa: E402

CHIMERAX = "/Applications/ChimeraX-1.11.1.app/Contents/MacOS/ChimeraX"
BASE = Path(__file__).resolve().parents[2] / "working/structures/fig02b_examples"
RNA_RES = {"A", "U", "G", "C"}
CANON = "ACUGCCGUAUAGGCAG"
SHORT = {"RFD3_LiMPNN": "RFD3+LiMPNN", "BG_LiMPNN": "BG+LiMPNN",
         "BG_BIF": "BG+BIF", "BD1_LiMPNN": "BD1+LiMPNN", "HD_LiMPNN": "HD+LiMPNN"}

# Head binding dropped on request 2026-11-05. The reference RNA frame is still
# panel a's design, so dropping a panel does not move the hairpin in the others.
# Groove-edge atoms, standard Watson-Crick geometry. The surface patch is
# coloured from THESE atoms, so the highlight is chemistry, not a hand-drawn
# blob: ChimeraX colours each surface vertex from its nearest atom.
#   major groove: purine N7/O6/N6, pyrimidine N4/O4, plus C5/C6 edge carbons
#   minor groove: purine N3/N2/C2, pyrimidine O2, plus the sugar O2'
GROOVE = {
    "major": "@N7,O6,N6,N4,O4,C5,C6",
    "minor": "@N3,N2,C2,O2,O2'",
}

# (tag, design, which groove this panel is about)
PANELS = [
    ("a_minorgroove", "csy4free_RFD3_LiMPNN_len80_8b13398e_design50_seq0", "minor"),
    ("c_majorgroove", "csy4free_RFD3_LiMPNN_len80_383063df_design6_seq0", "major"),
]


PAIRS = [(2, 16), (3, 15), (4, 14), (5, 13), (6, 12), (7, 11)]
RING = {"A": ["N9", "C8", "N7", "C5", "C6", "N1", "C2", "N3", "C4"],
        "G": ["N9", "C8", "N7", "C5", "C6", "N1", "C2", "N3", "C4"],
        "U": ["N1", "C2", "N3", "C4", "C5", "C6"],
        "C": ["N1", "C2", "N3", "C4", "C5", "C6"]}


def _atoms(path):
    """Model 1, hydrogens dropped. Accepts .cif and .pdb: the ordered-hit
    structures (working/structures/hits32) are PDB, everything else is CIF."""
    if str(path).endswith(".pdb"):
        a = pdbio.PDBFile.read(str(path)).get_structure(model=1)
    else:
        a = pdbx.get_structure(pdbx.CIFFile.read(str(path)), model=1)
    return a[a.element != "H"]


def canonical_orientation(path):
    """Commands that put the hairpin upright: helix axis vertical on screen,
    OPEN END (A1/G16) at the bottom, loop at the top, viewed side-on.

    Derived from the geometry, not guessed: the helix axis is the first
    principal direction of the six Watson-Crick pair midpoints, signed so the
    loop is positive. The camera's default is -Z view with +Y up, so moving the
    RNA centroid to the origin and rotating the axis onto +Y gives the side
    view directly, and a later `turn y` then spins about the hairpin axis,
    which keeps the bottom at the bottom in every frame.
    """
    a = _atoms(path)
    rc = [str(c) for c in np.unique(a.chain_id)
          if set(map(str, np.unique(a[a.chain_id == c].res_name))) <= RNA_RES][0]
    rna = a[a.chain_id == rc]
    ids, names = struc.get_residues(rna)
    off = "".join(map(str, names)).find(CANON)
    cen = {}
    for k in range(1, 17):
        rid = ids[off + k - 1]
        sel = rna[(rna.res_id == rid) & np.isin(rna.atom_name, RING[CANON[k - 1]])]
        cen[k] = sel.coord.mean(0).astype(float)
    mids = np.array([(cen[i] + cen[j]) / 2 for i, j in PAIRS])
    centre = np.array([cen[k] for k in range(1, 17)]).mean(0)
    axis = np.linalg.svd(mids - mids.mean(0))[2][0]
    if np.dot(cen[9] - centre, axis) < 0:        # sign so the loop is "up"
        axis = -axis
    # face the flipped-out U10 toward the camera: reproducible, and it is the
    # nucleotide every binder competes for
    f = cen[10] - centre
    f -= np.dot(f, axis) * axis
    f /= np.linalg.norm(f)
    ey, ez = axis, f
    ex = np.cross(ey, ez)
    R = np.vstack([ex, ey, ez])                  # scene -> canonical
    ang = np.degrees(np.arccos(np.clip((np.trace(R) - 1) / 2, -1, 1)))
    if ang < 1e-6:
        v = np.array([0.0, 1.0, 0.0])
    else:
        v = np.array([R[2, 1] - R[1, 2], R[0, 2] - R[2, 0], R[1, 0] - R[0, 1]])
        v /= np.linalg.norm(v)
    return centre, v, ang


def probe(path):
    a = _atoms(path)
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
    return pro_c[0], rna_c[0], int(ids[off]), int(ids[off + 15])


def stack_of(did):
    for k, v in SHORT.items():
        if f"_{k}_" in did:
            return v
    raise ValueError(did)


def build(out_dir, width, frames):
    # Design only: the re-prediction overlay was dropped on request, so each
    # panel is one model and the binder reads without a grey ghost behind it.
    files, meta = [], []
    for tag, did, groove in PANELS:
        p = BASE / "cifs" / f"{did}.cif"
        files.append(p)
        meta.append((tag, did, groove) + probe(p))

    L = [f"open {p}" for p in files]
    L += ["set bgColor white", "hide atoms", "cartoon"]

    ref_pro, ref_rna, ref0, ref1 = meta[0][3], meta[0][4], meta[0][5], meta[0][6]
    for i, m in enumerate(meta, start=1):
        if i == 1:
            continue
        _, _, _, pro, rna, a0, a1 = m
        L.append(f"align #{i}/{rna}:{a0}-{a1}@C1' toAtoms "
                 f"#1/{ref_rna}:{ref0}-{ref1}@C1' move structures")

    # one shared RNA: draw the reference copy only, hide every other RNA
    for i, m in enumerate(meta, start=1):
        rna = m[4]
        if i == 1:
            L += [f"surface #1/{ref_rna}:{ref0}-{ref1}",
                  f"color #1/{ref_rna} #d7dcdf target s",
                  f"hide #1/{ref_rna} cartoons"]
        else:
            L.append(f"hide #{i}/{rna} cartoons")
            L.append(f"hide #{i}/{rna} atoms")

    for i, m in enumerate(meta, start=1):
        L.append(f"color #{i}/{m[3]} {palette.stack_color(stack_of(m[1]))}")

    centre, v, ang = canonical_orientation(files[0])
    allm = f"#1-{len(files)}"
    for ax, d in zip("xyz", -centre):
        L.append(f"move {ax} {d:.4f} models {allm}")
    L.append(f"turn {v[0]:.6f},{v[1]:.6f},{v[2]:.6f} {ang:.4f} center 0,0,0 "
             f"models {allm}")

    L += ["lighting soft",
          f"graphics silhouettes true width 1.5 color {palette.NEUTRAL['700']}",
          "view", "zoom 1.0"]

    all_binders = [f"#{i}/{m[3]}" for i, m in enumerate(meta, start=1)]
    ref_rna_spec = f"#1/{ref_rna}:{ref0}-{ref1}"

    for f in range(frames):
        if f:
            L.append("turn y 45")
        for i, m in enumerate(meta, start=1):
            tag, groove = m[0], m[2]
            L.append("hide " + " ".join(all_binders) + " cartoons")
            L.append(f"show #{i}/{m[3]} cartoons")
            # retint the shared surface for the groove THIS panel is about
            L.append(f"color {ref_rna_spec} #d7dcdf target s")
            L.append(f"color {ref_rna_spec}{GROOVE[groove]} "
                     f"{palette.THEME['300']} target s")
            L.append(f"save {out_dir}/{tag}_y{45*f:03d}.png width {width} "
                     "supersample 3 transparentBackground true")
    L.append("exit")
    return "\n".join(L)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("out_dir")
    ap.add_argument("--width", type=int, default=1800)
    ap.add_argument("--frames", type=int, default=8)
    a = ap.parse_args()
    out = Path(a.out_dir).resolve()
    out.mkdir(parents=True, exist_ok=True)
    script = out / "render_shared_frame.cxc"
    script.write_text(build(out, a.width, a.frames) + "\n")
    print("script ->", script)
    r = subprocess.run([CHIMERAX, "--exit", "--cmd", f"open {script}"],
                       capture_output=True, text=True)
    for ln in (r.stdout + r.stderr).strip().splitlines()[-20:]:
        print("  cx|", ln)
    pngs = [p for p in out.glob("*.png") if not p.name.startswith("._")]
    print(f"wrote {len(pngs)} png (expect {a.frames*len(PANELS)})")


if __name__ == "__main__":
    main()
