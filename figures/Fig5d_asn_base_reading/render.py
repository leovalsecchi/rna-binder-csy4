"""Fig5d: asparagine reading a base pair.

New framing (2026-11-14): asparagine is the residue of interest. This panel
shows interface asparagines whose sidechain amide reads a BASE, not the
backbone, with the Watson-Crick partner drawn so the pair is visible.

WHICH DESIGNS, AND WHY NOT THE TOP THREE
The three highest by flow FC do not have this to show:

  twist 23  FC 1.56   0 interface asparagines
  twist 15  FC 1.54   0 interface asparagines
  twist 16  FC 1.47   1 interface asparagine (ASN24), reads backbone, not a base

So the panel uses the three binders whose asparagine actually reads a base and
whose contact survives a second structure prediction:

  twist 19  FC 1.34  ASN71 reads C5 (pair C5-G13) and G4; present in the
                     DESIGN model and in both refolds
  twist 12  FC 1.34  ASN7 OD1 to G12 O6 at 2.65 A (pair C6-G12), plus ASN4 to
                     G12 N7; ASN7 present in AF3 and ESM
  twist 07  FC 1.41  ASN54 ND2 to G7 O6 at 3.41 A; AF3 3.41, ESM 3.39

MODEL CHOICE: each panel is rendered on the model with the LOWEST
self-consistency RMSD among those in the file, not on AF3 by default.

  twist 19  design  ligand_RMSD 0 by definition (af3 8.30, b2 1.32, esm 2.93)
  twist 12  ESM     1.77 A   (af3 2.70, b2 22.73)
  twist 07  ESM     2.39 A   (af3 3.62, b2 6.72)

An earlier version rendered AF3 for all three. AF3 is the worst available model
for every one of them, and for twist 19 it sits 8.30 A from the design, far
outside the 3 A the gate asks for. Boltz-2 is closest for twist 19 but is not
in these PDB files, so the designed model is used: twist 19 is BG+BIF, whose
design carries sidechains. The cross-model check is in the README.

Interface = any binder residue with an atom within 5.0 A of the RNA.
Base atoms exclude sugar and phosphate, so "reads a base" is a real
distinction from backbone binding.

House style: code/common/CHIMERAX_RENDERING.md.
"""
import argparse
import subprocess
import sys
from pathlib import Path

import numpy as np
import biotite.structure as struc
import biotite.structure.io.pdb as pdbio


def rodrigues(axis, deg):
    """Rotation matrix for a right-handed turn of `deg` about `axis`."""
    a = np.asarray(axis, float)
    a = a / np.linalg.norm(a)
    t = np.radians(deg)
    K = np.array([[0, -a[2], a[1]], [a[2], 0, -a[0]], [-a[1], a[0], 0]])
    return np.eye(3) + np.sin(t) * K + (1 - np.cos(t)) * (K @ K)


def align_to_y(u):
    """Axis and angle of the rotation taking unit vector u onto +Y."""
    u = np.asarray(u, float) / np.linalg.norm(u)
    y = np.array([0.0, 1.0, 0.0])
    ax = np.cross(u, y)
    n = np.linalg.norm(ax)
    if n < 1e-9:
        return np.array([1.0, 0.0, 0.0]), (0.0 if u[1] > 0 else 180.0)
    return ax / n, float(np.degrees(np.arccos(np.clip(u @ y, -1, 1))))


def atom_xyz(path, model, chain, resid, name):
    cur = 0
    for L in open(path):
        if L.startswith("MODEL"):
            cur = int(L.split()[1])
        elif L.startswith(("ATOM", "HETATM")) and cur == model:
            if (L[21] == chain and int(L[22:26]) == resid
                    and L[12:16].strip() == name):
                return np.array([float(L[30:38]), float(L[38:46]),
                                 float(L[46:54])])
    raise KeyError(f"{path.name} model {model} {chain}:{resid}@{name}")

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "code"))
sys.path.insert(0, str(ROOT / "code" / "Fig2b_map_callouts"))
from common import palette                                        # noqa: E402
from render_shared_frame import canonical_orientation, stack_of   # noqa: E402

CHIMERAX = "/Applications/ChimeraX-1.11.1.app/Contents/MacOS/ChimeraX"
CAND = ROOT.parent / "LV_BinderDesign" / "LV_32Candidates"
RNA_RES = {"A", "U", "G", "C"}
CANON = "ACUGCCGUAUAGGCAG"
CXC = {"BG+BIF": "thesis_bg_bif", "BG+LiMPNN": "thesis_bg_limpnn",
       "RFD3+LiMPNN": "thesis_rfd3_limpnn",
       "BD1+LiMPNN": "thesis_boltzdesign_limpnn",
       "HD+LiMPNN": "thesis_hallu_limpnn"}
DESIGN, AF3, ESM = 1, 2, 3
ASN_C = palette.THEME["700"]        # the residue the panel is about
SWEEP = [0, 30, 60, 90, 120, 150]   # about the contact axis, so it stays in plane
# Chosen per panel from the sweep. Rotating about the contact axis never
# changes the drawn distance, so this is a pure legibility choice: least
# occlusion of the asparagine and the base it reads.
CHOSEN = {"19": 60, "12": 30, "07": 120}
PAIR_C = palette.NEUTRAL["600"]     # the base pair it reads

# twist, asparagines, (read nt, its WC partner), FC, label, measured contacts
# Contacts are drawn as explicit distances, not left to the H-bond geometry
# test: these are predicted models without hydrogens and the test misses two of
# the three. The values are the measured ones, printed on the panel.
PANELS = [
    ("19", DESIGN, [71], [(5, 13)], 1.34,
     "ASN71 reads C5 of the C5-G13 pair", [(71, "OD1", 5, "N4")]),
    ("12", ESM, [7], [(12, 6)], 1.34,
     "ASN7 reads G12 O6 of the C6-G12 pair", [(7, "ND2", 12, "O6")]),
    ("07", ESM, [54], [(7, 11)], 1.41,
     "ASN54 reads G7 O6", [(54, "ND2", 7, "O6")]),
]


def probe(path, model):
    a = pdbio.PDBFile.read(str(path)).get_structure(model=model)
    a = a[a.element != "H"]
    rc = [str(c) for c in np.unique(a.chain_id)
          if set(map(str, np.unique(a[a.chain_id == c].res_name))) <= RNA_RES][0]
    pc = [str(c) for c in np.unique(a.chain_id) if str(c) != rc][0]
    rna = a[a.chain_id == rc]
    ids, names = struc.get_residues(rna)
    off = "".join(map(str, names)).find(CANON)
    return pc, rc, int(ids[off]), int(ids[off + 15])


def build(out_dir, width, sweep=False):
    files = [next(CAND.glob(f"Q1_{t}_*.pdb")) for t, *_ in PANELS]
    MDL = [p[1] for p in PANELS]
    geo = [probe(p, mi) for p, mi in zip(files, MDL)]
    for p, g in zip(files, geo):
        print(f"  {p.name[:46]:<46} pro={g[0]} rna={g[1]} {g[2]}-{g[3]}")

    L = [f"open {ROOT / 'config' / 'thesis_colors.cxc'}"]
    L += [f"open {p}" for p in files]
    L += ["set bgColor white", "hide atoms", "hide cartoons", "hide surfaces"]

    rpro, rrna, r0, r1 = geo[0]
    ref = f"#1.{MDL[0]}/{rrna}:{r0}-{r1}"
    for i, (_, rna, a0, a1) in enumerate(geo, start=1):
        if i == 1:
            continue
        L.append(f"align #{i}.{MDL[i-1]}/{rna}:{a0}-{a1}@C1' toAtoms {ref}@C1' "
                 "move structures")
    centre, v, ang = canonical_orientation(files[0])
    globals()["centre"], globals()["v"], globals()["ang"] = centre, v, ang
    allm = f"#1-{len(files)}"
    for ax, dd in zip("xyz", -centre):
        L.append(f"move {ax} {dd:.4f} models {allm}")
    L.append(f"turn {v[0]:.6f},{v[1]:.6f},{v[2]:.6f} {ang:.4f} center 0,0,0 "
             f"models {allm}")
    L += ["lighting soft",
          f"graphics silhouettes true width 1.5 color {palette.NEUTRAL['700']}"]

    # ------------------------------------------------------------ overviews
    L += [f"surface {ref}", f"color {ref} {palette.NEUTRAL['300']} target s",
          f"hide #1.{MDL[0]}/{rrna} cartoons"]
    did = ["csy4free_" + p.name.split("csy4free_")[1][:-4] for p in files]
    for i in range(1, len(PANELS) + 1):
        L.append(f"color #{i}.{MDL[i-1]}/{geo[i-1][0]} {CXC[stack_of(did[i-1])]}")
    L += ["view", "zoom 0.85"]
    binders = " ".join(f"#{i}.{MDL[i-1]}/{geo[i-1][0]}"
                       for i in range(1, len(PANELS) + 1))
    for i, (tw, mi, asn, _p, _f, _l, _d) in enumerate(PANELS, start=1):
        pro = geo[i - 1][0]
        sel = ",".join(map(str, asn))
        L += [f"hide {binders} cartoons", f"hide #{i}.{mi} atoms",
              f"show #{i}.{mi}/{pro} cartoons",
              # the asparagines, on the overview too, so they can be located
              f"show #{i}.{mi}/{pro}:{sel} atoms",
              f"style #{i}.{mi}/{pro}:{sel} stick",
              f"size #{i}.{mi}/{pro}:{sel} stickRadius 0.22",
              f"color #{i}.{mi}/{pro}:{sel} {ASN_C} target a",
              f"save {out_dir}/t{tw}_overview.png width {width} "
              "supersample 3 transparentBackground true"]
    L += [f"hide {binders} cartoons", f"hide {ref} surfaces"]

    # ---------------------------------------------------------------- zooms
    # The camera is chosen, not inherited. For each panel the scene is turned
    # so the measured contact lies along screen Y, which means the drawn
    # distance is never foreshortened, and then swept about that same axis so
    # every frame keeps it at full length and only the surroundings rotate.
    Rc = rodrigues(v, ang)
    for i, (tw, mi, asn, pairs, _fc, _lab, dists) in enumerate(PANELS, start=1):
        pro, rna, a0, a1 = geo[i - 1]
        me = f"#{i}.{mi}"
        rspec, pspec = f"{me}/{rna}:{a0}-{a1}", f"{me}/{pro}"
        col = CXC[stack_of(did[i - 1])]
        sel = ",".join(map(str, asn))
        nts = sorted({n for p_ in pairs for n in p_})
        nt = ",".join(str(n + a0 - 1) for n in nts)

        ar, aa, rr, ra = dists[0]
        pa = Rc @ (atom_xyz(files[i - 1], mi, pro, ar, aa) - centre)
        pb = Rc @ (atom_xyz(files[i - 1], mi, rna, rr + a0 - 1, ra) - centre)
        axis, adeg = align_to_y(pb - pa)
        mid = (pa + pb) / 2
        cen = f"{mid[0]:.4f},{mid[1]:.4f},{mid[2]:.4f}"
        allm = f"#1-{len(PANELS)}"

        L += [f"hide {me} atoms", f"hide {me} cartoons",
              f"cartoon {rspec}", f"cartoon {pspec}",
              f"cartoon style {pspec} thickness 0.4",
              f"color {rspec} {palette.NEUTRAL['400']}",
              f"color {pspec} {col}",
              f"transparency {pspec} 80 target c",
              f"transparency {rspec} 35 target c",
              f"show {rspec} atoms", f"nucleotides {rspec} fill",
              f"color {rspec} {palette.NEUTRAL['400']} target ac",
              f"transparency {rspec} 70 target a",
              f"color {me}/{rna}:{nt} {palette.NEUTRAL['600']} target ac",
              f"transparency {me}/{rna}:{nt} 0 target ac",
              f"show {pspec}:{sel} atoms", f"style {pspec}:{sel} stick",
              f"size {pspec}:{sel} stickRadius 0.26",
              f"color {pspec}:{sel} {ASN_C} target a",
              f"transparency {pspec}:{sel} 0 target a"]
        L += [f"distance {pspec}:{q[0]}@{q[1]} {me}/{rna}:{q[2] + a0 - 1}@{q[3]} "
              f"color {ASN_C} radius 0.05 dashes 6" for q in dists]
        L.append(f"turn {axis[0]:.6f},{axis[1]:.6f},{axis[2]:.6f} {adeg:.4f} "
                 f"center {cen} models {allm}")
        angles = SWEEP if sweep else [CHOSEN[tw]]
        prev = 0
        for th in angles:
            if th != prev:
                L.append(f"turn y {th - prev} center {cen} models {allm}")
                prev = th
            suffix = f"_y{th:03d}" if sweep else ""
            L += [f"view {pspec}:{sel} {me}/{rna}:{nt}", "zoom 1.45",
                  f"save {out_dir}/t{tw}_asn{suffix}.png width {width} "
                  "supersample 3 transparentBackground true"]
        # put the scene back so the next panel starts from the shared frame
        L.append(f"turn y {-prev} center {cen} models {allm}")
        L.append(f"turn {axis[0]:.6f},{axis[1]:.6f},{axis[2]:.6f} {-adeg:.4f} "
                 f"center {cen} models {allm}")
        L += ["~distance", f"hide {me} atoms", f"hide {me} cartoons"]
    L.append("exit")
    return "\n".join(L)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("out_dir")
    ap.add_argument("--width", type=int, default=3000)
    ap.add_argument("--sweep", action="store_true",
                    help="render all SWEEP angles for picking")
    a = ap.parse_args()
    out = Path(a.out_dir).resolve()
    out.mkdir(parents=True, exist_ok=True)
    script = out / "render_asn.cxc"
    script.write_text(build(out, a.width, a.sweep) + "\n")
    r = subprocess.run([CHIMERAX, "--exit", "--cmd", f"open {script}"],
                       capture_output=True, text=True)
    for ln in (r.stdout + r.stderr).splitlines():
        if "rror" in ln:
            print("  cx|", ln)
    pngs = sorted(p.name for p in out.glob("*.png") if not p.name.startswith("._"))
    print(f"wrote {len(pngs)} png: {pngs}")


if __name__ == "__main__":
    main()
