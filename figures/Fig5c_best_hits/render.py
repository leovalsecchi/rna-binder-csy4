"""Figure 5 panel c: the three highest by flow fold change, plus twist 7's
asparagine on the right.

  twist 23  FC 1.56  BG+LiMPNN    ARG12, LYS23, TYR24 on G7 and the A11 phosphate
  twist 15  FC 1.54  BG+BIF       ARG67 into the major groove at G12 N7
  twist 16  FC 1.47  RFD3+LiMPNN  LYS42 on G7 O6 and the U10 phosphate
  twist 07  FC 1.41  BD1+LiMPNN   ASN54 on the A9 phosphate, reaching G7 O6

Each is rendered on the model with the LOWEST self-consistency RMSD among those
in the file, never AF3 by default:

  twist 23  AF3 1.55   (b2 10.89, esm 16.03)
  twist 15  ESM 2.08   (af3 2.41, b2 39.77)
  twist 16  ESM 1.83   (af3 7.86, b2 5.72)
  twist 07  ESM 2.39   (af3 3.62, b2 6.72)

DASHES ARE DRAWN ON TWO PANELS ONLY. Running ChimeraX `hbonds` on the chosen
models, most of the contacts worth showing either fail the geometry test or
pass at distances too short to be real:

  twist 23  ARG12 NH2-G7 N7 1.94 A and LYS23 NZ-A11 OP1 1.42 A pass the test
            but are overlap distances, not bonds. TYR24 OH-A11 OP2 2.59 A is
            the only clean one.
  twist 15  ARG67 to G12 N7 is 2.43 A and is NOT an H-bond by geometry.
  twist 16  only LYS42 NZ-U10 OP1 at 3.11 A passes.
  twist 07  ASN54 has no H-bond in this model.

So only two contacts are both geometrically valid and at a real bond distance,
and only those two are dashed: twist 23 TYR24 OH-A11 OP2 and twist 16 LYS42
NZ-U10 OP1. Twist 15 and twist 07 are shown without dashes. A reader will take
the absence of a dash as "no hydrogen bond here", which is what the test says,
but the figure legend must state that dashes mark contacts passing the ChimeraX
hbonds test so the asymmetry is not read as a rendering oversight.

House style: code/common/CHIMERAX_RENDERING.md.
"""
import argparse
import subprocess
import sys
from pathlib import Path

import numpy as np
import biotite.structure as struc
import biotite.structure.io.pdb as pdbio

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
RES_C = palette.THEME["700"]
NT_C = palette.NEUTRAL["600"]
SWEEP = [0, 30, 60, 90, 120, 150, 180, 210, 240, 270, 300, 330]

# twist, model, residues, nucleotides, chosen angle, FC, dashes
#
# Dashes are drawn ONLY for contacts that pass the ChimeraX `hbonds` geometry
# test AT A DISTANCE THAT IS PHYSICALLY A BOND. Two qualify across the four
# panels. Twist 23's ARG12-G7 N7 (1.94 A) and LYS23-A11 OP1 (1.42 A) pass the
# test but are overlap distances, so they are deliberately NOT dashed; twist 15
# and twist 07 have no passing contact at all.
PANELS = [
    ("23", AF3, [12, 23, 24], [7, 11], 60, 1.56, [(24, "OH", 11, "OP2")]),
    ("15", ESM, [67], [12], 60, 1.54, []),
    ("16", ESM, [42, 50], [7, 10], 90, 1.47, [(42, "NZ", 10, "OP1")]),
    ("07", ESM, [54, 53], [9, 7], 30, 1.41, []),
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


def build(out_dir, width, sweep=False, surface=False):
    files = [next(CAND.glob(f"Q1_{t}_*.pdb")) for t, *_ in PANELS]
    MDL = [p[1] for p in PANELS]
    geo = [probe(p, mi) for p, mi in zip(files, MDL)]
    did = ["csy4free_" + p.name.split("csy4free_")[1][:-4] for p in files]

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
    allm = f"#1-{len(PANELS)}"
    for ax, dd in zip("xyz", -centre):
        L.append(f"move {ax} {dd:.4f} models {allm}")
    L.append(f"turn {v[0]:.6f},{v[1]:.6f},{v[2]:.6f} {ang:.4f} center 0,0,0 "
             f"models {allm}")
    L += ["lighting soft",
          f"graphics silhouettes true width 1.5 color {palette.NEUTRAL['700']}",
          "view", "zoom 0.9", "distance style decimalPlaces 2"]

    panels = PANELS[:3] if surface else PANELS
    for i, (tw, mi, res, nts, chosen, _fc, dash) in enumerate(panels, start=1):
        pro, rna, a0, a1 = geo[i - 1]
        me = f"#{i}.{mi}"
        rspec, pspec = f"{me}/{rna}:{a0}-{a1}", f"{me}/{pro}"
        col = CXC[stack_of(did[i - 1])]
        sel = ",".join(map(str, res))
        nt = ",".join(str(n + a0 - 1) for n in nts)
        if surface:
            # Shape complementarity view: both molecular surfaces, nothing
            # else. Saved twice, the binder solid and then at 55%, so the RNA
            # groove it sits in can be read through it.
            L += [f"hide #1-{len(PANELS)} atoms",
                  f"hide #1-{len(PANELS)} cartoons",
                  f"hide #1-{len(PANELS)} surfaces", "~label",
                  f"surface {rspec}", f"surface {pspec}",
                  f"color {rspec} {palette.NEUTRAL['300']} target s",
                  f"color {pspec} {col} target s",
                  f"transparency {rspec} 0 target s",
                  f"transparency {pspec} 0 target s",
                  f"view {rspec} {pspec}", "zoom 0.95",
                  f"save {out_dir}/s{tw}_solid.png width {width} "
                  "supersample 3 transparentBackground true",
                  f"transparency {pspec} 55 target s",
                  f"save {out_dir}/s{tw}_open.png width {width} "
                  "supersample 3 transparentBackground true"]
            continue
        L += [f"hide #1-{len(PANELS)} atoms", f"hide #1-{len(PANELS)} cartoons",
              "~label",
              f"cartoon {rspec}", f"cartoon {pspec}",
              f"cartoon style {pspec} thickness 0.4",
              f"color {rspec} {palette.NEUTRAL['400']}",
              f"color {pspec} {col}",
              f"transparency {pspec} 65 target c",
              f"transparency {rspec} 25 target c",
              f"show {rspec} atoms", f"nucleotides {rspec} fill",
              f"color {rspec} {palette.NEUTRAL['400']} target ac",
              f"transparency {rspec} 60 target a",
              f"color {me}/{rna}:{nt} {NT_C} target ac",
              f"transparency {me}/{rna}:{nt} 0 target ac",
              f"label {me}/{rna}:{nt} residues height 1.7 color {NT_C}"]
        # 2026-11-18 (Leo): "remove all the highlighted amino acids". The
        # side-chain sticks and their residue labels are gone. The DASHES went
        # with them, not as a separate decision: every dash terminated on an
        # atom of one of those side chains, so leaving them would have drawn
        # H-bonds to a cartoon with nothing at the other end.
        #
        # Kept: the contacted NUCLEOTIDE's colour and label. It is not an
        # amino acid, and it still says which base each binder reads. Say so
        # if it should go too - it is one line.
        #
        # `dash`, RES_C and the hbonds provenance in this module's docstring
        # are now unused by the render but left in place: the contacts were
        # verified with ChimeraX `hbonds` (verify_hbonds.cxc / .log) and that
        # record should not be deleted just because the dashes are not drawn.
        angles = SWEEP if sweep else [chosen]
        prev = 0
        for th in angles:
            if th != prev:
                L.append(f"turn y {th - prev} center 0,0,0 models {allm}")
                prev = th
            suf = f"_y{th:03d}" if sweep else ""
            L += [f"view {pspec}:{sel} {me}/{rna}:{nt}", "zoom 1.15",
                  f"save {out_dir}/t{tw}{suf}.png width {width} "
                  "supersample 3 transparentBackground true"]
        # (the matching teardown for those dashes went with them; leaving
        # the loop here was a NameError on the first panel)
        L.append(f"turn y {-prev} center 0,0,0 models {allm}")
    L.append("exit")
    return "\n".join(L)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("out_dir")
    ap.add_argument("--width", type=int, default=2400)
    ap.add_argument("--sweep", action="store_true")
    ap.add_argument("--surface", action="store_true",
                    help="shape-complementarity view of the top three, both "
                         "surfaces, no sticks or labels")
    a = ap.parse_args()
    out = Path(a.out_dir).resolve()
    out.mkdir(parents=True, exist_ok=True)
    script = out / ("render_surface.cxc" if a.surface else "render_besthits.cxc")
    script.write_text(build(out, a.width, a.sweep, a.surface) + "\n")
    r = subprocess.run([CHIMERAX, "--exit", "--cmd", f"open {script}"],
                       capture_output=True, text=True)
    for ln in (r.stdout + r.stderr).splitlines():
        if "rror" in ln:
            print("  cx|", ln)
    print("wrote", len([p for p in out.glob('*.png')
                        if not p.name.startswith('._')]), "png")


if __name__ == "__main__":
    main()
