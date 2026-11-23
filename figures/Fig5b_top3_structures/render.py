"""Fig5b: the three LV2.33 candidates with the highest reporter fold change.

Ranked on the readout Leo asked for: mean %ON over the two hairpin targets
(posA, posB) divided by the mean over the two negatives (NegA, Neg_noScaf),
recomputed from LV_FlowCytometry/analysis/out/variants.csv.

  V15  BG+BIF          1.38x
  V16  RFD3+LiMPNN     1.35x
  V26  BD1+LiMPNN      1.34x

READ THE README BEFORE USING THIS PANEL. None of these is distinguishable from
the empirical null (all q = 1.0), the Csy4 positive control is 14x on the same
readout, and 12 of the 32 ordered candidates have no flow data at all. This is
a structure panel of the three top-ranked, not evidence that they bind.

House convention is code/common/CHIMERAX_RENDERING.md: binder in its stack
colour, RNA a light neutral, soft lighting, neutral-700 silhouettes width 1.5,
white bg, transparent PNG, supersample 3, `view` called once.

Two panel kinds per binder:
  overview  the DESIGN model, cartoon in stack colour, over the shared RNA
            surface, at three angles
  contact   the AF3 REFOLD, because two of the three designs are backbone only
            (320 atoms, no CB) so their interactions exist only in a refold.
            Nucleotides filled, contacting sidechains as sticks, cross-chain
            H-bonds as pseudobonds, everything else faded.

All nine submodels are aligned independently onto ONE reference RNA (V15's
design) by C1', so the hairpin is in the identical pose in every panel and only
the binder moves.
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
# The 3-model files (design + AF3 + ESM) live with the ordered candidates.
# working/structures/hits32/ carries the SAME REMARK header claiming three
# models but holds only model 1, so a refold read from there is silently the
# design. Checked 2026-11-13.
HITS = (ROOT.parent / "LV_BinderDesign" / "LV_32Candidates")
RNA_RES = {"A", "U", "G", "C"}
CANON = "ACUGCCGUAUAGGCAG"
CXC_COLOR = {"BG+BIF": "thesis_bg_bif", "RFD3+LiMPNN": "thesis_rfd3_limpnn",
             "BD1+LiMPNN": "thesis_boltzdesign_limpnn"}

DESIGN, AF3 = 1, 2        # submodel numbers inside each PDB

# (tag, design_id, fold change, [(zoom tag, binder residues, caption residues)])
# Residues were measured, not guessed: polar contacts <= 3.5 A and aromatic /
# base centroid pairs <= 5.0 A in the AF3 model. See README for the table.
PANELS = [
    ("V15", "csy4free_BG_BIF_len80_25a69495_design18_seq0", 1.380,
     [("hub", [67, 64, 63], [5, 7, 9, 10, 12], "ARG67 reads G7, U10 and G12"),
      ("stack", [60, 63], [7, 9, 10], "HIS60 stacks on G7")]),
    ("V16", "csy4free_RFD3_LiMPNN_len80_73d23e78_design49_seq4", 1.352,
     [("hub", [47, 46, 50], [1, 14, 15, 16],
       "ARG47 clamps A1, GLN46 and ARG50 read the helix end"),
      ("loop", [27, 34], [9, 10, 11], "HIS27 and GLU34 reach the loop at U10")]),
    ("V26", "csy4free_BD1_LiMPNN_len80_d4d5ef6f_design56_seq1", 1.335,
     [("stack", [37, 32], [8, 10], "TRP37 on U10 and TYR32 on U8"),
      ("hub", [20, 18, 50], [6, 7, 9, 10, 11],
       "ARG20 and ARG18 read C6 and G7")]),
]


def probe_model(path, model):
    """Chain roles and the 16-nt canonical range, for ONE model of a PDB.

    Needed per model because a BoltzDesign design carries an 18-nt hairpin
    while its refolds carry the 16-nt core; see CHIMERAX_RENDERING.md.
    """
    a = pdbio.PDBFile.read(str(path)).get_structure(model=model)
    a = a[a.element != "H"]
    rna_c = [str(c) for c in np.unique(a.chain_id)
             if set(map(str, np.unique(a[a.chain_id == c].res_name))) <= RNA_RES]
    pro_c = [str(c) for c in np.unique(a.chain_id) if str(c) not in rna_c]
    if len(rna_c) != 1 or len(pro_c) != 1:
        raise ValueError(f"{path.name} model {model}: rna={rna_c} prot={pro_c}")
    rna = a[a.chain_id == rna_c[0]]
    ids, names = struc.get_residues(rna)
    off = "".join(map(str, names)).find(CANON)
    if off < 0:
        raise ValueError(f"{path.name} model {model}: canonical frame not found")
    return pro_c[0], rna_c[0], int(ids[off]), int(ids[off + 15])


def build(out_dir, width):
    files = [next(HITS.glob(f"Q1_*_{did}.pdb")) for _, did, _, _ in PANELS]
    for p in files:
        if not p.exists():
            raise SystemExit(f"missing structure: {p}")
    geo = {}
    for i, p in enumerate(files, start=1):
        for m in (DESIGN, AF3):
            geo[(i, m)] = probe_model(p, m)
            print(f"  #{i}.{m} {p.name[:34]:<34} "
                  f"pro={geo[(i,m)][0]} rna={geo[(i,m)][1]} "
                  f"{geo[(i,m)][2]}-{geo[(i,m)][3]}")

    L = [f"open {ROOT / 'config' / 'thesis_colors.cxc'}"]
    L += [f"open {p}" for p in files]
    L += ["set bgColor white", "hide atoms", "hide cartoons", "hide surfaces"]

    # one shared frame: every submodel's RNA onto the reference, by C1'
    rpro, rrna, r0, r1 = geo[(1, DESIGN)]
    ref = f"#1.{DESIGN}/{rrna}:{r0}-{r1}"
    for (i, m), (_, rna, a0, a1) in geo.items():
        if (i, m) == (1, DESIGN):
            continue
        L.append(f"align #{i}.{m}/{rna}:{a0}-{a1}@C1' toAtoms {ref}@C1' "
                 "move structures")

    centre, v, ang = canonical_orientation(files[0])
    allm = "#1-3"
    for ax, dd in zip("xyz", -centre):
        L.append(f"move {ax} {dd:.4f} models {allm}")
    L.append(f"turn {v[0]:.6f},{v[1]:.6f},{v[2]:.6f} {ang:.4f} center 0,0,0 "
             f"models {allm}")
    L += ["lighting soft",
          f"graphics silhouettes true width 1.5 color {palette.NEUTRAL['700']}"]

    # ---------------------------------------------------------- overviews
    L += [f"surface {ref}", f"color {ref} {palette.NEUTRAL['300']} target s",
          f"hide #1.{DESIGN}/{rrna} cartoons"]
    for i, (_, did, _, _) in enumerate(PANELS, start=1):
        pro = geo[(i, DESIGN)][0]
        L.append(f"color #{i}.{DESIGN}/{pro} {CXC_COLOR[stack_of(did)]}")
    L += ["view", "zoom 0.85"]

    binders = " ".join(f"#{i}.{DESIGN}/{geo[(i,DESIGN)][0]}"
                       for i in range(1, len(PANELS) + 1))
    for f in range(3):
        if f:
            L.append("turn y 120")
        for i, (tag, _, _, _) in enumerate(PANELS, start=1):
            pro = geo[(i, DESIGN)][0]
            L.append(f"hide {binders} cartoons")
            L.append(f"show #{i}.{DESIGN}/{pro} cartoons")
            L.append(f"save {out_dir}/{tag}_overview_y{120*f:03d}.png "
                     f"width {width} supersample 3 transparentBackground true")
    L.append(f"hide {binders} cartoons")
    L.append(f"hide {ref} surfaces")

    # ----------------------------------------------------------- contacts
    # AF3 refold, its own RNA (already in the shared frame).
    for i, (tag, did, _, zooms) in enumerate(PANELS, start=1):
        pro, rna, a0, a1 = geo[(i, AF3)]
        me = f"#{i}.{AF3}"
        rspec, pspec = f"{me}/{rna}:{a0}-{a1}", f"{me}/{pro}"
        col = CXC_COLOR[stack_of(did)]
        for zt, res, nts, _ in zooms:
            sel = ",".join(str(r) for r in res)
            nt = ",".join(str(n + a0 - 1) for n in nts)   # 1-16 -> file numbering
            L += [f"hide {me} cartoons", f"hide {me} atoms",
                  f"cartoon {rspec}", f"cartoon {pspec}",
                  f"cartoon style {pspec} thickness 0.4",
                  f"color {rspec} {palette.NEUTRAL['400']}",
                  f"color {pspec} {col}",
                  f"transparency {pspec} 55 target c",
                  f"transparency {rspec} 35 target c",
                  f"show {rspec} atoms", f"nucleotides {rspec} fill",
                  f"color {rspec} {palette.NEUTRAL['400']} target a",
                  f"transparency {rspec} 45 target a",
                  # the nucleotides this cluster actually touches, full strength
                  f"color {me}/{rna}:{nt} {palette.NEUTRAL['600']} target ac",
                  f"transparency {me}/{rna}:{nt} 0 target ac",
                  f"style {me}/{rna}:{nt} stick",
                  # the binder sidechains, full strength
                  f"show {pspec}:{sel} atoms",
                  f"style {pspec}:{sel} stick",
                  f"size {pspec}:{sel} stickRadius 0.18",
                  f"color {pspec}:{sel} {col} target a",
                  f"transparency {pspec}:{sel} 0 target ac",
                  # relaxed criteria: AF3 models carry no hydrogens and some
                  # predicted polar contacts are tighter than ideal geometry
                  f"hbonds {pspec}:{sel} restrict {rspec} reveal true "
                  f"distSlop 0.5 angleSlop 25 "
                  f"color {palette.THEME['700']} radius 0.09 dashes 6",
                  f"view {pspec}:{sel} {me}/{rna}:{nt}", "zoom 1.25",
                  f"save {out_dir}/{tag}_contact_{zt}.png width {width} "
                  "supersample 3 transparentBackground true",
                  "~hbonds"]
        L.append(f"hide {me} atoms")
        L.append(f"hide {me} cartoons")
    L.append("exit")
    return "\n".join(L)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("out_dir")
    ap.add_argument("--width", type=int, default=3000)
    a = ap.parse_args()
    out = Path(a.out_dir).resolve()
    out.mkdir(parents=True, exist_ok=True)
    script = out / "render_top3.cxc"
    script.write_text(build(out, a.width) + "\n")
    print("script ->", script)
    r = subprocess.run([CHIMERAX, "--exit", "--cmd", f"open {script}"],
                       capture_output=True, text=True)
    log = (r.stdout + r.stderr).strip().splitlines()
    for ln in log:
        if "bond" in ln.lower() or "error" in ln.lower():
            print("  cx|", ln)
    for ln in log[-6:]:
        print("  cx|", ln)
    pngs = sorted(p.name for p in out.glob("*.png")
                  if not p.name.startswith("._"))
    print(f"wrote {len(pngs)} png "
          f"(expect {3*len(PANELS) + sum(len(z) for *_, z in PANELS)})")
    for n in pngs:
        print("   ", n)


if __name__ == "__main__":
    main()
