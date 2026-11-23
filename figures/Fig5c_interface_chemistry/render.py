"""Fig5c: the two interface-chemistry metrics from HANDOVER_METRICS.md, shown.

Metric 1, the acid-phosphate kill flag: a Asp/Glu carboxylate pressed against
the RNA backbone phosphate. Like charges repel, and the flag is binary:
0 of 15 binders carry one, 7 of 17 non-binders do (Fisher p = 0.0077).

Metric 2, the Csy4 direction: natural Csy4 reads RNA with arginine and not
lysine (Arg 0.127 / Lys 0.021 at its interface) while the library's interfaces
average 15.7% Lys. The top binder sits at Arg 0.50 / Lys 0.05.

Two structures, chosen as the clearest case of each:

  acid   csy4free_BG_BIF_len80_d503930b_design39_seq0    FC 0.84, flagged
         GLU4 carboxylate on A1's phosphate: OE2-OP2 1.45 A, OE1-OP1 2.03 A.
         Those are overlap distances, not contacts. ASP51 sits 4.24 A off C2.
  clean  csy4free_BG_LiMPNN_len80_e8cb617b_design26_seq2 FC 1.56, top binder
         No carboxylate within 5 A of any phosphate. ARG12 salt-bridges C6's
         OP2 at 3.31 A.

NOT claimed: that ARG12 is bidentate. It makes one NH1-OP2 contact under 3.6 A.
The bidentate geometry in the handover is the general rationale for preferring
Arg, not a measurement on this model.

All measurements are on the AF3 refold, which is where the metric is defined
and the only model with sidechains for the LiMPNN stacks.

House style: code/common/CHIMERAX_RENDERING.md.
`NEED:` the repulsive pair is drawn in the diverging warm end, which COLORS.md
rule 4 reserves for signed data. CHIMERAX_RENDERING.md already records this as
an open NEED for clash colouring; resolve by adding a synbio entry to
palette.json under rule 6.
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
CXC_COLOR = {"BG+BIF": "thesis_bg_bif", "BG+LiMPNN": "thesis_bg_limpnn",
             "RFD3+LiMPNN": "thesis_rfd3_limpnn",
             "BD1+LiMPNN": "thesis_boltzdesign_limpnn",
             "HD+LiMPNN": "thesis_hallu_limpnn"}
AF3 = 2
WARN = palette.DIVERGING[-1]          # NEED, see module docstring

PANELS = [
    ("acid", "d503930b_design39_seq0", 0.842,
     [4, 51], [1, 2], "repulsion",
     "GLU4 carboxylate on the A1 phosphate, 1.45 A"),
    ("clean", "e8cb617b_design26_seq2", 1.563,
     [12], [6], "basic",
     "ARG12 salt-bridges the C6 phosphate, no carboxylate within 5 A"),
]


def probe_model(path, model):
    a = pdbio.PDBFile.read(str(path)).get_structure(model=model)
    a = a[a.element != "H"]
    rna_c = [str(c) for c in np.unique(a.chain_id)
             if set(map(str, np.unique(a[a.chain_id == c].res_name))) <= RNA_RES][0]
    pro_c = [str(c) for c in np.unique(a.chain_id) if str(c) != rna_c][0]
    rna = a[a.chain_id == rna_c]
    ids, names = struc.get_residues(rna)
    off = "".join(map(str, names)).find(CANON)
    return pro_c, rna_c, int(ids[off]), int(ids[off + 15])


def build(out_dir, width):
    files = [next(CAND.glob(f"Q1_*{k}.pdb")) for _, k, *_ in PANELS]
    geo = [probe_model(p, AF3) for p in files]
    for p, g in zip(files, geo):
        print(f"  {p.name[:44]:<44} pro={g[0]} rna={g[1]} {g[2]}-{g[3]}")

    L = [f"open {ROOT / 'config' / 'thesis_colors.cxc'}"]
    L += [f"open {p}" for p in files]
    L += ["set bgColor white", "hide atoms", "hide cartoons", "hide surfaces"]

    rpro, rrna, r0, r1 = geo[0]
    ref = f"#1.{AF3}/{rrna}:{r0}-{r1}"
    for i, (_, rna, a0, a1) in enumerate(geo, start=1):
        if i == 1:
            continue
        L.append(f"align #{i}.{AF3}/{rna}:{a0}-{a1}@C1' toAtoms {ref}@C1' "
                 "move structures")
    centre, v, ang = canonical_orientation(files[0])
    allm = f"#1-{len(files)}"
    for ax, dd in zip("xyz", -centre):
        L.append(f"move {ax} {dd:.4f} models {allm}")
    L.append(f"turn {v[0]:.6f},{v[1]:.6f},{v[2]:.6f} {ang:.4f} center 0,0,0 "
             f"models {allm}")
    L += ["lighting soft",
          f"graphics silhouettes true width 1.5 color {palette.NEUTRAL['700']}"]

    # ------------------------------------------------------------ overviews
    L += [f"surface {ref}", f"color {ref} {palette.NEUTRAL['300']} target s",
          f"hide #1.{AF3}/{rrna} cartoons"]
    for i, (_, key, *_rest) in enumerate(PANELS, start=1):
        did = files[i - 1].name.split("csy4free_")[1][:-4]
        L.append(f"color #{i}.{AF3}/{geo[i-1][0]} "
                 f"{CXC_COLOR[stack_of('csy4free_' + did)]}")
    L += ["view", "zoom 0.85"]
    binders = " ".join(f"#{i}.{AF3}/{geo[i-1][0]}"
                       for i in range(1, len(PANELS) + 1))
    for i, (tag, *_r) in enumerate(PANELS, start=1):
        L += [f"hide {binders} cartoons",
              f"show #{i}.{AF3}/{geo[i-1][0]} cartoons",
              f"save {out_dir}/{tag}_overview.png width {width} "
              "supersample 3 transparentBackground true"]
    L += [f"hide {binders} cartoons", f"hide {ref} surfaces"]

    # -------------------------------------------------------------- zooms
    for i, (tag, key, fc, res, nts, kind, _) in enumerate(PANELS, start=1):
        pro, rna, a0, a1 = geo[i - 1]
        me = f"#{i}.{AF3}"
        rspec, pspec = f"{me}/{rna}:{a0}-{a1}", f"{me}/{pro}"
        did = files[i - 1].name.split("csy4free_")[1][:-4]
        col = CXC_COLOR[stack_of("csy4free_" + did)]
        sel = ",".join(map(str, res))
        nt = ",".join(str(n + a0 - 1) for n in nts)
        hl = WARN if kind == "repulsion" else palette.THEME["700"]
        L += [f"cartoon {rspec}", f"cartoon {pspec}",
              f"cartoon style {pspec} thickness 0.4",
              f"color {rspec} {palette.NEUTRAL['400']}",
              f"color {pspec} {col}",
              f"transparency {pspec} 80 target c",
              f"transparency {rspec} 30 target c",
              f"show {rspec} atoms", f"nucleotides {rspec} fill",
              f"color {rspec} {palette.NEUTRAL['400']} target a",
              f"transparency {rspec} 50 target a",
              # the nucleotides whose backbone is involved
              f"color {me}/{rna}:{nt} {palette.NEUTRAL['600']} target ac",
              f"transparency {me}/{rna}:{nt} 0 target ac",
              f"style {me}/{rna}:{nt} stick",
              f"show {me}/{rna}:{nt}@P,OP1,OP2,O5',O3' atoms",
              # the binder residues the metric is about
              f"show {pspec}:{sel} atoms", f"style {pspec}:{sel} stick",
              f"size {pspec}:{sel} stickRadius 0.24",
              f"color {pspec}:{sel} {hl} target a",
              f"transparency {pspec}:{sel} 0 target a"]
        if kind == "repulsion":
            L += [f"distance {pspec}:4@OE2 {me}/{rna}:{a0}@OP2 "
                  f"color {WARN} radius 0.05 dashes 6",
                  f"distance {pspec}:4@OE1 {me}/{rna}:{a0}@OP1 "
                  f"color {WARN} radius 0.05 dashes 6"]
        else:
            # the same measurement as the acid panel, opposite chemistry, so
            # the two zooms are read the same way
            L.append(f"distance {pspec}:12@NH1 {me}/{rna}:{a0+5}@OP2 "
                     f"color {palette.THEME['700']} radius 0.05 dashes 6")
        L += [f"view {pspec}:{sel} {me}/{rna}:{nt}", "zoom 0.95",
              f"save {out_dir}/{tag}_zoom.png width {width} "
              "supersample 3 transparentBackground true",
              "~distance", "~hbonds",
              f"hide {me} atoms", f"hide {me} cartoons"]
    L.append("exit")
    return "\n".join(L)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("out_dir")
    ap.add_argument("--width", type=int, default=3000)
    a = ap.parse_args()
    out = Path(a.out_dir).resolve()
    out.mkdir(parents=True, exist_ok=True)
    script = out / "render_chem.cxc"
    script.write_text(build(out, a.width) + "\n")
    r = subprocess.run([CHIMERAX, "--exit", "--cmd", f"open {script}"],
                       capture_output=True, text=True)
    for ln in (r.stdout + r.stderr).splitlines():
        if "rror" in ln:
            print("  cx|", ln)
    pngs = sorted(p.name for p in out.glob("*.png")
                  if not p.name.startswith("._"))
    print(f"wrote {len(pngs)} png: {pngs}")


if __name__ == "__main__":
    main()
