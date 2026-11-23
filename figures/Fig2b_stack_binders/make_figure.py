#!/usr/bin/env python3
"""Figure 3b: four binders, one per generator, in the Fig 2b shared RNA frame.

Same frame as render_shared_frame.py: every complex is aligned by C1' onto the
SAME reference hairpin (panel a's design) and only the reference RNA is drawn,
so the hairpin sits identically in every frame of every panel and only the
binder moves. Binder colour is the stack colour, so colour means generator.

Designs are backbone-only CIFs, so these are cartoon overviews, no sidechains.

Usage: python render_stack_panels.py <out_dir> [--width 1200] [--frames 12]
"""
import argparse
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from common import palette                                      # noqa: E402
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "Fig2b_map_callouts"))
from render_shared_frame import (CHIMERAX, BASE, canonical_orientation,
                                 probe, stack_of)               # noqa: E402

# the hairpin frame every Fig 2b panel already uses
REF = "csy4free_RFD3_LiMPNN_len80_8b13398e_design50_seq0"

# "library": the four 12k picks. All in csy4_binder_library_12k_v6_k30.csv
# (rfd3 n=5918, bgbif n=960, bd1 n=6112, hd n=502).
# "clusters": the largest backbone cluster per generator. The BG entry is
# BG+LiMPNN, not BG+BIF: both stacks carry a design called
# 5de20158_design12_seq0 but they are different poses (20.9 A max CA
# separation), and only the LiMPNN one reproduces the 37.5 / 26.2, 2H / 3E
# secondary structure the cluster table quotes.
SETS = {
    "library": [
        ("rfd3", "csy4free_RFD3_LiMPNN_len80_88b66226_design33_seq1"),
        ("bgbif", "csy4free_BG_BIF_len80_ef8ee7a5_design8_seq2"),
        ("bd1", "csy4free_BD1_LiMPNN_len80_8bca6080_design31_seq0"),
        ("hd", "csy4free_HD_LiMPNN_len80_1ef1a7dc_recycle11_sample1_pack1"),
    ],
    "clusters": [
        ("rfd3_cl3", "csy4free_RFD3_LiMPNN_len80_1f19b627_design7_seq0"),
        ("bg_cl9", "csy4free_BG_LiMPNN_len80_5de20158_design12_seq0"),
        ("bd1_cl0", "csy4free_BD1_LiMPNN_len80_f2bbc6ab_design16_seq0"),
        ("hd_cl1", "csy4free_HD_LiMPNN_len80_a083553a_recycle8_sample0_pack1"),
    ],
    # Second round of cluster representatives, 2026-11-07. Here BOTH BG stacks
    # carry 97eb1c8c_design25_seq0 with the IDENTICAL fold (37.5 / 27.5,
    # 2H / 3E in both), so unlike the first round the secondary structure in
    # the cluster table cannot say which is meant. They are docked 22.1 A apart
    # (max CA), so the panels differ. Both are rendered; pick one.
    "clusters2": [
        ("rfd3_cl3", "csy4free_RFD3_LiMPNN_len80_8eaafbb8_design14_seq0"),
        ("bg_cl9_limpnn", "csy4free_BG_LiMPNN_len80_97eb1c8c_design25_seq0"),
        ("bg_cl9_bif", "csy4free_BG_BIF_len80_97eb1c8c_design25_seq0"),
        ("bd1_cl0", "csy4free_BD1_LiMPNN_len80_0ba0a531_design20_seq0"),
        ("hd_cl1", "csy4free_HD_LiMPNN_len80_f70e1daf_recycle11_sample0_pack1"),
    ],
    # Third round, 2026-11-07, drawn from a candidate pool rather than raw
    # cluster size: BD1, RFD3, HD and BG+LiMPNN all pass every gate term, and
    # the HD pick is in the 12k. BG is BG+LiMPNN: the BG+BIF design of the same
    # name fails ligand_RMSD (13.1) and DockQ (0.172), so it is not a candidate.
    "clusters3": [
        ("rfd3", "csy4free_RFD3_LiMPNN_len80_8eaafbb8_design14_seq0"),
        ("bg", "csy4free_BG_LiMPNN_len80_81878aac_design15_seq0"),
        ("bd1", "csy4free_BD1_LiMPNN_len80_20a7e2a8_design7_seq0"),
        ("hd", "csy4free_HD_LiMPNN_len80_58be5280_recycle15_sample0_pack1"),
    ],
    # Eight map callouts, 2026-11-08, chosen by pick8.py on Euler: backbone
    # representatives that pass every gate term, scored on pose then interface,
    # two per generator forced >4 UMAP units apart. All eight are in the 12k.
    # Tag suffix is the fold cluster, so the tag says where on the map it sits.
    "map8": [
        ("rfd3_a", "csy4free_RFD3_LiMPNN_len80_4bcad7a8_design62_seq0"),
        ("rfd3_b", "csy4free_RFD3_LiMPNN_len80_7c0a8297_design17_seq0"),
        ("bg_a", "csy4free_BG_LiMPNN_len80_aebd763e_design11_seq0"),
        ("bg_b", "csy4free_BG_LiMPNN_len80_e1e86ec3_design9_seq0"),
        ("bd1_a", "csy4free_BD1_LiMPNN_len80_b6ccc04e_design28_seq0"),
        ("bd1_b", "csy4free_BD1_LiMPNN_len80_731e2185_design57_seq0"),
        ("hd_a", "csy4free_HD_LiMPNN_len80_98a6ba23_recycle9_sample0_pack1"),
        ("hd_b", "csy4free_HD_LiMPNN_len80_13538ba3_recycle15_sample0_pack1"),
    ],
}

# The 32 ordered candidates (LV2.27 Echo order, screened in LV2.33). Read from
# working/structures/hits32/index.tsv so the set can never drift from the order
# itself. These are PDB complexes WITH sidechains, unlike the backbone-only
# design CIFs; tag is the index number, which is the ON/OFF rank.
HITS32 = Path(__file__).resolve().parents[2] / "working/structures/hits32"


def _hits32():
    rows = [r.split("\t") for r in
            HITS32.joinpath("index.tsv").read_text().splitlines()[1:] if r.strip()]
    return [(f.split("_", 1)[0], f.split("_", 1)[1][:-4]) for f, _, _ in rows]


SETS["hits32"] = _hits32()
PANELS = SETS["clusters3"]            # the published set


def _find(design):
    cif = BASE / "cifs" / f"{design}.cif"
    return cif if cif.exists() else HITS32 / f"{design}.pdb"


def build(out_dir, width, frames, only=None):
    ids = [REF] + [d for _, d in PANELS]
    files = [_find(d) for d in ids]
    meta = [probe(p) for p in files]                 # (pro, rna, res0, res1)

    L = [f"open {p}" for p in files]
    L += ["set bgColor white", "hide atoms", "cartoon"]

    rpro, rrna, r0, r1 = meta[0]
    for i in range(2, len(files) + 1):
        _, rna, a0, a1 = meta[i - 1]
        L.append(f"align #{i}/{rna}:{a0}-{a1}@C1' toAtoms "
                 f"#1/{rrna}:{r0}-{r1}@C1' move structures")

    # one shared RNA surface: the reference copy only
    L += [f"surface #1/{rrna}:{r0}-{r1}",
          f"color #1/{rrna} #d7dcdf target s",
          f"hide #1/{rrna} cartoons",
          f"hide #1/{rpro} cartoons"]                # ref binder never drawn
    for i in range(2, len(files) + 1):
        rna = meta[i - 1][1]
        L += [f"hide #{i}/{rna} cartoons", f"hide #{i}/{rna} atoms"]

    for i in range(2, len(files) + 1):
        L.append(f"color #{i}/{meta[i-1][0]} "
                 f"{palette.stack_color(stack_of(ids[i-1]))}")

    centre, v, ang = canonical_orientation(files[0])
    allm = f"#1-{len(files)}"
    for ax, d in zip("xyz", -centre):
        L.append(f"move {ax} {d:.4f} models {allm}")
    L.append(f"turn {v[0]:.6f},{v[1]:.6f},{v[2]:.6f} {ang:.4f} center 0,0,0 "
             f"models {allm}")

    L += ["lighting soft",
          f"graphics silhouettes true width 1.5 color {palette.NEUTRAL['700']}",
          "view", "zoom 1.0"]

    binders = [f"#{i}/{meta[i-1][0]}" for i in range(2, len(files) + 1)]
    step = 360 // frames
    for f in range(frames):
        if f:
            L.append(f"turn y {step}")
        for k, (tag, _) in enumerate(PANELS):
            L.append("hide " + " ".join(binders) + " cartoons")
            L.append(f"show {binders[k]} cartoons")
            if only is not None and step * f != only:
                continue
            L.append(f"save {out_dir}/{tag}_y{step*f:03d}.png width {width} "
                     "supersample 3 transparentBackground true")
    L.append("exit")
    return "\n".join(L)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("out_dir")
    ap.add_argument("--width", type=int, default=1200)
    ap.add_argument("--frames", type=int, default=12)
    ap.add_argument("--set", dest="which", choices=sorted(SETS), default="clusters3")
    ap.add_argument("--angle", type=int, default=None,
                    help="save only this y angle (must be a multiple of 360/frames)")
    a = ap.parse_args()
    global PANELS
    PANELS = SETS[a.which]
    out = Path(a.out_dir).resolve()
    out.mkdir(parents=True, exist_ok=True)
    script = out / "render_stack_panels.cxc"
    script.write_text(build(out, a.width, a.frames, a.angle) + "\n")
    print("script ->", script)
    r = subprocess.run([CHIMERAX, "--exit", "--cmd", f"open {script}"],
                       capture_output=True, text=True)
    for ln in (r.stdout + r.stderr).strip().splitlines()[-20:]:
        print("  cx|", ln)
    pngs = [p for p in out.glob("*.png") if not p.name.startswith("._")]
    n = len(PANELS) * (1 if a.angle is not None else a.frames)
    print(f"wrote {len(pngs)} png (expect {n})")


if __name__ == "__main__":
    main()
