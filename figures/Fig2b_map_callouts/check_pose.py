#!/usr/bin/env python3
"""Does a refold reproduce the design's pose? Binder RMSD after RNA superposition.

The exact check that caught the fault: superpose on the RNA, then measure the
BINDER where it sits. A refold can fold the binder correctly and still dock it
somewhere else, and only this comparison sees that.

    python3 check_pose.py <design.cif> <model.cif> [<model2.cif> ...]

Prints, per model, the RNA fit quality, the in-place binder CA RMSD (the pose
number), and the binder's own fold RMSD after superposing the binder on itself
(the fold number) -- a design can be reproduced as a shape and misplaced as a
pose, and the two columns separate those. Protein is matched by residue order
with an equal-sequence assert, since refolding does not change the sequence; the
backbone-only design files have no side chains, so CA carries the comparison.
"""
import sys
from pathlib import Path

import numpy as np
import biotite.structure as struc
import biotite.structure.io.pdbx as pdbx

RNA_RES = {"A", "U", "G", "C"}
AA3 = {"ALA","ARG","ASN","ASP","CYS","GLN","GLU","GLY","HIS","ILE","LEU","LYS",
       "MET","PHE","PRO","SER","THR","TRP","TYR","VAL"}


def split(path):
    a = pdbx.get_structure(pdbx.CIFFile.read(str(path)), model=1)
    a = a[a.element != "H"]
    rc = [c for c in np.unique(a.chain_id)
          if set(map(str, np.unique(a[a.chain_id == c].res_name))) <= RNA_RES]
    rna = a[np.isin(a.chain_id, rc)]
    pro = a[~np.isin(a.chain_id, rc)]
    return rna, pro[np.isin(pro.res_name, list(AA3))]


def rna_index(r):
    ids, _ = struc.get_residues(r)
    rank = {rid: k for k, rid in enumerate(ids)}
    return {(rank[int(r.res_id[i])], str(r.atom_name[i]).strip()): i
            for i in range(r.array_length())}


def ca(pro):
    ids, names = struc.get_residues(pro)
    xyz = np.array([pro[(pro.res_id == r) & (pro.atom_name == "CA")].coord[0]
                    for r in ids])
    return ids, [str(n) for n in names], xyz


def report(design, model):
    rd, pd_ = split(design)
    rm, pm = split(model)
    a, b = rna_index(rd), rna_index(rm)
    keys = sorted(set(a) & set(b))
    if len(keys) < 30:
        return "only %d shared RNA atoms" % len(keys)
    # transform maps the MODEL onto the DESIGN's frame
    fitted, tf = struc.superimpose(rd[[a[k] for k in keys]],
                                   rm[[b[k] for k in keys]])
    rna_fit = float(struc.rmsd(rd[[a[k] for k in keys]], fitted))

    _, sd, ca_d = ca(pd_)
    _, sm, ca_m = ca(pm)
    if len(sd) != len(sm):
        return "protein length differs: design %d vs model %d" % (len(sd), len(sm))
    if sd != sm:
        n = sum(x != y for x, y in zip(sd, sm))
        return "protein sequence differs at %d of %d residues" % (n, len(sd))

    pose = np.linalg.norm(ca_d - tf.apply(ca_m), axis=1)
    pose_rmsd = float(np.sqrt((pose ** 2).mean()))
    # fold alone: best-fit the binder onto itself, ignoring where it sits
    fold = struc.superimpose(ca_d, ca_m)[0]
    fold_rmsd = float(np.sqrt((np.linalg.norm(ca_d - fold, axis=1) ** 2).mean()))
    verdict = "REPRODUCES" if pose_rmsd <= 3.0 else "DOES NOT reproduce"
    return ("RNA fit %.2f A / %d atoms | POSE %6.2f A (max %5.1f) | fold %.2f A "
            "| n=%d  -> %s" % (rna_fit, len(keys), pose_rmsd, pose.max(),
                               fold_rmsd, len(pose), verdict))


if __name__ == "__main__":
    design = Path(sys.argv[1])
    print("design: %s" % design.name)
    for m in sys.argv[2:]:
        print("  vs %-52s %s" % (Path(m).name, report(design, Path(m))))
