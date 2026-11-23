#!/usr/bin/env python3
"""Which pi-stacking binders stack on the TOP flipped-out nucleotides?

pi_stack_score in the store is a single number with no position, so it cannot
distinguish a stack on the flipped-out U10 at the top of the hairpin from one
on G16 at the bottom. This measures the position directly in each AF3 refold.

Flipped-out is computed, not assumed: a base is flipped out when its ring
centroid sits far off the helix axis fitted to the six Watson-Crick pair
midpoints. In the reference hairpin U10 is 13.0 A off-axis at the TOP and A1 is
13.6 A off-axis at the BOTTOM, while every stem base sits at 2-5 A.

A stack is a protein aromatic ring (PHE/TYR/TRP/HIS) whose centroid is within
5.0 A of a base ring centroid with interplanar angle under 30 degrees.

Usage: python find_loop_stacker.py <refold_dir>
"""
import sys
from pathlib import Path
import numpy as np
import biotite.structure as struc
import biotite.structure.io.pdbx as pdbx

RNA_RES = {"A", "U", "G", "C"}
CANON = "ACUGCCGUAUAGGCAG"
PAIRS = [(2, 16), (3, 15), (4, 14), (5, 13), (6, 12), (7, 11)]
TOP = {8, 9, 10}
RING = {"A": ["N9", "C8", "N7", "C5", "C6", "N1", "C2", "N3", "C4"],
        "G": ["N9", "C8", "N7", "C5", "C6", "N1", "C2", "N3", "C4"],
        "U": ["N1", "C2", "N3", "C4", "C5", "C6"],
        "C": ["N1", "C2", "N3", "C4", "C5", "C6"]}
AROM = {"PHE": ["CG", "CD1", "CD2", "CE1", "CE2", "CZ"],
        "TYR": ["CG", "CD1", "CD2", "CE1", "CE2", "CZ"],
        "TRP": ["CD2", "CE2", "CE3", "CZ2", "CZ3", "CH2"],
        "HIS": ["CG", "ND1", "CD2", "CE1", "NE2"]}


def normal(c):
    return np.linalg.svd(c - c.mean(0))[2][2]


def scan(path):
    a = pdbx.get_structure(pdbx.CIFFile.read(str(path)), model=1)
    a = a[a.element != "H"]
    rc = [str(c) for c in np.unique(a.chain_id)
          if set(map(str, np.unique(a[a.chain_id == c].res_name))) <= RNA_RES]
    rna = a[np.isin(a.chain_id, rc)]
    pro = a[~np.isin(a.chain_id, rc)]
    ids, names = struc.get_residues(rna)
    off = "".join(map(str, names)).find(CANON)
    if off < 0:
        raise ValueError("canonical frame not found")
    cen, nrm = {}, {}
    for k in range(1, 17):
        rid = ids[off + k - 1]
        sel = rna[(rna.res_id == rid) & np.isin(rna.atom_name, RING[CANON[k - 1]])]
        if len(sel) < 5:
            continue
        cen[k], nrm[k] = sel.coord.mean(0).astype(float), normal(sel.coord.astype(float))
    mids = np.array([(cen[i] + cen[j]) / 2 for i, j in PAIRS if i in cen and j in cen])
    centre = np.array(list(cen.values())).mean(0)
    axis = np.linalg.svd(mids - mids.mean(0))[2][0]
    if 9 in cen and np.dot(cen[9] - centre, axis) < 0:
        axis = -axis

    out = []
    for pid in np.unique(pro.res_id):
        pr = pro[pro.res_id == pid]
        rn = str(pr.res_name[0]).strip()
        if rn not in AROM:
            continue
        sel = pr[np.isin(pr.atom_name, AROM[rn])]
        if len(sel) < 5:
            continue
        ac, an = sel.coord.mean(0).astype(float), normal(sel.coord.astype(float))
        for k in cen:
            d = float(np.linalg.norm(ac - cen[k]))
            if d > 5.0:
                continue
            ang = float(np.degrees(np.arccos(min(1.0, abs(np.dot(an, nrm[k]))))))
            if ang >= 30:
                continue
            v = cen[k] - centre
            h = float(np.dot(v, axis))
            rad = float(np.linalg.norm(v - h * axis))
            out.append(dict(res=f"{rn}{pid}", nt=f"{CANON[k-1]}{k}", k=k,
                            d=round(d, 2), ang=round(ang, 1),
                            height=round(h, 1), radial=round(rad, 1),
                            top=k in TOP))
    return out


if __name__ == "__main__":
    d = Path(sys.argv[1])
    hits = []
    for f in sorted(p for p in d.glob("*.cif") if not p.name.startswith("._")):
        try:
            st = scan(f)
        except Exception as e:
            print(f"SKIP {f.name[:50]}: {type(e).__name__}: {e}")
            continue
        did = f.name.replace("_AF3_seed1_sample0.cif", "")
        top = [s for s in st if s["top"]]
        if st:
            tag = "TOP-LOOP" if top else "other"
            for s in (top or st):
                print(f"{tag:9s} {did[:46]:48s} {s['res']:>7s} . {s['nt']:>4s} "
                      f"d={s['d']:<5} ang={s['ang']:<5} height={s['height']:<6} "
                      f"radial={s['radial']}")
        if top:
            hits.append((did, top))
    print(f"\n{len(hits)} designs stack on the top loop (nt 8-10)")
