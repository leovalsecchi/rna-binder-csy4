#!/usr/bin/env python3
"""Describe what each example binder actually does to the RNA, from the
structure itself, so the figure caption states a measured feature rather than
an asserted one.

Reports per design: contact split by moiety, per-nucleotide profile, loop vs
stem, major/minor groove edge counts, and aromatic stacking pairs.

RNA frame verified from the structures: canonical 16-mer ACUGCCGUAUAGGCAG,
pairs (2,16)(3,15)(4,14)(5,13)(6,12)(7,11), so loop = nt 8,9,10.
Chain roles detected by content (RFD3 design files are flipped).
"""
import sys
from pathlib import Path
import numpy as np
import biotite.structure as struc
import biotite.structure.io.pdbx as pdbx

CUT = 4.5
RNA_RES = {"A", "U", "G", "C"}
LOOP, STEM = [8, 9, 10], [2, 3, 4, 5, 6, 7, 11, 12, 13, 14, 15, 16]
CANON = "ACUGCCGUAUAGGCAG"

PHOS = {"P", "OP1", "OP2", "OP3", "O1P", "O2P", "O3P", "O5'", "O3'"}
SUGAR = {"C1'", "C2'", "C3'", "C4'", "C5'", "O2'", "O4'"}
MAJOR = {"A": {"N7", "N6", "C5"}, "G": {"N7", "O6", "C5"},
         "U": {"O4", "C5", "C6"}, "C": {"N4", "C5", "C6"}}
MINOR = {"A": {"N3", "C2", "O2'"}, "G": {"N3", "N2", "C2", "O2'"},
         "U": {"O2", "O2'"}, "C": {"O2", "O2'"}}
AROM = {"PHE": ["CG", "CD1", "CD2", "CE1", "CE2", "CZ"],
        "TYR": ["CG", "CD1", "CD2", "CE1", "CE2", "CZ"],
        "TRP": ["CD2", "CE2", "CE3", "CZ2", "CZ3", "CH2"],
        "HIS": ["CG", "ND1", "CD2", "CE1", "NE2"]}
BASE_RING = {"A": ["N9", "C8", "N7", "C5", "C6", "N1", "C2", "N3", "C4"],
             "G": ["N9", "C8", "N7", "C5", "C6", "N1", "C2", "N3", "C4"],
             "U": ["N1", "C2", "N3", "C4", "C5", "C6"],
             "C": ["N1", "C2", "N3", "C4", "C5", "C6"]}


def moiety(n):
    return "phosphate" if n in PHOS else ("sugar" if n in SUGAR else "base")


def ring_normal(coords):
    c = coords - coords.mean(0)
    return np.linalg.svd(c)[2][2]


def analyse(path):
    a = pdbx.get_structure(pdbx.CIFFile.read(str(path)), model=1)
    a = a[a.element != "H"]
    rc = [c for c in np.unique(a.chain_id)
          if set(map(str, np.unique(a[a.chain_id == c].res_name))) <= RNA_RES]
    rna = a[np.isin(a.chain_id, rc)]
    pro = a[~np.isin(a.chain_id, rc)]

    ids, names = struc.get_residues(rna)
    seq = "".join(map(str, names))
    off = seq.find(CANON)
    if off < 0:
        raise ValueError(f"canonical frame not in {seq}")
    pos = {rid: i - off + 1 for i, rid in enumerate(ids)}

    # brute force: a few hundred atoms either side, so no kd-tree needed
    d2 = ((rna.coord[:, None, :].astype(float)
           - pro.coord[None, :, :].astype(float)) ** 2).sum(-1)
    cnt = (d2 <= CUT ** 2).sum(axis=1)

    moi = {"phosphate": 0, "sugar": 0, "base": 0}
    per = {i: 0 for i in range(1, 17)}
    mj = mn = 0
    for i, c in enumerate(cnt):
        if not c:
            continue
        p = pos.get(rna.res_id[i], 0)
        if not (1 <= p <= 16):
            continue
        rn, an = str(rna.res_name[i]).strip(), str(rna.atom_name[i]).strip()
        moi[moiety(an)] += c
        per[p] += c
        if an in MAJOR.get(rn, ()):
            mj += c
        elif an in MINOR.get(rn, ()):
            mn += c

    tot = sum(moi.values())
    # aromatic stacking: ring centroids < 5 A and interplanar angle < 30 deg
    stacks, near = [], []
    for rid in ids:
        p = pos.get(rid, 0)
        if not (1 <= p <= 16):
            continue
        rn = str(names[list(ids).index(rid)])
        bsel = rna[(rna.res_id == rid) & np.isin(rna.atom_name, BASE_RING[rn])]
        if len(bsel) < 5:
            continue
        bc, bnv = bsel.coord.mean(0), ring_normal(bsel.coord)
        for pid in np.unique(pro.res_id):
            pr = pro[pro.res_id == pid]
            rnm = str(pr.res_name[0]).strip()
            if rnm not in AROM:
                continue
            asel = pr[np.isin(pr.atom_name, AROM[rnm])]
            if len(asel) < 5:
                continue
            ac, anv = asel.coord.mean(0), ring_normal(asel.coord)
            d = np.linalg.norm(ac - bc)
            if d > 7.5:
                continue
            ang = np.degrees(np.arccos(min(1, abs(np.dot(bnv, anv)))))
            near.append((round(float(d), 2), round(float(ang), 1),
                         f"{rnm}{pid}", f"{CANON[p-1]}{p}"))
            # PLIP-standard face-to-face: BOTH the 5.5 A centroid distance and
            # the 30 deg interplanar angle. The angle alone admitted pairs out to
            # 7.5 A, which are not stacks (e.g. HIS55-G4 at 7.28 A).
            if d <= 5.5 and ang < 30:
                stacks.append((f"{rnm}{pid}", f"{rn}{p}", round(float(d), 2),
                               round(float(ang), 1)))

    return dict(seq=seq, off=off, total=tot, moi=moi, per=per, major=mj,
                minor=mn, stacks=stacks, near=sorted(near)[:3])


if __name__ == "__main__":
    for p in sys.argv[1:]:
        r = analyse(Path(p))
        t = max(r["total"], 1)
        loop = sum(r["per"][i] for i in LOOP)
        stem = sum(r["per"][i] for i in STEM)
        top = sorted(r["per"].items(), key=lambda kv: -kv[1])[:4]
        print(f"\n### {Path(p).stem}")
        print(f"  contacts {r['total']}   phosphate {100*r['moi']['phosphate']/t:.0f}% "
              f"sugar {100*r['moi']['sugar']/t:.0f}%  base {100*r['moi']['base']/t:.0f}%")
        print(f"  backbone(P+S) {100*(r['moi']['phosphate']+r['moi']['sugar'])/t:.0f}%"
              f"   loop(8-10) {100*loop/t:.0f}%   stem {100*stem/t:.0f}%")
        print(f"  groove edges: major {r['major']}  minor {r['minor']}  -> "
              + ("MAJOR" if r['major'] >= 2*r['minor'] else
                 "MINOR" if r['minor'] >= 2*r['major'] else "mixed"))
        print("  top nt: " + ", ".join(f"{CANON[i-1]}{i}={c}" for i, c in top))
        print(f"  aromatic stacks (<=5.5A, <30deg): {r['stacks'] if r['stacks'] else 'NONE'}")
        print(f"  nearest arom/base ring pairs (d, angle): {r['near'] if r['near'] else 'none within 7.5A'}")
        print(f"  frame offset {r['off']} (RNA len {len(r['seq'])})")
