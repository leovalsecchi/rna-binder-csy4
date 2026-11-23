#!/usr/bin/env python3
"""Atom-level major/minor groove classification for a protein-RNA complex.

The 48-dim epitope fingerprint bins RNA atoms as phosphate / sugar / base only,
so it cannot tell major groove from minor groove. This does, by which base-edge
atoms the protein actually touches.

Edge assignment follows the standard Watson-Crick geometry: in a duplex, the
major-groove edge of a purine is N7/O6/N6, of a pyrimidine N4/O4/C5/C6; the
minor-groove edge is N3/C2/N2/O2 plus the sugar O2'. Atoms on the
Watson-Crick face itself (N1/N3 pairing atoms) are excluded, since in a paired
stem they are buried in the pair and unavailable to the protein.

Chain identity is detected by CONTENT, never by chain id: RFD3 design files
have the chains the other way round from every other stack (see
code/common/CHIMERAX_RENDERING.md).

Usage:  python groove.py <complex.cif> [...]
"""
import sys
import numpy as np
import biotite.structure as struc
import biotite.structure.io.pdbx as pdbx

CUTOFF = 4.5                      # must match binder_map.RNA_CONTACT_CUTOFF
RNA_RES = {"A", "U", "G", "C"}
PURINE = {"A", "G"}

MAJOR = {"A": {"N7", "N6", "C5"}, "G": {"N7", "O6", "C5"},
         "U": {"O4", "C5", "C6"}, "C": {"N4", "C5", "C6"}}
MINOR = {"A": {"N3", "C2", "O2'"}, "G": {"N3", "N2", "C2", "O2'"},
         "U": {"O2", "O2'"}, "C": {"O2", "O2'"}}


def split_chains(atoms):
    """Return (protein, rna) by residue content, not by chain id."""
    rna_c, pro_c = [], []
    for c in np.unique(atoms.chain_id):
        names = set(map(str, np.unique(atoms[atoms.chain_id == c].res_name)))
        (rna_c if names <= RNA_RES else pro_c).append(c)
    if not rna_c or not pro_c:
        raise ValueError(f"could not split chains: rna={rna_c} prot={pro_c}")
    return (atoms[np.isin(atoms.chain_id, pro_c)],
            atoms[np.isin(atoms.chain_id, rna_c)])


def groove_counts(path):
    a = pdbx.get_structure(pdbx.CIFFile.read(path), model=1)
    a = a[a.element != "H"]
    pro, rna = split_chains(a)
    from scipy.spatial import cKDTree
    tree = cKDTree(pro.coord.astype(float))
    n = tree.query_ball_point(rna.coord.astype(float), r=CUTOFF)
    counts = np.fromiter((len(x) for x in n), dtype=int, count=len(n))

    major = minor = 0
    for i, c in enumerate(counts):
        if not c:
            continue
        res = str(rna.res_name[i]).strip()
        atom = str(rna.atom_name[i]).strip()
        if atom in MAJOR.get(res, ()):
            major += c
        elif atom in MINOR.get(res, ()):
            minor += c
    return major, minor, int(counts.sum())


if __name__ == "__main__":
    print(f"{'file':<58} {'major':>7} {'minor':>7} {'total':>7}  verdict")
    for p in sys.argv[1:]:
        try:
            mj, mn, tot = groove_counts(p)
        except Exception as e:
            print(f"{p.split('/')[-1]:<58} {type(e).__name__}: {e}")
            continue
        edge = mj + mn
        if edge < 8:
            v = "backbone//no groove reading"
        elif mj >= 2 * mn:
            v = "MAJOR groove"
        elif mn >= 2 * mj:
            v = "MINOR groove"
        else:
            v = "mixed"
        print(f"{p.split('/')[-1]:<58} {mj:>7} {mn:>7} {tot:>7}  {v}")
