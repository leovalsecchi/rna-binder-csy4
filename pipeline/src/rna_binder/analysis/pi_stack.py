#!/usr/bin/env python3
"""pi_stack.py — protein<->RNA aromatic pi-stacking descriptors (Biotite).

Part of the analyses/ metric pipeline (sibling of interface_geom.py). Detects and
scores pi-stacking between protein aromatic / planar side chains and RNA
nucleobases -- the dominant specific-recognition force many strong RNA-binders
exploit. Chain splitting reuses rna_binder.analysis.structural.split_chains
(RNA = the chain whose residues subset {A,C,G,U}), exactly like interface_geom.

pi systems considered
  RNA base (the pi partner on the RNA side) -- every nucleobase is aromatic; the
  ribose + phosphate have no pi system, so stacking is base-mediated:
    purines  A, G : bicyclic (imidazole + pyrimidine), large polarizable pi -> strong
    pyrimidines C, U : monocyclic six-membered ring, smaller pi             -> weaker
  Protein side chains:
    TRP indole (bicyclic, largest) > TYR phenol ~ PHE benzene > HIS imidazole
    ARG guanidinium -- NOT a ring, but a planar delocalised pi that stacks on
        nucleobases very frequently in RNA interfaces (kept, down-weighted).
  (LYS and other non-planar groups do cation-pi but not true face stacking, so
   they are excluded here -- they are already captured by the contact/H-bond
   metrics in interface_geom.py / hbonds_pyrosetta.py.)

Detection (per aromatic-side-chain <-> base pair)
  Each ring -> centroid c and unit plane normal n (SVD of the ring-atom coords).
  With v = c_nt - c_aa, d = |v|, inter-plane angle a = angle(n_aa, n_nt) in [0,90],
  and lateral offset = component of v perpendicular to n_aa:
    face-to-face (parallel) : d <= DIST_MAX AND a <= ANGLE_PAR AND offset <= OFFSET_MAX
    edge-to-face (T-shaped) : d <= DIST_MAX AND a >= ANGLE_PERP
  Face-to-face is the "real" stack; T-shaped is reported separately.

Strength (graded, so "strong" != merely "present")
  s = f_dist(d) * f_angle(a) * w_aa * w_nt   in (0, 1]
    f_dist  = exp(-((d - D0)/DW)^2)   Gaussian peaked at the ideal spacing D0 = 3.7 A
    f_angle = cos(a)^2                1 at parallel, 0 at perpendicular
    w_aa, w_nt = identity weights (ring size / polarizability; see W_AA / W_NT).
  All thresholds and weights are module constants -- tune them in one place.

Per-design columns produced (PISTACK_COLUMNS)
  pi_stack_score        sum of face-to-face pair strengths   (PRIMARY, higher=better)
  pi_stack_count        number of face-to-face stacks (strength >= PAIR_MIN)
  pi_stack_best         strongest single stack strength (is there one dominant stack)
  pi_stack_res_aa       distinct protein aromatic residues that stack
  pi_stack_res_nt       distinct RNA bases that stack
  pi_stack_count_tshape edge-to-face (T-shaped) stacks (secondary)
  pi_stack_<AA>         face-to-face stack count per residue type (TRP/TYR/PHE/HIS/ARG)
"""
from __future__ import annotations

from pathlib import Path

import numpy as np

from rna_binder.analysis.structural import load_structure, split_chains

# ---- ring-atom sets (standard PDB / AF3 / Boltz atom names) -----------------
PURINE_RING = {"N1", "C2", "N3", "C4", "C5", "C6", "N7", "C8", "N9"}       # A, G
PYRIMIDINE_RING = {"N1", "C2", "N3", "C4", "C5", "C6"}                      # C, U
BASE_RING = {"A": PURINE_RING, "G": PURINE_RING,
             "C": PYRIMIDINE_RING, "U": PYRIMIDINE_RING}
AA_RING = {
    "PHE": {"CG", "CD1", "CD2", "CE1", "CE2", "CZ"},                        # benzene
    "TYR": {"CG", "CD1", "CD2", "CE1", "CE2", "CZ"},                        # phenol
    "TRP": {"CG", "CD1", "NE1", "CE2", "CD2", "CE3", "CZ2", "CZ3", "CH2"},  # indole
    "HIS": {"CG", "ND1", "CD2", "CE1", "NE2"},                             # imidazole
    "ARG": {"NE", "CZ", "NH1", "NH2"},                                     # guanidinium
}
AA_TYPES = ["TRP", "TYR", "PHE", "HIS", "ARG"]

# ---- geometry thresholds ---------------------------------------------------
DIST_MAX = 5.5      # A, centroid-centroid cutoff for any stack (offset stacks reach ~5)
D0 = 3.7            # A, ideal stacking separation (Gaussian peak)
DW = 0.9            # A, Gaussian width
ANGLE_PAR = 30.0    # deg, inter-plane angle <= this => face-to-face (parallel)
ANGLE_PERP = 60.0   # deg, inter-plane angle >= this => edge-to-face (T-shaped)
OFFSET_MAX = 3.0    # A, max lateral (in-plane) offset for a face-to-face stack
PAIR_MIN = 0.10     # min strength for a pair to count in pi_stack_count

# ---- identity weights (heuristic: ring size / polarizability) --------------
W_AA = {"TRP": 1.00, "TYR": 0.85, "PHE": 0.80, "HIS": 0.60, "ARG": 0.70}
W_NT = {"G": 1.00, "A": 0.90, "C": 0.70, "U": 0.65}

PISTACK_COLUMNS = (
    ["pi_stack_score", "pi_stack_count", "pi_stack_best",
     "pi_stack_res_aa", "pi_stack_res_nt", "pi_stack_count_tshape"]
    + [f"pi_stack_{aa}" for aa in AA_TYPES]
)


def _ring_plane(coords: np.ndarray):
    """Centroid + unit plane normal (SVD) for a set of >=3 ring-atom coords."""
    c = coords.mean(axis=0)
    # smallest right-singular vector of the centered coords = best-fit plane normal
    _, _, vt = np.linalg.svd(coords - c, full_matrices=False)
    return c, vt[-1]


def _rings(atoms, ring_map):
    """List of (res_key, res_name, centroid, normal) for every residue in `atoms`
    whose ring atoms (per ring_map, keyed by res_name) are present (>=3 atoms)."""
    chain, resid = atoms.chain_id, atoms.res_id
    resname, aname, coord = atoms.res_name, atoms.atom_name, atoms.coord
    groups: dict = {}
    for i in range(atoms.array_length()):
        groups.setdefault((chain[i], int(resid[i])), []).append(i)
    rings = []
    for key, idxs in groups.items():
        rn = str(resname[idxs[0]])
        want = ring_map.get(rn)
        if not want:
            continue
        sel = [i for i in idxs if str(aname[i]) in want]
        if len(sel) < 3:
            continue
        c, n = _ring_plane(coord[sel])
        rings.append((key, rn, c, n))
    return rings


def _plane_angle(n1: np.ndarray, n2: np.ndarray) -> float:
    """Angle between two plane normals, folded to [0, 90] (sign of normal irrelevant)."""
    d = abs(float(np.dot(n1, n2)))
    return float(np.degrees(np.arccos(min(1.0, max(0.0, d)))))


def pi_stack(atoms) -> dict:
    """Aromatic pi-stacking descriptors for one protein<->RNA complex. Never raises
    upward -- returns {'pi_stack_error': ...} on failure."""
    protein, rna = split_chains(atoms)
    if protein is None or rna is None:
        return {"pi_stack_error": "could not split protein/RNA chains"}

    aa_rings = _rings(protein, AA_RING)
    nt_rings = _rings(rna, BASE_RING)

    score = best = 0.0
    count = tshape = 0
    res_aa, res_nt = set(), set()
    per_aa = {aa: 0 for aa in AA_TYPES}

    for ak, arn, ac, an in aa_rings:
        for nk, nrn, nc, nn in nt_rings:
            v = nc - ac
            d = float(np.linalg.norm(v))
            if d > DIST_MAX:
                continue
            a = _plane_angle(an, nn)
            axial = abs(float(np.dot(v, an)))               # separation along aa normal
            lateral = float(np.sqrt(max(0.0, d * d - axial * axial)))
            if a <= ANGLE_PAR and lateral <= OFFSET_MAX:    # face-to-face
                f_dist = float(np.exp(-((d - D0) / DW) ** 2))
                f_ang = float(np.cos(np.radians(a)) ** 2)
                s = f_dist * f_ang * W_AA.get(arn, 0.7) * W_NT.get(nrn, 0.7)
                if s >= PAIR_MIN:
                    score += s
                    best = max(best, s)
                    count += 1
                    res_aa.add(ak)
                    res_nt.add(nk)
                    per_aa[arn] = per_aa.get(arn, 0) + 1
            elif a >= ANGLE_PERP:                            # edge-to-face (T-shaped)
                tshape += 1

    out = dict(
        pi_stack_score=round(score, 4),
        pi_stack_count=int(count),
        pi_stack_best=round(best, 4),
        pi_stack_res_aa=int(len(res_aa)),
        pi_stack_res_nt=int(len(res_nt)),
        pi_stack_count_tshape=int(tshape),
    )
    for aa in AA_TYPES:
        out[f"pi_stack_{aa}"] = int(per_aa[aa])
    return out


def pi_stack_from_path(path, **_) -> dict:
    try:
        atoms = load_structure(Path(path))
    except Exception as e:
        return {"pi_stack_error": f"load: {type(e).__name__}: {e}"}
    return pi_stack(atoms)


if __name__ == "__main__":
    import sys
    for p in sys.argv[1:]:
        print(p)
        for k, v in pi_stack_from_path(p).items():
            print(f"  {k:22s} {v}")
