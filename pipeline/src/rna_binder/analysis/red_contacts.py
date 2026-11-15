"""Unfavorable ("red") protein–RNA interface contacts (biotite port).

A batch/quantitative port of the RED annotations in the ChimeraX preset
``AutoRNAproteincontactscore.py`` (the "RNAcon" view). Where that preset draws
the unfavorable interactions for one open structure, this computes them as
per-structure COUNTS so they can be tabulated across many designs and compared
to natural RBP–RNA complexes. Companion to ``tip_contact.py`` (favorable
contacts) and the per-nucleotide counts in ``structural.py``.

The three interface red categories (identical geometry/cutoffs to the preset):

    clash        : protein heavy atom vs RNA heavy atom whose van-der-Waals
                   spheres overlap by >= 0.6 A (a 0.4 A allowance is subtracted
                   for polar–polar pairs that could H-bond) — steric overlap.
                   [preset: ``clashes (bndiface) restrict (rnaiface)``]
    acid_phos    : Asp/Glu sidechain carboxylate O within 4.0 A (centre-to-
                   centre) of an RNA phosphate O (OP1/OP2/O1P/O2P) — like-charge
                   electrostatic repulsion (negative on negative).
                   [preset: ``contacts (:asp,glu & sc & @O*) restrict (rnaphos) 4.0``]
    hphob_polar  : nonpolar sidechain C (A,V,L,I,M,F,W,P,C — NOT Tyr) within
                   3.6 A of a polar RNA backbone O (phosphate + ribose
                   O2'/O3'/O4'/O5') — greasy atom buried on a polar/charged
                   backbone (desolvation / poor complementarity).
                   [preset: ``contacts (hydrophobic sc @C*) restrict (rnapolarbackbone) 3.6``]

These are geometric HEURISTIC flags, not energies (as the preset itself notes).
Counts are cutoff-sensitive; use them for relative distributions / ranking, and
normalise by interface size (the ``*_per_ifaceres`` fields) before comparing
structures with different interface sizes. Comparing predicted design models to
experimental natural structures carries a structure-source confound — read the
counts as a QC/benchmark signal, not ground truth.

Public entry point: ``compute_red_contacts(path, prot_chains=None,
rna_chains=None) -> dict`` returns a flat ``red_*`` dict (never raises; failures
become a ``red_contacts_error`` string, matching ``structural.compute_all_structural``).
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import Optional, Sequence

import numpy as np

from .structural import load_structure, split_chains

# --- cutoffs (match AutoRNAproteincontactscore.py) ---------------------------
ACID_PHOS_DIST = 4.0        # A, Asp/Glu sidechain O  <-> RNA phosphate O
HPHOB_POLAR_DIST = 3.6      # A, hydrophobic sidechain C <-> RNA polar backbone O
INTERFACE_CUTOFF = 5.0      # A, protein residue counted as interface (matches structural.py)

# ChimeraX `clashes` defaults
CLASH_OVERLAP = 0.6         # A, VDW overlap to call a clash
CLASH_HBOND_ALLOWANCE = 0.4 # A, subtracted from overlap for polar–polar (H-bond-capable) pairs

# --- residue / atom sets (match the preset selections) -----------------------
ACID_RES = {"ASP", "GLU"}
# preset red hydrophobic set EXCLUDES Tyr (its OH is polar); MSE≈Met added for crystal structures
RED_HYDROPHOBIC_RES = {"ALA", "VAL", "LEU", "ILE", "MET", "PHE", "TRP", "PRO", "CYS", "MSE"}
PROT_BACKBONE = {"N", "CA", "C", "O", "OXT"}

RNA_PHOS_ATOMS = {"OP1", "OP2", "O1P", "O2P"}
RNA_POLAR_BACKBONE_O = {"OP1", "OP2", "O1P", "O2P", "O2'", "O3'", "O4'", "O5'"}

# Bondi van-der-Waals radii (A); ChimeraX's defaults are effectively these.
VDW_RADII = {"C": 1.70, "N": 1.55, "O": 1.52, "P": 1.80, "S": 1.80, "SE": 1.90,
             "F": 1.47, "CL": 1.75, "BR": 1.85, "MG": 1.73, "NA": 2.27,
             "ZN": 1.39, "K": 2.75, "H": 1.20}
_VDW_DEFAULT = 1.70
_MAX_VDW = 1.90            # largest radius we care about for the clash query window


def _vdw(elem: str) -> float:
    return VDW_RADII.get(elem.upper(), _VDW_DEFAULT)


def _records(chain) -> list:
    """Heavy-atom records for a biotite chain/array: chain_id/seqid/resname/atomname/elem/xyz."""
    if chain is None or chain.array_length() == 0:
        return []
    heavy = chain[chain.element != "H"]
    out = []
    for i in range(heavy.array_length()):
        out.append({
            "chain": str(heavy.chain_id[i]),
            "seqid": int(heavy.res_id[i]),
            "resname": str(heavy.res_name[i]).strip().upper(),
            "atomname": str(heavy.atom_name[i]).replace("*", "'"),
            "elem": (str(heavy.element[i]).strip() or str(heavy.atom_name[i])[0]).upper(),
            "xyz": heavy.coord[i].astype(float),
        })
    return out


def score_red(prot_records: list, rna_records: list) -> dict:
    """Core red-contact scorer over pre-extracted protein/RNA heavy-atom records.

    Raises ValueError on empty input (caught by the wrapper)."""
    if not prot_records:
        raise ValueError("no protein heavy atoms")
    if not rna_records:
        raise ValueError("no RNA heavy atoms")

    rna_xyz = np.array([a["xyz"] for a in rna_records])
    prot_xyz = np.array([a["xyz"] for a in prot_records])

    from scipy.spatial import cKDTree
    tree = cKDTree(rna_xyz)
    # one neighbour query at the largest window (interface cutoff) covers every category
    query_r = max(INTERFACE_CUTOFF, ACID_PHOS_DIST, HPHOB_POLAR_DIST, 2 * _MAX_VDW)
    neigh = tree.query_ball_point(prot_xyz, r=query_r)

    n_clash = n_acid = n_hphob = 0
    clash_res, acid_res, hphob_res = set(), set(), set()
    iface_prot_res, iface_rna_res = set(), set()
    n_iface_heavy = 0

    for i, pa in enumerate(prot_records):
        p_bb = pa["atomname"] in PROT_BACKBONE
        p_elem = pa["elem"]
        p_acid_o = (pa["resname"] in ACID_RES) and (p_elem == "O") and (not p_bb)
        p_hphob_c = (pa["resname"] in RED_HYDROPHOBIC_RES) and (p_elem == "C") and (not p_bb)
        p_polar = p_elem in ("N", "O")
        p_vdw = _vdw(p_elem)
        p_res = (pa["chain"], pa["seqid"])
        pxyz = pa["xyz"]

        for j in neigh[i]:
            ra = rna_records[j]
            d = float(np.linalg.norm(pxyz - ra["xyz"]))
            if d > INTERFACE_CUTOFF:
                continue
            # interface bookkeeping (any protein heavy atom within 5 A of RNA)
            iface_prot_res.add(p_res)
            iface_rna_res.add((ra["chain"], ra["seqid"]))
            n_iface_heavy += 1

            r_an, r_elem = ra["atomname"], ra["elem"]

            # 1) steric clash (VDW overlap)
            overlap = p_vdw + _vdw(r_elem) - d
            allow = CLASH_HBOND_ALLOWANCE if (p_polar and r_elem in ("N", "O")) else 0.0
            if overlap - allow >= CLASH_OVERLAP:
                n_clash += 1
                clash_res.add(p_res)

            # 2) acidic sidechain O vs phosphate O (like-charge repulsion)
            if p_acid_o and r_an in RNA_PHOS_ATOMS and d <= ACID_PHOS_DIST:
                n_acid += 1
                acid_res.add(p_res)

            # 3) hydrophobic sidechain C vs polar RNA backbone O
            if p_hphob_c and r_an in RNA_POLAR_BACKBONE_O and d <= HPHOB_POLAR_DIST:
                n_hphob += 1
                hphob_res.add(p_res)

    n_if = len(iface_prot_res)
    has_sc = any(a["atomname"] not in PROT_BACKBONE for a in prot_records)

    def per_if(n):
        return round(n / n_if, 4) if n_if else None

    red_total = n_clash + n_acid + n_hphob
    return {
        "red_n_clashes": n_clash,
        "red_n_acid_phos_pairs": n_acid,
        "red_n_hphob_polar_pairs": n_hphob,
        "red_n_total_pairs": red_total,
        "red_n_clash_res": len(clash_res),
        "red_n_acid_phos_res": len(acid_res),
        "red_n_hphob_polar_res": len(hphob_res),
        "red_n_interface_res_prot": n_if,
        "red_n_interface_res_rna": len(iface_rna_res),
        "red_n_interface_heavy_contacts": n_iface_heavy,
        "red_clashes_per_ifaceres": per_if(n_clash),
        "red_acid_phos_per_ifaceres": per_if(n_acid),
        "red_hphob_polar_per_ifaceres": per_if(n_hphob),
        "red_total_per_ifaceres": per_if(red_total),
        "red_has_sidechains": bool(has_sc),
    }


def _select_chains(atoms, prot_chains: Optional[Sequence[str]],
                   rna_chains: Optional[Sequence[str]]):
    """Return (prot_records, rna_records). If explicit chain lists are given
    (natural multi-chain complexes), subset by chain_id; else fall back to
    ``split_chains`` (the single-binder + single-RNA design case)."""
    if prot_chains or rna_chains:
        pset, rset = set(prot_chains or ()), set(rna_chains or ())
        prot = atoms[np.isin(atoms.chain_id, list(pset))] if pset else None
        rna = atoms[np.isin(atoms.chain_id, list(rset))] if rset else None
        return _records(prot), _records(rna)
    prot, rna = split_chains(atoms)
    return _records(prot), _records(rna)


def compute_red_contacts(path: Optional[Path],
                         prot_chains: Optional[Sequence[str]] = None,
                         rna_chains: Optional[Sequence[str]] = None) -> dict:
    """Red-contact metrics for one structure. Never raises: any failure yields a
    ``red_contacts_error`` string so the row is kept with NaN red metrics."""
    if path is None:
        return {"red_contacts_error": "no path"}
    try:
        atoms = load_structure(Path(path))
        prot_rec, rna_rec = _select_chains(atoms, prot_chains, rna_chains)
        if not prot_rec:
            return {"red_contacts_error": "no protein chain"}
        if not rna_rec:
            return {"red_contacts_error": "no RNA chain"}
        return score_red(prot_rec, rna_rec)
    except Exception as e:  # noqa: BLE001 - keep the row, report the error
        return {"red_contacts_error": f"{type(e).__name__}: {e}"}


RED_COLUMNS = (
    "red_n_clashes", "red_n_acid_phos_pairs", "red_n_hphob_polar_pairs", "red_n_total_pairs",
    "red_n_clash_res", "red_n_acid_phos_res", "red_n_hphob_polar_res",
    "red_n_interface_res_prot", "red_n_interface_res_rna", "red_n_interface_heavy_contacts",
    "red_clashes_per_ifaceres", "red_acid_phos_per_ifaceres", "red_hphob_polar_per_ifaceres",
    "red_total_per_ifaceres", "red_has_sidechains",
)


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("structure", help=".pdb / .cif / .cif.gz (RNA + protein)")
    ap.add_argument("--prot-chains", nargs="+", default=None,
                    help="explicit protein chain ids (natural complexes); else auto-detect")
    ap.add_argument("--rna-chains", nargs="+", default=None,
                    help="explicit RNA chain ids (natural complexes); else auto-detect")
    args = ap.parse_args()
    out = compute_red_contacts(Path(args.structure), args.prot_chains, args.rna_chains)
    out["file"] = args.structure
    print(json.dumps(out, indent=2))


if __name__ == "__main__":
    main()
