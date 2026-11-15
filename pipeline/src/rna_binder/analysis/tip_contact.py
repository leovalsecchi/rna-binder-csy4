"""Tip-weighted protein–RNA contact score (biotite port).

A geometric scorer for de novo RNA-binder designs, complementary to the
per-nucleotide `rna_contact_counts` in ``structural.py``. It scores favorable
protein–RNA CONTACTS by how close they sit to the loop *tip* of the target
hairpin versus the open-end/stem-bottom (the scaffold-adjacent end), and
reports the fraction of contacts that rely on the stem bottom — unreliable
once a scaffold is attached there.

Physical model
--------------
Coaxial stem axis via PCA on the RNA C1' atoms, oriented from the open end
("bottom", scaffold junction) toward the loop apex ("tip"). Each RNA atom gets
a normalized axial coordinate ``t``::

    t = 0   at the bottom anchor (open end / scaffold junction)
    t = 1   at the loop anchor   (apex)

Favorable protein–RNA contacts are detected in four categories, mirroring the
original ``tip_weighted_contact_score.py``:

    hbond    : polar heavy atom (N/O) <-> polar RNA atom (N/O), close range
    basic_p  : Lys/Arg/His sidechain N <-> RNA phosphate O (salt bridge)
    base     : aromatic/Arg protein atom <-> RNA base-ring atom
    vdw      : hydrophobic protein C <-> RNA C (packing)

Each favorable contact is weighted by ``weight(t) = t ** gamma``, so contacts
near the tip (t~1) count fully and contacts near the bottom (t~0) count close
to zero. ``gamma > 1`` sharpens the tip preference; ``gamma < 1`` softens it.

Bottom/loop residues are **auto-detected** from the RNA chain (no manual seqid
picks): the design outputs number the RNA differently per stack (RFD3 1-16,
HD/BD 1-18, BoltzGen 2-17) so absolute seqids are not portable. For a hairpin
the 5' and 3' termini both sit at the open (scaffold) end, so:

    bottom = the first END_N and last END_N nucleotides in 5'->3' order
    loop   = the LOOP_N nucleotides whose C1' project farthest along the axis
             from the bottom anchor (the apex)

This is validated against the fixed csy4 target geometry (termini cluster at
one axial end, the sequence-middle nucleotides at the far/apex end).

Run on the DESIGN-NATIVE structure (generator output), where the RNA is the
fixed docking target, so contact geometry reflects the as-designed interface.
Note some stacks emit backbone-only design structures (e.g. RFD3+LiMPNN threads
a sequence but no sidechain coords); ``tip_has_sidechains`` flags that, since the
sidechain-dependent categories (basic_p/base/vdw) collapse without sidechains.

Public entry point: ``compute_tip_weighted_contact(path) -> dict`` returns a
flat, column-friendly dict (``tip_*`` keys) for one design, following the
``structural.compute_all_structural`` convention (never raises; failures become
a ``tip_contact_error`` string).
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import Iterable, Optional

import numpy as np

from .structural import load_structure, split_chains

# --- tunables (module-level, mirroring structural.py's constants style) ------
GAMMA = 1.0            # tip-preference sharpness: weight = t ** GAMMA
BOTTOM_ZONE = 0.3      # t below this = "bottom zone" for the bottom_fraction diagnostic
BOTTOM_FRAC_CUTOFF = 0.5  # tip_pass = bottom_fraction <= this (hard-filter knob)
END_N = 2              # nucleotides at each terminus taken as the open-end "bottom"
LOOP_N = 3             # apex nucleotides (max axial projection) taken as the "loop"

HBOND_DIST = 3.5       # Å, polar–polar
PHOS_DIST = 4.0        # Å, basic sidechain N – phosphate O
BASE_DIST = 4.2        # Å, aromatic/Arg – base ring
VDW_DIST = 4.0         # Å, hydrophobic C – RNA C
_MAX_CUT = max(HBOND_DIST, PHOS_DIST, BASE_DIST, VDW_DIST)

BASIC_RES = {"LYS", "ARG", "HIS"}
AROMATIC_RES = {"PHE", "TYR", "TRP", "HIS", "ARG"}
HYDROPHOBIC_RES = {"ALA", "VAL", "LEU", "ILE", "MET", "PHE", "TRP", "TYR", "PRO", "CYS"}
PROT_BACKBONE = {"N", "CA", "C", "O", "OXT"}

RNA_PHOS_ATOMS = {"OP1", "OP2", "O1P", "O2P"}
RNA_BACKBONE_ATOMS = {
    "P", "OP1", "OP2", "O1P", "O2P", "O5'", "C5'", "C4'", "O4'",
    "C3'", "O3'", "C2'", "O2'", "C1'",
}

_CONTACT_TYPES = ("hbond", "basic_p", "base", "vdw")


# ---------------------------------------------------------------------------
# atom records (from a biotite AtomArray chain), mirroring the gemmi original
# ---------------------------------------------------------------------------

def _records(chain) -> list:
    """Heavy-atom records for a biotite chain: dicts of seqid/resname/atomname/elem/xyz."""
    if chain is None:
        return []
    heavy = chain[chain.element != "H"]
    out = []
    for i in range(heavy.array_length()):
        out.append({
            "seqid": int(heavy.res_id[i]),
            "resname": str(heavy.res_name[i]).strip(),
            "atomname": str(heavy.atom_name[i]).replace("*", "'"),
            "elem": str(heavy.element[i]).strip() or str(heavy.atom_name[i])[0],
            "xyz": heavy.coord[i].astype(float),
        })
    return out


def _c1prime(rna_records, residues: Optional[Iterable[int]] = None):
    resset = set(residues) if residues is not None else None
    pts = [a["xyz"] for a in rna_records
           if a["atomname"] == "C1'" and (resset is None or a["seqid"] in resset)]
    return np.asarray(pts, float)


def _c1prime_by_res(rna_records):
    """{seqid: C1' xyz} in ascending-seqid (5'->3') order."""
    d = {}
    for a in rna_records:
        if a["atomname"] == "C1'" and a["seqid"] not in d:
            d[a["seqid"]] = a["xyz"]
    return dict(sorted(d.items()))


# ---------------------------------------------------------------------------
# stem axis with auto-detected bottom (termini) and loop (apex)
# ---------------------------------------------------------------------------

def stem_axis_auto(rna_records, bottom_residues=None, loop_residues=None,
                   end_n=END_N, loop_n=LOOP_N):
    """PCA stem axis oriented bottom -> loop, with anchors auto-detected from the
    RNA hairpin unless explicit ``bottom_residues``/``loop_residues`` are given.

    Returns (axis_unit, bottom_anchor, loop_anchor, bottom_ids, loop_ids)."""
    by_res = _c1prime_by_res(rna_records)
    if len(by_res) < 4:
        raise ValueError("Too few C1' atoms to estimate a stem axis.")
    ids = list(by_res.keys())
    allc1 = np.asarray([by_res[i] for i in ids], float)
    centroid = allc1.mean(0)
    _, _, vt = np.linalg.svd(allc1 - centroid)
    axis = vt[0] / np.linalg.norm(vt[0])

    # bottom: explicit picks, else the first/last end_n nucleotides (both termini
    # of a hairpin sit at the open, scaffold-adjacent end).
    if bottom_residues:
        bottom_ids = [i for i in ids if i in set(bottom_residues)]
    else:
        k = min(end_n, len(ids) // 2)
        bottom_ids = ids[:k] + ids[-k:]
    if not bottom_ids:
        raise ValueError("No C1' atoms found for the bottom residues.")
    bottom_anchor = np.asarray([by_res[i] for i in bottom_ids], float).mean(0)

    # orient axis so projections increase going from the bottom toward the loop
    if np.dot(centroid - bottom_anchor, axis) < 0:
        axis = -axis

    # loop: explicit picks, else the loop_n nucleotides projecting farthest from
    # the bottom anchor along the axis (the apex).
    if loop_residues:
        loop_ids = [i for i in ids if i in set(loop_residues)]
        if not loop_ids:
            raise ValueError("No C1' atoms found for the loop residues.")
    else:
        proj = {i: float(np.dot(by_res[i] - bottom_anchor, axis)) for i in ids}
        loop_ids = sorted(proj, key=proj.get)[-min(loop_n, len(ids)):]
    loop_anchor = np.asarray([by_res[i] for i in loop_ids], float).mean(0)
    return axis, bottom_anchor, loop_anchor, bottom_ids, sorted(loop_ids)


def axial_coord(xyz, axis, bottom_anchor, stem_length):
    """Normalized axial position t in [0, 1] (clipped): 0=bottom, 1=loop tip."""
    t = float(np.dot(xyz - bottom_anchor, axis)) / stem_length
    return float(np.clip(t, 0.0, 1.0))


# ---------------------------------------------------------------------------
# favorable contact detection (faithful port of the original four categories)
# ---------------------------------------------------------------------------

def find_favorable_contacts(prot_atoms, rna_atoms,
                            hbond_dist=HBOND_DIST, phos_dist=PHOS_DIST,
                            base_dist=BASE_DIST, vdw_dist=VDW_DIST):
    """List of contact dicts across the four categories (hbond/basic_p/base/vdw).

    A given atom pair may satisfy more than one category; each satisfied category
    is recorded separately (they represent different physical interactions)."""
    contacts = []
    if not prot_atoms or not rna_atoms:
        return contacts
    max_cut = max(hbond_dist, phos_dist, base_dist, vdw_dist)

    rna_xyz = np.array([a["xyz"] for a in rna_atoms])
    prot_xyz = np.array([a["xyz"] for a in prot_atoms])
    try:
        from scipy.spatial import cKDTree
        pairs = cKDTree(rna_xyz).query_ball_point(prot_xyz, r=max_cut)

        def rna_near(i):
            return pairs[i]
    except ImportError:  # pragma: no cover - scipy is present in the runtime env
        def rna_near(i):
            return range(len(rna_atoms))

    for i, pa in enumerate(prot_atoms):
        p_res, p_an, p_elem = pa["resname"], pa["atomname"], pa["elem"]
        p_backbone = p_an in PROT_BACKBONE
        p_is_polar = p_elem in ("N", "O")
        p_basic_sc_n = (p_res in BASIC_RES) and p_elem == "N" and not p_backbone
        p_arom_sc = (p_res in AROMATIC_RES) and not p_backbone
        p_hphob_c = (p_res in HYDROPHOBIC_RES) and p_elem == "C" and not p_backbone
        pxyz = pa["xyz"]

        for j in rna_near(i):
            ra = rna_atoms[j]
            d = float(np.linalg.norm(pxyz - ra["xyz"]))
            if d > max_cut:
                continue
            r_elem, r_an = ra["elem"], ra["atomname"]
            r_is_polar = r_elem in ("N", "O")
            r_is_phos = r_an in RNA_PHOS_ATOMS
            r_is_base = r_an not in RNA_BACKBONE_ATOMS

            if p_is_polar and r_is_polar and d <= hbond_dist:
                contacts.append(_c("hbond", d, pa, ra))
            if p_basic_sc_n and r_is_phos and d <= phos_dist:
                contacts.append(_c("basic_p", d, pa, ra))
            if p_arom_sc and r_is_base and p_elem in ("C", "N") and d <= base_dist:
                contacts.append(_c("base", d, pa, ra))
            if p_hphob_c and r_elem == "C" and d <= vdw_dist:
                contacts.append(_c("vdw", d, pa, ra))
    return contacts


def _c(ctype, d, pa, ra):
    return dict(type=ctype, dist=round(d, 2),
                prot_seqid=pa["seqid"], prot_resname=pa["resname"], prot_atom=pa["atomname"],
                rna_seqid=ra["seqid"], rna_resname=ra["resname"], rna_atom=ra["atomname"])


# ---------------------------------------------------------------------------
# scoring
# ---------------------------------------------------------------------------

def score_atoms(prot_atoms, rna_atoms, gamma=GAMMA, bottom_zone=BOTTOM_ZONE,
                bottom_frac_cutoff=BOTTOM_FRAC_CUTOFF,
                bottom_residues=None, loop_residues=None,
                keep_contacts=False):
    """Core scorer over pre-extracted protein/RNA atom records. Returns a summary
    dict; raises ValueError on degenerate geometry (caught by the wrapper)."""
    axis, bottom_anchor, loop_anchor, bottom_ids, loop_ids = stem_axis_auto(
        rna_atoms, bottom_residues, loop_residues)
    stem_length = float(np.dot(loop_anchor - bottom_anchor, axis))
    if stem_length <= 0:
        raise ValueError("Non-positive stem length; check bottom/loop picks.")

    has_sidechains = any(a["atomname"] not in PROT_BACKBONE for a in prot_atoms)

    contacts = find_favorable_contacts(prot_atoms, rna_atoms)
    rna_by_key = {(a["seqid"], a["atomname"]): a for a in rna_atoms}
    for c in contacts:
        ra = rna_by_key[(c["rna_seqid"], c["rna_atom"])]
        t = axial_coord(ra["xyz"], axis, bottom_anchor, stem_length)
        c["t"] = round(t, 3)
        c["weight"] = round(t ** gamma, 4)

    n_fav = len(contacts)
    tip_weighted_score = round(sum(c["weight"] for c in contacts), 4)
    mean_weight = round(tip_weighted_score / n_fav, 4) if n_fav else 0.0
    n_bottom = sum(1 for c in contacts if c["t"] < bottom_zone)
    bottom_fraction = round(n_bottom / n_fav, 4) if n_fav else None

    by_type = {t: 0.0 for t in _CONTACT_TYPES}
    by_type_n = {t: 0 for t in _CONTACT_TYPES}
    for c in contacts:
        by_type[c["type"]] += c["weight"]
        by_type_n[c["type"]] += 1

    out = {
        "tip_stem_length_A": round(stem_length, 2),
        "tip_n_favorable_contacts": n_fav,
        "tip_weighted_score": tip_weighted_score,
        "tip_mean_contact_weight": mean_weight,
        "tip_bottom_fraction": bottom_fraction,
        "tip_n_bottom_zone_contacts": n_bottom,
        "tip_has_sidechains": bool(has_sidechains),
        "tip_pass": (bottom_fraction is None) or (bottom_fraction <= bottom_frac_cutoff),
        "tip_bottom_ids": bottom_ids,
        "tip_loop_ids": loop_ids,
    }
    for t in _CONTACT_TYPES:
        out[f"tip_{t}_wsum"] = round(by_type[t], 4)
        out[f"tip_{t}_n"] = by_type_n[t]
    if keep_contacts:
        out["contacts"] = contacts
    return out


def compute_tip_weighted_contact(designed_path: Optional[Path],
                                 gamma=GAMMA, bottom_zone=BOTTOM_ZONE,
                                 bottom_frac_cutoff=BOTTOM_FRAC_CUTOFF,
                                 bottom_residues=None, loop_residues=None) -> dict:
    """Tip-weighted contact metrics for one DESIGN-NATIVE structure.

    Never raises: a missing path or any failure yields a ``tip_contact_error``
    string so the row is kept with NaN tip metrics (matching structural.py)."""
    if designed_path is None:
        return {"tip_contact_error": "no designed_cif_path"}
    try:
        atoms = load_structure(Path(designed_path))
        prot, rna = split_chains(atoms)
        if prot is None:
            return {"tip_contact_error": "no protein chain"}
        if rna is None:
            return {"tip_contact_error": "no RNA chain"}
        return score_atoms(_records(prot), _records(rna),
                           gamma=gamma, bottom_zone=bottom_zone,
                           bottom_frac_cutoff=bottom_frac_cutoff,
                           bottom_residues=bottom_residues, loop_residues=loop_residues)
    except Exception as e:  # noqa: BLE001 - keep the row, report the error
        return {"tip_contact_error": f"{type(e).__name__}: {e}"}


# columns this module contributes (for callers that want the schema up front)
TIP_COLUMNS = (
    "tip_stem_length_A", "tip_n_favorable_contacts", "tip_weighted_score",
    "tip_mean_contact_weight", "tip_bottom_fraction", "tip_n_bottom_zone_contacts",
    "tip_has_sidechains", "tip_pass",
) + tuple(f"tip_{t}_wsum" for t in _CONTACT_TYPES) \
  + tuple(f"tip_{t}_n" for t in _CONTACT_TYPES)


# ---------------------------------------------------------------------------
# CLI — single design, JSON out (parity with the original standalone script)
# ---------------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("structure", help=".pdb / .cif / .cif.gz design (RNA + protein)")
    ap.add_argument("--bottom", type=int, nargs="+", default=None,
                    help="RNA seqids at the open end (default: auto-detect termini)")
    ap.add_argument("--loop", type=int, nargs="+", default=None,
                    help="RNA seqids in the loop (default: auto-detect apex)")
    ap.add_argument("--gamma", type=float, default=GAMMA)
    ap.add_argument("--bottom-zone", type=float, default=BOTTOM_ZONE)
    ap.add_argument("--bottom-frac-cutoff", type=float, default=BOTTOM_FRAC_CUTOFF)
    ap.add_argument("--contacts", action="store_true",
                    help="include the per-contact list in the JSON output")
    args = ap.parse_args()

    atoms = load_structure(Path(args.structure))
    prot, rna = split_chains(atoms)
    if prot is None or rna is None:
        sys.exit("Could not find both a protein and an RNA chain.")
    result = score_atoms(
        _records(prot), _records(rna), gamma=args.gamma,
        bottom_zone=args.bottom_zone, bottom_frac_cutoff=args.bottom_frac_cutoff,
        bottom_residues=args.bottom, loop_residues=args.loop,
        keep_contacts=args.contacts)
    result["file"] = args.structure
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
