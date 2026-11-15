#!/usr/bin/env python3
"""interface_geom.py — exact geometric protein<->RNA interface descriptors (Biotite).

Part of the analyses/ metric pipeline. Owns the two metrics that are exact by
construction (no numerical/definition approximation):

  * molecular contacts   -- heavy-atom protein<->RNA atom-pair counts
  * pairwise AA<->{backbone,base} table -- per-amino-acid contact counts

The two dSASA metrics (total + hydrophobic/NPSA) are NOT computed here: they are
produced field-standard by Rosetta InterfaceAnalyzer (interface_rosetta.py), and
H-bonds by hbonds_pyrosetta.py. A Biotite dSASA is kept available as an optional
cross-check (interface_geom(..., include_sasa=True)) but is OFF by default so the
canonical Rosetta numbers are authoritative.

Chain splitting reuses rna_binder.analysis.structural.split_chains (RNA = the
chain whose residues subset {A,C,G,U}; the other = protein), so chain labels need
not be A/B.

Per-design columns produced (CONTACT_COLUMNS):

  Molecular contacts  (protein heavy atom <-> RNA heavy atom)
    n_heavy_contacts_4A      atom-atom pairs within 4.0 A (primary)
    n_heavy_contacts_4p5A    atom-atom pairs within 4.5 A (matches rna_contact_nt)
    n_contact_res_prot       distinct protein residues in contact
    n_contact_res_rna        distinct RNA residues in contact

  Base vs backbone partition of the 4.0 A contacts
    n_contacts_base / n_contacts_backbone / n_contacts_phosphate / n_contacts_sugar
    frac_base_contacts       base / (base + backbone)  -> reads bases vs grips cage

  Pairwise AA <-> {backbone, base} table (metric 3), per design
    aa_base_<AA>  and  aa_bb_<AA>  for each of the 20 amino acids (40 columns).
    e.g. aa_bb_ARG = Arg heavy-atom contacts to the RNA backbone, aa_base_ARG to
    bases. The per-AA columns sum exactly to n_contacts_base / n_contacts_backbone.
    aggregate_aa_matrix() rolls these up across designs into the campaign 20 x
    {backbone,base} frequency table for the propensity heatmap artifact.

Optional Biotite cross-check dSASA (include_sasa=True): dsasa_total_biotite,
dsasa_apolar_biotite, dsasa_polar_biotite, dsasa_apolar_frac_biotite.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import biotite.structure as struc

from rna_binder.analysis.structural import load_structure, split_chains

# ---- geometry / classification constants -----------------------------------
CONTACT_CUTOFF = 4.0       # A, primary heavy-atom contact + base/backbone split
CONTACT_CUTOFF_WIDE = 4.5  # A, second count for continuity with rna_contact_nt
SASA_POINTS = 500          # Shrake-Rupley sphere points (optional cross-check only)
APOLAR_ELEMENTS = ["C", "S"]   # hydrophobic atoms (cross-check); polar = N, O, P, ...

# van der Waals radii (Bondi) for the interface steric-clash metric. A clash = a
# protein<->RNA heavy-atom pair whose vdW spheres overlap by more than CLASH_TOL Å
# (MolProbity's 0.4 Å threshold, heavy-atom version — no hydrogens on AF3/Boltz CIFs).
VDW = {"C": 1.70, "N": 1.55, "O": 1.52, "P": 1.80, "S": 1.80, "H": 1.20,
       "F": 1.47, "CL": 1.75, "BR": 1.85, "MG": 1.73, "NA": 2.27, "K": 2.75, "ZN": 1.39}
VDW_DEFAULT = 1.70
CLASH_TOL = 0.4                # Å of overlap beyond which a contact counts as a clash

# 20 canonical amino acids (3-letter, as in AF3/Boltz CIFs) -- rows of the
# pairwise AA x {backbone, base} contact table (metric 3).
AA3 = ["ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE",
       "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL"]

# RNA atom-name partition (primes as written in AF3/Boltz CIFs, e.g. C1')
PHOSPHATE_ATOMS = {"P", "OP1", "OP2", "OP3", "O1P", "O2P", "O3P"}
SUGAR_ATOMS = {"C1'", "C2'", "C3'", "C4'", "C5'", "O2'", "O3'", "O4'", "O5'",
               "C1*", "C2*", "C3*", "C4*", "C5*", "O2*", "O3*", "O4*", "O5*"}
# anything else on a nucleotide = base

CONTACT_COLUMNS = (
    ["n_heavy_contacts_4A", "n_heavy_contacts_4p5A",
     "n_contact_res_prot", "n_contact_res_rna",
     "n_contacts_base", "n_contacts_backbone", "n_contacts_phosphate",
     "n_contacts_sugar", "frac_base_contacts",
     "n_interface_clashes", "clash_overlap_sum", "worst_clash"]
    + [f"aa_base_{aa}" for aa in AA3]      # metric 3: AA <-> RNA bases
    + [f"aa_bb_{aa}" for aa in AA3]        # metric 3: AA <-> RNA backbone
)
# --- base-edge partition for the groove split --------------------------------
# Every base atom is assigned to EXACTLY ONE edge so the parts sum to the base
# total (same additivity check as aa_base_* / aa_bb_*). Edges per the standard
# Seeman-Rosenberg-Rich definition. Atoms chemically shared between the Watson-
# Crick edge and a groove face (A:N6, G:O6, G:N2, C:N4, U:O4, C/U:O2) are assigned
# to the GROOVE face; the ring nitrogen (purine N1, pyrimidine N3) stays WC.
# That choice is nearly free: a WC-PAIRED base has ~0 free SASA on its WC atoms,
# so they bury nothing either way -- it only bites on UNPAIRED bases (the loop and
# the ss 3' tail), which is exactly where a protein really can read that face.
# The glycosidic N (purine N9, pyrimidine N1) goes on the minor face -- the side
# the sugar sits on -- and is near-fully occluded regardless.
#
# CAVEAT: major/minor groove is only strictly defined along a duplex. In loops and
# single-stranded regions these are base-EDGE labels, a chemistry proxy, not a
# geometric groove assignment.
BASE_EDGES = {
    "A": {"major": {"N6", "N7", "C5", "C6", "C8"},
          "minor": {"N3", "C2", "C4", "N9"},
          "wc":    {"N1"}},
    "G": {"major": {"O6", "N7", "C5", "C6", "C8"},
          "minor": {"N2", "N3", "C2", "C4", "N9"},
          "wc":    {"N1"}},
    "C": {"major": {"N4", "C4", "C5", "C6"},
          "minor": {"O2", "C2", "N1"},
          "wc":    {"N3"}},
    "U": {"major": {"O4", "C4", "C5", "C6"},
          "minor": {"O2", "C2", "N1"},
          "wc":    {"N3"}},
}
O2PRIME_ATOMS = {"O2'", "O2*"}   # 2'-OH: counted inside `sugar`, also reported alone

# RNA-SIDE ONLY. dsasa_total (Rosetta) counts BOTH sides of the interface, so
# dsasa_rna_total is roughly half of it -- different denominator, do not mix.
SASA_MOIETY_COLUMNS = [
    "dsasa_rna_total", "dsasa_rna_backbone", "dsasa_rna_phosphate",
    "dsasa_rna_sugar", "dsasa_rna_o2prime", "dsasa_rna_base",
    "dsasa_rna_base_major", "dsasa_rna_base_minor", "dsasa_rna_base_wc",
    "dsasa_rna_base_unassigned", "frac_dsasa_base", "frac_dsasa_major",
]

SASA_BIOTITE_COLUMNS = (["dsasa_total_biotite", "dsasa_apolar_biotite",
                         "dsasa_polar_biotite", "dsasa_apolar_frac_biotite"]
                        + SASA_MOIETY_COLUMNS)


def _rna_class(atom_name: str) -> str:
    if atom_name in PHOSPHATE_ATOMS:
        return "phosphate"
    if atom_name in SUGAR_ATOMS:
        return "sugar"
    return "base"


def _buried(sasa_free, sasa_in_complex, mask=None) -> float:
    if mask is None:
        return float(np.nansum(sasa_free) - np.nansum(sasa_in_complex))
    return float(np.nansum(sasa_free[mask]) - np.nansum(sasa_in_complex[mask]))


def _base_type(res_name):
    """Map a residue name onto one of A/G/C/U for the edge partition (None = other)."""
    r = str(res_name).strip().upper()
    if r in BASE_EDGES:
        return r
    if len(r) == 2 and r[0] == "D" and r[1] in BASE_EDGES:
        return r[1]
    if r in ("T", "DT"):
        return "U"          # thymine shares the uracil edges (5-Me sits on major)
    return None


def _rna_moiety_masks(rna) -> dict:
    """Boolean atom masks over `rna`. phosphate|sugar|base partition the array;
    base_major|base_minor|base_wc|base_unassigned partition the base atoms."""
    names = np.asarray(rna.atom_name, dtype=object)
    resn = np.asarray(rna.res_name, dtype=object)
    phos = np.array([n in PHOSPHATE_ATOMS for n in names], dtype=bool)
    sug = np.array([n in SUGAR_ATOMS for n in names], dtype=bool)
    base = ~(phos | sug)
    o2p = np.array([n in O2PRIME_ATOMS for n in names], dtype=bool)
    major = np.zeros(len(names), dtype=bool)
    minor = np.zeros(len(names), dtype=bool)
    wc = np.zeros(len(names), dtype=bool)
    for i in np.flatnonzero(base):
        bt = _base_type(resn[i])
        if bt is None:
            continue
        edges, n = BASE_EDGES[bt], names[i]
        if n in edges["major"]:
            major[i] = True
        elif n in edges["minor"]:
            minor[i] = True
        elif n in edges["wc"]:
            wc[i] = True
    return {"phosphate": phos, "sugar": sug, "base": base, "o2prime": o2p,
            "base_major": major, "base_minor": minor, "base_wc": wc,
            "base_unassigned": base & ~(major | minor | wc)}


def _sasa_rna_moieties(sasa_r_free, sasa_r_bound, rna) -> dict:
    """RNA-side buried SASA split by moiety (phosphate/sugar/base) and, within the
    bases, by groove-facing edge. Reuses SASA arrays already computed by
    _sasa_biotite -- no extra Shrake-Rupley passes."""
    out = {f"dsasa_rna_{k}": round(_buried(sasa_r_free, sasa_r_bound, m), 2)
           for k, m in _rna_moiety_masks(rna).items()}
    out["dsasa_rna_total"] = round(_buried(sasa_r_free, sasa_r_bound), 2)
    out["dsasa_rna_backbone"] = round(out["dsasa_rna_phosphate"]
                                      + out["dsasa_rna_sugar"], 2)
    tot, bse = out["dsasa_rna_total"], out["dsasa_rna_base"]
    out["frac_dsasa_base"] = round(bse / tot, 4) if tot > 0 else None
    out["frac_dsasa_major"] = (round(out["dsasa_rna_base_major"] / bse, 4)
                               if bse > 0 else None)
    return out


def _sasa_biotite(atoms, protein, rna) -> dict:
    """Optional heavy-atom Biotite dSASA cross-check (Rosetta is authoritative)."""
    try:
        sasa_c = struc.sasa(atoms, point_number=SASA_POINTS)
        sasa_p = struc.sasa(protein, point_number=SASA_POINTS)
        sasa_r = struc.sasa(rna, point_number=SASA_POINTS)
        prot_in_c = np.isin(atoms.chain_id, np.unique(protein.chain_id))
        rna_in_c = np.isin(atoms.chain_id, np.unique(rna.chain_id))
        ap_p = np.isin(protein.element, APOLAR_ELEMENTS)
        ap_r = np.isin(rna.element, APOLAR_ELEMENTS)
        ap_c = np.isin(atoms.element, APOLAR_ELEMENTS)
        total = _buried(sasa_p, sasa_c[prot_in_c]) + _buried(sasa_r, sasa_c[rna_in_c])
        apolar = (_buried(sasa_p[ap_p], sasa_c[prot_in_c & ap_c])
                  + _buried(sasa_r[ap_r], sasa_c[rna_in_c & ap_c]))
        out = dict(
            dsasa_total_biotite=round(total, 2),
            dsasa_apolar_biotite=round(apolar, 2),
            dsasa_polar_biotite=round(total - apolar, 2),
            dsasa_apolar_frac_biotite=round(apolar / total, 4) if total > 0 else None,
        )
        out.update(_sasa_rna_moieties(sasa_r, sasa_c[rna_in_c], rna))
        return out
    except Exception as e:
        return {"dsasa_biotite_error": f"{type(e).__name__}: {e}"}


def interface_geom(atoms, include_sasa: bool = False) -> dict:
    """Exact contact + pairwise-AA metrics for one complex (+ optional Biotite dSASA
    cross-check). Never raises upward -- returns {'interface_error': ...} on failure."""
    protein, rna = split_chains(atoms)
    if protein is None or rna is None:
        return {"interface_error": "could not split protein/RNA chains"}

    out: dict = {}
    if include_sasa:
        out.update(_sasa_biotite(atoms, protein, rna))

    p_heavy = protein[protein.element != "H"]
    r_heavy = rna[rna.element != "H"]
    if p_heavy.array_length() == 0 or r_heavy.array_length() == 0:
        out["contacts_error"] = "no heavy atoms on one side"
        return out

    cell = struc.CellList(r_heavy, cell_size=CONTACT_CUTOFF_WIDE)
    n4 = n45 = 0
    n_base = n_bb = n_phos = n_sug = 0
    prot_res_hit, rna_res_hit = set(), set()
    aa_base = {aa: 0 for aa in AA3}    # metric 3: AA <-> bases
    aa_bb = {aa: 0 for aa in AA3}      # metric 3: AA <-> backbone
    # steric-clash accumulators (vdW-overlap on the same neighbour search)
    p_vdw = np.array([VDW.get(str(e), VDW_DEFAULT) for e in p_heavy.element])
    r_vdw = np.array([VDW.get(str(e), VDW_DEFAULT) for e in r_heavy.element])
    n_clash = 0
    clash_sum = 0.0
    worst_clash = 0.0
    for i in range(p_heavy.array_length()):
        idx = cell.get_atoms(p_heavy.coord[i], radius=CONTACT_CUTOFF_WIDE)
        if len(idx) == 0:
            continue
        d = np.linalg.norm(r_heavy.coord[idx] - p_heavy.coord[i], axis=1)
        ov = (p_vdw[i] + r_vdw[idx]) - d          # vdW overlap for each neighbour
        cm = ov > CLASH_TOL
        if cm.any():
            ovc = ov[cm]
            n_clash += int(cm.sum())
            clash_sum += float(ovc.sum())
            worst_clash = max(worst_clash, float(ovc.max()))
        within45 = idx[d <= CONTACT_CUTOFF_WIDE]
        within4 = idx[d <= CONTACT_CUTOFF]
        n45 += len(within45)
        n4 += len(within4)
        if len(within45):
            prot_res_hit.add(p_heavy.res_id[i])
        aa = str(p_heavy.res_name[i])           # protein partner identity (metric 3)
        for j in within4:                        # base/backbone split at 4.0 A
            cls = _rna_class(str(r_heavy.atom_name[j]))
            if cls == "base":
                n_base += 1
                if aa in aa_base:
                    aa_base[aa] += 1
            else:
                n_bb += 1
                if cls == "phosphate":
                    n_phos += 1
                else:
                    n_sug += 1
                if aa in aa_bb:
                    aa_bb[aa] += 1
            rna_res_hit.add(r_heavy.res_id[j])

    out.update(
        n_heavy_contacts_4A=int(n4),
        n_heavy_contacts_4p5A=int(n45),
        n_contact_res_prot=int(len(prot_res_hit)),
        n_contact_res_rna=int(len(rna_res_hit)),
        n_contacts_base=int(n_base),
        n_contacts_backbone=int(n_bb),
        n_contacts_phosphate=int(n_phos),
        n_contacts_sugar=int(n_sug),
        frac_base_contacts=round(n_base / (n_base + n_bb), 4) if (n_base + n_bb) else None,
        n_interface_clashes=int(n_clash),
        clash_overlap_sum=round(clash_sum, 3),
        worst_clash=round(worst_clash, 3),
    )
    for aa in AA3:
        out[f"aa_base_{aa}"] = int(aa_base[aa])
        out[f"aa_bb_{aa}"] = int(aa_bb[aa])
    return out


def interface_geom_from_path(path, include_sasa: bool = False) -> dict:
    try:
        atoms = load_structure(Path(path))
    except Exception as e:
        return {"interface_error": f"load: {type(e).__name__}: {e}"}
    return interface_geom(atoms, include_sasa=include_sasa)


def aggregate_aa_matrix(rows):
    """Roll per-design aa_base_*/aa_bb_* counts into the campaign 20 x {backbone,base}
    frequency table (metric 3 heatmap). `rows` = iterable of per-design metric dicts.
    Returns {aa: {'base': freq, 'backbone': freq}} (each channel sums to 1) plus
    '_totals'."""
    base = {aa: 0 for aa in AA3}
    bb = {aa: 0 for aa in AA3}
    for r in rows:
        for aa in AA3:
            base[aa] += int(r.get(f"aa_base_{aa}", 0) or 0)
            bb[aa] += int(r.get(f"aa_bb_{aa}", 0) or 0)
    tb, tbb = sum(base.values()), sum(bb.values())
    freq = {aa: {"base": (base[aa] / tb if tb else 0.0),
                 "backbone": (bb[aa] / tbb if tbb else 0.0)} for aa in AA3}
    freq["_totals"] = {"base": tb, "backbone": tbb}
    return freq


if __name__ == "__main__":
    import sys
    for p in sys.argv[1:]:
        print(p)
        for k, v in interface_geom_from_path(p, include_sasa=True).items():
            print(f"  {k:26s} {v}")
