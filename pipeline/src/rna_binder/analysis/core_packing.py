"""Binder core packing metrics (biotite port of the ChimeraX ``Bbindercore`` preset).

A batch/quantitative port of the ``rnacon``/``Bbindercore.cxc`` core view: where that
preset colours the buried core of the folded binder yellow (hydrophobic / good) vs red
(polar / clashy / bad), this computes the same picture as per-design NUMBERS so it can be
tabulated across designs. Companion to ``red_contacts.py`` (interface reds) and
``tip_contact.py`` (favourable interface contacts) — but this metric is about the FOLD of
the binder itself (its monomer core), independent of the RNA interface.

Two families of columns:

  FAITHFUL PORT (reproduces the preset exactly)
    core = protein residues whose SASA (binder ALONE, RNA dropped) < 30 Å²  [preset ::area<30]
    split hydrophobic {A V L I M F W Y P C}  (yellow/good) vs polar/charged (red/bad)
      - core_hydrophobic_frac : buried hydrophobic ÷ all buried core   (the "yellow fraction")
      - core_good_pack_contacts : C–C sidechain pairs ≤3.8 Å among hydrophobic core (yellow)
      - core_bad_polar_contacts : buried polar N/O/S sidechain pairs ≤3.5 Å (red)
      - core_clashes            : intra-core sidechain VdW overlaps ≥0.6 Å (red)

  IMPROVED — geometric packing QUALITY (density, not composition)
    per core residue: Cβ coordination = # other protein Cβ within 10 Å (higher = tighter)
      - core_packing_score  : mean over core residues of min(cn/cn_ref, 1) ∈ [0,1]  (higher=better)
      - core_underpacked_frac : fraction of core residues with cn < underpack_cn  (void proxy)
      - core_mean_cb_cn / core_mean_heavy_cn : raw coordination numbers

SASA is measured on the isolated binder monomer (the preset tiles models apart for exactly
this), so "core" means intrinsically buried fold residues, NOT interface burial. cn_ref /
cn_radius / underpack_cn are first-pass constants to be calibrated on real data + naturals.

Public entry point: ``compute_core_packing(path, ...) -> dict`` returns a flat ``core_*``
dict; never raises (failures become ``core_packing_error``), matching ``compute_red_contacts``.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import Optional

import numpy as np

from .structural import load_structure, split_chains

# --- tunables (defaults; core cutoff matches Bbindercore.cxc ::area<30.0) -----
SASA_CORE_CUTOFF = 30.0     # Å², per-residue SASA of the isolated binder -> "buried core"
PACK_DIST = 3.8             # Å, hydrophobic-core C–C packing contact (preset)
POLAR_DIST = 3.5            # Å, buried polar N/O/S contact (preset)
CN_RADIUS = 8.0             # Å, Cβ coordination shell for packing density (first+second shell)
# core_packing_score linearly maps per-core-residue Cβ coordination onto [0,1] between CN_LO
# (under-packed) and CN_HI (well-packed). CALIBRATED on a stratified 100k sample (see
# calibrate_core_packing.py); placeholders until then.
# Calibrated on real AF3 refold cores (per-core-residue Cβ CN@8Å: p5=8, p50=12, p95=15):
CN_LO = 8.0                 # coordination at/below which a core residue is fully "under-packed" (~p5)
CN_HI = 15.0                # coordination at/above which a core residue is fully "well-packed" (~p95)
UNDERPACK_CN = 9.0          # Cβ coordination below this = under-packed core residue (void proxy, ~p10)
HEAVY_CN_DIST = 5.0         # Å, inter-residue heavy-atom contact shell

# ChimeraX `clashes` defaults (same as red_contacts.py)
CLASH_OVERLAP = 0.6
CLASH_HBOND_ALLOWANCE = 0.4

# hydrophobic core set — Bbindercore INCLUDES Tyr (unlike the red metric); MSE≈Met for crystals
CORE_HYDROPHOBIC_RES = {"ALA", "VAL", "LEU", "ILE", "MET", "PHE", "TRP", "TYR", "PRO", "CYS", "MSE"}
BACKBONE = {"N", "CA", "C", "O", "OXT"}

VDW_RADII = {"C": 1.70, "N": 1.55, "O": 1.52, "P": 1.80, "S": 1.80, "SE": 1.90, "H": 1.20}
_VDW_DEFAULT = 1.70


def _vdw(elem: str) -> float:
    return VDW_RADII.get(elem.upper(), _VDW_DEFAULT)


def score_core(protein) -> dict:
    """Core-packing metrics over a biotite protein AtomArray (binder alone, no RNA).

    Raises ValueError on unusable input (caught by the wrapper)."""
    import biotite.structure as struc
    from scipy.spatial import cKDTree

    protein = protein[protein.element != "H"]          # heavy atoms (cifs are already heavy-only)
    if protein.array_length() == 0:
        raise ValueError("no protein heavy atoms")

    # per-residue SASA of the ISOLATED binder (matches the preset's model-tiling isolation)
    sasa = struc.sasa(protein, point_number=200)
    starts = struc.get_residue_starts(protein, add_exclusive_stop=True)
    n_res = len(starts) - 1
    if n_res == 0:
        raise ValueError("no residues")

    res_sasa = np.zeros(n_res)
    res_name = np.empty(n_res, dtype=object)
    rep = np.full((n_res, 3), np.nan)                  # Cβ (or CA) per residue for coordination
    names = protein.atom_name
    coords = protein.coord
    for k in range(n_res):
        s, e = starts[k], starts[k + 1]
        res_sasa[k] = np.nansum(sasa[s:e])
        res_name[k] = str(protein.res_name[s]).strip().upper()
        seg = names[s:e]
        cb = np.where(seg == "CB")[0]
        ca = np.where(seg == "CA")[0]
        if len(cb):
            rep[k] = coords[s + cb[0]]
        elif len(ca):
            rep[k] = coords[s + ca[0]]

    core = res_sasa < SASA_CORE_CUTOFF
    core_n = int(core.sum())
    if core_n == 0:
        # a fully-exposed "fold" (e.g. extended/unpacked) — valid, just an empty core
        return {**{c: (0 if c.startswith(("core_n", "core_good", "core_bad", "core_clash")) else None)
                   for c in CORE_COLUMNS},
                "core_n_res": 0, "core_frac": 0.0,
                "core_has_sidechains": bool(any(a not in BACKBONE for a in np.unique(names)))}

    hydroph = np.array([n in CORE_HYDROPHOBIC_RES for n in res_name])
    core_hyd = core & hydroph
    n_polar = int((core & ~hydroph).sum())

    # --- geometric packing density: Cβ coordination of core residues -----------
    valid = ~np.isnan(rep).any(axis=1)
    cn = np.full(n_res, np.nan)
    if valid.sum() >= 2:
        tree = cKDTree(rep[valid])
        vidx = np.where(valid)[0]
        for k in np.where(core & valid)[0]:
            cn[k] = len(tree.query_ball_point(rep[k], CN_RADIUS)) - 1   # exclude self
    core_cn = cn[core & valid]
    core_cn = core_cn[~np.isnan(core_cn)]

    # inter-residue heavy-atom contact number per core residue
    atom_res = np.repeat(np.arange(n_res), np.diff(starts))
    atree = cKDTree(coords)
    heavy_cn = []
    for k in np.where(core)[0]:
        s, e = starts[k], starts[k + 1]
        seen = 0
        for ai in range(s, e):
            for j in atree.query_ball_point(coords[ai], HEAVY_CN_DIST):
                if atom_res[j] != k:
                    seen += 1
        heavy_cn.append(seen)
    heavy_cn = np.array(heavy_cn, dtype=float)

    # --- contact hygiene: good C–C pack / bad polar / clashes (inter-residue) ---
    is_sc = np.array([a not in BACKBONE for a in names])
    elem = np.array([str(x).upper() for x in protein.element])
    core_atom = core[atom_res]
    hyd_atom = core_hyd[atom_res]
    polar_res_atom = (core & ~hydroph)[atom_res]

    def _pairs(mask, radius):
        idx = np.where(mask)[0]
        if len(idx) < 2:
            return np.empty((0, 2), dtype=int)
        pr = cKDTree(coords[idx]).query_pairs(radius, output_type="ndarray")
        return idx[pr] if len(pr) else np.empty((0, 2), dtype=int)

    # GOOD: hydrophobic-core sidechain carbons, C–C ≤3.8, inter-residue
    good_m = hyd_atom & is_sc & (elem == "C")
    gp = _pairs(good_m, PACK_DIST)
    core_good = int(np.sum(atom_res[gp[:, 0]] != atom_res[gp[:, 1]])) if len(gp) else 0

    # BAD polar: core sidechain N/O/S ≤3.5, inter-residue, at least one on a polar core residue
    polar_m = core_atom & is_sc & np.isin(elem, ["N", "O", "S"])
    bp = _pairs(polar_m, POLAR_DIST)
    core_bad = 0
    if len(bp):
        diff = atom_res[bp[:, 0]] != atom_res[bp[:, 1]]
        one_polar = polar_res_atom[bp[:, 0]] | polar_res_atom[bp[:, 1]]
        core_bad = int(np.sum(diff & one_polar))

    # CLASHES: core sidechain heavy atoms, VdW overlap ≥0.6 (0.4 allowance polar–polar), inter-res
    clash_m = core_atom & is_sc
    cp = _pairs(clash_m, 2 * max(VDW_RADII.values()))
    core_clash = 0
    for a, b in cp:
        if atom_res[a] == atom_res[b]:
            continue
        d = float(np.linalg.norm(coords[a] - coords[b]))
        overlap = _vdw(elem[a]) + _vdw(elem[b]) - d
        allow = CLASH_HBOND_ALLOWANCE if (elem[a] in ("N", "O") and elem[b] in ("N", "O")) else 0.0
        if overlap - allow >= CLASH_OVERLAP:
            core_clash += 1

    def _f(x):
        return float(x) if x is not None and np.isfinite(x) else None

    mean_cn = _f(np.mean(core_cn)) if len(core_cn) else None
    return {
        "core_n_res": core_n,
        "core_frac": round(core_n / n_res, 4),
        "core_hydrophobic_frac": round(int(core_hyd.sum()) / core_n, 4),
        "core_n_buried_polar": n_polar,
        "core_buried_polar_frac": round(n_polar / core_n, 4),
        "core_good_pack_contacts": core_good,
        "core_bad_polar_contacts": core_bad,
        "core_clashes": core_clash,
        "core_mean_cb_cn": round(mean_cn, 3) if mean_cn is not None else None,
        "core_mean_heavy_cn": round(float(np.mean(heavy_cn)), 3) if len(heavy_cn) else None,
        "core_underpacked_frac": round(float(np.mean(core_cn < UNDERPACK_CN)), 4) if len(core_cn) else None,
        "core_packing_score": round(float(np.mean(np.clip((core_cn - CN_LO) / (CN_HI - CN_LO), 0.0, 1.0))), 4) if len(core_cn) else None,
        "core_has_sidechains": bool(is_sc.any()),
    }


def compute_core_packing(path: Optional[Path]) -> dict:
    """Core-packing metrics for one structure (binder chain auto-detected, RNA dropped).
    Never raises: any failure yields a ``core_packing_error`` string so the row is kept."""
    if path is None:
        return {"core_packing_error": "no path"}
    try:
        atoms = load_structure(Path(path))
        protein, _rna = split_chains(atoms)
        if protein is None or protein.array_length() == 0:
            return {"core_packing_error": "no protein chain"}
        return score_core(protein)
    except Exception as e:  # noqa: BLE001 — keep the row, report the error
        return {"core_packing_error": f"{type(e).__name__}: {e}"}


CORE_COLUMNS = (
    "core_n_res", "core_frac", "core_hydrophobic_frac", "core_n_buried_polar",
    "core_buried_polar_frac", "core_good_pack_contacts", "core_bad_polar_contacts",
    "core_clashes", "core_mean_cb_cn", "core_mean_heavy_cn", "core_underpacked_frac",
    "core_packing_score", "core_has_sidechains",
)


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("structure", help=".pdb / .cif / .cif.gz (binder + RNA, or binder alone)")
    args = ap.parse_args()
    out = compute_core_packing(Path(args.structure))
    out["file"] = args.structure
    print(json.dumps(out, indent=2))


if __name__ == "__main__":
    main()
