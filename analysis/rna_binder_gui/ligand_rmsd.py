#!/usr/bin/env python
"""Ligand (RNA) RMSD vs designed — protein-frame self-consistency for binder designs.

The companion of `ca_rmsd_vs_designed`. That metric asks "did the refold fold the
binder backbone as designed?". This asks "did the refold keep the RNA sitting on
the binder *where the design put it*?" — i.e. complex-level folding consistency.

Method (standard CAPRI ligand-RMSD, binder = receptor):
  1. Superpose the refold's protein Cα onto the DESIGNED protein Cα (Kabsch).
  2. Apply that one transform to the whole refold (no second superposition).
  3. ligand_rmsd = RMSD over RNA backbone atoms, refold-vs-designed, under that
     protein-frame superposition.  <-- the headline number.

Plus an interpretable decomposition of *why* the RNA is off:
  - rot_deg / centroid_shift : the rigid-body misplacement of the RNA relative
    to the binder (orientation error + lever-arm shift).
  - rna_internal_rmsd        : RNA's own conformational drift (best-fit RNA→RNA),
    i.e. error that is NOT explained by rigid misplacement.

Chains are auto-detected by residue type and matched chain-letter-agnostically
(the designed file and the refold often label protein/RNA with swapped letters),
by (residue, atom) with an order-based fallback.

Single pair:
    ligand_rmsd.py DESIGNED.cif[.gz] REFOLD.cif
Batch over a run's AF3 refolds:
    ligand_rmsd.py --batch RUN_DIR --out out.csv
"""
from __future__ import annotations

import argparse
import csv
import json
import sys
from dataclasses import dataclass, asdict
from pathlib import Path

import gemmi
import numpy as np

# RNA backbone atom subsets. Default to the phosphate/ribose backbone — robust to
# base-flip noise, which is irrelevant to "where does the RNA sit".
_BACKBONE = {"P", "O5'", "C5'", "C4'", "C3'", "O3'", "C1'"}
_ATOM_SETS = {"backbone": _BACKBONE, "c1p": {"C1'", "P"}}  # "all" handled specially


@dataclass
class Result:
    protein_ca_rmsd: float        # sanity check; should match ca_rmsd_vs_designed
    ligand_rmsd: float            # headline: RNA RMSD under protein superposition
    rna_internal_rmsd: float      # RNA conformational drift (rigid-body removed)
    rot_deg: float                # rigid rotation of RNA vs binder frame
    centroid_shift: float         # rigid translation of RNA centroid (Angstrom)
    n_ca: int                     # protein Cα atoms matched
    n_rna_atoms: int              # RNA backbone atoms matched
    match_mode: str               # "<protein>/<rna>" e.g. "order/resid"


# --------------------------------------------------------------------------- #
# Structure parsing
# --------------------------------------------------------------------------- #
def _kind(resname: str) -> str:
    info = gemmi.find_tabulated_residue(resname)
    if info is not None:
        if info.is_amino_acid():
            return "protein"
        if info.is_nucleic_acid():
            return "na"
    return "other"


def _atoms(structure_path: str):
    """Return (protein_ca, rna_atoms) for the first model, chain-agnostic.

    protein_ca : list of (res_id, coord) for every Cα, in file order.
    rna_atoms  : list of (res_id, atom_name, coord) for every nucleic-acid atom.

    Chain identity is intentionally dropped: designed and refold often use
    different chain letters (designed RNA=A/protein=B, refold the other way), so
    we match on molecule type + residue instead. Picks the first altloc and the
    first model only. gemmi reads .cif, .cif.gz and .pdb transparently.
    """
    st = gemmi.read_structure(structure_path)
    if len(st) == 0:
        raise ValueError(f"{structure_path}: no models")
    model = st[0]
    prot_ca, rna = [], []
    for chain in model:
        for res in chain:
            kind = _kind(res.name)
            if kind == "other":
                continue
            for atom in res:
                # gemmi uses '\x00' for "no altloc"; keep those plus the primary
                # alternate conformer ('A'), drop B/C/... duplicates.
                if atom.altloc not in ("\x00", "", "A"):
                    continue
                # Some mmCIF files keep literal quoting on primed atom names
                # (e.g. '"C1\'"' instead of C1'); strip a matched outer quote pair
                # without eating the real prime.
                aname = atom.name.strip()
                if len(aname) >= 2 and aname[0] == aname[-1] and aname[0] in ('"', "'"):
                    aname = aname[1:-1]
                c = np.array([atom.pos.x, atom.pos.y, atom.pos.z])
                if kind == "protein":
                    if aname == "CA":
                        prot_ca.append((res.seqid.num, c))
                else:  # nucleic acid
                    rna.append((res.seqid.num, aname, c))
    return prot_ca, rna


# --------------------------------------------------------------------------- #
# Geometry
# --------------------------------------------------------------------------- #
def _kabsch(P: np.ndarray, Q: np.ndarray):
    """Rigid transform (R, t) minimizing ||R@P + t - Q||; rows are points."""
    Pc, Qc = P.mean(0), Q.mean(0)
    H = (P - Pc).T @ (Q - Qc)
    U, _, Vt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    R = Vt.T @ np.diag([1.0, 1.0, d]) @ U.T
    return R, Qc - R @ Pc


def _apply(R, t, X):
    return X @ R.T + t


def _rmsd(A, B):
    return float(np.sqrt(((A - B) ** 2).sum(axis=1).mean()))


def _rot_angle_deg(R):
    cos = (np.trace(R) - 1.0) / 2.0
    return float(np.degrees(np.arccos(np.clip(cos, -1.0, 1.0))))


# --------------------------------------------------------------------------- #
# Correspondence (chain-agnostic: match within molecule type)
# --------------------------------------------------------------------------- #
def _pair_protein(des, ref):
    """Pair Cα between designed/refold. Prefer residue-id match; fall back to order."""
    dd = {rid: c for rid, c in des}
    dr = {rid: c for rid, c in ref}
    keys = sorted(set(dd) & set(dr))
    # residue-id match — only trusted if ids are unique (single chain) and overlap well
    if (len(dd) == len(des) and len(dr) == len(ref) and keys
            and len(keys) >= 0.5 * min(len(dd), len(dr))):
        return (np.array([dd[k] for k in keys]),
                np.array([dr[k] for k in keys]), "resid")
    if len(des) == len(ref) and des:
        return (np.array([c for _, c in des]),
                np.array([c for _, c in ref]), "order")
    if keys:
        return (np.array([dd[k] for k in keys]),
                np.array([dr[k] for k in keys]), "resid?")
    raise ValueError("no protein Cα correspondence")


def _pair_rna(des, ref, atom_ok):
    """Pair RNA atoms by (res_id, atom_name); fall back to order."""
    fd = [(rid, a, c) for rid, a, c in des if atom_ok(a)]
    fr = [(rid, a, c) for rid, a, c in ref if atom_ok(a)]
    dd = {(rid, a): c for rid, a, c in fd}
    dr = {(rid, a): c for rid, a, c in fr}
    keys = sorted(set(dd) & set(dr))
    if (len(dd) == len(fd) and len(dr) == len(fr) and keys
            and len(keys) >= 0.5 * min(len(dd), len(dr))):
        return (np.array([dd[k] for k in keys]),
                np.array([dr[k] for k in keys]), "resid")
    if len(fd) == len(fr) and fd:
        return (np.array([c for _, _, c in fd]),
                np.array([c for _, _, c in fr]), "order")
    if keys:
        return (np.array([dd[k] for k in keys]),
                np.array([dr[k] for k in keys]), "resid?")
    raise ValueError("no RNA correspondence")


def compute(designed: str, refold: str, rna_atoms: str = "backbone") -> Result:
    prot_des, rna_des = _atoms(designed)
    prot_ref, rna_ref = _atoms(refold)
    atom_ok = (lambda n: True) if rna_atoms == "all" else (lambda n: n in _ATOM_SETS[rna_atoms])

    P_des, P_ref, pm = _pair_protein(prot_des, prot_ref)
    L_des, L_ref, rm = _pair_rna(rna_des, rna_ref, atom_ok)

    # superpose on protein, score the RNA under that one transform
    R, t = _kabsch(P_ref, P_des)
    protein_ca_rmsd = _rmsd(_apply(R, t, P_ref), P_des)
    L_ref_pf = _apply(R, t, L_ref)             # RNA in the binder's frame
    ligand_rmsd = _rmsd(L_ref_pf, L_des)

    # decompose RNA error into rigid misplacement vs internal drift
    Rr, tr = _kabsch(L_ref_pf, L_des)
    rna_internal_rmsd = _rmsd(_apply(Rr, tr, L_ref_pf), L_des)
    rot_deg = _rot_angle_deg(Rr)
    centroid_shift = float(np.linalg.norm(L_ref_pf.mean(0) - L_des.mean(0)))

    return Result(round(protein_ca_rmsd, 3), round(ligand_rmsd, 3),
                  round(rna_internal_rmsd, 3), round(rot_deg, 2),
                  round(centroid_shift, 3), len(P_des), len(L_des), f"{pm}/{rm}")


# --------------------------------------------------------------------------- #
# Batch over a run directory's AF3 refolds
# --------------------------------------------------------------------------- #
_CSV_COLS = ["design", "refolder", "msa", "ligand_rmsd", "rot_deg", "centroid_shift",
             "rna_internal_rmsd", "protein_ca_rmsd", "n_ca", "n_rna_atoms",
             "match_mode", "error"]


def _af3_pairs(run: Path):
    """Yield (design_name, msa, model_cif) by CONSTRUCTING paths from designs/.

    Avoids a recursive glob over the network mount (pathologically slow on sshfs).
    Layout: refold_af3/<msa>/outputs/<DESIGN>/<lc>/<lc>_model.cif, where <lc> is
    the lowercased design name (AF3 job-name convention).
    """
    designs = sorted((run / "designs").glob("*.cif.gz"))
    af3 = run / "refold_af3"
    msa_dirs = [d.name for d in sorted(af3.iterdir()) if d.is_dir()] if af3.exists() else []
    for dz in designs:
        design = dz.name[:-len(".cif.gz")]
        lc = design.lower()
        for msa in msa_dirs:
            out = af3 / msa / "outputs" / design
            mc = out / lc / f"{lc}_model.cif"
            if mc.exists():
                yield design, msa, str(mc)
            else:  # one cheap single-dir glob as fallback if naming differs
                alt = sorted(out.glob("*/*_model.cif"))
                yield design, msa, (str(alt[0]) if alt else None)


def batch(run_dir: str, out_csv: str, rna_atoms: str = "backbone", progress: bool = True):
    run = Path(run_dir)
    designs_dir = run / "designs"
    pairs = list(_af3_pairs(run))
    rows = []
    for i, (design, msa, mc) in enumerate(pairs, 1):
        designed = designs_dir / f"{design}.cif.gz"
        row = {c: None for c in _CSV_COLS}
        row.update(design=design, refolder="af3", msa=msa, error="")
        if mc is None:
            row["error"] = "refold model.cif missing"
        elif not designed.exists():
            row["error"] = "designed file missing"
        else:
            try:
                r = compute(str(designed), mc, rna_atoms)
                for k in ("ligand_rmsd", "rot_deg", "centroid_shift", "rna_internal_rmsd",
                          "protein_ca_rmsd", "n_ca", "n_rna_atoms", "match_mode"):
                    row[k] = getattr(r, k)
            except Exception as e:  # one bad file must not kill the batch
                row["error"] = repr(e)
        rows.append(row)
        if progress:
            print(f"[{i}/{len(pairs)}] {design} {msa}: lig={row['ligand_rmsd']} "
                  f"rot={row['rot_deg']} mode={row['match_mode']} {row['error']}", flush=True)
    with open(out_csv, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=_CSV_COLS)
        w.writeheader()
        w.writerows(rows)
    ok = sum(1 for r in rows if not r["error"])
    print(f"\nWrote {len(rows)} rows ({ok} ok, {len(rows)-ok} errors) to {out_csv}", flush=True)


# --------------------------------------------------------------------------- #
def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("designed", nargs="?", help="designed complex (PDB/mmCIF[.gz]) — reference pose")
    ap.add_argument("refold", nargs="?", help="refold prediction (PDB/mmCIF) — AF3 or Boltz-2")
    ap.add_argument("--batch", metavar="RUN_DIR", help="compute over all AF3 refolds in a run dir")
    ap.add_argument("--out", default="ligand_rmsd.csv", help="batch output CSV (default: ligand_rmsd.csv)")
    ap.add_argument("--rna-atoms", choices=["backbone", "c1p", "all"], default="backbone")
    ap.add_argument("--json", action="store_true", help="single-pair: emit one JSON object")
    args = ap.parse_args(argv)

    if args.batch:
        batch(args.batch, args.out, args.rna_atoms)
        return
    if not (args.designed and args.refold):
        ap.error("provide DESIGNED and REFOLD, or --batch RUN_DIR")

    res = compute(args.designed, args.refold, args.rna_atoms)
    if args.json:
        print(json.dumps(asdict(res)))
        return
    print(f"  protein_ca_rmsd    {res.protein_ca_rmsd:7.3f} A   ({res.n_ca} Cα, match={res.match_mode})")
    print(f"  ligand_rmsd        {res.ligand_rmsd:7.3f} A   ({res.n_rna_atoms} RNA atoms)  <- headline")
    print(f"  rna_internal_rmsd  {res.rna_internal_rmsd:7.3f} A")
    print(f"  rot_deg            {res.rot_deg:7.2f} deg")
    print(f"  centroid_shift     {res.centroid_shift:7.3f} A")


if __name__ == "__main__":
    sys.exit(main())
