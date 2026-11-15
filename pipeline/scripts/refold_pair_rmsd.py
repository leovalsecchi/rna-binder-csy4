#!/usr/bin/env python3
"""refold_pair_rmsd.py — PAIRWISE agreement between the three refold backends.

For every design it loads the af3 / b2 / esm refold CIFs and, for each of the three
backend pairs (af3-b2, af3-esm, b2-esm), reports three RMSDs:

  rmsd_prot_<pair>  binder Ca RMSD, binder-on-binder Kabsch fit   -> do they agree on the FOLD?
  rmsd_rna_<pair>   RNA heavy-atom RMSD, RNA-on-RNA fit           -> do they agree on the TARGET?
  rmsd_pose_<pair>  binder Ca RMSD after superposing on the RNA   -> do they agree on the POSE?
                    (no second fit; the RNA transform is applied to the binder)

pose is the discriminating one: two refolders can each build a good helix bundle
(low rmsd_prot) yet dock it on opposite faces of the hairpin (high rmsd_pose).

Atoms are matched by (res_id, atom_name), so unequal atom counts degrade to the
common subset instead of returning None.

Work items are ARCHIVE TRIPLES, not designs: every design in a packed refold unit is
scored while that unit's three archives are already expanded in node-local scratch
(archive_read LRU), so each archive is extracted once per chunk.

Runs INSIDE boltzgen.sif (biotite). PYTHONPATH must include code/src and code/scripts.

    apptainer exec --bind /cluster/work/wicky --bind /tmp --env MPLCONFIGDIR=/tmp \\
        /cluster/work/wicky/containers/boltzgen.sif \\
        env PYTHONPATH=code/src:code/scripts python3 code/scripts/refold_pair_rmsd.py \\
        --csv chunk_00000.csv --output out.csv --jobs 8
"""
from __future__ import annotations

import argparse
import csv
import os
import sys
from collections import defaultdict
from multiprocessing import Pool

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "src"))

import biotite.structure as struc

from archive_read import resolve
from rna_binder.analysis.structural import load_structure, split_chains

csv.field_size_limit(10 ** 7)

# backend key -> path column in the input CSV
BACKENDS = {
    "af3": "af3_msaoff_cif_path",
    "b2": "b2_cif_path",
    "esm": "esm_cif_path",
}
PAIRS = [("af3", "b2"), ("af3", "esm"), ("b2", "esm")]
KINDS = ("prot", "rna", "pose")

RMSD_COLUMNS = [f"rmsd_{k}_{a}_{b}" for a, b in PAIRS for k in KINDS] + ["pair_rmsd_error"]

_MIN_ATOMS = 3


def _key_match(a, b):
    """Indices into (a, b) of atoms shared by (res_id, atom_name), in a's order."""
    pos = {}
    for i, (r, n) in enumerate(zip(a.res_id, a.atom_name)):
        pos.setdefault((int(r), str(n)), i)
    ia, ib = [], []
    for j, (r, n) in enumerate(zip(b.res_id, b.atom_name)):
        i = pos.get((int(r), str(n)))
        if i is not None:
            ia.append(i)
            ib.append(j)
    return np.asarray(ia, dtype=int), np.asarray(ib, dtype=int)


def _apply(transform, atoms):
    """biotite >=0.40 returns an AffineTransformation; older returns raw params."""
    if hasattr(transform, "apply"):
        return transform.apply(atoms)
    return struc.superimpose_apply(atoms, transform)


def _parts(path):
    """(protein Ca AtomArray, RNA heavy-atom AtomArray) for one refold CIF."""
    atoms = load_structure(path)
    atoms = atoms[atoms.element != "H"]
    prot, rna = split_chains(atoms)
    ca = prot[prot.atom_name == "CA"] if prot is not None else None
    return ca, rna


def _pair_rmsds(pa, pb):
    """Three RMSDs for one backend pair. pa/pb are (ca, rna) tuples."""
    out = {}
    ca_a, rna_a = pa
    ca_b, rna_b = pb

    # --- fold agreement: binder on binder
    if ca_a is not None and ca_b is not None:
        ia, ib = _key_match(ca_a, ca_b)
        if len(ia) >= _MIN_ATOMS:
            fit, _ = struc.superimpose(ca_a[ia], ca_b[ib])
            out["prot"] = float(struc.rmsd(ca_a[ia], fit))

    # --- target agreement: RNA on RNA (and the transform reused for the pose)
    rna_tr = None
    if rna_a is not None and rna_b is not None:
        ja, jb = _key_match(rna_a, rna_b)
        if len(ja) >= _MIN_ATOMS:
            fit, rna_tr = struc.superimpose(rna_a[ja], rna_b[jb])
            out["rna"] = float(struc.rmsd(rna_a[ja], fit))

    # --- pose agreement: RNA-frame superposition, binder RMSD with NO second fit
    if rna_tr is not None and ca_a is not None and ca_b is not None:
        ia, ib = _key_match(ca_a, ca_b)
        if len(ia) >= _MIN_ATOMS:
            moved = _apply(rna_tr, ca_b[ib])
            out["pose"] = float(struc.rmsd(ca_a[ia], moved))
    return out


def _score_group(item):
    """One archive triple: resolve its 3 archives once, score every design in it."""
    _triple, rows = item
    results = []
    for did, refs in rows:
        rec = {}
        try:
            parts = {}
            for be in BACKENDS:
                ref = refs.get(be)
                parts[be] = _parts(resolve(ref)) if ref else (None, None)
            for a, b in PAIRS:
                for kind, val in _pair_rmsds(parts[a], parts[b]).items():
                    rec[f"rmsd_{kind}_{a}_{b}"] = val
        except Exception as e:                       # one bad design must not kill the group
            rec["pair_rmsd_error"] = f"{type(e).__name__}: {e}"
        results.append((did, rec))
    return results


def _archive(ref):
    return ref.split(".tar.zst!")[0] if ".tar.zst!" in ref else ref


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--csv", required=True, help="input with design_id + the 3 refold path cols")
    ap.add_argument("--output", required=True, help="sidecar: design_id + RMSD_COLUMNS")
    ap.add_argument("--jobs", type=int, default=min(8, os.cpu_count() or 1))
    ap.add_argument("--limit", type=int, default=0, help="score only the first N designs (debug)")
    ap.add_argument("--progress", action="store_true")
    args = ap.parse_args()

    with open(args.csv, newline="") as f:
        rows = list(csv.DictReader(f))
    if args.limit:
        rows = rows[:args.limit]

    # group by archive triple so each archive is expanded once per worker
    groups = defaultdict(list)
    n_skip = 0
    for r in rows:
        refs = {be: (r.get(col) or "").strip() for be, col in BACKENDS.items()}
        if not all(refs.values()):
            n_skip += 1
            continue
        groups[tuple(_archive(refs[be]) for be in BACKENDS)].append((r["design_id"], refs))
    items = sorted(groups.items(), key=lambda kv: kv[0])
    print(f"{len(rows)} rows, {len(items)} archive triples, {n_skip} missing a refold, "
          f"jobs={args.jobs}", file=sys.stderr)

    scored, n_err = {}, 0
    if args.jobs > 1 and len(items) > 1:
        with Pool(args.jobs) as pool:
            for i, res in enumerate(pool.imap_unordered(_score_group, items)):
                for did, rec in res:
                    scored[did] = rec
                    n_err += int("pair_rmsd_error" in rec)
                if args.progress:
                    print(f"  triple {i + 1}/{len(items)} ({len(scored)} designs)", file=sys.stderr)
    else:
        for i, item in enumerate(items):
            for did, rec in _score_group(item):
                scored[did] = rec
                n_err += int("pair_rmsd_error" in rec)
            if args.progress:
                print(f"  triple {i + 1}/{len(items)} ({len(scored)} designs)", file=sys.stderr)

    tmp = args.output + ".tmp"
    with open(tmp, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["design_id"] + RMSD_COLUMNS, extrasaction="ignore")
        w.writeheader()
        for r in rows:
            did = r["design_id"]
            if did in scored:
                w.writerow({"design_id": did, **scored[did]})
    os.replace(tmp, args.output)
    print(f"wrote {args.output}: {len(scored)} designs, {n_err} errors", file=sys.stderr)


if __name__ == "__main__":
    main()
