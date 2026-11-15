#!/usr/bin/env python3
"""Recompute protein–RNA interface H-bond counts with PyRosetta.

Runs INSIDE pyrosetta.sif. The structural aggregator (aggregate_metrics.py
--structural) runs in boltzgen.sif, which has no usable hydrogen-adder
(reduce / OpenBabel / pdbfixer / PyMOL all absent, RDKit mis-protonates
macromolecules). AF3 and Boltz CIFs are heavy-atom only, so biotite's
struc.hbond — which needs explicit D–H…A triplets — cannot be used there.
Rosetta builds full-atom hydrogens and uses orientation-dependent H-bond
detection (distance AND angle), including nucleic acids, so it gives a real
count instead of the heavy-atom distance proxy.

This pass reads a structural-metrics CSV (the aggregate_metrics.py --structural
output), recomputes `n_h_bonds_interface` for every row from its CIF path
(`--path-col`, default `abs_path`), and writes the CSV back in place (or to
`--output`). An interface H-bond is one where one partner residue is protein
and the other is a nucleic acid (RNA/DNA), regardless of chain IDs.

Stdlib-only (csv, multiprocessing) so it does not depend on pandas being in the
container. PyRosetta is initialised once per worker process.

Usage (inside pyrosetta.sif):
    apptainer exec --bind /cluster/work/wicky --bind /tmp --env MPLCONFIGDIR=/tmp \\
        /cluster/work/wicky/containers/pyrosetta.sif \\
        python3 code/scripts/hbonds_pyrosetta.py \\
        --csv analyses/midrun_v2_all/aggregated_metrics_structural.csv --progress
"""
from __future__ import annotations

import argparse
import csv
import os
import sys
from multiprocessing import Pool

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from archive_read import resolve   # archive-aware path resolution (.tar.zst!member)

PATH_COL_DEFAULT = "abs_path"
HBOND_COL = "n_h_bonds_interface"

# Per-worker PyRosetta scorefunction, set in the pool initializer.
_SF = None


def _worker_init():
    global _SF
    import pyrosetta
    pyrosetta.init("-mute all -ignore_unrecognized_res true -ignore_zero_occupancy false")
    _SF = pyrosetta.get_fa_scorefxn()


def _count_one(path):
    """Return (path, count) for one structure, or (path, None) on any failure.

    Counts H-bonds where one residue is protein and the other is a nucleic
    acid — i.e. the protein–RNA interface — over the whole complex.
    """
    import pyrosetta
    from pyrosetta.rosetta.core.scoring.hbonds import HBondSet, fill_hbond_set
    try:
        pose = pyrosetta.pose_from_file(resolve(path))
        _SF(pose)                       # populate energies the HBondSet needs
        pose.update_residue_neighbors()
        hb = HBondSet()
        fill_hbond_set(pose, False, hb)
        n = 0
        for i in range(1, hb.nhbonds() + 1):
            b = hb.hbond(i)
            rd = pose.residue(b.don_res())
            ra = pose.residue(b.acc_res())
            d_na = rd.is_RNA() or rd.is_DNA()
            a_na = ra.is_RNA() or ra.is_DNA()
            if (rd.is_protein() and a_na) or (ra.is_protein() and d_na):
                n += 1
        return path, n
    except Exception:
        return path, None


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--csv", required=True,
                    help="structural-metrics CSV to update (in place by default)")
    ap.add_argument("--output", default=None,
                    help="write here instead of overwriting --csv")
    ap.add_argument("--path-col", default=PATH_COL_DEFAULT,
                    help=f"column holding CIF paths (default: {PATH_COL_DEFAULT})")
    ap.add_argument("--jobs", type=int, default=min(8, os.cpu_count() or 1),
                    help="parallel worker processes (default: min(8, ncpu))")
    ap.add_argument("--progress", action="store_true",
                    help="print progress to stderr")
    args = ap.parse_args()

    with open(args.csv, newline="") as f:
        reader = csv.DictReader(f)
        rows = list(reader)
        fieldnames = list(reader.fieldnames or [])
    if args.path_col not in fieldnames:
        sys.exit(f"--path-col '{args.path_col}' not in CSV columns: {fieldnames}")
    if HBOND_COL not in fieldnames:
        fieldnames.append(HBOND_COL)

    # One Rosetta run per distinct structure, even if several CSV rows share it.
    paths = sorted({r.get(args.path_col, "") for r in rows if r.get(args.path_col)})
    print(f"{len(rows)} rows, {len(paths)} unique structures, jobs={args.jobs}",
          file=sys.stderr)

    counts: dict[str, object] = {}
    n_err = 0
    if args.jobs > 1:
        with Pool(args.jobs, initializer=_worker_init) as pool:
            for i, (path, n) in enumerate(
                    pool.imap_unordered(_count_one, paths, chunksize=8)):
                counts[path] = n
                n_err += int(n is None)
                if args.progress and i % 500 == 0:
                    print(f"  {i}/{len(paths)} structures", file=sys.stderr)
    else:
        _worker_init()
        for i, path in enumerate(paths):
            _, n = _count_one(path)
            counts[path] = n
            n_err += int(n is None)
            if args.progress and i % 500 == 0:
                print(f"  {i}/{len(paths)} structures", file=sys.stderr)

    n_filled = 0
    for r in rows:
        c = counts.get(r.get(args.path_col, ""))
        r[HBOND_COL] = "" if c is None else c
        n_filled += int(c is not None)

    # Atomic write so a crash can't corrupt the 15k-row CSV.
    out = args.output or args.csv
    tmp = out + ".tmp"
    with open(tmp, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()
        w.writerows(rows)
    os.replace(tmp, out)
    print(f"Wrote {HBOND_COL} for {n_filled}/{len(rows)} rows "
          f"({n_err} structures failed) -> {out}")


if __name__ == "__main__":
    main()
