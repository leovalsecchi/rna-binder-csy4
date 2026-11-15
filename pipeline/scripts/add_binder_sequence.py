#!/usr/bin/env python3
"""Backfill the `binder_sequence` column into a structural-metrics CSV.

Runs INSIDE boltzgen.sif (biotite). For each row's CIF (`--path-col`, default
`abs_path`) it extracts the protein-chain 1-letter sequence via
rna_binder.analysis.structural.protein_sequence and adds/overwrites the
`binder_sequence` column. This avoids a full `--structural` re-aggregation when
the sequence is the only thing missing.

Stdlib csv + multiprocessing; atomic write. `abs_path` is project-root-relative,
so run with cwd = the rna_binder project root (the sbatch cd's there).

Usage (inside boltzgen.sif):
    apptainer exec --bind /cluster/work/wicky --bind /tmp --env MPLCONFIGDIR=/tmp \\
        /cluster/work/wicky/containers/boltzgen.sif \\
        python3 code/scripts/add_binder_sequence.py \\
        --csv analyses/midrun_v2_all/aggregated_metrics_structural.csv --progress
"""
from __future__ import annotations

import argparse
import csv
import os
import sys
from multiprocessing import Pool

# Make `rna_binder` importable (code/src), like aggregate_metrics.py does.
sys.path.insert(0, os.path.join(
    os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "src"))

from rna_binder.analysis.structural import load_structure, protein_sequence  # noqa: E402

PATH_COL_DEFAULT = "abs_path"
SEQ_COL = "binder_sequence"


def _seq_one(path):
    """Return (path, sequence) for one CIF, or (path, "") on any failure."""
    try:
        return path, (protein_sequence(load_structure(path)) or "")
    except Exception:
        return path, ""


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--csv", required=True, help="CSV to update (in place by default)")
    ap.add_argument("--output", default=None, help="write here instead of overwriting")
    ap.add_argument("--path-col", default=PATH_COL_DEFAULT,
                    help=f"column holding CIF paths (default: {PATH_COL_DEFAULT})")
    ap.add_argument("--jobs", type=int, default=min(8, os.cpu_count() or 1),
                    help="parallel worker processes (default: min(8, ncpu))")
    ap.add_argument("--progress", action="store_true", help="print progress to stderr")
    args = ap.parse_args()

    with open(args.csv, newline="") as f:
        reader = csv.DictReader(f)
        rows = list(reader)
        fieldnames = list(reader.fieldnames or [])
    if args.path_col not in fieldnames:
        sys.exit(f"--path-col '{args.path_col}' not in CSV columns: {fieldnames}")
    if SEQ_COL not in fieldnames:
        fieldnames.append(SEQ_COL)

    paths = sorted({r.get(args.path_col, "") for r in rows if r.get(args.path_col)})
    print(f"{len(rows)} rows, {len(paths)} unique structures, jobs={args.jobs}",
          file=sys.stderr)

    seqs: dict[str, str] = {}
    n_err = 0
    if args.jobs > 1:
        with Pool(args.jobs) as pool:
            for i, (p, s) in enumerate(pool.imap_unordered(_seq_one, paths, chunksize=16)):
                seqs[p] = s
                n_err += int(not s)
                if args.progress and i % 1000 == 0:
                    print(f"  {i}/{len(paths)} structures", file=sys.stderr)
    else:
        for i, p in enumerate(paths):
            _, s = _seq_one(p)
            seqs[p] = s
            n_err += int(not s)
            if args.progress and i % 1000 == 0:
                print(f"  {i}/{len(paths)} structures", file=sys.stderr)

    n_filled = 0
    for r in rows:
        s = seqs.get(r.get(args.path_col, ""), "")
        r[SEQ_COL] = s
        n_filled += int(bool(s))

    out = args.output or args.csv
    tmp = out + ".tmp"
    with open(tmp, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()
        w.writerows(rows)
    os.replace(tmp, out)
    print(f"Wrote {SEQ_COL} for {n_filled}/{len(rows)} rows "
          f"({n_err} structures failed) -> {out}")


if __name__ == "__main__":
    main()
