#!/usr/bin/env python3
"""pi_stack_sidecar.py — fill the aromatic protein<->RNA pi-stacking columns in a
metrics CSV, in place. Runs INSIDE boltzgen.sif (Biotite).

Thin CSV driver over rna_binder.analysis.pi_stack.pi_stack_from_path. For every row
(matched by --path-col, default af3_msaoff_cif_path) it computes the pi-stacking
descriptors (see PISTACK_COLUMNS). One compute per distinct structure. Same pattern
as contacts_sidecar.py.

Usage (inside boltzgen.sif, from the repo root):
    apptainer exec --bind /cluster/work/wicky --bind /tmp --env MPLCONFIGDIR=/tmp \\
        /cluster/work/wicky/containers/boltzgen.sif \\
        env PYTHONPATH=code/src:code/scripts python3 \\
        code/scripts/pi_stack_sidecar.py --csv analyses/<name>/sample10k_metrics.csv --progress
"""
from __future__ import annotations

import argparse
import csv
import os
import sys
from multiprocessing import Pool

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "src"))
from archive_read import resolve
from rna_binder.analysis.pi_stack import pi_stack_from_path, PISTACK_COLUMNS

csv.field_size_limit(10 ** 7)
PATH_COL_DEFAULT = "af3_msaoff_cif_path"


def _compute_one(ref):
    """Return (ref, metrics dict). Never raises -- pi_stack_from_path traps."""
    try:
        local = resolve(ref)
    except Exception as e:
        return ref, {"pi_stack_error": f"resolve: {type(e).__name__}: {e}"}
    return ref, pi_stack_from_path(local)


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--csv", required=True)
    ap.add_argument("--output", default=None)
    ap.add_argument("--path-col", default=PATH_COL_DEFAULT)
    ap.add_argument("--col-suffix", default="",
                    help="suffix appended to every output column (e.g. _b2) for per-backend runs")
    ap.add_argument("--jobs", type=int, default=min(16, os.cpu_count() or 1))
    ap.add_argument("--progress", action="store_true")
    args = ap.parse_args()
    out_cols = [c + args.col_suffix for c in PISTACK_COLUMNS]

    with open(args.csv, newline="") as f:
        reader = csv.DictReader(f)
        rows = list(reader)
        fieldnames = list(reader.fieldnames or [])
    if args.path_col not in fieldnames:
        sys.exit(f"--path-col '{args.path_col}' not in CSV columns")
    for c in out_cols:
        if c not in fieldnames:
            fieldnames.append(c)

    refs = sorted({r.get(args.path_col, "") for r in rows if r.get(args.path_col)})
    print(f"{len(rows)} rows, {len(refs)} unique structures, jobs={args.jobs}",
          file=sys.stderr)

    res = {}
    n_err = 0
    if args.jobs > 1:
        with Pool(args.jobs) as pool:
            for i, (ref, m) in enumerate(pool.imap_unordered(_compute_one, refs, chunksize=8)):
                res[ref] = m
                n_err += int("pi_stack_error" in m)
                if args.progress and i % 500 == 0:
                    print(f"  {i}/{len(refs)} structures", file=sys.stderr)
    else:
        for i, ref in enumerate(refs):
            _, m = _compute_one(ref)
            res[ref] = m
            n_err += int("pi_stack_error" in m)
            if args.progress and i % 500 == 0:
                print(f"  {i}/{len(refs)} structures", file=sys.stderr)

    n_filled = 0
    for r in rows:
        m = res.get(r.get(args.path_col, ""), {})
        ok = "pi_stack_score" in m
        n_filled += int(ok)
        for c in PISTACK_COLUMNS:
            v = m.get(c)
            r[c + args.col_suffix] = "" if v is None else v

    out = args.output or args.csv
    tmp = out + ".tmp"
    with open(tmp, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()
        w.writerows(rows)
    os.replace(tmp, out)
    print(f"Wrote {len(out_cols)} pi-stacking cols (suffix='{args.col_suffix}') for "
          f"{n_filled}/{len(rows)} rows ({n_err} structures failed) -> {out}", file=sys.stderr)


if __name__ == "__main__":
    main()
