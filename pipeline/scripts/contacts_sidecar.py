#!/usr/bin/env python3
"""contacts_sidecar.py — fill the exact molecular-contact + pairwise AA<->{backbone,
base} columns in a metrics CSV, in place. Runs INSIDE boltzgen.sif (Biotite).

Thin CSV driver over rna_binder.analysis.interface_geom.interface_geom_from_path.
For every row (matched by --path-col, default af3_msaoff_cif_path) it computes the
heavy-atom contact counts and the 20 x {backbone,base} pairwise table (see
interface_geom for the column list). One compute per distinct structure.

Stdlib-only IO (csv, multiprocessing); the numeric work is numpy/biotite inside
interface_geom. Requires PYTHONPATH to include code/src (for the rna_binder
package) and code/scripts (for archive_read) -- the caller sets this.

Usage (inside boltzgen.sif, from the repo root):
    apptainer exec --bind /cluster/work/wicky --bind /tmp --env MPLCONFIGDIR=/tmp \\
        /cluster/work/wicky/containers/boltzgen.sif \\
        env PYTHONPATH=code/src:code/scripts python3 \\
        code/scripts/contacts_sidecar.py --csv analyses/<name>/sample10k_metrics.csv --progress
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
from rna_binder.analysis.interface_geom import (
    interface_geom_from_path, CONTACT_COLUMNS, SASA_BIOTITE_COLUMNS)

csv.field_size_limit(10 ** 7)
PATH_COL_DEFAULT = "af3_msaoff_cif_path"


WITH_SASA = False       # set from --sasa in main(); read by _compute_one workers


def _compute_one(ref):
    """Return (ref, metrics dict). Never raises -- interface_geom_from_path traps."""
    try:
        local = resolve(ref)
    except Exception as e:
        return ref, {"interface_error": f"resolve: {type(e).__name__}: {e}"}
    return ref, interface_geom_from_path(local, include_sasa=WITH_SASA)


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--csv", required=True)
    ap.add_argument("--output", default=None)
    ap.add_argument("--path-col", default=PATH_COL_DEFAULT)
    ap.add_argument("--jobs", type=int, default=min(16, os.cpu_count() or 1))
    ap.add_argument("--progress", action="store_true")
    ap.add_argument("--sasa", action="store_true",
                    help="also compute the Biotite dSASA cross-check AND the RNA "
                         "moiety/groove decomposition (dsasa_rna_*, frac_dsasa_*). "
                         "Adds one Shrake-Rupley pass per structure.")
    args = ap.parse_args()

    with open(args.csv, newline="") as f:
        reader = csv.DictReader(f)
        rows = list(reader)
        fieldnames = list(reader.fieldnames or [])
    if args.path_col not in fieldnames:
        sys.exit(f"--path-col '{args.path_col}' not in CSV columns")
    global WITH_SASA, COLUMNS
    WITH_SASA = args.sasa
    COLUMNS = CONTACT_COLUMNS + (SASA_BIOTITE_COLUMNS if args.sasa else [])
    for c in COLUMNS:
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
                n_err += int("interface_error" in m or "contacts_error" in m)
                if args.progress and i % 500 == 0:
                    print(f"  {i}/{len(refs)} structures", file=sys.stderr)
    else:
        for i, ref in enumerate(refs):
            _, m = _compute_one(ref)
            res[ref] = m
            n_err += int("interface_error" in m or "contacts_error" in m)
            if args.progress and i % 500 == 0:
                print(f"  {i}/{len(refs)} structures", file=sys.stderr)

    n_filled = 0
    for r in rows:
        m = res.get(r.get(args.path_col, ""), {})
        ok = "n_heavy_contacts_4A" in m
        n_filled += int(ok)
        for c in COLUMNS:
            v = m.get(c)
            r[c] = "" if v is None else v

    out = args.output or args.csv
    tmp = out + ".tmp"
    with open(tmp, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()
        w.writerows(rows)
    os.replace(tmp, out)
    print(f"Wrote {len(COLUMNS)} contact/pairwise cols for {n_filled}/{len(rows)} "
          f"rows ({n_err} structures failed) -> {out}", file=sys.stderr)


if __name__ == "__main__":
    main()
