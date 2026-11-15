#!/usr/bin/env python3
"""Fill DockQ metrics into a structural-metrics CSV (separate analysis pass).

Mirrors hbonds_pyrosetta.py. The structural aggregator (aggregate_metrics.py
--structural) runs in boltzgen.sif, which can't host DockQ (no compiler, no
prebuilt wheel). So this pass runs afterward in a DockQ-capable venv and fills
the DockQ columns for every row from its model/native CIF paths — exactly how
the H-bond column is filled by the PyRosetta pass.

DockQ (Mirabello & Wallner, Bioinformatics 2024) reports CAPRI metrics for the
refold (model) vs the designed complex (native): LRMSD (ligand RMSD), iRMSD,
fnat, and the combined DockQ score (0-1). Protein<->protein / RNA<->RNA chain
mapping is built by composition, since designed and refold files often label
chains with swapped letters.

Columns written: dockq, dockq_lrmsd, dockq_irmsd, dockq_fnat.
Reference (native) may be .cif / .pdb / .cif.gz; .gz is decompressed to a temp.

Usage (in a venv with DockQ + gemmi):
    <venv>/bin/python code/scripts/dockq_metrics.py \\
        --csv analyses/<run>/aggregated_metrics_structural.csv --jobs 8 --progress
"""
from __future__ import annotations

import argparse
import csv
import gzip
import os
import sys
import tempfile
from multiprocessing import Pool

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from archive_read import resolve   # archive-aware path resolution (.tar.zst!member)

MODEL_COL_DEFAULT = "abs_path"          # refold cif (the prediction)
NATIVE_COL_DEFAULT = "designed_source"  # designed complex (the reference)
DOCKQ_COLS = ["dockq", "dockq_lrmsd", "dockq_irmsd", "dockq_fnat"]


def _kind(resname: str) -> str:
    import gemmi
    info = gemmi.find_tabulated_residue(resname)
    if info is not None:
        if info.is_amino_acid():
            return "protein"
        if info.is_nucleic_acid():
            return "na"
    return "other"


def _chain_types(path: str) -> dict:
    """chain name -> 'protein'|'na' by residue composition (gemmi reads .gz too)."""
    import gemmi
    st = gemmi.read_structure(str(path))
    out = {}
    for c in st[0]:
        kinds = {_kind(r.name) for r in c} - {"other"}
        if "protein" in kinds:
            out[c.name] = "protein"
        elif "na" in kinds:
            out[c.name] = "na"
    return out


def _chain_map(model_path: str, native_path: str) -> dict:
    """{model_chain: native_chain} pairing protein<->protein and RNA<->RNA."""
    mt = _chain_types(model_path)
    n_by_type = {t: ch for ch, t in _chain_types(native_path).items()}
    return {mc: n_by_type[t] for mc, t in mt.items() if t in n_by_type}


def _dockq_one(pair):
    """pair = (model_path, native_path) -> (pair, metrics dict | None)."""
    model_path, native_path = pair
    model_path = resolve(model_path)
    native_path = resolve(native_path)
    tmp_path = None
    try:
        from DockQ.DockQ import load_PDB, run_on_all_native_interfaces
        native = native_path
        if native_path.endswith(".gz"):
            with gzip.open(native_path, "rt") as fh:
                content = fh.read()
            tf = tempfile.NamedTemporaryFile(mode="w", suffix=".cif", delete=False)
            tf.write(content)
            tf.close()
            tmp_path = tf.name
            native = tmp_path
        cmap = _chain_map(model_path, native)
        if not cmap:
            return pair, None
        m = load_PDB(model_path)
        n = load_PDB(native)
        res, _ = run_on_all_native_interfaces(m, n, chain_map=cmap)
        if not res:
            return pair, None
        v = next(iter(res.values()))
        return pair, {
            "dockq": round(float(v["DockQ"]), 4),
            "dockq_lrmsd": round(float(v["LRMSD"]), 3),
            "dockq_irmsd": round(float(v["iRMSD"]), 3),
            "dockq_fnat": round(float(v["fnat"]), 4),
        }
    except Exception:
        return pair, None
    finally:
        if tmp_path is not None:
            try:
                os.unlink(tmp_path)
            except OSError:
                pass


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--csv", required=True,
                    help="structural-metrics CSV to update (in place by default)")
    ap.add_argument("--output", default=None, help="write here instead of overwriting --csv")
    ap.add_argument("--model-col", default=MODEL_COL_DEFAULT,
                    help=f"column with the refold/model cif path (default: {MODEL_COL_DEFAULT})")
    ap.add_argument("--native-col", default=NATIVE_COL_DEFAULT,
                    help=f"column with the designed/native path (default: {NATIVE_COL_DEFAULT})")
    ap.add_argument("--jobs", type=int, default=min(8, os.cpu_count() or 1),
                    help="parallel worker processes (default: min(8, ncpu))")
    ap.add_argument("--progress", action="store_true", help="print progress to stderr")
    args = ap.parse_args()

    with open(args.csv, newline="") as f:
        reader = csv.DictReader(f)
        rows = list(reader)
        fieldnames = list(reader.fieldnames or [])
    for col in (args.model_col, args.native_col):
        if col not in fieldnames:
            sys.exit(f"column '{col}' not in CSV columns: {fieldnames}")
    for col in DOCKQ_COLS:
        if col not in fieldnames:
            fieldnames.append(col)

    # One DockQ run per distinct (model, native) pair, even if rows repeat it.
    pairs = sorted({(r.get(args.model_col, ""), r.get(args.native_col, ""))
                    for r in rows
                    if r.get(args.model_col) and r.get(args.native_col)})
    print(f"{len(rows)} rows, {len(pairs)} unique (model,native) pairs, jobs={args.jobs}",
          file=sys.stderr)

    results = {}
    n_err = 0
    if args.jobs > 1:
        with Pool(args.jobs) as pool:
            for i, (pair, d) in enumerate(
                    pool.imap_unordered(_dockq_one, pairs, chunksize=4)):
                results[pair] = d
                n_err += int(d is None)
                if args.progress and i % 200 == 0:
                    print(f"  {i}/{len(pairs)} pairs", file=sys.stderr)
    else:
        for i, pair in enumerate(pairs):
            _, d = _dockq_one(pair)
            results[pair] = d
            n_err += int(d is None)
            if args.progress and i % 200 == 0:
                print(f"  {i}/{len(pairs)} pairs", file=sys.stderr)

    n_filled = 0
    for r in rows:
        d = results.get((r.get(args.model_col, ""), r.get(args.native_col, "")))
        for col in DOCKQ_COLS:
            r[col] = "" if not d else d[col]
        n_filled += int(bool(d))

    # Atomic write so a crash can't corrupt the CSV.
    out = args.output or args.csv
    tmp = out + ".tmp"
    with open(tmp, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()
        w.writerows(rows)
    os.replace(tmp, out)
    print(f"Wrote DockQ cols for {n_filled}/{len(rows)} rows "
          f"({n_err} pairs failed) -> {out}")


if __name__ == "__main__":
    main()
