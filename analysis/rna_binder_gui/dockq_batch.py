#!/usr/bin/env python
"""Batch DockQ over a run's AF3 refolds vs designed references (citable metric).

DockQ (Mirabello & Wallner, Bioinformatics 2024) reports the CAPRI metrics:
    LRMSD   ligand RMSD  (superpose receptor/binder, RMSD over the RNA) — the
            headline; this is what our in-house ligand_rmsd reproduces.
    iRMSD   interface RMSD (local, interface residues only)
    fnat    fraction of designed interface contacts recovered
    DockQ   combined 0-1 quality score

Reference = the designed complex, so this is a *self-consistency* measure. The
protein↔protein / RNA↔RNA chain mapping is built by composition, so it is correct
regardless of how each file labels chains (designed: RNA=A/protein=B; refold the
other way).

    dockq_batch.py RUN_DIR --out dockq.csv [--limit N]
"""
import argparse
import csv
import gzip
import sys
import tempfile
from pathlib import Path

import gemmi
from DockQ.DockQ import load_PDB, run_on_all_native_interfaces

sys.path.insert(0, str(Path(__file__).resolve().parent))
from ligand_rmsd import _af3_pairs, _kind  # reuse enumeration + residue classifier

COLS = ["design", "refolder", "msa", "DockQ", "LRMSD", "iRMSD", "fnat", "error"]


def _chain_types(path):
    """Map chain name -> 'protein'|'na' by residue composition."""
    st = gemmi.read_structure(str(path))
    out = {}
    for c in st[0]:
        kinds = {_kind(r.name) for r in c} - {"other"}
        if "protein" in kinds:
            out[c.name] = "protein"
        elif "na" in kinds:
            out[c.name] = "na"
    return out


def _chain_map(model_path, native_path):
    """{model_chain: native_chain} pairing protein↔protein and RNA↔RNA."""
    mt = _chain_types(model_path)
    n_by_type = {t: ch for ch, t in _chain_types(native_path).items()}
    return {mc: n_by_type[t] for mc, t in mt.items() if t in n_by_type}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("run_dir")
    ap.add_argument("--out", default="dockq_midrun_v2.csv")
    ap.add_argument("--limit", type=int, default=0, help="process only first N pairs (testing)")
    args = ap.parse_args()

    run = Path(args.run_dir)
    designs = run / "designs"
    pairs = list(_af3_pairs(run))
    if args.limit:
        pairs = pairs[:args.limit]

    tmp = Path(tempfile.mkdtemp(prefix="dockq_"))
    native_cache = {}  # design -> (decompressed_path, loaded_structure)
    rows = []
    for i, (design, msa, mc) in enumerate(pairs, 1):
        row = {c: None for c in COLS}
        row.update(design=design, refolder="af3", msa=msa, error="")
        try:
            if mc is None:
                raise FileNotFoundError("refold model.cif missing")
            if design not in native_cache:
                dpath = tmp / f"{design}.cif"
                with gzip.open(designs / f"{design}.cif.gz", "rt") as f:
                    dpath.write_text(f.read())
                native_cache[design] = (dpath, load_PDB(str(dpath)))
            dpath, native = native_cache[design]
            model = load_PDB(mc)
            res, _ = run_on_all_native_interfaces(model, native, chain_map=_chain_map(mc, dpath))
            if res:
                v = next(iter(res.values()))
                row.update(DockQ=round(v["DockQ"], 4), LRMSD=round(v["LRMSD"], 3),
                           iRMSD=round(v["iRMSD"], 3), fnat=round(v["fnat"], 4))
            else:
                row["error"] = "no interface found"
        except Exception as e:  # one bad pair must not kill the batch
            row["error"] = repr(e)
        rows.append(row)
        print(f"[{i}/{len(pairs)}] {design} {msa}: DockQ={row['DockQ']} "
              f"LRMSD={row['LRMSD']} {row['error']}", flush=True)

    with open(args.out, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=COLS)
        w.writeheader()
        w.writerows(rows)
    ok = sum(1 for r in rows if not r["error"])
    print(f"\nWrote {len(rows)} rows ({ok} ok, {len(rows) - ok} errors) to {args.out}", flush=True)


if __name__ == "__main__":
    main()
