#!/usr/bin/env python3
"""Build protein-ONLY ESMFold2 input YAMLs from an apo chunk TSV (id, sequence, ...).
Protein-only (no RNA entry) => esmfold2_infer.py folds the binder alone. Node-local out dir."""
import argparse
import csv
import os

ap = argparse.ArgumentParser()
ap.add_argument("manifest")
ap.add_argument("outdir")
a = ap.parse_args()
os.makedirs(a.outdir, exist_ok=True)

n = 0
with open(a.manifest) as fh:
    for row in csv.DictReader(fh, delimiter="\t"):
        seq = (row.get("sequence") or "").strip()
        if len(seq) < 5:
            continue
        # minimal Boltz-style YAML, protein chain A only (no rna) — hand-written to avoid a
        # yaml dependency here; esmfold2_infer.parse_yaml_to_spi reads id/sequence only.
        with open(os.path.join(a.outdir, f"{row['id']}.yaml"), "w") as out:
            out.write("version: 1\nsequences:\n")
            out.write("- protein:\n")
            out.write(f"    id: A\n    sequence: {seq}\n")
        n += 1
print(f"wrote {n} protein-only YAMLs -> {a.outdir}")
