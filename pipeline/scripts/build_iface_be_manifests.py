#!/usr/bin/env python3
"""Build slim per-backend interface manifests (design_id, <backend>_cif_path) from the
master for the B2/ESM interface recompute. The producers match rows by --path-col, so each
manifest keeps design_id + that backend's CIF ref column.

    python3 build_iface_be_manifests.py <master.csv> <outdir>
"""
import csv
import os
import sys

csv.field_size_limit(10 ** 7)

SRC, OUT = sys.argv[1], sys.argv[2]
os.makedirs(OUT, exist_ok=True)
rows = list(csv.DictReader(open(SRC)))
for be, col in [("b2", "b2_cif_path"), ("esm", "esm_cif_path")]:
    if not rows or col not in rows[0]:
        sys.exit(f"missing column {col}")
    recs = [(r["design_id"], r[col]) for r in rows if r.get(col)]
    out = os.path.join(OUT, f"iface_manifest_{be}.csv")
    with open(out, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["design_id", col])
        w.writerows(recs)
    print(f"{be}: {len(recs)} refs -> {out}")
