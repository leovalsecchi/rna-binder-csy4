#!/usr/bin/env python3
"""Build slim per-backend manifests (design_id, <backend>_cif_path) for the groove/moiety
dSASA pass. Sibling of build_iface_be_manifests.py, but covers af3 as well -- the master
has NO dsasa_rna_* columns yet, so all three backends must be computed as sidecars.

    python3 build_groove_manifests.py <master.csv> <outdir>
"""
import csv
import os
import sys

csv.field_size_limit(10 ** 7)

SRC, OUT = sys.argv[1], sys.argv[2]
os.makedirs(OUT, exist_ok=True)
rows = list(csv.DictReader(open(SRC)))
for be, col in [("af3", "af3_msaoff_cif_path"), ("b2", "b2_cif_path"),
                ("esm", "esm_cif_path")]:
    if not rows or col not in rows[0]:
        sys.exit(f"missing column {col}")
    recs = [(r["design_id"], r[col]) for r in rows if r.get(col)]
    out = os.path.join(OUT, f"groove_manifest_{be}.csv")
    with open(out, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["design_id", col])
        w.writerows(recs)
    print(f"{be}: {len(recs)} refs -> {out}")
