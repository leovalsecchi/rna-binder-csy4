#!/usr/bin/env python3
"""Build B2 + ESM interface-H-bond input CSVs (id, abs_path, stack) from
app_metrics_100k.csv, chunked. abs_path = <backend>_cif_path (an
'archive.tar.zst!member' ref that hbonds_pyrosetta.py resolves via archive_read).

    python3 build_hbond_be_manifests.py <app_metrics_100k.csv> <outdir> [N=20]
"""
import csv, os, sys
csv.field_size_limit(10 ** 7)

SRC, OUTDIR = sys.argv[1], sys.argv[2]
N = int(sys.argv[3]) if len(sys.argv) > 3 else 20
os.makedirs(OUTDIR, exist_ok=True)

rows = list(csv.DictReader(open(SRC)))
for be, col in [("b2", "b2_cif_path"), ("esm", "esm_cif_path")]:
    recs = [(r["design_id"], r[col], r.get("stack", "")) for r in rows if r.get(col)]
    for i in range(N):
        chunk = recs[i::N]
        with open(f"{OUTDIR}/chunk_{i:02d}_{be}_hb.csv", "w", newline="") as fh:
            w = csv.writer(fh)
            w.writerow(["id", "abs_path", "stack"])
            w.writerows(chunk)
    print(f"{be}: {len(recs)} refs -> {N} chunks in {OUTDIR}")
