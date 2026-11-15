#!/usr/bin/env python3
"""Build B2 + ESM red-scoring chunk manifests (id, ref, stack) from app_metrics_100k.csv.
ref = <backend>_cif_path (already an 'archive.tar.zst!member' ref that red_score_refs
handles). Round-robin split into N chunks per backend.

    python3 build_red_be_manifests.py <app_metrics_100k.csv> <outdir> [N=20]
"""
import csv, os, sys
csv.field_size_limit(10 ** 7)

SRC, OUTDIR = sys.argv[1], sys.argv[2]
N = int(sys.argv[3]) if len(sys.argv) > 3 else 20
os.makedirs(OUTDIR, exist_ok=True)

rows = list(csv.DictReader(open(SRC)))
for be, col in [("b2", "b2_cif_path"), ("esm", "esm_cif_path")]:
    if col not in rows[0]:
        sys.exit(f"missing column {col}")
    recs = [(r["design_id"], r[col], r.get("stack", "")) for r in rows if r.get(col)]
    for i in range(N):
        chunk = recs[i::N]
        with open(f"{OUTDIR}/chunk_{i:02d}_{be}.tsv", "w", newline="") as fh:
            w = csv.writer(fh, delimiter="\t")
            w.writerow(["id", "ref", "stack"])
            w.writerows(chunk)
    print(f"{be}: {len(recs)} refs -> {N} chunks in {OUTDIR}")
