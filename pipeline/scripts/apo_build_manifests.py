#!/usr/bin/env python3
"""Build apo-fold chunk manifests (id, sequence, esm_ref) from a source metrics CSV.
Round-robin split so every chunk spans all stacks. Skips rows missing sequence or esm ref."""
import argparse
import csv
import os

RB = "/cluster/work/wicky/lvalsecchi/rna_binder"

ap = argparse.ArgumentParser()
ap.add_argument("--src", default=f"{RB}/analyses/cross_campaign/app_metrics_100k.csv")
ap.add_argument("--out", default=f"{RB}/apo_chunks")
ap.add_argument("-n", "--nchunks", type=int, default=50)
a = ap.parse_args()

os.makedirs(a.out, exist_ok=True)
writers = [open(f"{a.out}/chunk_{i:02d}.tsv", "w", newline="") for i in range(a.nchunks)]
ws = [csv.writer(w, delimiter="\t") for w in writers]
for w in ws:
    w.writerow(["id", "sequence", "esm_ref", "stack"])

total = kept = 0
with open(a.src) as fh:
    for d in csv.DictReader(fh):
        total += 1
        seq = (d.get("sequence") or "").strip()
        esm = (d.get("esm_cif_path") or "").strip()
        if not seq or not esm:
            continue
        ws[kept % a.nchunks].writerow([d["design_id"], seq, esm, d.get("stack", "")])
        kept += 1
for w in writers:
    w.close()
print(f"rows={total} kept={kept} skipped={total-kept} -> {a.nchunks} chunks in {a.out} (~{kept//a.nchunks}/chunk)")
