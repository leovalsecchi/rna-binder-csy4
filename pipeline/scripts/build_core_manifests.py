#!/usr/bin/env python3
"""Build 20 core-backfill chunk manifests (id, ref, stack) from app_metrics_100k.csv.
Round-robin split so every chunk spans all stacks (balanced load + representative sample)."""
import csv
import os

RB = "/cluster/work/wicky/lvalsecchi/rna_binder"
SRC = f"{RB}/analyses/cross_campaign/app_metrics_100k.csv"
OUT = f"{RB}/core_chunks"
N = 20

os.makedirs(OUT, exist_ok=True)
writers = [open(f"{OUT}/chunk_{i:02d}.tsv", "w", newline="") for i in range(N)]
ws = [csv.writer(w, delimiter="\t") for w in writers]
for w in ws:
    w.writerow(["id", "ref", "stack"])

total = kept = 0
with open(SRC) as fh:
    for d in csv.DictReader(fh):
        total += 1
        ref = (d.get("af3_msaoff_cif_path") or "").strip()
        if not ref:
            continue
        ws[kept % N].writerow([d["design_id"], ref, d.get("stack", "")])
        kept += 1
for w in writers:
    w.close()
print(f"rows={total} kept={kept} -> {N} chunks in {OUT} (~{kept//N}/chunk)")
