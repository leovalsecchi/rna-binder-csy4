#!/usr/bin/env python3
"""Merge apo_* chunk CSVs onto app_metrics_100k.csv by design_id. Idempotent + backup.
Clone of core_integrate.py."""
import glob
import shutil
import sys
import pandas as pd

RB = "/cluster/work/wicky/lvalsecchi/rna_binder"
TARGET = f"{RB}/analyses/cross_campaign/app_metrics_100k.csv"

chunks = sorted(glob.glob(f"{RB}/apo_chunks/chunk_*_apo.csv"))
if not chunks:
    sys.exit("no apo chunk CSVs found")
apo = pd.concat([pd.read_csv(c) for c in chunks], ignore_index=True).drop_duplicates("id")
apo_cols = [c for c in apo.columns if c.startswith("apo_")]
n_scored = int(apo["apo_rmsd_vs_esm"].notna().sum()) if "apo_rmsd_vs_esm" in apo else 0
print(f"apo chunks: {len(chunks)} | rows: {len(apo)} | rmsd scored: {n_scored} | apo cols: {len(apo_cols)}")

apo = apo[["id"] + apo_cols].rename(columns={"id": "design_id"})
base = pd.read_csv(TARGET, low_memory=False)
base = base[[c for c in base.columns if not c.startswith("apo_")]]   # idempotent
merged = base.merge(apo, on="design_id", how="left")

shutil.copy(TARGET, TARGET + ".bak_pre_apo")
merged.to_csv(TARGET, index=False)
matched = int(merged["apo_rmsd_vs_esm"].notna().sum())
print(f"wrote {TARGET}: {len(merged)} rows x {len(merged.columns)} cols | apo matched: {matched}/{len(merged)}")
