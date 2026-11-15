#!/usr/bin/env python3
"""Merge the core_* chunk CSVs and join them onto app_metrics_100k.csv by design_id.
Backs up the target first. Idempotent: drops any pre-existing core_* columns before join.
Clone of red_integrate.py.
"""
import glob
import shutil
import sys
import pandas as pd

RB = "/cluster/work/wicky/lvalsecchi/rna_binder"
CC = f"{RB}/analyses/cross_campaign"
TARGET = f"{CC}/app_metrics_100k.csv"

chunks = sorted(glob.glob(f"{RB}/core_chunks/chunk_*_core.csv"))
if not chunks:
    sys.exit("no core chunk CSVs found")
core = pd.concat([pd.read_csv(c) for c in chunks], ignore_index=True)
core = core.drop_duplicates("id")
core_cols = [c for c in core.columns if c.startswith("core_")]
n_scored = int(core["core_packing_score"].notna().sum()) if "core_packing_score" in core else 0
print(f"core chunks: {len(chunks)} | rows: {len(core)} | scored: {n_scored} | core cols: {len(core_cols)}")

core = core[["id"] + core_cols].rename(columns={"id": "design_id"})

base = pd.read_csv(TARGET, low_memory=False)
base = base[[c for c in base.columns if not c.startswith("core_")]]   # idempotent
merged = base.merge(core, on="design_id", how="left")

shutil.copy(TARGET, TARGET + ".bak_pre_core")
merged.to_csv(TARGET, index=False)
matched = int(merged["core_packing_score"].notna().sum())
print(f"wrote {TARGET}: {len(merged)} rows x {len(merged.columns)} cols | core matched: {matched}/{len(merged)}")
