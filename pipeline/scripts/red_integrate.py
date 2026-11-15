#!/usr/bin/env python3
"""Merge the red_* chunk CSVs and join them onto app_metrics_100k.csv by design_id.
Backs up the target first. Idempotent: drops any pre-existing red_* columns before join.
"""
import glob
import shutil
import sys
import pandas as pd

RB = "/cluster/work/wicky/lvalsecchi/rna_binder"
CC = f"{RB}/analyses/cross_campaign"
TARGET = f"{CC}/app_metrics_100k.csv"

chunks = sorted(glob.glob(f"{RB}/red_chunks/chunk_*_red.csv"))
if not chunks:
    sys.exit("no red chunk CSVs found")
red = pd.concat([pd.read_csv(c) for c in chunks], ignore_index=True)
red = red.drop_duplicates("id")
red_cols = [c for c in red.columns if c.startswith("red_")]
n_scored = int(red["red_total_per_ifaceres"].notna().sum()) if "red_total_per_ifaceres" in red else 0
print(f"red chunks: {len(chunks)} | rows: {len(red)} | scored: {n_scored} | red cols: {len(red_cols)}")

red = red[["id"] + red_cols].rename(columns={"id": "design_id"})

base = pd.read_csv(TARGET, low_memory=False)
base = base[[c for c in base.columns if not c.startswith("red_")]]   # idempotent
merged = base.merge(red, on="design_id", how="left")

shutil.copy(TARGET, TARGET + ".bak_pre_red")
merged.to_csv(TARGET, index=False)
matched = int(merged["red_total_per_ifaceres"].notna().sum())
print(f"wrote {TARGET}: {len(merged)} rows x {len(merged.columns)} cols | red matched: {matched}/{len(merged)}")
