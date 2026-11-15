#!/usr/bin/env python3
"""Merge B2 + ESM red chunks onto app_metrics_100k.csv as red_*_b2 / red_*_esm, then
build the MAX-consensus for the 4 explorer per-ifaceres red metrics. red is a BAD
signal (lower=better), so the strict "clean only if ALL 3 refolders agree it's clean"
filter is the worst-case = max(af3,b2,esm) (any backend flagging a bad contact counts).
Keep bare red_*_per_ifaceres = max(af3,b2,esm); preserve the AF3 value as
red_*_per_ifaceres_af3. Backs up first; idempotent.
"""
import glob, shutil, sys
import pandas as pd

RB = "/cluster/work/wicky/lvalsecchi/rna_binder"
TARGET = f"{RB}/analyses/cross_campaign/app_metrics_100k.csv"
CHUNKDIR = f"{RB}/red_chunks_be"

CONS = ["red_clashes_per_ifaceres", "red_acid_phos_per_ifaceres",
        "red_hphob_polar_per_ifaceres", "red_total_per_ifaceres"]


def load_be(be):
    chunks = sorted(glob.glob(f"{CHUNKDIR}/chunk_*_{be}_red.csv"))
    if not chunks:
        sys.exit(f"no {be} red chunks in {CHUNKDIR}")
    d = pd.concat([pd.read_csv(c) for c in chunks], ignore_index=True).drop_duplicates("id")
    red_cols = [c for c in d.columns if c.startswith("red_")]
    d = d[["id"] + red_cols].rename(columns={"id": "design_id"})
    d = d.rename(columns={c: f"{c}_{be}" for c in red_cols})
    scored = int(d[f"red_total_per_ifaceres_{be}"].notna().sum())
    print(f"{be}: {len(chunks)} chunks, {len(d)} rows, {scored} scored")
    return d


base = pd.read_csv(TARGET, low_memory=False)
# idempotent: drop any prior _af3/_b2/_esm red cols
base = base[[c for c in base.columns
             if not (c.startswith("red_") and c.rsplit("_", 1)[-1] in ("af3", "b2", "esm"))]]
b2, esm = load_be("b2"), load_be("esm")
m = base.merge(b2, on="design_id", how="left").merge(esm, on="design_id", how="left")

# consensus: current bare col IS AF3 -> save as _af3, rewrite bare = max(af3,b2,esm)
# (worst-case across refolders: a design is clean only if ALL 3 agree it's clean)
for c in CONS:
    m[f"{c}_af3"] = m[c]
    m[c] = m[[f"{c}_af3", f"{c}_b2", f"{c}_esm"]].max(axis=1)  # skipna=True -> max of available

shutil.copy(TARGET, TARGET + ".bak_pre_red_consensus")
m.to_csv(TARGET, index=False)
n = int(m["red_total_per_ifaceres"].notna().sum())
agree = int((m[[f"{c}_af3" for c in CONS[:1]]].notna().all(axis=1) &
             m[["red_total_per_ifaceres_b2", "red_total_per_ifaceres_esm"]].notna().all(axis=1)).sum())
print(f"wrote {TARGET}: {len(m)} rows x {len(m.columns)} cols | consensus red_total non-null: {n} | all-3-scored: {agree}")
