#!/usr/bin/env python3
"""Collect B2 + ESM interface-H-bond counts from the hbond chunk CSVs, join onto
app_metrics_100k.csv as n_h_bonds_interface_b2 / _esm, then build the MIN-consensus:
H-bonds are a GOOD signal (higher=better), so "clean only if all 3 refolders agree it
H-bonds" = min(af3,b2,esm). Preserve the AF3 value as n_h_bonds_interface_af3.
Backs up first; idempotent.
"""
import glob, shutil, sys
import pandas as pd

RB = "/cluster/work/wicky/lvalsecchi/rna_binder"
TARGET = f"{RB}/analyses/cross_campaign/app_metrics_100k.csv"
CHUNKDIR = f"{RB}/hbond_chunks_be"
COL = "n_h_bonds_interface"


def load_be(be):
    chunks = sorted(glob.glob(f"{CHUNKDIR}/chunk_*_{be}_hb.csv"))
    if not chunks:
        sys.exit(f"no {be} hbond chunks")
    d = pd.concat([pd.read_csv(c) for c in chunks], ignore_index=True).drop_duplicates("id")
    if COL not in d.columns:
        sys.exit(f"{be} chunks missing {COL} (arrays not finished?)")
    d = d[["id", COL]].rename(columns={"id": "design_id", COL: f"{COL}_{be}"})
    print(f"{be}: {len(chunks)} chunks, {len(d)} rows, {int(d[f'{COL}_{be}'].notna().sum())} scored")
    return d


base = pd.read_csv(TARGET, low_memory=False)
base = base[[c for c in base.columns if c not in (f"{COL}_af3", f"{COL}_b2", f"{COL}_esm")]]
b2, esm = load_be("b2"), load_be("esm")
m = base.merge(b2, on="design_id", how="left").merge(esm, on="design_id", how="left")

m[f"{COL}_af3"] = m[COL]                                   # current bare col IS AF3
m[COL] = m[[f"{COL}_af3", f"{COL}_b2", f"{COL}_esm"]].min(axis=1)  # min-consensus

shutil.copy(TARGET, TARGET + ".bak_pre_hb_consensus")
m.to_csv(TARGET, index=False)
med_af3 = m[f"{COL}_af3"].median()
med_con = m[COL].median()
print(f"wrote {TARGET}: {len(m)} rows x {len(m.columns)} cols | "
      f"{COL} AF3 median={med_af3} -> MIN-consensus median={med_con}")
