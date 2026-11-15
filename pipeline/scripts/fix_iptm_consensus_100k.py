#!/usr/bin/env python3
"""fix_iptm_consensus_100k.py — give the 100k table honest ipTM consensus columns.

`consensus_iptm` in app_metrics_100k.csv is min(AF3, Boltz-2) — ESM is excluded (verified
100.00% over 99,998 rows), so it is NOT the project's min-of-3 consensus and it is not an
adaptive best-2-of-3 either. `mean_iptm` is likewise mean(AF3, Boltz-2).

Adds both real aggregates and makes the mislabelled one correct, preserving the originals:

    iptm_min3        min(af3, b2, esm)          strict consensus, all three must agree
    iptm_best2       median of the three        = min of the BEST two, drops the worst
    consensus_iptm   <- iptm_min3               (old value kept as consensus_iptm_af3b2)
    mean_iptm        untouched                  (old value; mean_iptm3 added alongside)

Both new columns are gated on all three backends being present, so a design missing a
refold gets NaN rather than a quietly-2-refolder number.
"""
import argparse
import os
import shutil

import numpy as np
import pandas as pd

BE = ["af3_msaoff_iptm", "b2_iptm", "esm_iptm"]

ap = argparse.ArgumentParser()
ap.add_argument("--table", required=True)
ap.add_argument("--backup-dir", required=True)
a = ap.parse_args()

df = pd.read_csv(a.table, low_memory=False)
print(f"read {os.path.basename(a.table)}: {len(df):,} rows x {len(df.columns)} cols", flush=True)

vals = df[BE].apply(pd.to_numeric, errors="coerce")
n_ok = vals.notna().sum(axis=1)
all3 = n_ok.eq(3)
print(f"  designs with all three ipTMs: {int(all3.sum()):,} ({100*all3.mean():.2f}%)")

df["iptm_min3"] = vals.min(axis=1).where(all3)
df["iptm_best2"] = vals.median(axis=1).where(all3)      # median of 3 == min of the best two
df["mean_iptm3"] = vals.mean(axis=1).where(all3)

if "consensus_iptm" in df.columns:
    old = pd.to_numeric(df["consensus_iptm"], errors="coerce")
    ab = vals[["af3_msaoff_iptm", "b2_iptm"]].min(axis=1)
    k = old.notna() & ab.notna()
    print(f"  old consensus_iptm == min(AF3,B2): {100*np.isclose(old[k], ab[k], atol=1e-6).mean():.2f}%"
          f"  -> preserving as consensus_iptm_af3b2")
    df["consensus_iptm_af3b2"] = df["consensus_iptm"]
    df["consensus_iptm"] = df["iptm_min3"]

for c in ("iptm_min3", "iptm_best2", "mean_iptm3"):
    s = df[c]
    print(f"  {c:12s} n={int(s.notna().sum()):>7,}  median {s.median():.4f}")

os.makedirs(a.backup_dir, exist_ok=True)
bak = os.path.join(a.backup_dir, os.path.basename(a.table))
tmp = a.table + ".new"
df.to_csv(tmp, index=False)
if sum(1 for _ in open(tmp)) != len(df) + 1:
    raise SystemExit("ABORT: row count mismatch in the replacement; original untouched")
shutil.move(a.table, bak)
shutil.move(tmp, a.table)
print(f"\nwrote {a.table} ({len(df.columns)} cols); original -> {bak}", flush=True)
