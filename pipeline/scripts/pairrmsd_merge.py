#!/usr/bin/env python3
"""pairrmsd_merge.py — join the per-chunk refolder-agreement RMSD sidecars with the
ipTM/stack columns of the master table into ONE slim analysis table for plotting.

  design_id, stack, campaign, af3_msaoff_iptm, b2_iptm, esm_iptm, consensus_iptm,
  rmsd_{prot,rna,pose}_{af3_b2,af3_esm,b2_esm}

Only the ipTM columns are read from the 6 GB master (usecols), so this stays light.
"""
import argparse
import glob
import os

import pandas as pd

ID = "design_id"
META = ["stack", "campaign", "af3_msaoff_iptm", "b2_iptm", "esm_iptm", "consensus_iptm"]
PAIRS = ["af3_b2", "af3_esm", "b2_esm"]
RMSD = [f"rmsd_{k}_{p}" for p in PAIRS for k in ("prot", "rna", "pose")]

ap = argparse.ArgumentParser()
ap.add_argument("--master", required=True)
ap.add_argument("--sidecar-dir", required=True)
ap.add_argument("--out", required=True)
a = ap.parse_args()

files = sorted(glob.glob(os.path.join(a.sidecar_dir, "pairrmsd_*.csv")))
print(f"{len(files)} sidecars", flush=True)
sc = pd.concat([pd.read_csv(f, usecols=[ID] + RMSD, low_memory=False) for f in files],
               ignore_index=True).drop_duplicates(ID)
print(f"sidecar rows: {len(sc)}", flush=True)

m = pd.read_csv(a.master, usecols=[ID] + META, low_memory=False)
print(f"master rows: {len(m)}", flush=True)

df = m.merge(sc, on=ID, how="left")
for c in RMSD:
    df[c] = pd.to_numeric(df[c], errors="coerce").astype("float32")

n = df[RMSD].notna().any(axis=1).sum()
print(f"merged: {len(df)} rows, {n} with at least one pair RMSD "
      f"({100.0 * n / max(len(df), 1):.1f}%)", flush=True)
for c in RMSD:
    print(f"  {c:22s} n={int(df[c].notna().sum()):>9d}  median={df[c].median():.3f}", flush=True)

if a.out.endswith(".parquet"):
    df.to_parquet(a.out, index=False)
else:
    df.to_csv(a.out, index=False)
print(f"wrote {a.out}", flush=True)
