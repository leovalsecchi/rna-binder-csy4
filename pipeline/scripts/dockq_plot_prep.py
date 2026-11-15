#!/usr/bin/env python3
"""dockq_plot_prep.py — pull plot-ready extracts out of the (rewritten) master tables.

The masters are 6.8 GB / 2.3 GB, so plotting locally means shipping a small extract. This
emits three, all on the RNA-aligned convention that is now the default in `dockq_self*`
with the stock values preserved in `dockq_std_self*`:

  <out>_sample.csv.gz   stratified per-stack sample, raw rows -> violin / scatter shapes
  <out>_stats.csv       EXACT per (campaign, stack, backend) summary over ALL rows, so the
                        reported medians and pass rates are never sample estimates
  <out>_seed2k.csv      the rows for the seed self-consistency designs, whatever is asked

Usage: dockq_plot_prep.py --master M --campaign v2v7 --out DIR/v2v7 [--ids ids.txt]
"""
import argparse
import os

import numpy as np
import pandas as pd

BE = ["af3", "b2", "esm"]
GATE = 0.4
WANT = (["design_id", "stack", "campaign", "length", "mean_iptm", "consensus_iptm",
         "af3_msaoff_iptm", "b2_iptm", "esm_iptm", "dockq_consensus_n",
         "dockq_self", "dockq_std_self"] +
        [f"dockq_self_{b}" for b in BE] + [f"dockq_std_self_{b}" for b in BE])

ap = argparse.ArgumentParser()
ap.add_argument("--master", required=True)
ap.add_argument("--campaign", required=True)
ap.add_argument("--out", required=True, help="output prefix")
ap.add_argument("--ids", help="newline-separated design_ids to extract in full")
ap.add_argument("--per-stack", type=int, default=60000, help="sample rows per stack")
ap.add_argument("--seed", type=int, default=0)
a = ap.parse_args()

head = pd.read_csv(a.master, nrows=0).columns.tolist()
use = [c for c in WANT if c in head]
print(f"reading {os.path.basename(a.master)}: {len(use)} of {len(WANT)} wanted cols",
      flush=True)
df = pd.read_csv(a.master, usecols=use, low_memory=False)
if "campaign" not in df.columns:
    df["campaign"] = a.campaign
print(f"  {len(df):,} rows", flush=True)

# ---- exact stats over ALL rows -------------------------------------------------------
rows = []
for (camp, st), g in df.groupby(["campaign", "stack"], dropna=False):
    for be in BE + ["consensus"]:
        for conv, pre in (("rna", "dockq_self"), ("std", "dockq_std_self")):
            col = pre if be == "consensus" else f"{pre}_{be}"
            if col not in g.columns:
                continue
            s = g[col].dropna()
            if s.empty:
                continue
            rows.append({
                "campaign": camp, "stack": st, "backend": be, "convention": conv,
                "n": len(s), "median": s.median(), "mean": s.mean(),
                "q1": s.quantile(.25), "q3": s.quantile(.75),
                "p90": s.quantile(.90), "p99": s.quantile(.99), "max": s.max(),
                "pass_frac": float((s >= GATE).mean()), "n_pass": int((s >= GATE).sum()),
            })
pd.DataFrame(rows).to_csv(f"{a.out}_stats.csv", index=False)
print(f"wrote {a.out}_stats.csv ({len(rows)} groups)", flush=True)

# ---- stratified sample for distribution shapes ---------------------------------------
rng = np.random.RandomState(a.seed)
parts = [g if len(g) <= a.per_stack
         else g.iloc[rng.choice(len(g), a.per_stack, replace=False)]
         for _, g in df.groupby("stack", dropna=False)]
samp = pd.concat(parts, ignore_index=True)
samp.to_csv(f"{a.out}_sample.csv.gz", index=False, compression="gzip")
print(f"wrote {a.out}_sample.csv.gz ({len(samp):,} rows)", flush=True)

# ---- the seed-experiment designs, in full --------------------------------------------
if a.ids and os.path.exists(a.ids):
    ids = {x.strip() for x in open(a.ids) if x.strip()}
    sub = df[df.design_id.isin(ids)]
    sub.to_csv(f"{a.out}_seed2k.csv", index=False)
    print(f"wrote {a.out}_seed2k.csv ({len(sub):,} of {len(ids):,} ids)", flush=True)
