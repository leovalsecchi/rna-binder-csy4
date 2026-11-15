#!/usr/bin/env python3
"""balance_stacks.py — downsample a per-design metrics CSV so EVERY stack has the same
number of rows (= the smallest stack's count), by RANDOMLY discarding rows. This only
reduces rows in a NEW csv — the design files on disk are never touched. Deterministic
(seed-fixed) so it's reproducible.

  python3 balance_stacks.py -i sample10k_metrics.csv -o basil_balanced_metrics.csv [--seed 0]
"""
import argparse
import pandas as pd

ap = argparse.ArgumentParser()
ap.add_argument("--input", "-i", required=True)
ap.add_argument("--output", "-o", required=True)
ap.add_argument("--seed", type=int, default=0)
ap.add_argument("--stack-col", default="stack")
a = ap.parse_args()

df = pd.read_csv(a.input, low_memory=False)
counts = df.groupby(a.stack_col).size().sort_values()
n = int(counts.min())
print("per-stack counts (before):")
print(counts.to_string())
print(f"\nbalancing every stack to n = {n} (smallest: {counts.idxmin()}), seed={a.seed}")

parts = [g.sample(n=n, random_state=a.seed) for _, g in df.groupby(a.stack_col)]
bal = (pd.concat(parts)
       .sort_values([a.stack_col, "design_id"] if "design_id" in df.columns else [a.stack_col])
       .reset_index(drop=True))
bal.to_csv(a.output, index=False)
print(f"\nwrote {len(bal)} rows ({n}/stack x {len(counts)} stacks) -> {a.output}")
print("per-stack counts (after):")
print(bal.groupby(a.stack_col).size().to_string())
