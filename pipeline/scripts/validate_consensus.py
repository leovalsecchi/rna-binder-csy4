#!/usr/bin/env python3
"""Assert the consensus invariants on a backfilled master.
    python3 validate_consensus.py <master.csv>
Exits non-zero on any violation."""
import sys

import pandas as pd

from interface_consensus_integrate import MIN_COLS, MAX_COLS

m = pd.read_csv(sys.argv[1], low_memory=False)
bad = 0
for col in MIN_COLS + MAX_COLS:
    trio = [f"{col}_af3", f"{col}_b2", f"{col}_esm"]
    if col not in m.columns or not all(t in m.columns for t in trio):
        print(f"MISSING cols for {col}")
        bad += 1
        continue
    vals = m[trio].apply(pd.to_numeric, errors="coerce")
    want = vals.min(axis=1) if col in MIN_COLS else vals.max(axis=1)
    got = pd.to_numeric(m[col], errors="coerce")
    mism = int((~((got - want).abs() < 1e-6) & want.notna()).sum())
    if mism:
        print(f"VIOLATION {col}: {mism} rows bare != {'min' if col in MIN_COLS else 'max'}(trio)")
        bad += 1
print(f"consensus_n_structures value counts:\n"
      f"{m['consensus_n_structures'].value_counts().sort_index()}")
print("OK" if not bad else f"{bad} PROBLEMS")
sys.exit(1 if bad else 0)
