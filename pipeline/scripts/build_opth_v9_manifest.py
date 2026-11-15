#!/usr/bin/env python3
"""build_opth_v9_manifest.py — per-backend optH manifest for midrun_v9, mirroring
build_opth_full_manifest.py (which covers v2v7 + v8 only). Reads the v9 clean table
directly; cols match what opth_parallel.sh expects: design_id, stack, <path_col>."""
import sys
import pandas as pd

RB = "/cluster/work/wicky/lvalsecchi/rna_binder/"
CC = RB + "analyses/cross_campaign/"
BE = sys.argv[1] if len(sys.argv) > 1 else "af3"
PATHCOL = {"af3": "af3_msaoff_cif_path", "b2": "b2_cif_path",
           "esm": "esm_cif_path"}[BE]
COLS = ["design_id", "stack", PATHCOL]

parts, n = [], 0
for ch in pd.read_csv(RB + "analyses/midrun_v9/collected_clean.csv", usecols=COLS,
                      chunksize=500_000, low_memory=False):
    parts.append(ch)
    n += len(ch)
df = pd.concat(parts, ignore_index=True)
before = len(df)
df = df.drop_duplicates(subset="design_id")
df = df[df[PATHCOL].notna() & (df[PATHCOL] != "")]
out = CC + f"opth_manifest_v9_{BE}.csv"
df.to_csv(out, index=False)
print(f"v9 rows read {n}, deduped {before}->{len(df)} -> {out}", flush=True)
