#!/usr/bin/env python3
"""Add the footn_<tag>_<backend> columns to an existing binder_filter store.

ADDITIVE: writes new .npy files beside the existing ones and leaves every current column
untouched, so the store does not need rebuilding and nothing already working can regress.

The store's columns are POSITIONAL — row i of every .npy is the same design — so the merged
CSV has to be mapped onto the store's design_id order, not appended in file order. Any design
in the store with no row in the CSV stays NaN, which the filter treats as FAIL (the store's
"NaN fails closed" rule), so a partial merge can never silently pass designs.

usage: build_footn_cols.py <store_dir> <foot_counts_merged.csv> [--dry-run]
"""
import csv
import os
import sys

import numpy as np

TAGS = ("005", "010", "015", "020", "025", "030", "034", "040", "050", "075", "100")
BACKENDS = ("af3", "esm")          # b2 is excluded from this tool by policy, at every layer


def main():
    store_dir, csv_path = sys.argv[1], sys.argv[2]
    dry = "--dry-run" in sys.argv

    ids = np.load(os.path.join(store_dir, "design_id.npy"), allow_pickle=False)
    n = ids.shape[0]
    print(f"store {store_dir}: {n:,} rows", flush=True)
    pos = {}
    for i, b in enumerate(ids):
        pos[b.decode() if isinstance(b, bytes) else str(b)] = i
    print(f"indexed {len(pos):,} design_ids", flush=True)
    if len(pos) != n:
        print(f"*** {n - len(pos):,} DUPLICATE design_ids in store — aborting"); sys.exit(1)

    cols = {f"footn_{t}_{be}": np.full(n, np.nan, dtype=np.float32)
            for t in TAGS for be in BACKENDS}
    src = {f"footn_{t}_{be}": f"foot_tp_f{t}__{be}" for t in TAGS for be in BACKENDS}

    seen = miss = 0
    with open(csv_path, newline="") as fh:
        r = csv.DictReader(fh)
        missing_src = [s for s in src.values() if s not in r.fieldnames]
        if missing_src:
            print(f"*** merged CSV lacks columns: {missing_src[:4]} — aborting"); sys.exit(1)
        for row in r:
            i = pos.get(row["design_id"])
            if i is None:
                miss += 1
                continue
            for out_name, in_name in src.items():
                v = row[in_name]
                if v != "":
                    cols[out_name][i] = float(v)
            seen += 1
            if seen % 1_000_000 == 0:
                print(f"  {seen:,} rows mapped", flush=True)

    print(f"mapped {seen:,} rows; {miss:,} CSV rows not in store")
    ref = cols[f"footn_010_af3"]
    cov = np.isfinite(ref).mean() * 100
    print(f"coverage of footn_010_af3: {cov:.2f}% of store rows")
    if cov < 99.0:
        print(f"*** coverage below 99% — check the merge is complete before shipping")

    if dry:
        print("dry run — nothing written"); return
    for name, arr in cols.items():
        np.save(os.path.join(store_dir, name + ".npy"), arr)
    print(f"wrote {len(cols)} columns to {store_dir}")


if __name__ == "__main__":
    main()
