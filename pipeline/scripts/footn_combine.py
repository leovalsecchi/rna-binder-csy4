#!/usr/bin/env python3
"""Combine footn shards into store columns.
usage: footn_combine.py <store_dir> <shard_dir> [col_prefix]"""
import glob
import os
import sys
import time

import numpy as np
import pandas as pd

TAGS = ("005", "010", "015", "020", "025", "030", "034", "040", "050", "075", "100")


def main():
    store_dir, shard_dir = sys.argv[1], sys.argv[2]
    PFX = sys.argv[3] if len(sys.argv) > 3 else "footn"
    t0 = time.time()
    sids = np.load(os.path.join(store_dir, "design_id.npy"), allow_pickle=False)
    n = sids.shape[0]
    pos = pd.Series(np.arange(n, dtype=np.int64),
                    index=pd.Index([b.decode() for b in sids]))
    if pos.index.has_duplicates:
        sys.exit("*** duplicate design_ids in store")
    print(f"indexed {n:,} store rows in {time.time()-t0:.0f}s", flush=True)

    cols = {f"{PFX}_{t}_{be}": np.full(n, np.nan, dtype=np.float32)
            for t in TAGS for be in ("af3", "esm")}
    total = 0
    for f in sorted(glob.glob(f"{shard_dir}/shard_*.npz")):
        z = np.load(f)
        idx = pos.reindex([b.decode() for b in z["design_id"]]).values
        ok = ~np.isnan(idx)
        i = idx[ok].astype(np.int64)
        for t in TAGS:
            cols[f"{PFX}_{t}_af3"][i] = z[f"af3_{t}"][ok]
            cols[f"{PFX}_{t}_esm"][i] = z[f"esm_{t}"][ok]
        total += i.size
    print(f"mapped {total:,} rows from shards ({time.time()-t0:.0f}s)")
    for be in ("af3", "esm"):
        cov = np.isfinite(cols[f"{PFX}_010_{be}"]).mean() * 100
        print(f"coverage {PFX}_010_{be}: {cov:.2f}%")
        if cov < 99.0:
            sys.exit("*** coverage under 99% — refusing to write")
    for name, arr in cols.items():
        np.save(os.path.join(store_dir, name + ".npy"), arr)
    print(f"wrote {len(cols)} columns in {time.time()-t0:.0f}s")


if __name__ == "__main__":
    main()
