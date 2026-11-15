#!/usr/bin/env python3
"""Build the footn_* store columns DIRECTLY from the sidecar chunk files.

Skips the merged CSV entirely. The previous route (csv.DictReader -> csv.DictWriter over 178
columns) ran at ~1,300 rows/s because it materialised two Python dicts per row; pandas' C
parser reading only the 12 columns we need is ~50x faster, and writing .npy instead of a 7 GB
CSV removes the other half of the cost.

Chunk pairs are independent, so this is trivially parallel too — but with the right reader a
single core is enough, so the parallel machinery is not worth its complexity here.

usage: footn_direct.py <store_dir> [--dry-run]
"""
import glob
import os
import sys
import time

import numpy as np
import pandas as pd

RB = "/cluster/work/wicky/lvalsecchi/rna_binder"
SRCS = [f"{RB}/analyses/cross_campaign/foot_clash/out2",
        f"{RB}/analyses/cross_campaign/foot_clash_v9/out2"]
TAGS = ("005", "010", "015", "020", "025", "030", "034", "040", "050", "075", "100")
BACKENDS = ("af3", "esm")
USE = ["design_id", "foot_status"] + [f"foot_tp_f{t}" for t in TAGS]


def main():
    store_dir = sys.argv[1]
    dry = "--dry-run" in sys.argv
    t0 = time.time()

    ids = np.load(os.path.join(store_dir, "design_id.npy"), allow_pickle=False)
    n = ids.shape[0]
    pos = pd.Series(np.arange(n, dtype=np.int64),
                    index=pd.Index([b.decode() for b in ids]))
    if pos.index.has_duplicates:
        sys.exit("*** duplicate design_ids in store — aborting")
    print(f"store {n:,} rows indexed in {time.time()-t0:.0f}s", flush=True)

    cols = {f"footn_{t}_{be}": np.full(n, np.nan, dtype=np.float32)
            for t in TAGS for be in BACKENDS}

    for be in BACKENDS:
        files = []
        for d in SRCS:
            files += sorted(glob.glob(f"{d}/foot_{be}_*.csv"))
        got = 0
        for k, f in enumerate(files):
            df = pd.read_csv(f, usecols=USE, low_memory=False)
            df = df[df["foot_status"] == "ok"]
            if df.empty:
                continue
            idx = pos.reindex(df["design_id"].values).values
            ok = ~np.isnan(idx)
            idx = idx[ok].astype(np.int64)
            for t in TAGS:
                cols[f"footn_{t}_{be}"][idx] = df[f"foot_tp_f{t}"].values[ok]
            got += idx.size
            if (k + 1) % 500 == 0:
                print(f"  {be} {k+1}/{len(files)} files, {got:,} rows, "
                      f"{time.time()-t0:.0f}s", flush=True)
        print(f"{be}: {got:,} rows mapped ({time.time()-t0:.0f}s)", flush=True)

    for be in BACKENDS:
        cov = np.isfinite(cols[f"footn_010_{be}"]).mean() * 100
        print(f"coverage footn_010_{be}: {cov:.2f}%")
        if cov < 99.0:
            print("*** coverage under 99% — refusing to write"); sys.exit(1)

    if dry:
        print("dry run, nothing written"); return
    for name, arr in cols.items():
        np.save(os.path.join(store_dir, name + ".npy"), arr)
    print(f"wrote {len(cols)} columns in {time.time()-t0:.0f}s total")


if __name__ == "__main__":
    main()
