#!/usr/bin/env python3
"""One shard of the footn column build: read a slice of chunk files, emit a compact .npz.

The cost here is PER FILE (~0.41 s of pandas setup + Lustre open latency on a 2500-row file),
not per row, which is exactly the case parallelism fixes: 5,460 files split N ways finishes in
1/N the time. Each shard emits only its own rows, so the partials stay small; the combine step
does the design_id -> store-row mapping once, vectorised.

af3 and esm chunk files with the same index derive from the same input chunk, so their rows are
aligned 1:1. That is ASSERTED here, not assumed.

usage: footn_shard.py <shard_idx> <n_shards> <outdir>
"""
import glob
import os
import sys

import numpy as np
import pandas as pd

RB = "/cluster/work/wicky/lvalsecchi/rna_binder"
SRCS = [f"{RB}/analyses/cross_campaign/foot_clash/out2",
        f"{RB}/analyses/cross_campaign/foot_clash_v9/out2"]
TAGS = ("005", "010", "015", "020", "025", "030", "034", "040", "050", "075", "100")
FAM = os.environ.get("FOOTN_FAMILY", "f")   # f = foot-side atoms, b = BINDER-side atoms
assert FAM in ("f", "b")
USE = ["design_id", "foot_status"] + [f"foot_tp_{FAM}{t}" for t in TAGS]


def file_list():
    fs = []
    for d in SRCS:
        fs += sorted(glob.glob(f"{d}/foot_af3_*.csv"))
    return fs


def main():
    shard, nsh, outdir = int(sys.argv[1]), int(sys.argv[2]), sys.argv[3]
    os.makedirs(outdir, exist_ok=True)
    fs = file_list()
    mine = fs[shard::nsh]                      # strided, so shards are I/O-balanced
    ids, af, es = [], {t: [] for t in TAGS}, {t: [] for t in TAGS}
    nbad = 0
    for f in mine:
        g = f.replace("foot_af3_", "foot_esm_")
        if not os.path.exists(g):
            continue
        a = pd.read_csv(f, usecols=USE, low_memory=False)
        b = pd.read_csv(g, usecols=USE, low_memory=False)
        if len(a) != len(b) or not (a["design_id"].values == b["design_id"].values).all():
            nbad += 1
            continue
        ok = (a["foot_status"] == "ok").values & (b["foot_status"] == "ok").values
        if not ok.any():
            continue
        ids.append(a["design_id"].values[ok])
        for t in TAGS:
            af[t].append(a[f"foot_tp_{FAM}{t}"].values[ok].astype(np.int16))
            es[t].append(b[f"foot_tp_{FAM}{t}"].values[ok].astype(np.int16))
    if not ids:
        print(f"shard {shard}: nothing"); return
    out = {"design_id": np.concatenate(ids).astype("S64")}
    for t in TAGS:
        out[f"af3_{t}"] = np.concatenate(af[t])
        out[f"esm_{t}"] = np.concatenate(es[t])
    np.savez_compressed(os.path.join(outdir, f"shard_{shard:03d}.npz"), **out)
    print(f"shard {shard}: {len(mine)} files, {out['design_id'].size:,} rows, "
          f"{nbad} misaligned pairs")


if __name__ == "__main__":
    main()
