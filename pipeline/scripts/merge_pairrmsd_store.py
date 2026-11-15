#!/usr/bin/env python3
"""merge_pairrmsd_store.py — add AF3<->ESM refolder-agreement columns to the filter store.

Writes two new .npy columns next to the existing ones, joined on design_id:

    pose_disagree  <- rmsd_pose_af3_esm   how far apart the two refolds DOCK (RNA-aligned,
                                          no second fit) — "do they agree where it binds?"
    fold_disagree  <- rmsd_prot_af3_esm   binder Ca RMSD after a binder-on-binder fit —
                                          "do they agree what it looks like?"

Sources are two files with the SAME schema but different vintages:
  v2-v7  analyses/cross_campaign/pairrmsd_v2v7.csv.gz   (the 2026-07-30 run)
  v8     analyses/cross_campaign/pairrmsd_chunks_v8/*   (sidecars, this run)

B2 columns exist in both sources and are deliberately NOT read: the tool has no b2 column and
must not gain one. This is the same rule the store builder enforces with _assert_no_consensus.

A design present in the store but absent from both sources stays NaN, which the filter treats
as failing — so the coverage report at the end is the number that matters, not a nicety.
"""
from __future__ import annotations

import argparse
import csv
import glob
import gzip
import os
import sys

import numpy as np

# (store column, source column). Nothing here may name b2.
COLS = [("pose_disagree", "rmsd_pose_af3_esm"),
        ("fold_disagree", "rmsd_prot_af3_esm")]

ap = argparse.ArgumentParser()
ap.add_argument("--store", required=True)
ap.add_argument("--v2v7", required=True, help="pairrmsd_v2v7.csv.gz")
ap.add_argument("--v8-dir", required=True, help="dir of pairrmsd_*.csv sidecars")
ap.add_argument("--dry-run", action="store_true")
a = ap.parse_args()

for _, src in COLS:
    toks = src.lower().split("_")
    if "b2" in toks:
        sys.exit(f"REFUSING: {src} names b2; the tool must never read a Boltz-2 value")

did = np.load(os.path.join(a.store, "design_id.npy"))
n = did.shape[0]
print(f"store: {n:,} rows", flush=True)
pos = {d: i for i, d in enumerate(did.tolist())}
del did

out = {k: np.full(n, np.nan, dtype=np.float32) for k, _ in COLS}
filled = {k: 0 for k, _ in COLS}
unmatched = 0


def ingest(reader):
    """Fill `out` from one source. Returns rows matched into the store."""
    global unmatched
    got = 0
    for row in reader:
        i = pos.get(row["design_id"].encode())
        if i is None:
            unmatched += 1
            continue
        got += 1
        for key, src in COLS:
            v = row.get(src, "")
            if v not in ("", "NA", "nan", "NaN"):
                try:
                    out[key][i] = float(v)
                    filled[key] += 1
                except ValueError:
                    pass
    return got


with gzip.open(a.v2v7, "rt") as f:
    print(f"{os.path.basename(a.v2v7)}: matched {ingest(csv.DictReader(f)):,}", flush=True)

files = sorted(glob.glob(os.path.join(a.v8_dir, "pairrmsd_*.csv")))
print(f"v8 sidecars: {len(files)}", flush=True)
tot = 0
for j, fn in enumerate(files):
    with open(fn) as f:
        tot += ingest(csv.DictReader(f))
    if (j + 1) % 150 == 0:
        print(f"  ...{j+1}/{len(files)} sidecars, {tot:,} matched", flush=True)
print(f"v8 sidecars: matched {tot:,}", flush=True)

print(f"\nunmatched source rows (not in store): {unmatched:,}")
for key, src in COLS:
    v = out[key]
    ok = np.isfinite(v)
    print(f"\n{key}  <- {src}")
    print(f"   coverage {ok.sum():,} / {n:,}  ({100*ok.sum()/n:.2f}%)   missing {n-ok.sum():,}")
    if ok.any():
        q = np.percentile(v[ok], [5, 25, 50, 75, 95])
        print("   " + "  ".join(f"p{p}={x:.1f}" for p, x in zip((5, 25, 50, 75, 95), q)) + " A")

if a.dry_run:
    print("\ndry run — nothing written")
    sys.exit(0)

for key, _ in COLS:
    p = os.path.join(a.store, key + ".npy")
    np.save(p, out[key])
    print(f"wrote {p}  ({os.path.getsize(p)/1e6:.1f} MB)")
