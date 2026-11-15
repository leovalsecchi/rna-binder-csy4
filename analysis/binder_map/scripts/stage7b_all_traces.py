#!/usr/bin/env python3
"""Extract the CA trace of EVERY backbone, not just cluster medoids.

Clicking a point should show that design's structure. Showing its cluster's
medoid instead is a different molecule, which is worse than showing nothing.

Storage: coordinates are centred per structure and quantised to int16 at 0.01 A,
which is far finer than the backbone accuracy of any of these models and about
a quarter the size of float32 JSON. Written as one flat binary in node_id order
so the browser can slice it without parsing anything.

Usage:  python scripts/stage7b_all_traces.py --layer design
"""
from __future__ import annotations

import argparse
import gzip
import tarfile
from pathlib import Path

import numpy as np
import pandas as pd

FM = Path("/cluster/work/wicky/lvalsecchi/rna_binder/analyses/foldseek_map")
N_RES = 80
# 0.1 A per unit. These are backbone models whose real accuracy is on the order
# of an Angstrom, so 0.1 A is already finer than the data justifies, and the
# coarser grid compresses far better than 0.01 A did.
SCALE = 10.0


def ca_from_pdb(text: str) -> np.ndarray:
    out = []
    for line in text.splitlines():
        if line.startswith("ATOM") and line[12:16].strip() == "CA":
            out.append((float(line[30:38]), float(line[38:46]), float(line[46:54])))
    return np.asarray(out, dtype=np.float64)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--layer", required=True, choices=["design", "af3", "esm"])
    ap.add_argument("--results", type=Path, default=FM / "results")
    args = ap.parse_args()

    out = args.results / args.layer
    nodes = pd.read_csv(out / "dedup" / "nodes.tsv", sep="\t")
    want = {name: i for i, name in enumerate(nodes["struct_name"])}
    n = len(nodes)
    print(f"{args.layer}: {n:,} backbones to extract")

    coords = np.zeros((n, N_RES, 3), dtype=np.int16)
    centres = np.zeros((n, 3), dtype=np.float32)
    found = np.zeros(n, dtype=bool)

    for shard in sorted((out / "shards").glob("shard_*.tar.gz")):
        with tarfile.open(shard, "r:gz") as tar:
            for member in tar:
                name = member.name.removesuffix(".pdb")
                idx = want.get(name)
                if idx is None or found[idx]:
                    continue
                fh = tar.extractfile(member)
                if fh is None:
                    continue
                ca = ca_from_pdb(fh.read().decode())
                if len(ca) != N_RES:
                    continue
                centre = ca.mean(axis=0)
                q = np.clip(np.round((ca - centre) * SCALE), -32767, 32767)
                # Delta-encode along the chain: consecutive CA atoms sit ~3.8 A
                # apart, so deltas occupy a tiny fraction of the range absolute
                # coordinates need and gzip finds far more structure in them.
                q[1:] -= q[:-1].copy()
                coords[idx] = q.astype(np.int16)
                centres[idx] = centre.astype(np.float32)
                found[idx] = True

    print(f"  extracted {int(found.sum()):,} / {n:,}")
    missing = int((~found).sum())
    if missing:
        print(f"  WARNING: {missing:,} backbones had no usable CA trace")

    dest = out / "traces"
    dest.mkdir(parents=True, exist_ok=True)
    blob = coords.tobytes()
    packed = gzip.compress(blob, 6)
    (dest / "ca_int16.bin.gz").write_bytes(packed)
    np.save(dest / "found.npy", found)
    print(f"  wrote {len(packed) / 1e6:.1f} MB gzipped "
          f"(from {len(blob) / 1e6:.1f} MB raw) -> {dest}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
