#!/usr/bin/env python3
"""Merge stage7c shard partials into the single blobs the page embeds."""
from __future__ import annotations

import argparse
import gzip
from pathlib import Path

import numpy as np
import pandas as pd

FM = Path("/cluster/work/wicky/lvalsecchi/rna_binder/analyses/foldseek_map")


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--layer", required=True)
    ap.add_argument("--results", type=Path, default=FM / "results")
    args = ap.parse_args()

    out = args.results / args.layer / "complex"
    parts = sorted((out / "parts").glob("shard_*.npz"))
    if not parts:
        print("no partials found")
        return 1
    n = len(pd.read_csv(args.results / args.layer / "dedup" / "nodes.tsv", sep="\t"))

    coords = sse = None
    found = np.zeros(n, dtype=bool)
    for f in parts:
        z = np.load(f)
        if coords is None:
            coords = np.zeros_like(z["coords"])
            sse = np.full_like(z["sse"], ord("C"))
        m = z["found"]
        coords[m] = z["coords"][m]
        sse[m] = z["sse"][m]
        found |= m
    print(f"{args.layer}: merged {len(parts)} partials, "
          f"{int(found.sum()):,}/{n:,} complexes")

    cpacked = gzip.compress(coords.tobytes(), 6)
    spacked = gzip.compress(sse.tobytes(), 6)
    (out / "complex_int16.bin.gz").write_bytes(cpacked)
    (out / "sse.bin.gz").write_bytes(spacked)
    print(f"  coords {len(cpacked)/1e6:.1f} MB gz (raw {coords.nbytes/1e6:.1f} MB)")
    print(f"  sse    {len(spacked)/1e6:.2f} MB gz")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
