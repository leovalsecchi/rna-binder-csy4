#!/usr/bin/env python3
"""Stage 6: UMAP layout per axis.

Usage:  python scripts/stage6_embed.py --layer esm --axis fold
"""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

from binder_map.embedding import embed

FM = Path("/cluster/work/wicky/lvalsecchi/rna_binder/analyses/foldseek_map")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--layer", required=True, choices=["design", "af3", "esm"])
    parser.add_argument("--axis", required=True, choices=["fold", "epitope"])
    parser.add_argument("--k", type=int, default=15)
    parser.add_argument("--results", type=Path, default=FM / "results")
    args = parser.parse_args()

    out = args.results / args.layer
    graph = out / "graph" / args.axis
    nodes = pd.read_parquet(graph / "nodes.parquet")
    layout = pd.read_parquet(graph / "knn_layout.parquet")
    n = len(nodes)
    print(f"{args.layer}/{args.axis}: embedding {n:,d} nodes "
          f"from {len(layout):,d} layout edges")

    xy = embed(n, layout[["i", "j"]].to_numpy(dtype=np.int64),
               layout["sim"].to_numpy(dtype=np.float32), k=args.k)

    dest = out / "embed" / args.axis
    dest.mkdir(parents=True, exist_ok=True)
    pd.DataFrame({
        "node_id": nodes["node_id"].values,
        "x": xy[:, 0], "y": xy[:, 1],
    }).to_parquet(dest / "embed.parquet", index=False)
    print(f"  wrote coordinates -> {dest}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
