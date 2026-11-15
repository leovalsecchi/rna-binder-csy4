#!/usr/bin/env python3
"""Stage 5: Leiden clusters and medoids, per axis.

Usage:  python scripts/stage5_cluster.py --layer esm --axis fold
"""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

from binder_map.clustering import RESOLUTIONS, cluster_all, medoids

FM = Path("/cluster/work/wicky/lvalsecchi/rna_binder/analyses/foldseek_map")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--layer", required=True, choices=["design", "af3", "esm"])
    parser.add_argument("--axis", required=True, choices=["fold", "epitope"])
    parser.add_argument("--results", type=Path, default=FM / "results")
    args = parser.parse_args()

    out = args.results / args.layer
    graph = out / "graph" / args.axis
    nodes = pd.read_parquet(graph / "nodes.parquet")
    edges = pd.read_parquet(graph / "edges_cluster.parquet")
    n = len(nodes)

    pairs = edges[["i", "j"]].to_numpy(dtype=np.int64)
    weights = edges["sim"].to_numpy(dtype=np.float64)
    print(f"{args.layer}/{args.axis}: {n:,d} nodes, {len(pairs):,d} edges")

    labels = cluster_all(n, pairs, weights)
    frame = pd.DataFrame({"node_id": nodes["node_id"].values})
    for name, arr in labels.items():
        frame[name] = arr
        sizes = np.bincount(arr)
        print(f"  {name}: {len(sizes):,d} clusters, largest {sizes.max():,d}, "
              f"singletons {int((sizes == 1).sum()):,d}")

    rows = []
    for name, arr in labels.items():
        chosen = medoids(arr, pairs, weights)
        sizes = np.bincount(arr)
        for cluster, node in sorted(chosen.items()):
            rows.append({
                "resolution": name,
                "cluster": int(cluster),
                "node_id": int(nodes["node_id"].values[node]),
                "struct_name": nodes["struct_name"].values[node],
                "size": int(sizes[cluster]),
            })

    dest = out / "clusters" / args.axis
    dest.mkdir(parents=True, exist_ok=True)
    frame.to_parquet(dest / "clusters.parquet", index=False)
    pd.DataFrame(rows).to_parquet(dest / "medoids.parquet", index=False)
    print(f"  wrote {len(rows):,d} medoids -> {dest}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
