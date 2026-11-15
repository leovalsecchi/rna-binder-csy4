#!/usr/bin/env python3
"""Stage 4: turn Foldseek hits and epitope fingerprints into graphs.

Both axes go through the same code. Each produces two graphs — a thresholded
clustering graph and an unthresholded layout graph — for the reasons in
binder_map.knn.

Usage:
    python scripts/stage4_graph.py --layer esm --tm-floor 0.5 --epi-floor 0.5
"""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

from binder_map import FINGERPRINT_DIM
from binder_map.fingerprint import FINGERPRINT_COLUMNS
from binder_map.knn import knn_from_pairs, knn_from_vectors, symmetrise_min

FM = Path("/cluster/work/wicky/lvalsecchi/rna_binder/analyses/foldseek_map")
# Held at float32 to match `sim`, which is float32 now that hits stream at that
# width. Under NumPy's NEP 50 promotion a bare Python float already adopts the
# array's dtype, so `sim >= 0.95` agrees with the old float64 path today — but a
# float64 threshold does NOT: float32(0.95) is 0.9499999880..., strictly below
# the float64 literal, so every hit scoring exactly 0.95 would drop out of the
# redundancy count. Pinning the width makes that independent of promotion rules.
REDUNDANCY_TM = np.float32(0.95)


HIT_COLUMNS = ["query", "target", "qtmscore", "ttmscore", "alntmscore",
               "lddt", "fident", "alnlen", "evalue", "bits"]
HIT_CHUNK = 5_000_000


def load_hits(hit_files, index):
    """Stream hit TSVs into compact (qi, tj, qtm, ttm) arrays.

    The previous implementation did `pd.concat([pd.read_csv(f) for f in files])`,
    which materialises every row with two Python-string columns. At 31,763 nodes
    that was 1.3 GB of TSV and fine; at ~508k nodes it is ~20 GB and ~150M rows,
    which exhausts the node before the first edge is built.

    Mapping names to node ids per chunk and keeping int32/float32 costs ~12 bytes
    per directed hit, so the same data sits in ~1.8 GB.
    """
    qs, ts, qtms, ttms = [], [], [], []
    unknown = 0
    for path in hit_files:
        if Path(path).stat().st_size == 0:
            continue
        reader = pd.read_csv(
            path, sep="\t", names=HIT_COLUMNS, header=None,
            usecols=["query", "target", "qtmscore", "ttmscore"],
            dtype={"qtmscore": "float32", "ttmscore": "float32"},
            chunksize=HIT_CHUNK,
        )
        for chunk in reader:
            qi = index.reindex(chunk["query"].str.removesuffix(".pdb")).to_numpy()
            tj = index.reindex(chunk["target"].str.removesuffix(".pdb")).to_numpy()
            known = ~(np.isnan(qi) | np.isnan(tj))
            unknown += int((~known).sum())
            qs.append(qi[known].astype(np.int32))
            ts.append(tj[known].astype(np.int32))
            qtms.append(chunk["qtmscore"].to_numpy()[known])
            ttms.append(chunk["ttmscore"].to_numpy()[known])
    if not qs:
        empty_i = np.empty(0, np.int32)
        empty_f = np.empty(0, np.float32)
        return empty_i, empty_i.copy(), empty_f, empty_f.copy(), unknown
    return (np.concatenate(qs), np.concatenate(ts),
            np.concatenate(qtms), np.concatenate(ttms), unknown)


def write_axis(out_dir: Path, axis: str, nodes: pd.DataFrame, res, redundancy):
    dest = out_dir / "graph" / axis
    dest.mkdir(parents=True, exist_ok=True)
    pd.DataFrame({
        "i": res.cluster_edges[:, 0] if len(res.cluster_edges) else np.array([], np.int32),
        "j": res.cluster_edges[:, 1] if len(res.cluster_edges) else np.array([], np.int32),
        "sim": res.cluster_weights,
    }).to_parquet(dest / "edges_cluster.parquet", index=False)
    pd.DataFrame({
        "i": res.layout_edges[:, 0] if len(res.layout_edges) else np.array([], np.int32),
        "j": res.layout_edges[:, 1] if len(res.layout_edges) else np.array([], np.int32),
        "sim": res.layout_weights,
    }).to_parquet(dest / "knn_layout.parquet", index=False)
    pd.DataFrame({
        "node_id": nodes["node_id"].values,
        "struct_name": nodes["struct_name"].values,
        "degree": res.degree,
        "redundancy": redundancy,
        "is_isolate": res.is_isolate,
    }).to_parquet(dest / "nodes.parquet", index=False)
    print(f"  {axis}: {len(res.cluster_edges):,d} cluster edges, "
          f"{len(res.layout_edges):,d} layout edges, "
          f"{int(res.is_isolate.sum()):,d} isolates")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--layer", required=True, choices=["design", "af3", "esm"])
    parser.add_argument("--tm-floor", type=float, default=0.5)
    parser.add_argument("--epi-floor", type=float, default=0.5)
    parser.add_argument("--k", type=int, default=15)
    parser.add_argument("--mutual", action="store_true",
                        help="keep an edge only if both endpoints rank the "
                             "other in their own top-k")
    parser.add_argument("--results", type=Path, default=FM / "results")
    args = parser.parse_args()

    out = args.results / args.layer
    nodes = pd.read_csv(out / "dedup" / "nodes.tsv", sep="\t")
    index = pd.Series(nodes["node_id"].values, index=nodes["struct_name"].values)
    n = len(nodes)
    print(f"{args.layer}: {n:,d} nodes")

    # ---------- fold axis ----------
    hit_files = sorted((out / "hits").glob("hits_*.tsv"))
    if not hit_files:
        print("FATAL: no hit files; did stage 3 run?")
        return 1
    qi, tj, qtm, ttm, unknown = load_hits(hit_files, index)
    print(f"  {len(qi):,d} usable hits from {len(hit_files)} chunks")
    if unknown:
        print(f"  WARNING: {unknown:,d} hits reference unknown names")

    i, j, sim = symmetrise_min(qi, tj, qtm, ttm)
    print(f"  {len(i):,d} undirected pairs after min(qTM,tTM)")

    # Redundancy from the FULL hit list, before top-k truncation, else it
    # saturates at k. Still bounded above by --max-seqs (300).
    redundancy = np.zeros(n, dtype=np.int32)
    strong = sim >= REDUNDANCY_TM
    if strong.any():
        np.add.at(redundancy, i[strong], 1)
        np.add.at(redundancy, j[strong], 1)

    fold = knn_from_pairs(n, i, j, sim, k=args.k, floor=args.tm_floor,
                          mutual=args.mutual)
    if args.mutual:
        print(f"  fold: mutual-kNN filter on (edges must be reciprocal)")
    write_axis(out, "fold", nodes, fold, redundancy)

    # ---------- epitope axis ----------
    fp = pd.concat([pd.read_parquet(f) for f in sorted((out / "fp").glob("shard_*.parquet"))],
                   ignore_index=True)
    fp = fp.drop_duplicates("design_id").set_index("design_id")
    canonical = fp.reindex(nodes["struct_name"].values)
    # Restore design_id as a COLUMN before writing. It is the index here, and
    # `index=False` would drop it silently -- which it did, leaving stage 6b
    # unable to key fingerprints by design and quietly skipping the
    # cross-layer epitope agreement entirely.
    canonical.insert(0, "design_id", nodes["struct_name"].values)
    canonical.insert(0, "node_id", nodes["node_id"].values)
    canonical.to_parquet(out / "fp" / "canonical.parquet", index=False)

    vectors = canonical[FINGERPRINT_COLUMNS].to_numpy(dtype=np.float64)
    assert vectors.shape == (n, FINGERPRINT_DIM)

    # Zero-contact designs have no epitope. Excluding them is deliberate:
    # imputing one would invent the signal being measured.
    has_interface = canonical["contacts_total"].to_numpy() > 0
    print(f"  epitope: {int((~has_interface).sum()):,d} nodes with no interface, excluded")

    sub_index = np.flatnonzero(has_interface)
    epi_sub = knn_from_vectors(vectors[sub_index], k=args.k, floor=args.epi_floor,
                               mutual=args.mutual)

    # Lift the sub-graph node ids back into the full node space.
    def lift(edges):
        if len(edges) == 0:
            return np.empty((0, 2), np.int32)
        return sub_index[edges].astype(np.int32)

    from binder_map.knn import KNNResult
    degree = np.zeros(n, dtype=np.int32)
    degree[sub_index] = epi_sub.degree
    is_isolate = np.ones(n, dtype=bool)
    is_isolate[sub_index] = epi_sub.is_isolate
    epi = KNNResult(
        lift(epi_sub.cluster_edges), epi_sub.cluster_weights,
        lift(epi_sub.layout_edges), epi_sub.layout_weights,
        degree, is_isolate,
    )
    epi_redundancy = np.zeros(n, dtype=np.int32)
    write_axis(out, "epitope", nodes, epi, epi_redundancy)

    # Record which nodes carry an interface, for the explorer's greying-out.
    pd.DataFrame({
        "node_id": nodes["node_id"].values,
        "has_interface": has_interface,
        "contacts_total": canonical["contacts_total"].values,
    }).to_parquet(out / "graph" / "epitope" / "interface.parquet", index=False)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
