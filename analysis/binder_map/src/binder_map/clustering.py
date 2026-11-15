"""Leiden community detection at three resolutions, plus medoid election."""
from __future__ import annotations

import igraph as ig
import leidenalg
import numpy as np

from binder_map import SEED

RESOLUTIONS = {"cluster_r05": 0.5, "cluster_r10": 1.0, "cluster_r20": 2.0,
               "cluster_r40": 4.0, "cluster_r80": 8.0}


def leiden(n: int, edges, weights, resolution: float, seed: int = SEED) -> np.ndarray:
    """Cluster labels for all n nodes. Isolates become singleton clusters."""
    edges = np.asarray(edges).reshape(-1, 2)
    # igraph SILENTLY grows the graph to max(n, max_edge_id + 1) rather than
    # raising, so an out-of-range edge would return more labels than there are
    # nodes, with no error at the fault site. These labels are joined against
    # length-n node tables and published, so a silent length mismatch would
    # misalign every downstream row. Fail here instead.
    if len(edges) and int(edges.max()) >= n:
        raise ValueError(
            f"edge references node {int(edges.max())} but n={n}"
        )
    graph = ig.Graph(n=n, edges=[tuple(int(x) for x in e) for e in edges])
    graph.es["weight"] = [float(w) for w in np.asarray(weights).ravel()]

    partition = leidenalg.find_partition(
        graph,
        leidenalg.RBConfigurationVertexPartition,
        weights="weight",
        resolution_parameter=resolution,
        seed=seed,
    )
    return np.asarray(partition.membership, dtype=np.int32)


def cluster_all(n: int, edges, weights) -> dict[str, np.ndarray]:
    """Run Leiden at every configured resolution."""
    return {
        name: leiden(n, edges, weights, resolution)
        for name, resolution in RESOLUTIONS.items()
    }


def medoids(labels, edges, weights) -> dict[int, int]:
    """Cluster label -> node index with the highest summed within-cluster similarity.

    Computed over the sparse graph, so this is a within-cluster degree
    weighted by similarity rather than a true metric medoid. At this scale a
    dense medoid is not affordable, and the representative it picks is the
    best-connected member either way.
    """
    labels = np.asarray(labels)
    edges = np.asarray(edges).reshape(-1, 2)
    weights = np.asarray(weights, dtype=np.float64).ravel()

    score = np.zeros(len(labels), dtype=np.float64)
    if len(edges):
        same = labels[edges[:, 0]] == labels[edges[:, 1]]
        np.add.at(score, edges[same][:, 0], weights[same])
        np.add.at(score, edges[same][:, 1], weights[same])

    out: dict[int, int] = {}
    order = np.lexsort((np.arange(len(labels)), -score))
    for node in order:
        label = int(labels[node])
        if label not in out:
            out[label] = int(node)
    return out
