"""2D layout via UMAP over the precomputed layout kNN.

Screen distance is illustrative only. Cluster membership and similarity are
defined by the graph and the exact scores, never inferred from coordinates.
"""
from __future__ import annotations

import numpy as np
import umap
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components

from binder_map import SEED


def _to_dense_knn(n: int, edges, weights, k: int):
    """Build the (indices, distances) pair UMAP's precomputed_knn expects.

    Distance is 1 - similarity. Nodes with fewer than k observed neighbours
    are padded with self-references at maximum distance, which UMAP treats as
    non-informative rather than as a fabricated edge.
    """
    edges = np.asarray(edges).reshape(-1, 2)
    weights = np.asarray(weights, dtype=np.float32).ravel()

    neighbours: list[list[tuple[float, int]]] = [[] for _ in range(n)]
    for (a, b), w in zip(edges, weights):
        neighbours[int(a)].append((float(w), int(b)))
        neighbours[int(b)].append((float(w), int(a)))

    idx = np.zeros((n, k), dtype=np.int32)
    dist = np.ones((n, k), dtype=np.float32)
    for node, nbrs in enumerate(neighbours):
        nbrs.sort(key=lambda t: -t[0])
        seen: set[int] = set()
        row = []
        for w, other in nbrs:
            if other in seen or other == node:
                continue
            seen.add(other)
            row.append((other, 1.0 - w))
            if len(row) == k:
                break
        while len(row) < k:
            row.append((node, 1.0))
        idx[node] = [r[0] for r in row]
        dist[node] = [r[1] for r in row]
    return idx, dist


def _ring_positions(count: int, centre, radius: float) -> np.ndarray:
    """`count` points evenly spaced on a circle, in deterministic order."""
    angles = np.linspace(0.0, 2.0 * np.pi, count, endpoint=False)
    return np.stack([
        centre[0] + radius * np.cos(angles),
        centre[1] + radius * np.sin(angles),
    ], axis=1).astype(np.float32)


# UMAP needs at least this many points before an embedding is meaningful.
MIN_UMAP_NODES = 3


def embed(n: int, layout_edges, layout_weights, k: int = 15,
          seed: int = SEED) -> np.ndarray:
    """Return (n, 2) float32 coordinates, finite for every node.

    Nodes with NO real neighbour are removed BEFORE UMAP runs, then placed on
    a ring outside the embedded cloud.

    Removing them first is not tidiness, it is the fix for a hard crash.
    UMAP's default spectral init raises IndexError from inside
    `multi_component_layout` once three or more nodes are mutually
    disconnected — so passing them through and cleaning up afterwards fails
    before the cleanup can run. Switching to `init="random"` would also avoid
    the crash but degrades the layout for the connected majority, which is
    the part anyone actually reads. Excluding them keeps spectral init for
    the real data and makes the degenerate case structurally impossible.

    Such nodes are reachable: `symmetrise_min` discards self-hits, so a
    structure whose only search hit is itself arrives with zero edges. On an
    exploratory map that design is an outlier - the most interesting thing a
    point can be - so it is shown, not dropped.
    """
    edges = np.asarray(layout_edges).reshape(-1, 2)
    weights = np.asarray(layout_weights, dtype=np.float32).ravel()

    connected = np.zeros(n, dtype=bool)
    if len(edges):
        connected[edges.ravel()] = True
    keep = np.flatnonzero(connected)

    xy = np.zeros((n, 2), dtype=np.float32)

    if len(keep) >= MIN_UMAP_NODES:
        # Compact the connected subgraph to 0..len(keep)-1. Every edge is
        # between connected nodes by construction, so no index can be -1.
        remap = np.full(n, -1, dtype=np.int64)
        remap[keep] = np.arange(len(keep))
        sub_edges = remap[edges]
        assert sub_edges.min() >= 0

        k_eff = max(1, min(k, len(keep) - 1))
        idx, dist = _to_dense_knn(len(keep), sub_edges, weights, k_eff)

        # UMAP's spectral init crashes with an IndexError inside
        # multi_component_layout when the graph splits into several
        # disconnected components. Dropping zero-edge nodes is NOT enough:
        # a component can be several nodes joined to each other but not to
        # the main body, and every one of those nodes has edges. Observed at
        # 71,863 nodes (af3/esm) while 36,907 (design) survived, so this is
        # size-dependent and must be decided from the graph, not assumed.
        #
        # Spectral gives a better global layout, so keep it when the graph is
        # connected and fall back only when it is not.
        n_sub = len(keep)
        adj = coo_matrix(
            (np.ones(len(sub_edges), dtype=np.int8),
             (sub_edges[:, 0], sub_edges[:, 1])),
            shape=(n_sub, n_sub),
        )
        n_components = connected_components(adj, directed=False,
                                            return_labels=False)
        init = "spectral" if n_components == 1 else "random"
        if n_components > 1:
            print(f"  layout graph has {n_components} components -> init=random "
                  f"(spectral init crashes on multi-component graphs)")

        reducer = umap.UMAP(
            n_neighbors=k_eff,
            min_dist=0.1,
            n_components=2,
            metric="precomputed",
            precomputed_knn=(idx, dist),
            init=init,
            random_state=seed,
            verbose=False,
        )
        # UMAP still wants an X of the right length when a precomputed knn is given.
        sub = reducer.fit_transform(np.zeros((len(keep), 1), dtype=np.float32))
        xy[keep] = np.asarray(sub, dtype=np.float32)
    else:
        # Too few connected nodes to embed meaningfully: everything goes on
        # the ring rather than pretending a layout exists.
        keep = np.array([], dtype=np.int64)

    lonely = np.setdiff1d(np.arange(n), keep, assume_unique=False)
    if len(lonely):
        if len(keep):
            centre = xy[keep].mean(axis=0)
            # True Euclidean farthest point, not the axis-aligned extent: a
            # Chebyshev bound understates the radius by up to sqrt(2) for a
            # point sitting off-axis, which could put the ring inside the cloud.
            radius = float(np.linalg.norm(xy[keep] - centre, axis=1).max()) * 1.15
        else:
            centre, radius = np.zeros(2, dtype=np.float32), 1.0
        if not np.isfinite(radius) or radius <= 0.0:
            # Explicit check, not `or 1.0`: NaN is truthy, so `nan or 1.0`
            # would keep the NaN and silently defeat the fallback.
            radius = 1.0
        xy[lonely] = _ring_positions(len(lonely), centre, radius)

    if not np.isfinite(xy).all():
        raise ValueError("embedding produced non-finite coordinates")
    return xy
