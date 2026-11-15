"""kNN graphs for both axes.

Two graphs come out of every call, and they are deliberately different:

  cluster graph  thresholded at `floor`, then top-k. A node with no neighbour
                 above the floor is a genuine isolate and becomes its own
                 singleton cluster. Absence of evidence stays absence of
                 evidence.
  layout graph   raw top-k with no floor, so every node has k neighbours.
                 UMAP's precomputed_knn needs a full k per node, and padding
                 the thresholded graph would fabricate exactly the edges the
                 floor exists to reject.

A point can therefore sit inside a visual cluster while being a singleton by
cluster id. That is informative, not a bug.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np


@dataclass(frozen=True)
class KNNResult:
    cluster_edges: np.ndarray    # (m, 2) int32, i < j
    cluster_weights: np.ndarray  # (m,) float32
    layout_edges: np.ndarray     # (p, 2) int32
    layout_weights: np.ndarray   # (p,) float32
    degree: np.ndarray           # (n,) int32, in the cluster graph
    is_isolate: np.ndarray       # (n,) bool


def symmetrise_min(qi, qj, qtm, ttm):
    """Collapse directed Foldseek hits into undirected min(qTM, tTM) pairs.

    Conservative: a pair only scores high if it is high in both
    normalisations. Where both directions were observed, the pair keeps the
    max of the two observed minima.
    """
    qi = np.asarray(qi)
    qj = np.asarray(qj)
    sim = np.minimum(np.asarray(qtm, dtype=np.float64),
                     np.asarray(ttm, dtype=np.float64))

    keep = qi != qj
    a = np.minimum(qi[keep], qj[keep])
    b = np.maximum(qi[keep], qj[keep])
    sim = sim[keep]

    if len(a) == 0:
        return (np.empty(0, np.int32), np.empty(0, np.int32),
                np.empty(0, np.float32))

    order = np.lexsort((b, a))
    a, b, sim = a[order], b[order], sim[order]
    new = np.ones(len(a), dtype=bool)
    new[1:] = (a[1:] != a[:-1]) | (b[1:] != b[:-1])
    group = np.cumsum(new) - 1
    best = np.zeros(group[-1] + 1, dtype=np.float64)
    np.maximum.at(best, group, sim)
    return (a[new].astype(np.int32), b[new].astype(np.int32),
            best.astype(np.float32))


def _topk_per_node(n, i, j, sim, k):
    """Keep the k strongest incident edges per node, as an undirected set."""
    if len(i) == 0:
        return np.empty((0, 2), np.int32), np.empty(0, np.float32)
    src = np.concatenate([i, j])
    dst = np.concatenate([j, i])
    val = np.concatenate([sim, sim])

    # Dedup identical directed (src, dst) pairs before ranking. Without this,
    # a *mutual* top-k pick (A has B in A's own top-k, and B independently
    # has A in B's own top-k) produces two rows for the same neighbour in
    # A's incident run: A's own forward entry, plus B's mirrored entry. That
    # duplicate silently eats a second rank slot, evicting a genuinely
    # distinct, lower-similarity neighbour from the top-k window entirely
    # (fewer than k neighbours survive for A, and a real edge can be dropped
    # from the cluster graph even though it is above the floor).
    order = np.lexsort((-val, dst, src))
    src, dst, val = src[order], dst[order], val[order]
    first = np.ones(len(src), dtype=bool)
    first[1:] = (src[1:] != src[:-1]) | (dst[1:] != dst[:-1])
    src, dst, val = src[first], dst[first], val[first]

    order = np.lexsort((-val, src))
    src, dst, val = src[order], dst[order], val[order]

    # Rank each edge within its source node's run: the arrays are already
    # sorted by (src, -val), so position-within-run IS the similarity rank.
    starts = np.flatnonzero(np.r_[True, src[1:] != src[:-1]])
    run_lengths = np.diff(np.r_[starts, len(src)])
    rank = np.arange(len(src)) - np.repeat(starts, run_lengths)
    keep = rank < k

    a = np.minimum(src[keep], dst[keep])
    b = np.maximum(src[keep], dst[keep])
    v = val[keep]
    order = np.lexsort((b, a))
    a, b, v = a[order], b[order], v[order]
    new = np.ones(len(a), dtype=bool)
    new[1:] = (a[1:] != a[:-1]) | (b[1:] != b[:-1])
    return (np.stack([a[new], b[new]], axis=1).astype(np.int32),
            v[new].astype(np.float32))


def _mutual_filter(n, edges, weights, i, j, sim, k):
    """Keep only edges where BOTH endpoints rank the other in their own top-k.

    A one-directional edge means A considers B one of its closest neighbours
    while B does not reciprocate - typically a hub attaching a loose structure
    to a dense family. Requiring reciprocity trades a little recall for
    markedly tighter clusters.
    """
    if len(edges) == 0:
        return edges, weights
    src = np.concatenate([i, j])
    dst = np.concatenate([j, i])
    val = np.concatenate([sim, sim])
    order = np.lexsort((-val, dst, src))
    src, dst, val = src[order], dst[order], val[order]
    first = np.ones(len(src), dtype=bool)
    first[1:] = (src[1:] != src[:-1]) | (dst[1:] != dst[:-1])
    src, dst = src[first], dst[first]
    val = val[first]
    order = np.lexsort((-val, src))
    src, dst = src[order], dst[order]
    starts = np.flatnonzero(np.r_[True, src[1:] != src[:-1]])
    runs = np.diff(np.r_[starts, len(src)])
    rank = np.arange(len(src)) - np.repeat(starts, runs)
    keep = rank < k
    directed = set(zip(src[keep].tolist(), dst[keep].tolist()))

    mask = np.array([(int(a), int(b)) in directed and (int(b), int(a)) in directed
                     for a, b in edges])
    return edges[mask], weights[mask]


def _assemble(n, i, j, sim, k, floor, mutual: bool = False) -> KNNResult:
    layout_edges, layout_weights = _topk_per_node(n, i, j, sim, k)

    above = sim >= floor
    cluster_edges, cluster_weights = _topk_per_node(
        n, i[above], j[above], sim[above], k
    )
    if mutual:
        cluster_edges, cluster_weights = _mutual_filter(
            n, cluster_edges, cluster_weights, i[above], j[above], sim[above], k
        )

    degree = np.zeros(n, dtype=np.int32)
    if len(cluster_edges):
        np.add.at(degree, cluster_edges.ravel(), 1)
    return KNNResult(cluster_edges, cluster_weights,
                     layout_edges, layout_weights,
                     degree, degree == 0)


def knn_from_pairs(n, i, j, sim, k, floor, mutual: bool = False) -> KNNResult:
    """Build both graphs from an existing sparse similarity list."""
    return _assemble(n, np.asarray(i), np.asarray(j),
                     np.asarray(sim, dtype=np.float32), k, floor, mutual)


def knn_from_vectors(vectors, k, floor, block: int = 1024,
                     mutual: bool = False) -> KNNResult:
    """Exact cosine kNN over L1-normalised rows, by blocked matrix product.

    `block` rows are scored against all n at once, so the temporary is
    block x n float32. At n = 508k a block of 4096 allocates 8.3 GB, which is
    why the default is 1024 (2 GB). The product stays exact; only the blocking
    changes.
    """
    v = np.asarray(vectors, dtype=np.float64)
    n = len(v)

    l1 = np.abs(v).sum(axis=1, keepdims=True)
    l1[l1 == 0] = 1.0
    v = v / l1
    norm = np.linalg.norm(v, axis=1, keepdims=True)
    norm[norm == 0] = 1.0
    unit = (v / norm).astype(np.float32)

    kk = min(k, max(n - 1, 1))
    all_i, all_j, all_s = [], [], []
    for start in range(0, n, block):
        stop = min(start + block, n)
        sims = unit[start:stop] @ unit.T
        rows = np.arange(start, stop)
        sims[np.arange(stop - start), rows] = -np.inf   # drop self
        idx = np.argpartition(-sims, kth=kk - 1, axis=1)[:, :kk]
        take = np.take_along_axis(sims, idx, axis=1)
        all_i.append(np.repeat(rows, kk))
        all_j.append(idx.ravel())
        all_s.append(take.ravel())

    i = np.concatenate(all_i)
    j = np.concatenate(all_j)
    s = np.concatenate(all_s).astype(np.float32)
    good = np.isfinite(s)
    return _assemble(n, i[good], j[good], s[good], k, floor, mutual)
