import numpy as np
import pytest

from binder_map.knn import knn_from_pairs, knn_from_vectors, symmetrise_min


def test_symmetrise_takes_min_of_the_two_normalisations():
    qi = np.array([0, 1])
    qj = np.array([1, 0])
    qtm = np.array([0.9, 0.4])
    ttm = np.array([0.5, 0.8])
    i, j, sim = symmetrise_min(qi, qj, qtm, ttm)
    assert len(i) == 1
    # direction 0->1 gives min(0.9,0.5)=0.5; 1->0 gives min(0.4,0.8)=0.4
    # the pair keeps the max of the observed mins
    assert sim[0] == pytest.approx(0.5)


def test_symmetrise_drops_self_hits():
    i, j, sim = symmetrise_min(
        np.array([0, 0]), np.array([0, 1]),
        np.array([1.0, 0.6]), np.array([1.0, 0.6]),
    )
    assert len(i) == 1
    assert (i[0], j[0]) == (0, 1)


def test_floor_excludes_weak_edges_from_clustering_but_not_layout():
    # 3 nodes: 0-1 strong, 0-2 weak
    i = np.array([0, 0])
    j = np.array([1, 2])
    sim = np.array([0.9, 0.1])
    res = knn_from_pairs(3, i, j, sim, k=2, floor=0.5)
    cluster_pairs = {tuple(sorted(e)) for e in res.cluster_edges}
    layout_pairs = {tuple(sorted(e)) for e in res.layout_edges}
    assert (0, 1) in cluster_pairs
    assert (0, 2) not in cluster_pairs   # below floor
    assert (0, 2) in layout_pairs        # layout keeps unthresholded top-k


def test_node_with_no_edge_above_floor_is_an_isolate():
    i = np.array([0])
    j = np.array([1])
    sim = np.array([0.2])
    res = knn_from_pairs(3, i, j, sim, k=2, floor=0.5)
    assert res.is_isolate[0]
    assert res.is_isolate[1]
    assert res.is_isolate[2]   # node 2 has no hits at all
    assert len(res.cluster_edges) == 0


def test_knn_from_vectors_finds_identical_vectors():
    v = np.array([
        [1.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],   # identical to 0
        [0.0, 1.0, 0.0],   # orthogonal
    ])
    res = knn_from_vectors(v, k=1, floor=0.5)
    pairs = {tuple(sorted(e)) for e in res.cluster_edges}
    assert (0, 1) in pairs
    assert (0, 2) not in pairs


def test_knn_from_vectors_is_scale_invariant():
    """L1 normalisation: a promiscuous binder and a sparse one that touch the
    same place must be neighbours."""
    v = np.array([
        [2.0, 1.0, 0.0],
        [20.0, 10.0, 0.0],   # same profile, 10x the contacts
        [0.0, 0.0, 5.0],
    ])
    res = knn_from_vectors(v, k=1, floor=0.9)
    pairs = {tuple(sorted(e)) for e in res.cluster_edges}
    assert (0, 1) in pairs


def test_knn_from_vectors_blocking_matches_unblocked():
    rng = np.random.default_rng(0)
    v = rng.random((200, 12))
    a = knn_from_vectors(v, k=5, floor=0.0, block=1000)
    b = knn_from_vectors(v, k=5, floor=0.0, block=32)
    assert {tuple(sorted(e)) for e in a.cluster_edges} == \
           {tuple(sorted(e)) for e in b.cluster_edges}


def test_layout_graph_gives_every_node_k_neighbours():
    rng = np.random.default_rng(1)
    v = rng.random((50, 8))
    res = knn_from_vectors(v, k=5, floor=0.99)
    counts = np.bincount(res.layout_edges.ravel(), minlength=50)
    assert (counts >= 5).all()
