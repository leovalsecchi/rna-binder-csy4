import numpy as np

from binder_map.clustering import RESOLUTIONS, cluster_all, leiden, medoids


def _two_triangles():
    """Two disjoint triangles: nodes 0,1,2 and 3,4,5."""
    edges = np.array([[0, 1], [0, 2], [1, 2], [3, 4], [3, 5], [4, 5]])
    weights = np.ones(len(edges), dtype=np.float32)
    return 6, edges, weights


def test_leiden_separates_disconnected_components():
    n, edges, weights = _two_triangles()
    labels = leiden(n, edges, weights, resolution=1.0)
    assert len(labels) == n
    assert len(set(labels[:3])) == 1
    assert len(set(labels[3:])) == 1
    assert labels[0] != labels[3]


def test_leiden_is_deterministic():
    n, edges, weights = _two_triangles()
    a = leiden(n, edges, weights, resolution=1.0)
    b = leiden(n, edges, weights, resolution=1.0)
    assert np.array_equal(a, b)


def test_isolates_get_their_own_cluster():
    edges = np.array([[0, 1]])
    weights = np.ones(1, dtype=np.float32)
    labels = leiden(3, edges, weights, resolution=1.0)
    assert labels[2] != labels[0]
    assert labels[2] != labels[1]


def test_cluster_all_returns_three_resolutions():
    n, edges, weights = _two_triangles()
    out = cluster_all(n, edges, weights)
    assert set(out) == set(RESOLUTIONS)
    for arr in out.values():
        assert arr.shape == (n,)


def test_out_of_range_edge_raises_instead_of_growing_the_graph():
    """igraph would silently grow the graph and return n+1 labels.

    These labels are joined against length-n node tables and published, so a
    silent length mismatch would misalign every downstream row.
    """
    import pytest

    with pytest.raises(ValueError, match="references node"):
        leiden(3, np.array([[0, 1], [1, 5]]), np.ones(2, dtype=np.float32),
               resolution=1.0)


def test_labels_always_have_length_n():
    n, edges, weights = _two_triangles()
    for arr in cluster_all(n, edges, weights).values():
        assert arr.shape == (n,)


def test_medoid_is_the_best_connected_member():
    # star: node 0 connected to 1, 2, 3; the others only to 0
    edges = np.array([[0, 1], [0, 2], [0, 3]])
    weights = np.array([1.0, 1.0, 1.0], dtype=np.float32)
    labels = np.zeros(4, dtype=np.int32)
    assert medoids(labels, edges, weights)[0] == 0


def test_medoid_exists_for_every_cluster():
    n, edges, weights = _two_triangles()
    labels = leiden(n, edges, weights, resolution=1.0)
    m = medoids(labels, edges, weights)
    assert set(m) == set(labels.tolist())
    for label, node in m.items():
        assert labels[node] == label
