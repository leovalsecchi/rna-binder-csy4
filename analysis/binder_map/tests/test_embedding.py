import numpy as np
import pytest

from binder_map.embedding import embed
from binder_map.knn import knn_from_vectors


def _two_blobs(seed=0, n=60):
    rng = np.random.default_rng(seed)
    a = rng.normal(0.0, 0.05, size=(n // 2, 6)) + np.array([1, 0, 0, 0, 0, 0])
    b = rng.normal(0.0, 0.05, size=(n // 2, 6)) + np.array([0, 0, 0, 0, 0, 1])
    return np.abs(np.vstack([a, b]))


def test_embed_returns_2d_coordinates():
    v = _two_blobs()
    res = knn_from_vectors(v, k=5, floor=0.5)
    xy = embed(len(v), res.layout_edges, res.layout_weights, k=5)
    assert xy.shape == (len(v), 2)
    assert np.isfinite(xy).all()


def test_embed_is_deterministic():
    v = _two_blobs()
    res = knn_from_vectors(v, k=5, floor=0.5)
    a = embed(len(v), res.layout_edges, res.layout_weights, k=5)
    b = embed(len(v), res.layout_edges, res.layout_weights, k=5)
    np.testing.assert_allclose(a, b)


def test_node_with_no_neighbours_gets_finite_coordinates():
    """A fully disconnected node must land somewhere, not become NaN.

    UMAP returns NaN for such a row. Reachable in production because
    symmetrise_min drops self-hits, so a structure whose only search hit is
    itself arrives with zero edges. A NaN would silently drop the point from
    the published map, and an outlier is exactly what an exploratory map most
    needs to show.
    """
    # nodes 0..3 form a connected block; node 4 has no edges at all.
    edges = np.array([[0, 1], [1, 2], [2, 3], [0, 3]])
    weights = np.ones(len(edges), dtype=np.float32)
    xy = embed(5, edges, weights, k=2)
    assert xy.shape == (5, 2)
    assert np.isfinite(xy).all(), "disconnected node produced NaN"
    # and it must sit outside the connected cloud, not on top of it
    centre = xy[:4].mean(axis=0)
    spread = float(np.linalg.norm(xy[:4] - centre, axis=1).max())
    assert float(np.linalg.norm(xy[4] - centre)) >= spread


@pytest.mark.parametrize("n_lonely", [1, 3, 5, 12])
def test_many_disconnected_nodes_do_not_crash_umap(n_lonely):
    """UMAP's spectral init raises IndexError at >=3 disconnected components.

    An earlier fix let those nodes reach UMAP and cleaned up afterwards, which
    could never work: the crash happens inside the library, before any
    cleanup. Testing only ONE lonely node hid it completely, so this is
    parametrized across the threshold.
    """
    v = _two_blobs(n=40)
    res = knn_from_vectors(v, k=5, floor=0.5)
    n = len(v) + n_lonely
    xy = embed(n, res.layout_edges, res.layout_weights, k=5)
    assert xy.shape == (n, 2)
    assert np.isfinite(xy).all()
    # every lonely node distinct: they must not stack on one point
    ring = xy[len(v):]
    assert len({tuple(np.round(p, 4)) for p in ring}) == n_lonely


def test_embed_with_no_edges_at_all_is_finite():
    """A completely edgeless graph must still yield usable coordinates."""
    xy = embed(6, np.empty((0, 2), dtype=np.int64),
               np.empty(0, dtype=np.float32), k=5)
    assert xy.shape == (6, 2)
    assert np.isfinite(xy).all()
    assert len({tuple(np.round(p, 4)) for p in xy}) == 6


def test_embed_keeps_distinct_groups_apart():
    """Two well-separated groups should not be interleaved in the layout."""
    v = _two_blobs()
    res = knn_from_vectors(v, k=5, floor=0.5)
    xy = embed(len(v), res.layout_edges, res.layout_weights, k=5)
    half = len(v) // 2
    ca, cb = xy[:half].mean(axis=0), xy[half:].mean(axis=0)
    spread = max(xy[:half].std(axis=0).mean(), xy[half:].std(axis=0).mean())
    assert np.linalg.norm(ca - cb) > spread


def test_multi_component_graph_does_not_crash():
    """UMAP's spectral init raises IndexError on multi-component graphs.

    Distinct from the zero-edge case: every node here HAS edges, they just
    form several islands with no path between them. Observed in production at
    71,863 nodes while 36,907 survived, so the failure is size-dependent and
    the init must be chosen from the graph's actual component count.
    """
    # five well-separated islands of four nodes each
    edges, weights = [], []
    for island in range(5):
        base = island * 4
        for a in range(4):
            for b in range(a + 1, 4):
                edges.append([base + a, base + b])
                weights.append(1.0)
    edges = np.array(edges)
    weights = np.array(weights, dtype=np.float32)

    xy = embed(20, edges, weights, k=3)
    assert xy.shape == (20, 2)
    assert np.isfinite(xy).all()
