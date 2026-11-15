### Task 9: UMAP embedding from a precomputed kNN

**Files:**
- Create: `src/binder_map/embedding.py`
- Create: `tests/test_embedding.py`

**Interfaces:**
- Consumes: `binder_map.SEED`.
- Produces: `embed(n: int, layout_edges: np.ndarray, layout_weights: np.ndarray, k: int = 15, seed: int = SEED) -> np.ndarray` — `(n, 2)` float32 coordinates.

- [ ] **Step 1: Write the failing test**

Create `tests/test_embedding.py`:

```python
import numpy as np

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


def test_embed_keeps_distinct_groups_apart():
    """Two well-separated groups should not be interleaved in the layout."""
    v = _two_blobs()
    res = knn_from_vectors(v, k=5, floor=0.5)
    xy = embed(len(v), res.layout_edges, res.layout_weights, k=5)
    half = len(v) // 2
    ca, cb = xy[:half].mean(axis=0), xy[half:].mean(axis=0)
    spread = max(xy[:half].std(axis=0).mean(), xy[half:].std(axis=0).mean())
    assert np.linalg.norm(ca - cb) > spread
```

- [ ] **Step 2: Run it to confirm it fails**

```bash
cd /Users/leo/Documents/binder_map
/Users/leo/.venvs/binder_map/bin/python -m pytest tests/test_embedding.py -v
```

Expected: FAIL — `ModuleNotFoundError: No module named 'binder_map.embedding'`

- [ ] **Step 3: Implement `src/binder_map/embedding.py`**

```python
"""2D layout via UMAP over the precomputed layout kNN.

Screen distance is illustrative only. Cluster membership and similarity are
defined by the graph and the exact scores, never inferred from coordinates.
"""
from __future__ import annotations

import numpy as np
import umap

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


def embed(n: int, layout_edges, layout_weights, k: int = 15,
          seed: int = SEED) -> np.ndarray:
    """Return (n, 2) float32 UMAP coordinates."""
    idx, dist = _to_dense_knn(n, layout_edges, layout_weights, k)
    reducer = umap.UMAP(
        n_neighbors=k,
        min_dist=0.1,
        n_components=2,
        metric="precomputed",
        precomputed_knn=(idx, dist),
        random_state=seed,
        verbose=False,
    )
    # UMAP still wants an X of the right length when a precomputed knn is given.
    xy = reducer.fit_transform(np.zeros((n, 1), dtype=np.float32))
    return np.asarray(xy, dtype=np.float32)
```

- [ ] **Step 4: Run the tests**

```bash
cd /Users/leo/Documents/binder_map
/Users/leo/.venvs/binder_map/bin/python -m pytest tests/test_embedding.py -v
```

Expected: 3 passed.

**API pre-verified 2026-08-23** against the installed `umap-learn` 0.5.12: `UMAP.__init__` accepts `precomputed_knn`, and `fit_transform(np.zeros((n, 1)))` works both with `metric="precomputed"` and with the default metric. So the call in Step 3 is known-good — if it fails, suspect the shape or dtype of the `(indices, distances)` pair, not the API. Never fall back to a dense distance matrix: at 72k nodes that is ~20 GB and will not fit.

- [ ] **Step 5: Commit**

```bash
cd /Users/leo/Documents/binder_map
git add src/binder_map/embedding.py tests/test_embedding.py
git commit -m "Add UMAP embedding over the precomputed layout kNN"
```

---

