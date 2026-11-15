# Task 8 report: Leiden clustering and medoid election

## What was implemented

Created two new files, per the brief, verbatim (no deviation):

- `src/binder_map/clustering.py`
  - `RESOLUTIONS: dict[str, float]` = `{"cluster_r05": 0.5, "cluster_r10": 1.0, "cluster_r20": 2.0}`
  - `leiden(n, edges, weights, resolution, seed=SEED) -> np.ndarray` — builds an
    `igraph.Graph` on `n` vertices from the edge list, attaches edge weights,
    and runs `leidenalg.find_partition` with `RBConfigurationVertexPartition`
    at the given resolution and seed. Returns `(n,)` int32 labels. Vertices
    with no incident edges are left untouched by igraph/leidenalg and end up
    each in their own singleton community, so isolates get their own label
    for free — no special-casing needed.
  - `cluster_all(n, edges, weights) -> dict[str, np.ndarray]` — runs `leiden`
    once per entry in `RESOLUTIONS`, keyed by resolution name.
  - `medoids(labels, edges, weights) -> dict[int, int]` — for each cluster,
    picks the node with the highest summed edge weight to same-cluster
    neighbours (a similarity-weighted within-cluster degree over the sparse
    graph, not a true metric medoid — this is intentional per the brief's
    note, and the docstring says so explicitly). Ties are broken
    deterministically by lowest node index via `np.lexsort`.

- `tests/test_clustering.py` — the six test cases from the brief, copied
  verbatim: disconnected-component separation, run-to-run determinism,
  isolate singleton clusters, `cluster_all` resolution coverage, medoid on a
  star graph, and medoid existence/self-consistency for every cluster
  produced by `leiden`.

## TDD sequence and actual command output

### Step 1/2: write the test file, confirm it fails for the right reason

```
cd /Users/leo/Documents/binder_map
/Users/leo/.venvs/binder_map/bin/python -m pytest tests/test_clustering.py -v
```

Output:
```
============================= test session starts ==============================
platform darwin -- Python 3.10.21, pytest-9.1.1, pluggy-1.6.0
rootdir: /Users/leo/Documents/binder_map
configfile: pyproject.toml
collected 0 items / 1 error

==================================== ERRORS ====================================
__________________ ERROR collecting tests/test_clustering.py ___________________
ImportError while importing test module '/Users/leo/Documents/binder_map/tests/test_clustering.py'.
Hint: make sure your test modules/packages have valid Python names.
Traceback:
/opt/homebrew/Cellar/python@3.10/3.10.21/Frameworks/Python.framework/Versions/3.10/lib/python3.10/importlib/__init__.py:126: in import_module
    return _bootstrap._gcd_import(name[level:], package, level)
tests/test_clustering.py:3: in <module>
    from binder_map.clustering import RESOLUTIONS, cluster_all, leiden, medoids
E   ModuleNotFoundError: No module named 'binder_map.clustering'
=========================== short test summary info ============================
ERROR tests/test_clustering.py
!!!!!!!!!!!!!!!!!!!! Interrupted: 1 error during collection !!!!!!!!!!!!!!!!!!!!
=============================== 1 error in 0.22s ===============================
```

Matches the brief's expected failure exactly (`ModuleNotFoundError`).

### Step 3/4: implement, run again

```
cd /Users/leo/Documents/binder_map
/Users/leo/.venvs/binder_map/bin/python -m pytest tests/test_clustering.py -v
```

Output:
```
============================= test session starts ==============================
platform darwin -- Python 3.10.21, pytest-9.1.1, pluggy-1.6.0
rootdir: /Users/leo/Documents/binder_map
configfile: pyproject.toml
collected 6 items

tests/test_clustering.py ......                                          [100%]

============================== 6 passed in 0.25s ===============================
```

### Full suite regression check

```
/Users/leo/.venvs/binder_map/bin/python -m pytest -q
```

Result: 71 passed, 1 failed (`tests/test_geometry.py::test_missing_frame_atoms_raise_a_clear_error`).
That failure belongs to `src/binder_map/geometry.py` / `tests/test_geometry.py`,
both of which show as modified-but-uncommitted in `git status` at the start of
this task (another agent is actively working on that file per the task
instructions). It is unrelated to clustering and was not touched here.

Re-running with that one test deselected, to confirm everything else
(66 pre-existing + 6 new) is green:

```
/Users/leo/.venvs/binder_map/bin/python -m pytest --deselect tests/test_geometry.py::test_missing_frame_atoms_raise_a_clear_error
```

Output tail:
```
72 passed, 1 deselected, 86 warnings in 1.67s
```

## Deviation from the brief

None. The test file and implementation were used verbatim as given in
`.superpowers/sdd/task-8-brief.md`.

## Self-review findings

- **Isolate handling**: verified directly — `leiden(3, [[0,1]], [1.0], resolution=1.0)`
  returns `[0, 0, 1]`: node 2 (no edges) gets its own singleton label distinct
  from the {0,1} pair. igraph/leidenalg naturally give disconnected/unlisted
  vertices their own community with no special-casing required in the
  wrapper, so there is no bug here.

- **Determinism, including across process boundaries**: the brief's
  requirement is stronger than what the two in-test-process calls check
  (`test_leiden_is_deterministic` only proves determinism within one
  interpreter run). Since cluster ids land in a published artifact and must
  be reproducible on a fresh run, I additionally checked determinism across
  two *separate* Python process invocations on a synthetic 500-node /
  ~2000-edge random weighted graph (seed=42, resolution=1.0):
  ```
  run 1: labels saved to disk
  run 2 (fresh process): np.array_equal(run1_labels, run2_labels) -> True
  ```
  Confirmed byte-identical across process restarts, not just within one. No
  bug found; `leidenalg.find_partition(..., seed=seed)` seeds its RNG
  explicitly rather than relying on process entropy, which is why this holds.

- **`medoids()` tie-breaking**: `np.lexsort((np.arange(len(labels)), -score))`
  sorts primarily by `-score` (highest summed similarity first) and, for
  exact ties, by ascending node index. This makes the "best-connected member"
  choice deterministic even when multiple nodes tie on score, which matters
  for the same byte-identical-artifact requirement as `leiden`. No bug found.

- **`medoids()` approximation**: per the task instructions, I did not
  "upgrade" this to a true metric medoid — the docstring already documents
  that it's a similarity-weighted within-cluster degree over the sparse
  graph, which is the deliberate, affordable approximation at ~72k-node
  scale. Confirmed this reads correctly and matches the star-graph test
  (`test_medoid_is_the_best_connected_member`), where node 0 — connected to
  all three others — wins over any leaf.

- **No genuine bug found** in the brief-supplied implementation for this
  task; nothing was silently worked around or weakened.

## Files touched

- `/Users/leo/Documents/binder_map/src/binder_map/clustering.py` (new)
- `/Users/leo/Documents/binder_map/tests/test_clustering.py` (new)

Only these two files were staged and committed
(`git add src/binder_map/clustering.py tests/test_clustering.py`), per the
shared-checkout instruction — `src/binder_map/geometry.py` and
`tests/test_geometry.py`, modified by another agent, were left untouched and
unstaged.

## Commit

```
eb8b2656066ef56ccaeb0f1fc505ace3e705cb1a Add Leiden clustering at three resolutions with medoid election
 src/binder_map/clustering.py | 61 ++++++++++++++++++++++++++++++++++++++++++++
 tests/test_clustering.py     | 59 ++++++++++++++++++++++++++++++++++++++++++
 2 files changed, 120 insertions(+)
```

---

## Fix: guard against out-of-range edge indices

**Problem:** igraph silently grows the graph to `max(n, max_edge_id + 1)` when an edge references a node index >= n, returning labels with length > n. These misaligned labels would corrupt downstream tables when joined against the n-node input.

**Test sequence:**

### Step 1: Baseline (6 tests pass)
```
/Users/leo/.venvs/binder_map/bin/python -m pytest tests/test_clustering.py -v
```
Output: `6 passed in 0.11s`

### Step 2: Add two new tests (first one FAILS, proving the bug is real)
Added:
- `test_out_of_range_edge_raises_instead_of_growing_the_graph`
- `test_labels_always_have_length_n`

```
/Users/leo/.venvs/binder_map/bin/python -m pytest tests/test_clustering.py -v
```
Output:
```
__________ test_out_of_range_edge_raises_instead_of_growing_the_graph __________
>       with pytest.raises(ValueError, match="references node"):
E       Failed: DID NOT RAISE ValueError
tests/test_clustering.py:53: Failed
=========================== 1 failed, 7 passed in 0.13s ==========================
```
✓ **Test genuinely failed before the fix** — igraph silently grew the graph and returned 6 labels for n=3.

### Step 3: Apply the guard to `leiden()` (lines 2064-2072 from plan)
Added check before graph construction:
```python
if len(edges) and int(edges.max()) >= n:
    raise ValueError(
        f"edge references node {int(edges.max())} but n={n}"
    )
```

### Step 4: Re-run clustering tests (all 8 pass)
```
/Users/leo/.venvs/binder_map/bin/python -m pytest tests/test_clustering.py -v
```
Output: `8 passed in 0.08s`

### Step 5: Full suite (all 79 pass, no failures outside clustering)
```
/Users/leo/.venvs/binder_map/bin/python -m pytest -q
```
Output: `79 passed, 102 warnings in 9.00s`

### Step 6: Commit
```
e7d0f48 Guard against out-of-range edge indices in leiden clustering
 src/binder_map/clustering.py | 9 ++++++++-
 tests/test_clustering.py     | 19 +++++++++++++++++++
 2 files changed, 28 insertions(+)
```

**Summary:** Guard prevents silent length mismatch. New test `test_out_of_range_edge_raises_instead_of_growing_the_graph` genuinely failed before fix (no exception was raised), confirming the bug was real. After guard, all 8 clustering + 79 full suite tests pass. Commit SHA: `e7d0f48`.
