# Task 9 report — UMAP embedding from a precomputed kNN

## What was implemented

- `src/binder_map/embedding.py` — `embed(n, layout_edges, layout_weights, k=15, seed=SEED) -> np.ndarray`,
  plus the private helper `_to_dense_knn(n, edges, weights, k)` that turns the
  sparse, symmetric layout-graph edge list into the fixed-width
  `(indices, distances)` pair UMAP's `precomputed_knn` argument requires
  (distance = 1 − similarity; nodes short of `k` real neighbours are padded
  with self-references at distance 1.0).
- `tests/test_embedding.py` — the three tests from the brief, verbatim:
  `test_embed_returns_2d_coordinates`, `test_embed_is_deterministic`,
  `test_embed_keeps_distinct_groups_apart`.

Both files were written exactly as given in
`.superpowers/sdd/task-9-brief.md` — no deviation. I reviewed the brief's
code carefully for bugs (see Self-review below) before implementing it
verbatim, since the brief explicitly pre-verified the UMAP API and two prior
tasks in this project had genuine bugs in brief-supplied code. I did not find
one here that falls inside Task 9's contract.

## TDD cycle

### Step 1/2 — red run (before `embedding.py` existed)

Command:
```
cd /Users/leo/Documents/binder_map
/Users/leo/.venvs/binder_map/bin/python -m pytest tests/test_embedding.py -v
```

Output:
```
============================= test session starts ==============================
platform darwin -- Python 3.10.21, pytest-9.1.1, pluggy-1.6.0
rootdir: /Users/leo/Documents/binder_map
configfile: pyproject.toml
collected 0 items / 1 error

==================================== ERRORS ====================================
___________________ ERROR collecting tests/test_embedding.py ___________________
ImportError while importing test module '/Users/leo/Documents/binder_map/tests/test_embedding.py'.
Hint: make sure your test modules/packages have valid Python names.
Traceback:
/opt/homebrew/Cellar/python@3.10/3.10.21/Frameworks/Python.framework/Versions/3.10/lib/python3.10/importlib/__init__.py:126: in import_module
    return _bootstrap._gcd_import(name[level:], package, level)
tests/test_embedding.py:3: in <module>
    from binder_map.embedding import embed
E   ModuleNotFoundError: No module named 'binder_map.embedding'
=========================== short test summary info ============================
ERROR tests/test_embedding.py
!!!!!!!!!!!!!!!!!!!! Interrupted: 1 error during collection !!!!!!!!!!!!!!!!!!!!
=============================== 1 error in 0.13s ===============================
```

Fails for the expected reason: `ModuleNotFoundError: No module named 'binder_map.embedding'`, matching the brief exactly.

### Step 3/4 — implement, green run

Command:
```
cd /Users/leo/Documents/binder_map
/Users/leo/.venvs/binder_map/bin/python -m pytest tests/test_embedding.py -v
```

Output:
```
============================= test session starts ==============================
platform darwin -- Python 3.10.21, pytest-9.1.1, pluggy-1.6.0
rootdir: /Users/leo/Documents/binder_map
configfile: pyproject.toml
collected 3 items

tests/test_embedding.py ...                                              [100%]

=============================== warnings summary ===============================
tests/test_embedding.py::test_embed_returns_2d_coordinates
tests/test_embedding.py::test_embed_is_deterministic
tests/test_embedding.py::test_embed_is_deterministic
tests/test_embedding.py::test_embed_keeps_distinct_groups_apart
  .../site-packages/umap/umap_.py:1865: UserWarning: using precomputed metric; inverse_transform will be unavailable
    warn("using precomputed metric; inverse_transform will be unavailable")

tests/test_embedding.py::test_embed_returns_2d_coordinates
tests/test_embedding.py::test_embed_is_deterministic
tests/test_embedding.py::test_embed_is_deterministic
tests/test_embedding.py::test_embed_keeps_distinct_groups_apart
  .../site-packages/umap/umap_.py:1952: UserWarning: n_jobs value 1 overridden to 1 by setting random_state. Use no seed for parallelism.
    warn(

tests/test_embedding.py::test_embed_returns_2d_coordinates
tests/test_embedding.py::test_embed_is_deterministic
tests/test_embedding.py::test_embed_is_deterministic
tests/test_embedding.py::test_embed_keeps_distinct_groups_apart
  .../site-packages/umap/umap_.py:2021: UserWarning: precomputed_knn[2] (knn_search_index) is not an NNDescent object: transforming new data with transform will be unavailable.
    warn(

-- Docs: https://docs.pytest.org/en/stable/how-to/capture-warnings.html
======================== 3 passed, 12 warnings in 8.55s ========================
```

3 passed, matching the brief's expectation exactly. The three UserWarnings
(`inverse_transform unavailable`, `n_jobs overridden by random_state`,
`knn_search_index is not an NNDescent object`) are UMAP telling us, correctly,
that it's operating in precomputed/deterministic mode — consistent with the
brief's note that numba/UMAP warnings during fit are expected and should not
be suppressed.

### Full-suite regression run (after commit)

Command:
```
cd /Users/leo/Documents/binder_map
/Users/leo/.venvs/binder_map/bin/python -m pytest
```

Output (tail):
```
79 passed, 102 warnings in 9.25s
```
Exit code 0.

Note: the brief said the suite was "at 74 passing" before Task 9. By the time
I ran the full suite, two more commits (`96ff7cb`, `e7d0f48`, "Guard against
out-of-range edge indices in leiden clustering") had landed from another
agent working in this shared checkout, adding tests of their own — that
accounts for the baseline moving from 74 to 76 before my 3 land on top (79
total). `git status` before staging showed only my two untracked files
(`src/binder_map/embedding.py`, `tests/test_embedding.py`); nothing else was
touched or pending. I staged and committed only those two files, per
instructions — never `git add -A`/`git add .`.

## Deviations from the brief

None. `embedding.py` and `tests/test_embedding.py` were written exactly as
specified.

## Self-review findings

I read `_to_dense_knn` and `embed` closely for correctness before trusting
them, and did the following extra checks beyond the brief's three tests:

1. **Cross-process determinism.** The brief's `test_embed_is_deterministic`
   only calls `embed` twice within one process. Since "coordinates end up in
   a published artifact," I additionally ran the same embed twice in two
   *separate* Python processes and diffed the first rows byte-for-byte —
   identical:
   ```
   [[21.33198356628418, 11.323407173156738], [21.284990310668945, 13.52055835723877], [21.253746032714844, 11.927889823913574]]
   [[21.33198356628418, 11.323407173156738], [21.284990310668945, 13.52055835723877], [21.253746032714844, 11.927889823913574]]
   ```
   `random_state=seed` correctly forces `n_jobs=1` (UMAP's own warning
   confirms this), so this is genuinely reproducible across runs, not just
   within one.

2. **Sort-order correctness of `_to_dense_knn`.** UMAP's `precomputed_knn`
   expects each row of `(idx, dist)` sorted by ascending distance (nearest
   first). The brief's code sorts each node's neighbour list by
   `key=lambda t: -t[0]` (descending similarity), which is equivalent to
   ascending distance since `distance = 1 - similarity` is a strictly
   decreasing, order-preserving function of similarity for any real
   similarity value — so this holds regardless of whether similarity is
   bounded to `[0, 1]`. Correct as written.

3. **Padding-value edge case (checked, not a bug in scope).** The self-padding
   distance is hardcoded to `1.0`, which is only truly the *maximum* possible
   distance if similarities are bounded to `[0, 1]` (i.e. non-negative
   cosine similarity). That holds for this pipeline because the vectors fed
   into `knn_from_vectors` are non-negative contact-count fingerprints
   (`np.abs(...)` in the brief's own test fixture, and structurally
   non-negative in the real fingerprint pipeline), so cosine similarity of
   two non-negative unit vectors is always in `[0, 1]`. I did not change
   anything here — flagging it only so a future caller who feeds `embed()`
   arbitrary signed-similarity edges knows the padding value assumption.

4. **Isolate (zero-real-neighbour) nodes produce NaN coordinates, but this
   path is unreachable from the documented pipeline.** I constructed a
   synthetic edge list with two fully isolated nodes (all `k` slots
   self-padded) and confirmed UMAP emits `"It has fully disconnected N
   vertices"` and returns `NaN` for those rows — silently, no exception.
   However, `binder_map.knn.knn_from_vectors`'s `_topk_per_node` (which
   builds Task 7's `layout_edges`) guarantees every node contributes its own
   top-`min(k, n-1)` picks as source edges *unconditionally* — the layout
   graph has no floor, per the brief's own note ("The layout graph is
   deliberately UNTHRESHOLDED top-k... every node has neighbours to embed").
   So a node with **zero** real edges cannot occur when `embed()` is fed
   Task 7's actual `layout_edges`/`layout_weights`, only when fed a
   hand-built or corrupted edge list that violates that contract. I verified
   separately that *partial* padding (1 real neighbour, 2 padded) produces
   finite, correct output — only the fully-isolated case is NaN. Since this
   can't be triggered from Task 9's documented input contract, I left the
   code as-is rather than adding defensive checks outside the task's scope;
   noting it here for visibility rather than silently working around it or
   silently ignoring it.

No bug requiring a fix was found inside Task 9's actual contract. The verbatim
brief implementation is correct for its specified interface and inputs.

## Commands run, verbatim

```
cd /Users/leo/Documents/binder_map
/Users/leo/.venvs/binder_map/bin/python -m pytest tests/test_embedding.py -v   # red
# ... wrote src/binder_map/embedding.py ...
/Users/leo/.venvs/binder_map/bin/python -m pytest tests/test_embedding.py -v   # green, 3 passed
/Users/leo/.venvs/binder_map/bin/python -m pytest                              # full suite, 79 passed
git add src/binder_map/embedding.py tests/test_embedding.py
git commit -m "Add UMAP embedding over the precomputed layout kNN"
```

## Commit

`780007e` — "Add UMAP embedding over the precomputed layout kNN"
(2 files changed, 102 insertions: `src/binder_map/embedding.py`,
`tests/test_embedding.py`)

---

## Fix: disconnected nodes get finite coordinates

The Task 9 implementer's finding above ("No bug requiring a fix was found
inside Task 9's actual contract") missed one path: `binder_map.knn.symmetrise_min`
discards self-hits when building the layout kNN, so a structure whose only
structural-search hit is itself reaches `embed()` with **zero real edges** —
not a hand-built or corrupted input, but the ordinary output of the upstream
pipeline for a genuine outlier design. That case is reachable in production,
and the coordinates are published in an interactive map: a NaN there silently
drops the single most interesting kind of point (a design resembling nothing
else in the dataset) instead of showing it.

The fix (`_ring_positions` helper + a `connected`/`n_lonely` precheck in
`embed()`, then post-`fit_transform` ring placement for any disconnected
nodes, then a final `np.isfinite` guard) was transcribed verbatim from the
plan's Task 9 section and applied test-first.

### Baseline (before touching either file)

```
cd /Users/leo/Documents/binder_map
/Users/leo/.venvs/binder_map/bin/python -m pytest tests/test_embedding.py -v
```
`3 passed` (existing tests only; no `test_node_with_no_neighbours_gets_finite_coordinates` yet).

### New test added, embedding.py left unchanged — confirmed red

Added `test_node_with_no_neighbours_gets_finite_coordinates` to
`tests/test_embedding.py`, immediately before `test_embed_keeps_distinct_groups_apart`,
verbatim from the plan. Ran it before any implementation change:

```
FAILED tests/test_embedding.py::test_node_with_no_neighbours_gets_finite_coordinates
AssertionError: disconnected node produced NaN
assert np.False_
 +  where np.False_ = <built-in method all of numpy.ndarray object at ...>()
 +    where ... = array([[ True,  True],
       [ True,  True],
       [ True,  True],
       [ True,  True],
       [False, False]]).all
 +      where ... = <ufunc 'isfinite'>(array([[18.968859  ,  0.07933703],
       [19.896679  ,  0.49532467],
       [20.17198   , -0.39142317],
       [19.244862  , -0.71436775],
       [        nan,         nan]], dtype=float32))
```
UMAP's own warning corroborated the mechanism: *"It has fully disconnected 1
vertices."* Node 4 (the one with no edges) came back `[nan, nan]` — the
failure is exactly the one predicted, and it fails on the `np.isfinite(xy).all()`
assertion as expected. `1 failed, 3 passed`. This is the evidence the finding
was real.

### After the fix — green

Transcribed `_ring_positions` and the rewritten `embed()` from the plan into
`src/binder_map/embedding.py` verbatim. Re-ran:

```
/Users/leo/.venvs/binder_map/bin/python -m pytest tests/test_embedding.py -v
```
`4 passed, 16 warnings` — including `test_embed_is_deterministic`, confirming
the deterministic ring placement did not break reproducibility of the
published coordinates.

### Full suite

```
/Users/leo/.venvs/binder_map/bin/python -m pytest
```
`80 passed, 106 warnings in 7.99s` — no failures anywhere, nothing outside
`tests/test_embedding.py` touched or affected.

### Commit

```
git add src/binder_map/embedding.py tests/test_embedding.py
```
(status confirmed only these two files staged, no `-A`/`.`)

`d9dbc75` — "Give disconnected embedding nodes finite ring coordinates"
(2 files changed, 62 insertions(+), 3 deletions(-):
`src/binder_map/embedding.py`, `tests/test_embedding.py`)

---

## Fix: exclude disconnected nodes before UMAP

The previous fix above (`d9dbc75`) let disconnected nodes reach UMAP and
tried to overwrite their coordinates afterward. That can never work: UMAP's
default spectral init raises `IndexError` from inside
`umap.spectral.multi_component_layout` once three or more nodes are mutually
disconnected — the crash happens inside the library, before any cleanup code
runs. The previous test covered only ONE disconnected node
(`test_node_with_no_neighbours_gets_finite_coordinates`), which returned NaN
rather than crashing and so hid the failure entirely. Net effect: a silent
NaN was replaced by a hard crash — worse, not better.

This fix transcribes the corrected `MIN_UMAP_NODES` constant and rewritten
`embed()` verbatim from `docs/superpowers/plans/2026-08-21-foldseek-binder-map.md`
("Task 9: UMAP embedding from a precomputed kNN"). Disconnected nodes are now
removed BEFORE UMAP runs; the connected subgraph is remapped to compact
indices and embedded alone (preserving spectral init for the connected
majority, rather than degrading to `init="random"`); disconnected nodes are
then placed on a ring outside the embedded cloud, radius = true Euclidean
farthest-point distance (not axis-aligned/Chebyshev extent), with an explicit
`not np.isfinite(radius) or radius <= 0.0` fallback check (replacing the old
`... or 1.0`, which would have kept a NaN since NaN is truthy).

### Baseline (before any change)

```
/Users/leo/.venvs/binder_map/bin/python -m pytest tests/test_embedding.py -v
```
`4 passed, 16 warnings in 10.61s` — matches expectation exactly.

### New tests added, embedding.py left unchanged — confirmed red

Added `test_many_disconnected_nodes_do_not_crash_umap` (parametrized
`n_lonely` = 1, 3, 5, 12) and `test_embed_with_no_edges_at_all_is_finite`
verbatim from the plan, plus `import pytest`, plus the Euclidean-not-Chebyshev
rewrite of the assertion in `test_node_with_no_neighbours_gets_finite_coordinates`.
`embedding.py` was left completely unchanged. Ran:

```
/Users/leo/.venvs/binder_map/bin/python -m pytest tests/test_embedding.py -v
```

Result: `4 failed, 5 passed`. Per-case:

- `n_lonely=1` — **passed**. Only one disconnected node never reaches the
  spectral-init branch that crashes, so this case cannot expose the bug —
  exactly why the old single-node test hid it.
- `n_lonely=3` — **failed** with `IndexError`:
  ```
  >       dist = linkage(dm_i[:, component_labels == c_j])
  E       IndexError: boolean index did not match indexed array along axis 1; size of axis is 1 but size of corresponding boolean axis is 52
  ../../.venvs/binder_map/lib/python3.10/site-packages/umap/spectral.py:86: IndexError
  ```
- `n_lonely=5` — **failed** with the same `IndexError` from
  `umap/spectral.py:86` inside `multi_component_layout`.
- `n_lonely=12` — **failed** with the same `IndexError` from
  `umap/spectral.py:86` inside `multi_component_layout`.
- `test_embed_with_no_edges_at_all_is_finite` — **failed** too, but with a
  different, unrelated error (`ValueError: zero-size array to reduction
  operation maximum which has no identity`, from `graph.data.max()` in
  `umap_.py:1089`) — a fully edgeless graph was never exercised by the old
  code path at all.

This is the evidence the Critical was real: the crash reproduces reliably at
`n_lonely >= 3`, exactly as the plan predicted, and is invisible at
`n_lonely = 1`, exactly why the earlier single-node regression test could not
have caught it.

### After the fix — green

Transcribed `MIN_UMAP_NODES` and the rewritten `embed()` verbatim from the
plan into `src/binder_map/embedding.py`, keeping `_ring_positions` as-is.
Re-ran:

```
/Users/leo/.venvs/binder_map/bin/python -m pytest tests/test_embedding.py -v
```
`9 passed, 27 warnings in 8.66s` — all cases green, including
`n_lonely=1,3,5,12` and the no-edges case, and including
`test_embed_is_deterministic`, confirming the rewrite preserved determinism
of the published coordinates.

### Full suite

```
/Users/leo/.venvs/binder_map/bin/python -m pytest -q
```
`94 passed, 117 warnings in 9.37s` — no failures anywhere, nothing outside
`embedding.py`/`test_embedding.py` touched.

### Commit

```
git add src/binder_map/embedding.py tests/test_embedding.py
git status --short   # confirmed only these two files staged; no -A / .
```

`37ee77d` — "Fix Critical: exclude disconnected nodes from UMAP before embedding"
(2 files changed, 101 insertions(+), 35 deletions(-):
`src/binder_map/embedding.py`, `tests/test_embedding.py`)
