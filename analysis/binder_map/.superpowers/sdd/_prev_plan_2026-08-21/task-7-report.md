# Task 7 report: kNN graph construction for both axes

Commit: `51182f56277df81ffdab3160bfddd03707efad0c` on branch `foldseek-binder-map`
("Add axis-agnostic kNN with separate cluster and layout graphs")

## What was implemented

- `src/binder_map/knn.py` — pure-numpy module (only `from dataclasses import
  dataclass` and `import numpy as np`; no `binder_map` imports, as required):
  - `KNNResult` frozen dataclass: `cluster_edges`, `cluster_weights`,
    `layout_edges`, `layout_weights`, `degree`, `is_isolate`.
  - `symmetrise_min(qi, qj, qtm, ttm)` — collapses directed Foldseek hits to
    undirected `min(qTM, tTM)` pairs, dropping self-hits, keeping the max of
    the two observed minima where both directions were seen.
  - `knn_from_pairs(n, i, j, sim, k, floor)` — builds both graphs from an
    existing sparse similarity list (fold axis).
  - `knn_from_vectors(vectors, k, floor, block=4096)` — exact cosine kNN over
    L1-then-L2-normalised rows via a blocked matrix product (epitope axis).
  - Private helpers `_topk_per_node` (per-node top-k as an undirected edge
    set) and `_assemble` (builds the floor-thresholded cluster graph and the
    unthresholded layout graph from the same directed candidate list).
- `tests/test_knn.py` — the 8 tests specified verbatim in the brief.

## TDD cycle, with actual commands and output

### Step 0: baseline (before touching anything)

```
$ /Users/leo/.venvs/binder_map/bin/python -m pytest -q
```
```
................................................                         [100%]
```
48 dots = 48 passing, matching the context note ("suite is at 48 passing").

### Step 1/2: write the test file, watch it fail (RED)

`tests/test_knn.py` was written exactly as given in the brief (Step 1 code
block, verbatim).

```
$ /Users/leo/.venvs/binder_map/bin/python -m pytest tests/test_knn.py -v
```
```
ERROR collecting tests/test_knn.py
ImportError while importing test module '.../tests/test_knn.py'.
tests/test_knn.py:4: in <module>
    from binder_map.knn import knn_from_pairs, knn_from_vectors, symmetrise_min
E   ModuleNotFoundError: No module named 'binder_map.knn'
=========================== short test summary info ============================
ERROR tests/test_knn.py
1 error in 0.16s
```
Failed for the expected reason (module doesn't exist yet).

### Step 3: implement `src/binder_map/knn.py` (first pass, verbatim from brief)

Wrote the module exactly as given in the brief's Step 3 code block, then ran
the tests.

```
$ /Users/leo/.venvs/binder_map/bin/python -m pytest tests/test_knn.py -v
```
```
tests/test_knn.py .......F                                               [100%]

=================================== FAILURES ===================================
_______________ test_layout_graph_gives_every_node_k_neighbours ________________
    res = knn_from_vectors(v, k=5, floor=0.99)
    counts = np.bincount(res.layout_edges.ravel(), minlength=50)
>   assert (counts >= 5).all()
E   assert np.False_
E    +  where ... = array([ 4, 5, 3, 3, 4, 11, 3, 4, 5, 3, 3, 4, 6, 4, 6, 5, 4,
       4, 4, 5, 7, 9, 3, 4, 3, 7, 3, 5, 5, 4, 4, 4, 5, 3, 4, 4, 3, 6, 8, 4, 4,
       3, 3, 6, 5, 6, 8, 7, 4, 5]) >= 5.all
tests/test_knn.py:91: AssertionError
========================= 1 failed, 7 passed in 0.09s ==========================
```
7 of 8 passed. `test_layout_graph_gives_every_node_k_neighbours` failed — this
is a real bug in the brief's given implementation, not a test problem (see
Deviation below).

### Root-cause debugging

Traced node 0 of the failing case (`rng=default_rng(1)`, 50×8 random
vectors, `k=5`): `knn_from_vectors`'s own directed top-k computation was
verified correct in isolation — node 0's true top-5 neighbours are
`[46, 18, 30, 12, 4]` with cosine sims `[0.970, 0.965, 0.943, 0.913, 0.901]`
(cross-checked against `np.argsort` on the raw similarity matrix). But the
final deduped `layout_edges` for node 0 were only
`[18, 24, 30, 46]` — neighbours 12 and 4 were missing, and neighbour 24
(sim 0.891, weaker than every one of node 0's true top-5) had leaked in.

Instrumenting `_topk_per_node`'s rank computation showed why: node 46 is a
**mutual** top-k pick — it is in node 0's own top-5 *and* node 0 is
independently in node 46's own top-5. Building `src=[i,j], dst=[j,i]`
mirrors every directed edge, so this one relationship produces **two** rows
in node 0's incident run (its own forward entry, plus node 46's mirrored
entry), both with rank 0 and 1. That duplicate silently eats two rank slots
for one neighbour, so the `rank < k` filter pushes node 0's genuinely
distinct 4th and 5th neighbours (12, 4) out of the window entirely — while a
weaker mirrored edge from some other node's pick of node 0 (24, at rank 9)
can survive elsewhere in the reshuffled ranking. Net effect: nodes with any
mutual top-k pick end up with **fewer than k** neighbours, and — more
seriously for the cluster graph — a genuinely above-floor edge can be
displaced by a duplicate and dropped from `cluster_edges` even though
nothing about it was weak.

### Deviation from the brief

**Fixed a bug in `_topk_per_node`.** Added a dedup pass over the directed
`(src, dst)` pairs (keeping the max weight per pair) *before* computing the
per-node rank, so a mutual pick occupies exactly one rank slot instead of
two:

```python
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
```

Everything else in the module matches the brief verbatim (interfaces,
`KNNResult` fields, `symmetrise_min`, `knn_from_pairs`, `knn_from_vectors`,
docstrings). Only this one block inside `_topk_per_node` was added; no other
logic was changed.

Why this had to be fixed rather than the test weakened: the test pins an
explicit interface contract stated in the brief itself
("the layout algorithm downstream requires a full k neighbours per node").
The bug is not cosmetic — at the stated production scale (~72k nodes, up to
300 hits/query) mutual top-k picks are common (any strong structural or
epitope match is very likely reciprocal), so this would have silently
starved a meaningful fraction of nodes of their full k in the layout graph,
and could drop legitimately above-floor edges from the cluster graph too.

### Step 4 (repeated): run the tests after the fix (GREEN)

```
$ /Users/leo/.venvs/binder_map/bin/python -m pytest tests/test_knn.py -v
```
```
tests/test_knn.py ........                                               [100%]
============================== 8 passed in 0.12s ===============================
```
All 8 tests specified in the brief pass.

### Extra verification (self-review, not in the brief)

To make sure the fix generalises and isn't just papering over the one seed
in the test, ran a stress check across 20 random seeds with varying n
(20–120) and dimensionality (4–20), floor set high enough to force every
node's own top-k to be the deciding factor:

```python
for seed in range(20):
    ...
    res = knn_from_vectors(v, k=k, floor=0.999999)
    counts = np.bincount(res.layout_edges.ravel(), minlength=n)
    assert (counts >= k).all()
```
Output: `stress test done` — no failures across all 20 seeds (previously,
before the fix, this pattern of failure was reproducible on essentially any
seed with a large enough n/k that mutual picks were likely).

### Full suite after the change

```
$ /Users/leo/.venvs/binder_map/bin/python -m pytest -q
```
```
..................................................................       [100%]
```
66 dots. Baseline was 48 before this task started; by the time of this run
other agents working concurrently in this shared checkout had added 10 more
tests to pre-existing files (confirmed via `pytest --collect-only -q`, and
via `git status` showing only `src/binder_map/knn.py` and `tests/test_knn.py`
as untracked before staging — everything else was already committed by
other agents). 58 pre-existing + 8 new = 66, all green. No regressions.

## Self-review findings

1. **Real bug found and fixed** (see Deviation above) — mutual top-k picks
   double-consumed a rank slot in `_topk_per_node`, violating the "every
   node gets k neighbours" layout-graph guarantee and risking silent edge
   loss in the cluster graph. Fixed with a pre-rank dedup pass; verified with
   the brief's test plus an independent 20-seed stress check.
2. Manually traced all 8 test cases by hand against the implementation logic
   (lexsort ordering, rank-via-run-position, floor masking, isolate/degree
   bookkeeping) before running anything, to have an independent check beyond
   "the brief said so."
3. Confirmed the module has zero `binder_map`-internal imports (`grep` of
   the file's imports: only `dataclasses` and `numpy`), matching the "pure
   numeric, standalone" requirement.
4. Benign observation, not fixed: `knn_from_vectors` occasionally triggers a
   `RuntimeWarning: invalid value encountered in matmul` from the blocked
   `unit[start:stop] @ unit.T` product on some random inputs (e.g. seed 19,
   n=78, d=10). Checked directly — the resulting `sims`/`layout_weights`/
   `cluster_weights` contain no NaN and node degrees are still correct; this
   looks like a transient floating-point flag from the BLAS backend (a known
   benign quirk on macOS Accelerate) rather than an actual invalid value
   reaching the output. Not in scope of the brief's tests; flagging in case
   it recurs at the real 72k-node scale on a different BLAS backend (Euler),
   where it would be worth a quick `np.isnan` sanity check on real data if
   ever observed there.
5. `_topk_per_node`'s `n` parameter is unused inside the function body (also
   true of the brief's original code) — harmless, since it's a private
   helper and callers pass it for symmetry with `_assemble`'s signature. Not
   changed, to keep the diff minimal and matching the brief's interface.
6. Ran `git status --short -- src/ tests/` before staging to confirm only
   the two intended files were untracked, then staged them by explicit path
   (`git add src/binder_map/knn.py tests/test_knn.py`) — no `git add -A`/`.`
   used, per the shared-checkout instruction.
