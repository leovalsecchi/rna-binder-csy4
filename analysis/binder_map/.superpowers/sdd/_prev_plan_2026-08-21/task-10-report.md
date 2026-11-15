# Task 10 report — Comparison metrics (`compare.py`)

## What was implemented

Created `src/binder_map/compare.py` and `tests/test_compare.py` exactly per the
brief, verbatim (no deviation from the code or test cases given in
`task-10-brief.md`). Four functions:

- `partition_scores(a, b) -> dict[str, float]` — AMI + ARI between two
  clusterings of the same nodes, keys `ami`, `ari`, `n`. Raises `ValueError`
  on length mismatch.
- `contingency(a, b) -> (matrix, labels_a, labels_b)` — dense co-membership
  count matrix between two partitions, via `np.unique(..., return_inverse=True)`
  + `np.add.at`.
- `fingerprint_agreement(fp_a, fp_b) -> np.ndarray` — row-wise cosine
  similarity; rows where either side has zero norm are forced to `0.0` via an
  explicit `good` mask (`na > 0 & nb > 0`), never computed through a
  divide-by-zero-then-patch path.
- `stack_profiles(per_nt, stacks) -> pd.DataFrame` — mean per-nucleotide
  profile per stack, columns `nt01..nt{RNA_LEN}` (16, from `binder_map.RNA_LEN`),
  grouped/sorted by stack label.

Imports are limited to `numpy`, `pandas`, `sklearn.metrics`, and
`binder_map.RNA_LEN` — no other `binder_map` submodule is imported, keeping
this the last pure-library piece as instructed.

## TDD cycle

### Red

```
cd /Users/leo/Documents/binder_map
/Users/leo/.venvs/binder_map/bin/python -m pytest tests/test_compare.py -v
```

Output (abridged):

```
ERROR collecting tests/test_compare.py
ImportError while importing test module '.../tests/test_compare.py'.
E   ModuleNotFoundError: No module named 'binder_map.compare'
=========================== short test summary info ============================
ERROR tests/test_compare.py
!!!!!!!!!!!!!!!!!!!! Interrupted: 1 error during collection !!!!!!!!!!!!!!!!!!!!
=============================== 1 error in 0.14s ===============================
```

Failed for the expected reason (module not yet created).

### Green — module tests

```
cd /Users/leo/Documents/binder_map
/Users/leo/.venvs/binder_map/bin/python -m pytest tests/test_compare.py -v
```

```
tests/test_compare.py::test_identical_partitions_score_one PASSED
tests/test_compare.py::test_relabelled_partition_still_scores_one PASSED
tests/test_compare.py::test_independent_partitions_score_near_zero PASSED
tests/test_compare.py::test_contingency_counts_co_membership PASSED
tests/test_compare.py::test_agreement_of_identical_fingerprints_is_one PASSED
tests/test_compare.py::test_agreement_of_orthogonal_fingerprints_is_zero PASSED
tests/test_compare.py::test_agreement_is_zero_when_a_row_has_no_contacts PASSED
tests/test_compare.py::test_agreement_ignores_magnitude PASSED
tests/test_compare.py::test_stack_profiles_averages_within_stack PASSED

============================== 9 passed in 1.60s ===============================
```

All 9 passed on the first implementation attempt — no edits were needed
between red and green.

### Green — full suite

```
cd /Users/leo/Documents/binder_map
/Users/leo/.venvs/binder_map/bin/python -m pytest -v
```

```
tests/test_backbone.py .......                                           [  7%]
tests/test_chains.py .......                                             [ 15%]
tests/test_clustering.py ........                                        [ 24%]
tests/test_compare.py .........                                          [ 34%]
tests/test_constants.py ...                                              [ 38%]
tests/test_embedding.py ....                                             [ 42%]
tests/test_fingerprint.py ..........                                     [ 53%]
tests/test_geometry.py .................                                 [ 73%]
tests/test_knn.py ........                                               [ 82%]
tests/test_rnaframe.py ..........                                        [ 93%]
tests/test_upstream_parity.py ......                                     [100%]

======================= 89 passed, 106 warnings in 8.59s =======================
```

Baseline before this task was 80 passed (confirmed by a pre-implementation
run: `/Users/leo/.venvs/binder_map/bin/python -m pytest --collect-only -q`
summed to 80 across `test_backbone` 7, `test_chains` 7, `test_clustering` 8,
`test_constants` 3, `test_embedding` 4, `test_fingerprint` 10, `test_geometry`
17, `test_knn` 8, `test_rnaframe` 10, `test_upstream_parity` 6). The brief's
own "67 passed" expectation is stale relative to this checkout — tasks 1–9
have since grown to 80 tests (extra `test_upstream_parity.py`, and more cases
than the brief anticipated in `clustering`/`embedding`/`geometry`/etc.). 80 + 9
new `compare` tests = 89, which is exactly what ran. This is not a deviation
in my code, just a note that the brief's step-5 count is out of date for this
branch.

No regressions: all 80 pre-existing tests remain green alongside the 9 new
ones.

## Deviations from the brief

None in the code or tests — both files were written verbatim as given. The
only discrepancy noted is the stale "67 passed" / per-module count in the
brief's Step 5 (see above), which is an artifact of the brief being written
before later tasks added more tests than it anticipated — not something I
changed.

## Self-review

I did not treat "the brief's code passed its own tests" as sufficient, given
this project's history of genuine bugs slipping through brief-supplied code.
I read each function against the two deliberate semantics called out in the
task (a low AMI/ARI is a real result and must never be nudged upward; an
all-zero fingerprint row must score 0.0 against everything, including another
all-zero row) and then probed beyond the pinned test cases:

- **`fingerprint_agreement`, both rows all-zero.** Not covered by the pinned
  tests (which only cover one-row-zero). Verified directly:
  `fingerprint_agreement([[0,0,0]], [[0,0,0]])` → `[0.]`. Correct — the `good`
  mask requires both norms `> 0`, so two absent epitopes never register as
  "agreeing," matching the brief's explicit warning not to manufacture
  agreement out of two absences.
- **`fingerprint_agreement`, mixed batch** (some zero rows, some not, to check
  no cross-row leakage from the masked assignment): rows `[[1,0],[0,0],[1,1]]`
  vs `[[1,0],[1,1],[1,1]]` → `[1., 0., 1.]` as expected — the masked
  vectorized write doesn't let a zero row's computation bleed into a good
  row's slot or vice versa.
- **`partition_scores`, non-numeric labels** (string cluster ids, matching how
  Leiden or stack names might realistically arrive) — works via sklearn's
  general label handling; not required by the interface but confirms no
  hidden assumption of integer labels.
- **`contingency` at realistic scale**: 5000 draws over up to 3000 labels each
  produced a 2465×2436 dense `int64` matrix (~48 MB), `sum() == 5000` as
  expected. This confirms the brief's own caveat: at fine Leiden resolutions
  on ~72k designs, a dense contingency matrix of a few thousand by a few
  thousand labels is on the order of tens of MB, not a problem on its own,
  but worth keeping in mind if `contingency` is ever called with the *finest*
  resolution's labels on the full 72k-design layer rather than on
  already-reduced cluster/stack labels — a naive n×n case (72k×72k) would be
  ~40 GB and must never happen. The function itself is fine; the risk is only
  in how a caller chooses `a`/`b`.
- **No deprecation/warning issues**: ran with
  `-W error::DeprecationWarning` against `partition_scores` — sklearn 1.7.2's
  `adjusted_mutual_info_score` default `average_method='arithmetic'` produces
  no warnings.

**No bugs found.** The brief's implementation is correct as given, matches
both deliberate semantics exactly, and needed no changes. This is unlike the
three prior tasks in this project where brief-supplied code had genuine bugs
— here the code was verified correct rather than assumed correct.

## Commit

```
git add src/binder_map/compare.py tests/test_compare.py
git commit -m "Add partition comparison, epitope agreement, and stack profiles"
```

Commit: `3602f56` on branch `foldseek-binder-map`. Only the two intended files
were staged (`git status --short` before staging showed exactly
`src/binder_map/compare.py` and `tests/test_compare.py` as untracked — no
`git add -A`/`git add .` was used, per the shared-checkout instruction).
