# Task 7 report: Stage 9 map comparison

## What was implemented

Exactly what the brief specified, verbatim, nothing more:

- `scripts/stage9_compare.py` — reads the three finished maps
  (`sample500k`, `passing60k`, `library12k`) and compares them:
  - `expand_members(members)` — explodes the comma-joined `design_ids` column
    of `dedup/members.tsv` into one row per `(struct_name, design_id)`.
  - `design_to_cluster(results, resolution)` — joins members -> nodes ->
    cluster assignment, returns a `design_id -> cluster` Series for one map.
  - `fold_key(design_id)` — collapses the dedup key's known over-splitting:
    strips `_seq<N>` for non-HD stacks, strips `_sample<S>_pack<P>` for HD.
  - `summarise(labels, name)` — per-map dict: node/design/cluster counts,
    singleton count, largest/median/p90 cluster size, `in_passing`/`in_12k`
    totals.
  - `crosswalk_scores(old, new)` — ARI/AMI between two cluster-label Series,
    restricted to the shared index (`n_shared, ari, ami`); returns NaN scores
    when there is no overlap.
  - `load_labels(results, build)` — per-node label table joining dedup output
    with `meta/nodes.parquet` (or, for the pre-stage-0 `passing60k` reference
    map, rebuilding membership by treating everything in it as `in_passing`
    and cross-checking `in_12k` against the library export).
  - `main()` — drives all three maps, writes `compare/map_summary.parquet`,
    `compare/cluster_table_<map>.parquet` (with `size` **and** `folds`
    columns side by side), and `compare/crosswalk.parquet` for the three
    pairwise comparisons.
- `tests/test_compare_maps.py` — the 11 tests specified in the brief, copied
  verbatim.

No other files were touched. I did not run `main()` against real result
trees — none of the three maps exist yet (confirmed: `results_full`,
`results_sw_k30_f50_mut`, `results_lib12k` all absent under
`/cluster/work/wicky/lvalsecchi/rna_binder/analyses/foldseek_map` on Euler).
Per the task instructions, I tested only against the brief's synthetic
fixtures and did not attempt to execute `main()`.

## TDD evidence

**RED** — `tests/test_compare_maps.py` written first, `scripts/stage9_compare.py`
did not exist yet.

Command:
```
ssh euler 'cd .../foldseek_map && PYTHONPATH=src:scripts .../binder_map_venv/bin/python -m pytest tests/test_compare_maps.py -v'
```
Output (tail):
```
ERROR collecting tests/test_compare_maps.py
...
E   ModuleNotFoundError: No module named 'stage9_compare'
=========================== short test summary info ============================
ERROR tests/test_compare_maps.py
!!!!!!!!!!!!!!!!!!!! Interrupted: 1 error during collection !!!!!!!!!!!!!!!!!!!!
============================== 1 error in 15.53s ===============================
```
This is the expected failure mode — the module genuinely doesn't exist yet,
not a fixture/path/import-order accident.

**GREEN** — after writing `scripts/stage9_compare.py` (brief's code, copied
verbatim) and syncing:
```
ssh euler 'cd .../foldseek_map && PYTHONPATH=src:scripts .../binder_map_venv/bin/python -m pytest tests/test_compare_maps.py -v'
```
Output:
```
collected 11 items
tests/test_compare_maps.py ...........                                   [100%]
============================= 11 passed in 21.08s ==============================
```
Matches the brief's "Expected: PASS, 11 passed" exactly.

## Discriminate-check evidence (the "test that cannot fail" trap)

Target: `fold_key`, and by extension the `folds` column of the cluster table
(`folds` is `nunique(fold_key(...))` over each cluster — there is no separate
test for the table itself in the brief's suite, since `main()` is untestable
without real data; the check on `fold_key` is the check that covers it, as
that function is the entire collapsing logic the `folds` column relies on).

Procedure: temporarily replaced `fold_key`'s body with the identity function
(`return design_id`, no `_seq`/`_sample_pack` stripping at all), synced to
Euler, and re-ran just the `fold_key`-related tests:

```
ssh euler '... pytest tests/test_compare_maps.py -v -k fold_key'
```
Result with the neutered `fold_key`:
```
tests/test_compare_maps.py FF.                                           [100%]
FAILED tests/test_compare_maps.py::test_fold_key_collapses_bg_bif_seq_variants
  AssertionError: assert 4 == 1
FAILED tests/test_compare_maps.py::test_fold_key_collapses_hd_sample_variants
  AssertionError: assert 'csy4free_HD_...sample8_pack1' == 'csy4free_HD_...sample0_pack1'
2 failed, 1 passed, 8 deselected in 19.37s
```
The two collapsing tests fail exactly as expected when collapsing is
disabled (`test_fold_key_keeps_different_hd_recycles_apart` still passes,
correctly, since identity never *incorrectly* collapses two different
recycles either). This confirms the tests are load-bearing, not
vacuously-passing. I then restored the real implementation from a backup
(`diff` confirmed byte-identical to the pre-edit version), re-synced, and
re-ran the full file: 11 passed again (same output as the GREEN run above).

## Shipped library CSV (`csy4_binder_library_12k_v6_k30.csv`)

Checked its header and a couple of rows:
```
n,design_id,sequence,length,stack,generator,inverse_folder,campaign,
backbone_id,fold_cluster,epitope_cluster,iptm_af3,iptm_esm,iptm_max,
sc_rmsd_af3,sc_rmsd_esm,dockq_af3,dockq_esm,dsasa_af3,dsasa_esm,
n_hbonds_af3,n_hbonds_esm
```
12,000 data rows, `fold_cluster` populated with plain integers (e.g. `136`,
`72`) keyed by `design_id`.

It is useful, but I did not build it into `stage9_compare.py` because the
brief doesn't call for it and doing so would be scope creep (YAGNI): the
brief's `load_labels`/`design_to_cluster` already re-derive `passing60k`'s
cluster labels directly from `results_sw_k30_f50_mut/design/clusters/fold/clusters.parquet`,
which is the authoritative source. What this CSV *does* give, for free, to
whoever eventually runs Stage 9 for real: a correctness check with zero extra
compute — join the CSV's `fold_cluster` (for the 12,000 library rows) against
the re-derived `design_to_cluster(results_sw_k30_f50_mut design)` on
`design_id`, and confirm they agree. If they don't, either the resolution
constant (`cluster_r80`) is wrong or the re-derivation logic has a bug — this
would be a cheap, independent sanity check worth running once real data
exists. I flag it here per the task's request rather than implementing it,
since it's not in the brief's required interface.

## Files changed

- `scripts/stage9_compare.py` (new, 189 lines)
- `tests/test_compare_maps.py` (new, 85 lines)

Both copied verbatim from the brief.

## Full suite

Baseline stated in the task: 147 passed, 12 skipped, 0 failed (159 collected).

Ran on Euler:
```
ssh euler 'cd .../foldseek_map && PYTHONPATH=src:scripts .../binder_map_venv/bin/python -m pytest -q -p no:warnings'
```
Exit code: 0. No `F` or `E` markers anywhere in the dot/letter output.

Note: the run's final numeric summary line ("N passed, M skipped in Ts")
did not render over the non-interactive ssh session in three separate
attempts (a pre-existing quirk of this environment, not something I could
attribute to the test change — same thing happened on a bare
`--collect-only` run). I do not want to assert pass/fail counts from a line I
never actually saw, so I verified independently instead:
- Counted every `.` (pass) and `s` (skip) character in the raw progress
  output: **158 dots, 12 s's, 170 total**, zero other characters (no `F`/`E`).
- Cross-checked with `pytest --collect-only -q`, summed per-file counts by
  hand: **170** collected, matching exactly.
- 158 = 147 (baseline) + 11 (this task's new tests). 12 skipped is unchanged
  (the pre-existing `passing.csv`/`upstream rna_binder` skips). 170 collected
  = 159 (baseline) + 11.

So: **158 passed, 12 skipped, 0 failed, 170 collected** — consistent with the
baseline plus exactly this task's 11 tests, no regressions elsewhere.

## Self-review

- **Completeness**: all four required interfaces implemented
  (`expand_members`, `design_to_cluster`, `summarise`, `crosswalk_scores`),
  plus the supporting `fold_key` and `load_labels`/`main` the brief specified.
  Nothing partially done.
- **Naming**: matches the brief's interface signatures exactly (function
  names, argument names, return shapes/keys).
- **No overbuilding**: implementation is the brief's code verbatim; I did not
  add the library-CSV cross-check, did not add CLI flags beyond `--out`, did
  not try to run `main()` against non-existent data trees.
- **Real-behaviour tests**: confirmed via the discriminate check above that
  `fold_key`'s tests fail when the collapsing logic is removed — this is not
  a vacuous test.
- **Pristine output**: `git status` after commit shows only the two intended
  files staged and committed; no `__pycache__`, no stray files. Confirmed
  with `git show --stat HEAD` — exactly `scripts/stage9_compare.py` (189
  insertions) and `tests/test_compare_maps.py` (85 insertions), 274 total,
  no deletions/modifications elsewhere.

## Concerns

- None functional. The only oddity encountered was environmental: two
  transient SSH connection resets during the authoring loop ("Connection
  closed... unexpected end of file" from rsync), which resolved on retry —
  not related to the code. Also noted (and worked around) that `pytest`'s
  final summary line doesn't print over this non-interactive ssh session;
  documented the independent verification method above so this doesn't
  recur as ambiguity for the next task.
- As instructed, `main()` was never executed against real trees since none
  of the three maps are built yet. This means the `load_labels`/`main`
  code path (joining `meta/nodes.parquet`, the `passing60k` reference-map
  reconstruction, and the parquet-writing) is exercised only by the brief's
  unit tests on synthetic fixtures — a smoke test against real data will be
  needed once `results_full`, `results_lib12k`, and the existing
  `results_sw_k30_f50_mut` are all in place with the expected `clusters_r80`
  resolution column present.
