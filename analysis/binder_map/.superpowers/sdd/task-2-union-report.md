# Task 2 — Union build driver (stages 2-5, fold axis)

**Status: COMPLETE.** Driver written TDD, 21 tests pass, mutation-checked, real
chain submitted and its dependencies verified with `scontrol`. Not committed —
git is refused in this environment.

---

## Files touched

| File | Change |
|---|---|
| `submit_union.sh` | **new** — the driver (stages 2-5 over `results_union`, fold axis) |
| `tests/test_union_submit.py` | **new** — 21 tests, `sbatch` stubbed on `PATH`, nothing submitted |
| `sync.sh` | **modified** — added `submit_union.sh` to the file push list. Without it the driver never reaches Euler and the tests cannot see it. |

`slurm/stage3_search.sbatch`, `scripts/stage2_builddb.sh`, `scripts/stage4_graph.py`,
`scripts/stage5_cluster.py` — **unmodified**, as required.

---

## What the driver does

No stage 1 or 1b (Task 1 already built the tree). Flat, four `sbatch` calls:

- **preflight** — refuses to submit unless `results_union/design/shards/shard_*.tar.gz`
  *and* `hashes/shard_*.tsv` both exist, and prints what it found. Fails closed.
- **s2** `normal.4h 04:00:00 8cpu 8G` — `RESULTS=results_union bash scripts/stage2_builddb.sh design`
  `&&` `stage3_split_queries.py --layer design --n-chunks 512 --results $FM/results_union`.
  No `--dependency`: its input is already on disk.
- **s3** `--array=0-511%128 --time=00:20:00 --mem-per-cpu=1G --export=ALL,ALN_TYPE=1,RESULTS=results_union slurm/stage3_search.sbatch design`.
  Submitted as the shipped sbatch, **not** re-wrapped, so its node-local staging of
  the target DB survives. Time/memory are the 2026-09-25 measurements (56 s, 927 MB peak).
- **s4** `normal.4h 04:00:00 16cpu 8G` — `stage4_graph.py --layer design --tm-floor 0.5 --epi-floor 0.6 --k 30 --mutual --results $FM/results_union`.
- **s5** `slurm/stage_single.sbatch scripts/stage5_cluster.py --layer design --axis fold --results $FM/results_union`. Once, fold only.

`set -euo pipefail` + `shopt -s inherit_errexit` + `|| exit 1` on all four captures.
The driver deliberately has **no** helper function, so no `$( )` here wraps more
than a single command — the guards stay anyway, documented as such in the header.

Two deliberate departures from `submit_maps.sh`:

- `FM="${FM:-/cluster/.../foldseek_map}"` — overridable *only* so the preflight can
  be tested against a fixture tree. Default is byte-identical to the hardcoded one.
- It does **not** write `logs/last_job_*.txt`. `submit_maps.sh` does; here a stubbed
  test run would clobber a real job id in it.

---

## Test run (real output)

RED, before `submit_union.sh` existed — all 21 fail, none vacuously:

```
$ ssh euler 'cd $FM && PYTHONPATH=src:scripts $PY -m pytest tests/test_union_submit.py -o addopts=' | tail -3
FAILED tests/test_union_submit.py::test_stage3_is_the_shipped_sbatch_not_a_wrap
FAILED tests/test_union_submit.py::test_stage4_carries_the_frozen_graph_parameters
============================== 21 failed in 3.75s ==============================
```

(The first RED run had 3 vacuous passes — `test_every_stage_targets_the_union_tree`,
`test_refuses_to_submit_when_the_hashes_are_missing` and
`test_a_failed_submission_stops_the_chain[1]` all passed over an empty row list.
Each got an extra assertion — `len(rows) == 4`, `"FATAL" in stderr`, and
`"stub sbatch: refusing submission" in stderr` — which is what took RED to 21/21.)

GREEN:

```
$ ssh euler 'cd $FM && PYTHONPATH=src:scripts $PY -m pytest tests/test_union_submit.py -o addopts= -v' | tail -1
============================== 21 passed in 3.79s ==============================
```

Whole suite:

```
$ ssh euler 'cd $FM && PYTHONPATH=src:scripts $PY -m pytest -o addopts=' | tail -1
========== 307 passed, 12 skipped, 112 warnings in 120.32s (0:02:00) ===========
```

**Reconciliation:** the handover figure was 241 passed / 12 skipped plus Task 1's
45 = 286. 286 + my 21 = **307 passed, 12 skipped**. Measured, exact match, no drift.

`bash -n submit_union.sh` → rc 0 (also asserted by `test_the_driver_parses`).

### Mutation checks (the tests are discriminating)

| Mutation applied to a copy on Euler | Result |
|---|---|
| `--axis fold` → `--axis epitope` | 1 failed (`test_stage5_runs_once_on_the_fold_axis`) |
| `--time=00:20:00` → `04:00:00`, and stage 4's `--dependency` deleted | 2 failed (`test_stage3_carries_the_measured_shape`, `test_the_chain_is_wired_s2_s3_s4_s5`) |
| `set -euo` → `set -uo`, `shopt` line deleted, every `) \|\| exit 1` → `)` | 6 failed — incl. all four `test_a_failed_submission_stops_the_chain[1..4]`, i.e. the driver really did carry on past a refused `sbatch` |

The file was restored after each and md5-verified identical to the local copy
(`6540918a08d3c983d1b31dfb7f49bff7`).

**Note on the brief's own test.** `test_submits_fold_axis_only` is weak alone: under
the axis mutation it still passed, because the echo label and the job name still
read "fold". `test_stage5_runs_once_on_the_fold_axis` (asserts `--axis fold` in the
stage-5 argv and that no submission carries `--axis epitope`) is what actually
holds that requirement.

---

## Submitted chain — job ids and verified dependencies

```
$ ssh euler 'cd $FM && bash submit_union.sh'
union tree: 99 shards, 99 hash tables in .../results_union/design
union s2 dedup+createdb+split   15387743
union s3 all-vs-all search      15387746
union s4 knn graph              15387748
union s5 fold clusters          15387750
chain: s2=15387743 -> s3=15387746 -> s4=15387748 -> s5 fold=15387750
```

```
$ ssh euler 'for j in 15387743 15387746 15387748 15387750; do scontrol show job $j | ...; done'
JobId=15387743 JobName=union_s2       JobState=RUNNING Reason=None       Dependency=(null)                        TimeLimit=04:00:00 Partition=normal.4h NumCPUs=8  MinMemoryCPU=8G
JobId=15387746 JobName=union_s3       JobState=PENDING Reason=Dependency Dependency=afterok:15387743(unfulfilled) TimeLimit=00:20:00 Partition=normal.4h NumCPUs=16 MinMemoryCPU=1G  ArrayTaskThrottle=128
JobId=15387748 JobName=union_s4       JobState=PENDING Reason=Dependency Dependency=afterok:15387746_*(unfulfilled) TimeLimit=04:00:00 Partition=normal.4h NumCPUs=16 MinMemoryCPU=8G
JobId=15387750 JobName=union_s5_fold  JobState=PENDING Reason=Dependency Dependency=afterok:15387748(unfulfilled) TimeLimit=04:00:00 Partition=normal.4h NumCPUs=16 MinMemoryCPU=4G
```

Also confirmed: `ArrayTaskId=0-511%128`, s3's `Command=.../slurm/stage3_search.sbatch`
(shipped file, node-local staging intact), s5's `Command=.../slurm/stage_single.sbatch`.

Every `Dependency=` points exactly where intended: s2 → none, s3 → s2, s4 → the whole
s3 array (`_*`), s5 → s4. Not waited on.

---

## Concerns

1. **Stage 4 builds BOTH axes, unavoidably.** `stage4_graph.py` has no `--axis`; it
   writes `graph/fold/` *and* `graph/epitope/` (plus `fp/canonical.parquet`) in one
   pass, and it reads `fp/shard_*.parquet` to do it. Fold-axis-only is enforced at
   **stage 5**, which is where the plan puts it. Two consequences: (a) `--epi-floor 0.6`
   is genuinely load-bearing, not decoration; (b) the later epitope axis needs only a
   stage-5 rerun, no new search. Sized as `results_full` (16 cpu × 8G = 128 GB), which
   cleared `knn_from_vectors` at 500k nodes; union is 526,467, +5%.
2. **Stage 2 is the long pole and the single point of failure.** Dedup of 533,596 hash
   rows + `foldseek createdb` over 99 shards inside 4 h. `results_full` did the same
   shape at 500k, but if it times out everything downstream goes
   `DependencyNeverSatisfied` and the chain must be resubmitted from s2.
3. **`FM` is now env-overridable.** An operator with a stale exported `FM` would submit
   against the wrong tree. Mitigated: the preflight refuses an incomplete tree and
   prints the path it resolved, which is the first line of output.
4. **Node-local scratch headroom is fine.** `results_full/design/db/uniq*` is 339 MB at
   500k nodes, so union ≈ 357 MB against the shipped `--tmp=10000` (10 GB). No change
   needed to `stage3_search.sbatch`.
5. **Not committed.** Git is refused here. `submit_union.sh`, `tests/test_union_submit.py`
   and `sync.sh` are staged in the working tree only.
