### Task 9: Hard gate between `full_s2` and `full_s3`

**Not part of the original 8-task plan.** Added 2026-09-24 after the Task 8 review,
by explicit user decision. Everything in the plan's Global Constraints still binds.

## Why

`submit_maps.sh` chains `full_s3` (a 128-task array, ~2,100 CPU-h) directly off
`full_s2` with `--dependency=afterok`. The plan's correctness gates (Steps 4 and 4b)
and its scaling probe are **manual steps a human is supposed to run in the window
between `full_s2` finishing and `full_s3` starting**, aborting with `scancel` if
they fail. That is a race, over a VPN that drops. Worse, Step 4b guards a failure
that is silent: if manifest `design_id`s and shard `.pdb` member stems disagree by
one character, every shard writes a well-formed EMPTY tar, nothing raises, and the
500k map comes out valid-looking and built on nothing.

This task converts that race into a Slurm dependency. `afterok` fails closed: if the
gate exits non-zero, crashes, or times out, `full_s3` never starts (it sits in
`DependencyNeverSatisfied` and Slurm reaps it). That property is the whole design —
**do not add exception handling that could let a failed check exit 0.**

## What to build

### 1. `scripts/gate_checks.py` (new)

Three correctness checks, as **pure functions plus a thin `main()`**, so they are
testable without a cluster. `main()` exits 0 only if all three pass.

- `check_hash_identity(shared_results, ref_results) -> (n_shared, n_mismatch)`
  Reads `<results>/design/hashes/shard_*.tsv` (columns `design_id`, `coord_hash`)
  from both trees, intersects on `design_id`, counts `coord_hash` mismatches.
  FAIL if `n_shared < 5000` (too few to be a real check) or `n_mismatch > 0`.
  Reference tree is `results_v6` — verified present on Euler with 64 hash shards.
  This proves the 2026-09-11 repack did not change what gets extracted.

- `check_id_format(manifest_tsv, shard_targz) -> (n_members, n_hits)`
  Reads `design_id` from the manifest TSV; opens the shard tar.gz and takes every
  `.pdb` member's stem (`m.name[:-4]`). FAIL if the shard has no `.pdb` members at
  all, or if any stem is not a manifest `design_id`.

- `check_split_nonempty(tree_results) -> n_members_first_shard`
  FAIL if a tree has no `shard_*.tar.gz`, or if its first shard holds 0 members.
  Run for both `results_full` and `results_lib12k`.

`main()` argument surface: `--shared`, `--ref`, `--manifest`, `--tree` (repeatable).
Print each check's numbers before its verdict — a gate that fails must say which
check failed and with what counts, because the operator reads only the log.

### 2. `scripts/probe_project.py` (new)

`project_cpu_h(t50, t150, n_full=500_000, n_probe_targets=50_000, n_probe_queries=2000, threads=16) -> (p, cpu_h)`

- `p = log(t150 / t50) / log(3)` — 3 because 150k/50k = 3.
- `wall_s = t50 * (n_full / n_probe_targets) ** p * (n_full / n_probe_queries)`
- `cpu_h = wall_s * threads / 3600`

This is the plan's own formula (`probe_scaling.sh` lines 31-33), now applied instead
of eyeballed. `main()` takes `--t50 --t150 --ceiling` (seconds, seconds, CPU-h),
prints `p` and the projection, and exits non-zero if `cpu_h > ceiling`.

**Fail closed on garbage input:** non-positive or non-finite `t50`/`t150` must exit
non-zero, not sail through. A probe arm that failed to run leaves an empty or absent
timing, and that must stop the chain, not release it.

Default ceiling **6000 CPU-h**, ~3x the plan's 2,100 estimate — a "something is badly
wrong" brake, not a budget. Overridable via `--ceiling`.

### 3. `slurm/gate_full.sbatch` (new)

One job, `--cpus-per-task=16 --time=04:00:00`, partition `normal.4h`. In order:

1. `PYTHONPATH=src $PY scripts/gate_checks.py --shared $FM/results_shared
   --ref $FM/results_v6 --manifest $FM/build_shared/manifests/design.tsv
   --tree $FM/results_full --tree $FM/results_lib12k` — abort on non-zero.
2. `module load stack/2024-06 gcc/12.2.0 openmpi/4.1.6 foldseek/8-ef4e960`, then the
   two probe arms, reusing `probe_scaling.sh`'s exact body: build sub-DBs with
   `createsubdb` over suffixes `"" _ss _ca` from `$FM/results_full/design/db/uniq*`,
   2000 queries, and time `foldseek search` with the **frozen** flags
   `-c 0.8 --cov-mode 0 -e 10 --max-seqs 300 --alignment-type 1 -a --threads 16`.
   Capture each arm's elapsed wall seconds into a shell variable.
3. `$PY scripts/probe_project.py --t50 <t50> --t150 <t150> --ceiling ${GATE_CPUH_CEILING:-6000}`

`set -e` so any step's failure fails the job. Work in `${TMPDIR:-/tmp}` with a
`trap ... EXIT` cleanup — **never expand anything under `results/`** (hard inode cap).

Time the arms with something you can parse reliably; `/usr/bin/time -v` writes to
stderr and its "Elapsed (wall clock)" is `h:mm:ss` or `m:ss`. Prefer capturing
`SECONDS` around each arm, or `date +%s` deltas — simpler to parse and no locale
dependence. Whatever you choose, if a timing cannot be parsed the job must fail.

### 4. `submit_maps.sh` (modify)

Add a 6th parameter to `submit_map()`: `GATE` (`0` or `1`). When `1`, submit
`slurm/gate_full.sbatch` with `--dependency=afterok:$J2` and make `J3` depend on
`afterok:$JGATE` instead of `afterok:$J2`. Call sites become:

```
J5_LIB=$(submit_map  lib12k results_lib12k 4   4  $J1B 0) || exit 1
J5_FULL=$(submit_map full   results_full   128 64 $J1B 1) || exit 1
```

Only the `full` map is gated — the lib12k chain is 4 chunks and cheap, and gating it
would make the gate depend on a tree it also validates. Echo the gate's job id in
`submit_map`'s stderr status line so the operator can find its log.

**Preserve the `|| exit 1` on both call sites.** They were added 2026-09-24 because
`submit_map` runs inside `$( )`, so its internal `|| exit 1` guards kill only the
subshell — without the outer guard a failed submission does not stop the script.
There is a sandbox repro of this; do not remove or "simplify" those guards.

### 5. `tests/test_gate.py` (new)

The gate's value is entirely in whether it fails when it should. Every test must be
**discriminating** — it must fail if the check's logic is inverted or removed. Build
real fixture files (small TSVs and real `tar.gz` archives written with `tarfile` in
`tmp_path`); do not mock the readers.

At minimum:
- hash identity: identical hashes → pass; one flipped `coord_hash` → fail; fewer
  than 5000 shared ids → fail even when every hash matches.
- id format: every stem a manifest id → pass; one stem with a trailing character →
  fail; shard with no `.pdb` members → fail.
- split non-empty: shard with members → pass; empty tar → fail; no shards → fail.
- projection: a (t50, t150) pair under the ceiling → exit 0; one over → non-zero;
  `t50=0`, negative, or NaN → non-zero.

## Verification

Run on Euler, from the laptop, after `bash sync.sh`:

```
ssh euler "cd $FM && PYTHONPATH=src:scripts $PY -m pytest tests/test_gate.py -v"
```

Then the whole suite, to prove nothing regressed:

```
ssh euler "cd $FM && PYTHONPATH=src:scripts $PY -m pytest -q"
```

Report both commands' actual output. The suite was last green at 158 passed /
12 skipped / 0 failed (170 collected) — if your number differs, reconcile it
explicitly rather than asserting a baseline you did not measure.

Also `bash -n slurm/gate_full.sbatch` and `bash -n submit_maps.sh`.

**Do not run `sbatch`.** Nothing in this task submits a cluster job.

## Constraints specific to this task

- **Do not touch stages 2, 3 or 5**, or any clustering parameter. The frozen set is
  `--alignment-type 1`, `-c 0.8 --cov-mode 0 -e 10 --max-seqs 300 -a`,
  `TM = min(qTM,tTM)`, `--tm-floor 0.5`, `--epi-floor 0.6`, `--k 30`, `--mutual`,
  Leiden fixed seed, `cluster_r80`.
- **Map 2 (`results_sw_k30_f50_mut`) is never written to.** `results_v6` is read-only
  here too — the gate compares against it, it does not modify it.
- **Never expand a nested archive to scratch.** Read tar members in memory.
- **Never write loose per-structure files** under `results/`.
- **No Boltz-2 / `b2` anywhere.**
