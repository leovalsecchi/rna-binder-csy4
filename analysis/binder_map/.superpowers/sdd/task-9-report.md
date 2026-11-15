# Task 9 report — hard gate between `full_s2` and `full_s3`

Date: 2026-09-24
Status: **DONE_WITH_CONCERNS** (concerns are operational, not defects — see the end)

## Files created

| Path | What |
|---|---|
| `/Users/leo/Documents/binder_map/scripts/gate_checks.py` | three correctness checks as pure functions + thin `main()` |
| `/Users/leo/Documents/binder_map/scripts/probe_project.py` | `project_cpu_h()` + `main()` with the CPU-h ceiling |
| `/Users/leo/Documents/binder_map/slurm/gate_full.sbatch` | the gate job: checks, two probe arms, projection |
| `/Users/leo/Documents/binder_map/tests/test_gate.py` | 45 tests covering both modules |

## Files modified

| Path | What |
|---|---|
| `/Users/leo/Documents/binder_map/submit_maps.sh` | `submit_map()` gained a 6th `GATE` parameter; `full` is gated, `lib12k` is not |

Nothing was committed — all changes are left in the working tree, as instructed.
`sync.sh` already covers `scripts/`, `slurm/`, `tests/` and `submit_maps.sh`, so
**no change to `sync.sh` was needed**; all four new files were pushed by the
existing script.

Stages 2, 3 and 5, every clustering parameter, `results_sw_k30_f50_mut` and
`results_v6` were all left untouched. `results_v6` is read from only.

---

## 1. `scripts/gate_checks.py`

Three pure functions over paths, plus a `main()` that only prints verdicts:

- `check_hash_identity(shared_results, ref_results) -> (n_shared, n_mismatch)`
  reads `<tree>/design/hashes/shard_*.tsv` (`usecols=["design_id","coord_hash"]`)
  from both trees, intersects on `design_id`, counts hash disagreements.
  `main()` fails on `n_shared < MIN_SHARED_HASHES` (5000) **or** `n_mismatch > 0`;
  both are reported independently so the log shows every reason it failed.
  A missing tree yields `(0, 0)` → falls under the floor → fails. No exception
  handling that could turn a bad tree into a pass.
- `check_id_format(manifest_tsv, shard_targz) -> (n_members, n_hits)`
  reads the tar in memory, takes each `.pdb` member's `name[:-4]` (the same stem
  rule `stage1b_split.split_tar` uses), and counts how many are manifest
  `design_id`s. Fails on `n_members == 0` or `n_hits != n_members`.
- `check_split_nonempty(tree_results) -> int`
  members of the first `shard_*.tar.gz` by sorted name, or the sentinel
  `NO_SHARDS = -1` when a tree has no shards at all. Fails on either `-1` or `0`,
  and the two print different messages so the operator can tell "the split ran
  and matched nothing" from "the split never ran".

Argument surface is exactly `--shared --ref --manifest --tree` (repeatable), with
`--tree` defaulting to nothing. One addition beyond the brief: **`main()` fails
when no `--tree` was given** ("refusing to pass a gate that checked no trees"),
so a typo in the sbatch cannot produce a green gate that checked two thirds of
what it claims. There is a test for it.

No `try/except` anywhere in the module. A malformed TSV or corrupt tar raises and
the job dies non-zero, which is the intended fail-closed behaviour.

## 2. `scripts/probe_project.py`

`project_cpu_h(t50, t150, n_full=500_000, n_probe_targets=50_000,
n_probe_queries=2000, threads=16) -> (p, cpu_h)`, implementing
`probe_scaling.sh` lines 31-33 verbatim:

```
p      = log(t150 / t50) / log(3)
wall_s = t50 * (n_full / n_probe_targets) ** p * (n_full / n_probe_queries)
cpu_h  = wall_s * threads / 3600
```

It raises `ValueError` on a non-positive or non-finite timing rather than
computing from garbage. `main()` takes `--t50 --t150 --ceiling` (default 6000
CPU-h), validates both timings via `valid_timing()` and returns 1 if either is
0, negative, NaN or inf; a non-numeric string dies at argparse with exit 2.
Every garbage path exits non-zero — none of them can reach exit 0.

## 3. `slurm/gate_full.sbatch`

`--partition=normal.4h --time=04:00:00 --cpus-per-task=16 --mem-per-cpu=4G`,
`set -euo pipefail`, with a comment saying in as many words that `set -e` is load
bearing and that `|| true` must never be added. Order:

1. `PYTHONPATH=src $PY scripts/gate_checks.py --shared … --ref $FM/results_v6
   --manifest $FM/build_shared/manifests/design.tsv --tree $FM/results_full
   --tree $FM/results_lib12k`
2. `module load stack/2024-06 gcc/12.2.0 openmpi/4.1.6 foldseek/8-ef4e960`, then
   the two probe arms.
3. `$PY scripts/probe_project.py --t50 … --t150 … --ceiling ${GATE_CPUH_CEILING:-6000}`

Probe arms reuse `probe_scaling.sh`'s body exactly: `createsubdb` over the
`"" _ss _ca` suffix triple from `$FM/results_full/design/db/uniq*`, `head -n 2000`
queries, and `foldseek search` with the frozen flags
`-c 0.8 --cov-mode 0 -e 10 --max-seqs 300 --alignment-type 1 -a --threads 16`.

Timing is `date +%s` deltas around each `search` (no `/usr/bin/time -v` parsing,
no locale dependence). Everything inside the arm function except the elapsed
integer is redirected to stderr, so the command substitution captures only the
number. After both arms the script rejects any timing that is empty,
non-numeric, or `<= 0` before handing it to the projector.

Two additional fail-closed guards:
- `uniq.lookup` must hold at least 150,000 entries, otherwise the 150k arm would
  silently get a smaller target set and produce a meaningless exponent.
- `command -v foldseek` after the module load.

All work is in `${TMPDIR:-/tmp}/gate_full_${SLURM_JOB_ID:-$$}` with
`trap 'rm -rf "$WORK"' EXIT`. Nothing under `results/` is expanded or written.

## 4. `submit_maps.sh`

`submit_map()` signature is now
`TAG RESULTS CHUNKS THROTTLE DEP GATE=${6:-0}`, with `local … JGATE S3DEP`.
When `GATE=1` it submits `slurm/gate_full.sbatch` with
`--dependency=afterok:$J2` and sets `S3DEP=$JGATE`; `J3` now uses
`--dependency=afterok:$S3DEP`. The status line became
`"$TAG: s2=$J2${JGATE:+ gate=$JGATE} s3=$J3 s4=$J4 s5=$J5" >&2`, so the gate's
job id is printed only for the gated map.

Call sites, verbatim from the brief, with both `|| exit 1` guards preserved (and
a comment recording why they exist, so a future reader does not "simplify" them):

```
J5_LIB=$(submit_map  lib12k results_lib12k 4   4  $J1B 0) || exit 1
J5_FULL=$(submit_map full   results_full   128 64 $J1B 1) || exit 1
```

### Sandbox proof the wiring is right

`sbatch` was stubbed with a file-backed counter (so each fake job gets a distinct
id even inside `$( )`) and the real `submit_map` was extracted from
`submit_maps.sh` with `sed` and sourced:

```
--- lib12k, GATE=0 ---
lib12k: s2=1001 s3=1002 s4=1003 s5=1004
gate submissions: 0
SBATCH[1002] --parsable --job-name=lib12k_s3 --dependency=afterok:1001 --array=0-3%4 --export=ALL,ALN_TYPE=1,RESULTS=results_lib12k slurm/stage3_search.sbatch design
--- full, GATE=1 ---
full: s2=1005 gate=1006 s3=1007 s4=1008 s5=1009
SBATCH[1005] --parsable --job-name=full_
SBATCH[1006] --parsable --job-name=full_gate --dependency=afterok:1005 slurm/gate_full.sbatch
SBATCH[1007] --parsable --job-name=full_s3 --dependency=afterok:1006 --array=0-127%64 --export=ALL,ALN_TYPE=1,RESULTS=results_full slurm/stage3_search.sbatch design
```

`lib12k` submits no gate and `s3` hangs off `s2` (1001). `full` submits the gate
(1006) off `s2` (1005), and `s3` (1007) hangs off **1006, not 1005**.

---

## Verification — commands and real output

### `bash -n` on the two shell files (laptop)

```
$ cd /Users/leo/Documents/binder_map && bash -n slurm/gate_full.sbatch; echo "gate_full.sbatch rc=$?"
gate_full.sbatch rc=0

$ cd /Users/leo/Documents/binder_map && bash -n submit_maps.sh; echo "submit_maps.sh rc=$?"
submit_maps.sh rc=0
```

### `bash sync.sh`

```
$ cd /Users/leo/Documents/binder_map && bash sync.sh
synced -> euler:/cluster/work/wicky/lvalsecchi/rna_binder/analyses/foldseek_map
```

### `pytest tests/test_gate.py -v` (Euler)

```
$ ssh euler 'cd $FM && PYTHONPATH=src:scripts $PY -m pytest tests/test_gate.py -v -p no:warnings'
============================= test session starts ==============================
platform linux -- Python 3.12.8, pytest-8.3.4, pluggy-1.5.0
rootdir: /cluster/work/wicky/lvalsecchi/rna_binder/analyses/foldseek_map
configfile: pyproject.toml
plugins: anyio-4.8.0, hypothesis-6.125.2
collected 45 items

tests/test_gate.py .............................................         [100%]

============================= 45 passed in 25.36s ==============================
```

(`-p no:warnings` only suppresses the ~40 lines of umap/numba `UserWarning`s the
suite already emitted before this task; it changes no result. The per-test `-v`
name list was trimmed here for length — 45 collected, 45 passed, 0 failed.)

### Whole suite (Euler)

The briefed command is `pytest -q`, but `pyproject.toml` already sets
`addopts = "-q"`, so adding another `-q` makes it `-qq` and pytest **suppresses
the count line**:

```
$ ssh euler 'cd $FM && PYTHONPATH=src:scripts $PY -m pytest -q -p no:warnings'
........................................................................ [ 33%]
........................................................................ [ 66%]
................................ssssss............................ssssss [100%]
```

Re-run without the redundant `-q` so the totals are actually printed:

```
$ ssh euler 'cd $FM && PYTHONPATH=src:scripts $PY -m pytest -p no:warnings'
........................................................................ [ 33%]
........................................................................ [ 66%]
................................ssssss............................ssssss [100%]
204 passed, 12 skipped in 70.84s (0:01:10)
```

### Baseline reconciliation — the brief's 158/12/170 is one test stale

I measured the baseline **before** creating any file, on the tree exactly as I
found it:

```
$ ssh euler 'cd $FM && PYTHONPATH=src:scripts $PY -m pytest -p no:warnings'
........................................................................ [ 42%]
...........................................................ssssss....... [ 84%]
.....................ssssss                                              [100%]
159 passed, 12 skipped in 88.81s (0:01:28)
```

That is **159 passed / 12 skipped / 0 failed = 171 collected**, i.e. one more
passing test than the brief's 158 / 12 / 0 = 170.

I confirmed the delta is not a laptop-vs-Euler skew: a dry-run rsync of all four
synced directories reported no differences at all *before* I pushed anything —

```
$ rsync -azn --itemize-changes --exclude='__pycache__' --exclude='*.pyc' tests/ euler:$FM/tests/
(no output)
$ ... src/   -> (no output)
$ ... scripts/ -> (no output)
$ ... slurm/ -> (no output)
```

— so laptop and Euler were byte-identical, and the 159 is a property of the
working tree as handed to me, not of my changes. Skips match the brief exactly
(12), so the difference is a single extra passing test that was added to the tree
after the brief's figure was recorded. I could not pin down which test, because
running git against this repo is refused in this session. **I am not asserting
158; the measured pre-change baseline is 159 passed / 12 skipped.**

The arithmetic then closes exactly:

```
baseline 171 collected (159 passed + 12 skipped)
      +   45 collected (tests/test_gate.py, all passing)
      =  216 collected  ->  204 passed + 12 skipped   ✓
```

204 − 159 = 45 = the full size of the new file, so **no pre-existing test changed
status**: nothing regressed, and none of the 12 skips moved.

### Discrimination — mutation check (Euler)

Tests that pass either way are worse than no gate, so each check's logic was
broken or inverted one at a time in a throwaway copy of the tree and the suite
re-run. Every mutation is caught:

```
[drop MIN_SHARED_HASHES floor]   1 failed, 44 passed
[drop coord_hash mismatch check] 1 failed, 44 passed
[drop id-format stem check]      1 failed, 44 passed
[drop empty-shard check]         1 failed, 44 passed
[drop empty-first-shard check]   1 failed, 44 passed
[drop no-shards check]           1 failed, 44 passed
[invert ceiling comparison]      4 failed, 41 passed
[drop ceiling check]             2 failed, 43 passed
[accept garbage timings]        11 failed, 34 passed
[wrong log base in exponent]     3 failed, 42 passed
```

10 mutations, 10 caught, 0 survivors. The scratch copy and the helper script were
deleted from Euler afterwards.

### Smoke test against real production data (read-only)

The unit tests use synthetic fixtures, so the readers were also pointed at real
trees on Euler — read-only, tars read in memory, nothing written:

```
$ ssh euler '... python -c "from gate_checks import check_hash_identity; ..."'
v6 vs v6 : (55706, 0)
v5 vs v6 : (55606, 0)

$ ssh euler '... python -c "from gate_checks import check_split_nonempty, check_id_format, first_shard; ..."'
first shard: results_v6/design/shards/shard_00.tar.gz
split_nonempty(results_v6): 500
id_format: (500, 500)
```

So: the hash reader handles the real 64-shard layout (`results_v6` has 55,706
design_ids, and all 55,606 it shares with `results_v5` agree on `coord_hash` —
0 mismatches); the shard reader finds 500 `.pdb` members in a real shard; and all
500 stems are `design_id`s in `build_v6/manifests/design.tsv`, confirming the
`name[:-4]` stem rule matches production naming (`<design_id>.pdb`, as
`stage2_builddb.sh`'s own FATAL message documents).

### Not run

`sbatch` was never invoked. Nothing in this task submitted a cluster job.

---

## Concerns

1. **The `n_shared >= 5000` floor is a false-positive risk that I could not
   measure.** `results_v6` holds 55,706 design_ids. The gate needs the 500k
   `results_shared` tree to share at least 5,000 of them. `build_shared/` does
   not exist on Euler yet, so the overlap is unknown — if the 500k sample was
   drawn from the 6.99M master without regard to the v6 passing set, the
   intersection could be small and the gate would block a perfectly good run.
   This fails in the safe direction (stage 3 never starts), but it would cost a
   cycle. **Check the overlap the moment `build_shared/manifests/design.tsv`
   exists, before relying on the gate:**
   ```
   ssh euler 'cd $FM && PYTHONPATH=src:scripts $PY -c "
   import pandas as pd, glob
   a=set(pd.read_csv(\"build_shared/manifests/design.tsv\",sep=\"\t\",usecols=[\"design_id\"])[\"design_id\"])
   b=set()
   for f in glob.glob(\"results_v6/design/hashes/shard_*.tsv\"):
       b|=set(pd.read_csv(f,sep=\"\t\",usecols=[\"design_id\"])[\"design_id\"])
   print(len(a&b))"'
   ```
   If that prints well under 5,000, the reference tree — not the threshold — is
   the thing to reconsider.

2. **Both probe arms now share one 4-hour job.** `probe_scaling.sh` gave each arm
   its own 4h job; the brief specifies `--time=04:00:00` on `normal.4h` for the
   gate, which must fit the hash checks *and* both arms. If the 150k arm is slow
   the gate times out. That is fail-closed by design (stage 3 sits in
   `DependencyNeverSatisfied`), so it is safe, but it costs a full stage-2 cycle
   to discover. `normal.4h` caps at 4h, so raising the limit would mean changing
   partition — out of scope for this task, flagged rather than done.

3. **The gate's end-to-end path has not been exercised against its real
   arguments,** because `results_shared`, `build_shared`, `results_full` and
   `results_lib12k` do not exist on Euler yet. Each individual check was smoke
   tested against `results_v6` / `build_v6`, which have the identical layout, and
   `main()` is covered end-to-end by the tests on fixture trees. The one thing
   genuinely unexercised is `foldseek createsubdb`/`search` inside
   `gate_full.sbatch`, which is copied verbatim from the already-proven
   `probe_scaling.sh`.

4. **Baseline figure differs from the brief** — measured 159 passed, brief says
   158. Reconciled in full above; the extra test predates my changes and I could
   not identify it without git.

---

# Task 9 review-fix report

Date: 2026-09-24
Status: **DONE** — 8 review findings fixed (CRITICAL 1, IMPORTANT 2-4, MINOR 5-8),
MINOR 9 left alone as instructed, plus the requested array widening 128/64 -> 512/256.

## Files touched

| Path | What |
|---|---|
| `/Users/leo/Documents/binder_map/slurm/gate_full.sbatch` | `shopt -s inherit_errexit`; `|| exit 1` on every command in `probe_arm`; sub-DB and `aln.index` output guards; 5 s per-arm timing floor; `#SBATCH --tmp=20000` |
| `/Users/leo/Documents/binder_map/scripts/gate_checks.py` | `MIN_SHARED_HASHES` 5000 -> 3000 with the measured margin in a comment; refuse `--shared` and `--ref` resolving to the same tree |
| `/Users/leo/Documents/binder_map/scripts/probe_project.py` | validate the CPU-h ceiling is finite and positive |
| `/Users/leo/Documents/binder_map/submit_maps.sh` | `case` rejecting any `GATE` that is not `0`/`1`; `full` call site now `512 256` |
| `/Users/leo/Documents/binder_map/tests/test_gate.py` | +8 tests (floor pin, same-tree refusal x2, ceiling validation x5) |
| `/Users/leo/Documents/binder_map/tests/test_gate_shell.py` | **new** — 15 tests driving the real shell functions against stub `foldseek` / `sbatch` |

Nothing committed; git was not run (refused in this session). Stages 2/3/5, every
frozen parameter, `results_v6` and `results_sw_k30_f50_mut` untouched. No `sbatch`.

## CRITICAL 1 — a failed probe arm exited 0 and released the array

Fixed both ways the review asked for:

- `shopt -s inherit_errexit` immediately after `set -euo pipefail`, with a comment
  saying why it is there. Requires bash >= 4.4 (Euler runs 5.1.16(1)); on an older
  shell the `shopt` itself fails and `set -e` kills the job — the closed direction.
- explicit `|| exit 1` on `mkdir`, both `head` calls (the `head | cut` pipeline
  included), all six `foldseek createsubdb` calls and the `foldseek search`, which
  is what the already-proven `slurm/stage3_search.sbatch:33-39` does.

Either fix alone is sufficient; both are in, as instructed.

## IMPORTANT 2 — the arm now checks it produced something

`probe_arm` refuses to return a timing unless:

- `sub_ss.index` and `q_ss.index` are non-empty after `createsubdb`
  (`stage3_search.sbatch:34`'s precedent), and
- `aln.index` is non-empty after `search`, and
- the elapsed time clears a 5 s floor (`local min_s=5`; the reference arms took
  189 s and 570 s, so anything instant did not search).

Evidence that `aln.index` is the right file to test: foldseek DBs are the
`<name>` / `<name>.index` / `<name>.dbtype` triple — visible in
`$FM/results/af3/db/` (`full`, `full.index`, `full.dbtype`, `full_ca.index`, ...) —
and the last step foldseek logs inside `search` is
`mvdb .../tmp/.../aln .../bm_s3_af3_00/aln` (`logs/bm_s3_af3_11554703_0.log:140`),
i.e. it moves exactly such a triple into place.

## IMPORTANT 3 — the shell logic now has tests, and they discriminate

`tests/test_gate_shell.py` extracts the **real** `probe_arm` and `submit_map`
functions from the shipped scripts (plus each script's own `set`/`shopt` preamble
and the real `J5_LIB=` / `J5_FULL=` call-site lines) and runs them against stub
executables on `PATH`. Taking the shell options from the file rather than
restating them is what keeps the tests honest: delete the `shopt` line and the
harness loses it too.

Before/after, both measured on Euler. "Before" is a scratch tree
(`~/gate_review_before/`) holding the **pre-fix** `slurm/gate_full.sbatch` and
`submit_maps.sh` (md5 `959f2c18...` / `d228e2ba...`, copied off Euler before
`sync.sh`) plus the new test file:

```
$ ssh euler 'cd $HOME/gate_review_before && $PY -m pytest tests/test_gate_shell.py -v -p no:warnings --no-header'
...
FAILED tests/test_gate_shell.py::test_probe_arm_aborts_when_foldseek_fails
FAILED tests/test_gate_shell.py::test_probe_arm_aborts_when_the_search_writes_no_alignment_db
FAILED tests/test_gate_shell.py::test_probe_arm_aborts_on_an_impossibly_fast_arm
FAILED tests/test_gate_shell.py::test_the_real_call_sites_are_what_the_script_ships
FAILED tests/test_gate_shell.py::test_submit_map_refuses_a_gate_flag_that_is_not_zero_or_one[true]
FAILED tests/test_gate_shell.py::test_submit_map_refuses_a_gate_flag_that_is_not_zero_or_one[yes]
FAILED tests/test_gate_shell.py::test_submit_map_refuses_a_gate_flag_that_is_not_zero_or_one[01]
FAILED tests/test_gate_shell.py::test_submit_map_refuses_a_gate_flag_that_is_not_zero_or_one[2]
FAILED tests/test_gate_shell.py::test_submit_map_refuses_a_gate_flag_that_is_not_zero_or_one[-1]
========================= 9 failed, 6 passed in 22.23s =========================
```

The CRITICAL-1 test reproduces the finding verbatim against the old script — a
`foldseek` that sleeps 1 s and exits 1 on every call leaves the gate holding a
timing and exiting 0:

```
>       assert proc.returncode != 0, f"gate would have released:\n{proc.stdout}{proc.stderr}"
E       AssertionError: gate would have released:
E         RELEASED t50=1
E         === probe arm N=50000 : 2000 queries against a 50000-target DB ===
E       assert 0 != 0
```

After the fix, same test file, against the shipped scripts:

```
$ ssh euler 'cd $FM && PYTHONPATH=src:scripts $PY -m pytest tests/test_gate_shell.py -o addopts= -v -p no:warnings --no-header'
collected 15 items

tests/test_gate_shell.py::test_the_extractor_finds_the_real_probe_arm PASSED [  6%]
tests/test_gate_shell.py::test_probe_arm_aborts_when_foldseek_fails PASSED [ 13%]
tests/test_gate_shell.py::test_probe_arm_aborts_when_the_search_writes_no_alignment_db PASSED [ 20%]
tests/test_gate_shell.py::test_probe_arm_aborts_on_an_impossibly_fast_arm PASSED [ 26%]
tests/test_gate_shell.py::test_probe_arm_returns_the_elapsed_seconds_when_the_search_really_ran PASSED [ 33%]
tests/test_gate_shell.py::test_the_real_call_sites_are_what_the_script_ships PASSED [ 40%]
tests/test_gate_shell.py::test_only_the_full_map_is_gated_and_its_s3_hangs_off_the_gate PASSED [ 46%]
tests/test_gate_shell.py::test_the_arrays_are_the_widths_the_call_sites_ask_for PASSED [ 53%]
tests/test_gate_shell.py::test_submit_map_refuses_a_gate_flag_that_is_not_zero_or_one[true] PASSED [ 60%]
tests/test_gate_shell.py::test_submit_map_refuses_a_gate_flag_that_is_not_zero_or_one[yes] PASSED [ 66%]
tests/test_gate_shell.py::test_submit_map_refuses_a_gate_flag_that_is_not_zero_or_one[01] PASSED [ 73%]
tests/test_gate_shell.py::test_submit_map_refuses_a_gate_flag_that_is_not_zero_or_one[2] PASSED [ 80%]
tests/test_gate_shell.py::test_submit_map_refuses_a_gate_flag_that_is_not_zero_or_one[-1] PASSED [ 86%]
tests/test_gate_shell.py::test_submit_map_accepts_both_gate_flags[0] PASSED [ 93%]
tests/test_gate_shell.py::test_submit_map_accepts_both_gate_flags[1] PASSED [100%]

============================= 15 passed in 15.87s ==============================
```

The four `probe_arm` tests are one positive control plus three failure modes, so
they cannot all pass for a trivial reason: `slow_ok` (stub sleeps 6 s and writes
`aln.index`) must reach `RELEASED t50=6` and exit 0, while `fail`, `no_output`
and `fast_ok` must not.

`test_only_the_full_map_is_gated_and_its_s3_hangs_off_the_gate` runs the real
call sites under a file-backed-counter `sbatch` stub and asserts `full_s3` carries
`--dependency=afterok:<gate id>` and **not** `afterok:<full_s2 id>`, that
`full_gate` depends on `full_s2` and runs `slurm/gate_full.sbatch`, that no
`lib12k_gate` is ever submitted, and that `lib12k_s3` still hangs off `lib12k_s2`.

## IMPORTANT 4 — the floor, re-measured

The reviewer's numbers are exact; I re-measured them myself, read-only, today:

```
$ ssh euler 'awk -F"\t" "FNR==1{next} FILENAME~/hashes/{v[\$1]=1; next} {m[\$2]=1} \
    END{print \"manifest_ids\", length(m); print \"v6_ids\", length(v); \
        n=0; for(k in m) if(k in v) n++; print \"intersection\", n}" \
    $FM/results_v6/design/hashes/shard_*.tsv $FM/build_shared/manifests/design.tsv'
manifest_ids 508207
v6_ids 55706
intersection 6377
```

`build_shared/` **does** exist on Euler (dated Sep 14) — the previous report's
concern 1 was wrong on that point — and `design_id` is column 2 of
`build_shared/manifests/design.tsv` (header:
`row_id design_id uuid8 stack campaign gen invfold outer_path inner_path member_path`).

`MIN_SHARED_HASHES` is now **3000**, with the measurement, its date and the
reasoning recorded beside the constant. `tests/test_gate.py` gains
`test_the_shared_id_floor_is_three_thousand`, which pins the value and carries the
same numbers; the existing floor tests are written against the constant and still
exercise both sides of the boundary.

## MINOR 5 / 6 / 7 / 8

- **5** `submit_map` now starts with `case "$GATE" in 0|1) ;; *) ... exit 1 ;; esac`.
  Five bad spellings are tested, and the test also asserts **nothing was submitted**
  before the refusal. An empty 6th argument is deliberately *not* rejected:
  `${6:-0}` substitutes the default for unset *or* empty, so `""` legitimately
  means ungated; that is noted at the test.
- **6** `probe_project.main()` rejects a ceiling that is not finite and positive
  before it does anything else, reusing the already-tested `valid_timing`
  predicate. `nan`, `inf`, `0` and `-5` are tested, plus a control that `1200`
  still passes. (`-inf` is not in the list: argparse reads it as a flag and dies
  with exit 2 — non-zero either way, but not through this branch.)
- **7** `gate_checks.main()` refuses when `--shared` and `--ref` resolve to the
  same tree, before comparing anything. Tested with the identical path and with a
  symlink, so the refusal is on the resolved path, not the spelling.
- **8** `#SBATCH --tmp=20000` on the gate, matching
  `slurm/stage1_extract_nested.sbatch`, with a comment naming the July 2026
  scratch-exhaustion incident.

## New change — 512 chunks / 256 concurrent

`J5_FULL=$(submit_map full   results_full   512 256 $J1B 1) || exit 1`.
`lib12k` stays at `4 4`. Both `|| exit 1` guards preserved verbatim and now
pinned by `test_the_real_call_sites_are_what_the_script_ships`.
`scripts/stage3_split_queries.py`, `slurm/stage3_search.sbatch` and
`scripts/stage4_graph.py` were **not** touched. The array geometry that reaches
`sbatch` is pinned by test (`--array=0-511%256`, `--n-chunks 512`). Comments in
`gate_full.sbatch`, `gate_checks.py`, `probe_project.py` and `submit_maps.sh`
that described the array as "128-task" now say 512.

## Verification — commands and real output

```
$ cd /Users/leo/Documents/binder_map && bash -n slurm/gate_full.sbatch; echo "rc=$?"
rc=0
$ cd /Users/leo/Documents/binder_map && bash -n submit_maps.sh; echo "rc=$?"
rc=0
$ ssh euler 'cd $FM && bash -n slurm/gate_full.sbatch; echo "gate_full.sbatch rc=$?"; bash -n submit_maps.sh; echo "submit_maps.sh rc=$?"'
gate_full.sbatch rc=0
submit_maps.sh rc=0

$ cd /Users/leo/Documents/binder_map && bash sync.sh
synced -> euler:/cluster/work/wicky/lvalsecchi/rna_binder/analyses/foldseek_map

$ ssh euler 'cd $FM && PYTHONPATH=src:scripts $PY -m pytest tests/test_gate.py -v -p no:warnings --no-header'
collected 53 items
tests/test_gate.py ..................................................... [100%]
============================= 53 passed in 28.71s ==============================

$ ssh euler 'cd $FM && PYTHONPATH=src:scripts $PY -m pytest -q -p no:warnings'
........................................................................ [ 30%]
........................................................................ [ 60%]
.......................................................ssssss........... [ 90%]
.................ssssss                                                  [100%]

$ ssh euler 'cd $FM && PYTHONPATH=src:scripts $PY -m pytest -p no:warnings'
........................................................................ [ 30%]
........................................................................ [ 60%]
.......................................................ssssss........... [ 90%]
.................ssssss                                                  [100%]
227 passed, 12 skipped in 85.10s (0:01:25)
```

Two quirks worth stating rather than hiding:

- `pyproject.toml` sets `addopts = "-q"`, so the briefed `pytest -q` becomes `-qq`
  and pytest suppresses the count line; the same run without the redundant `-q`
  prints the totals. Likewise `-v` on top of `-q` nets to normal verbosity, which
  is why the per-test names above needed `-o addopts=`.
- `-p no:warnings` only hides pre-existing umap/numba `UserWarning`s; it changes
  no result.

### Suite reconciliation

Baseline given in the review brief: 204 passed / 12 skipped / 0 failed = 216
collected. Measured now: **227 passed / 12 skipped / 0 failed = 239 collected**.

```
216 collected (204 passed + 12 skipped)      baseline
 +  8  tests/test_gate.py   (floor pin 1, same-tree 2, ceiling 5)
 + 15  tests/test_gate_shell.py
 = 239 collected -> 227 passed + 12 skipped   OK
```

227 - 204 = 23 = exactly the tests added, so no pre-existing test changed status
and none of the 12 skips moved. `tests/test_gate.py` went 45 -> 53.

### Not run

`sbatch` was never invoked. `foldseek` was never invoked (the shell tests use a
stub). `git` was not run. `results_v6`, `results_sw_k30_f50_mut` and
`build_shared/` were read only.

## Concerns

1. **The 5 s arm floor is a judgement call, not a measurement.** The reference
   arms took 189 s and 570 s, so 5 s has ~38x of margin, and a real arm can only
   get slower as the DB grows. It fails closed if it ever bites.
2. **`aln.index` is inferred, not observed.** No live `foldseek search` output
   tree survives to look at (stage 3 writes to node-local scratch and deletes it).
   The inference rests on the DB triple visible in `results/af3/db/` and on
   foldseek's own `mvdb ... aln` step in the stage-3 logs. If it is ever wrong the
   gate blocks a good run — costly, but the safe direction.
3. **The gate still shares one 4 h job with both probe arms** (the previous
   report's concern 2, unchanged). Chunking does not change total search work, so
   512 chunks do not affect the probe; the gate itself is still the same 4 h box.
4. **512/256 is not free of scheduler effects.** The width was verified only as
   far as chunk naming (`q%02d` widens identically past 99) and the
   `--array=0-511%256` string that reaches `sbatch`; Euler's `MaxArraySize` and
   the per-user backfill window are exercised by no test here, and nothing was
   submitted, so the first real evidence comes from the run itself.
5. **The design docs still say 128**
   (`docs/superpowers/specs/2026-09-14-campaign-wide-fold-map-design.md:223,278`
   and `docs/superpowers/plans/2026-09-14-campaign-wide-fold-map.md:1952,2062`).
   Left as the historical record of the plan; only executable code and its
   comments were updated.
