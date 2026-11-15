# Task 8 report: Scaling probe, then both runs

## What was implemented

Two new driver scripts at the repo root (synced to Euler's `foldseek_map/` by
`sync.sh`):

- **`probe_scaling.sh`** — submits two `sbatch` jobs (`probe_50000`,
  `probe_150000`) that each build a query/target sub-DB pair from the already-
  built `results_full/design/db/uniq*` Foldseek DB and time one
  `foldseek search` call (2,000 queries against a 50k- or 150k-target DB) with
  `/usr/bin/time -v`. Written **byte-identical** to the brief's code block
  (confirmed by diff — see below).

- **`submit_maps.sh`** — the seven-stage job chain: stage 1 (shared extraction,
  35-way array), stage 1b (split into `results_full` / `results_lib12k`
  trees), then stages 2 through 5 once per map via a `submit_map()` shell
  function, then stage 9 (three-way comparison). Matches the brief's code
  block with one deliberate correction plus two comment lines (see next
  section).

## Deviation from the brief's literal text, and why

The brief's own `submit_maps.sh` code block passes
`--map lib12k:build_lib12k:results_lib12k` (no census tag). The task
instructions explicitly called this out as a fact the brief doesn't carry:
`stage1b_split.py`'s `--map` argument takes an optional 4th `:census`
component, and **"The `lib12k` map MUST be passed with `:census`, or its
shortfall check silently never runs."** I verified this directly by reading
`scripts/stage1b_split.py` (lines 70–97): `should_warn_census()` only fires
when `is_census` is `True`, and `is_census` is only set when the spec's 4th
`:`-separated component equals the literal string `census`. Without the tag,
map 3 (a 9,759-node *complete census*, not a sample) could silently lose
designs during the split and nothing would report it.

I therefore changed that one line to
`--map lib12k:build_lib12k:results_lib12k:census` and added two comment lines
explaining why. This is the **only** functional deviation from the brief's
`submit_maps.sh`. `probe_scaling.sh` required no changes and is identical to
the brief.

Diff of `submit_maps.sh` against the brief's code block:
```
18a19,20
> # lib12k is a census (9,759 nodes, not a sample of a larger population) --
> # it MUST carry the :census tag or stage1b_split.py's shortfall check never runs.
25c27
<         --map lib12k:build_lib12k:results_lib12k'") || exit 1
---
>         --map lib12k:build_lib12k:results_lib12k:census'") || exit 1
```
`probe_scaling.sh` diff against the brief's code block: empty (identical).

## `bash -n` output

Local (`/Users/leo/Documents/binder_map`):
```
probe_scaling.sh: OK
submit_maps.sh: OK
```

On Euler, after `bash sync.sh`:
```
$ ssh euler 'cd /cluster/work/wicky/lvalsecchi/rna_binder/analyses/foldseek_map && bash -n probe_scaling.sh && bash -n submit_maps.sh && echo OK'
OK
```

## Job dependency chain, traced job by job

`submit_maps.sh` submits 11 `sbatch` calls total (1 + 1 + 4 + 4 + 1). Every
`--dependency=afterok:<id>` names a job captured earlier in the script via
`$(sbatch --parsable ...)`; none reference an undefined variable, and `set -u`
plus `|| exit 1` after every `sbatch` call means the script aborts immediately
if any submission fails, rather than continuing with an empty job id.

1. **`cwfm_s1`** (`J1`) — `slurm/stage1_extract_nested.sbatch`, array
   `0-34%35`, `BUILD=build_shared RESULTS=results_shared`. No dependency
   (first job).
2. **`cwfm_s1b`** (`J1B`) — `stage1b_split.py`, `--dependency=afterok:$J1`.
   Depends on (1).
3. **lib12k chain**, all branching off `J1B` (passed as `DEP` into
   `submit_map`):
   - `lib12k_s2` (`J2`) — `--dependency=afterok:$J1B`. Runs
     `stage2_builddb.sh design` then `stage3_split_queries.py --n-chunks 4`.
   - `lib12k_s3` (`J3`) — `--dependency=afterok:$J2`, array `0-3%4`,
     `stage3_search.sbatch design` with `ALN_TYPE=1`.
   - `lib12k_s4` (`J4`) — `--dependency=afterok:$J3`, `stage4_graph.py` with
     the frozen floors/k/mutual flags.
   - `lib12k_s5` (`J5`, returned as `J5_LIB`) — `--dependency=afterok:$J4`,
     `stage5_cluster.py --axis fold`.
4. **full chain**, also branching off `J1B` (independently of the lib12k
   chain — the two per-map chains do **not** depend on each other, only on
   the shared split):
   - `full_s2` (`J2`) — `--dependency=afterok:$J1B`. `stage2_builddb.sh
     design` then `stage3_split_queries.py --n-chunks 128`.
   - `full_s3` (`J3`) — `--dependency=afterok:$J2`, array `0-127%64` (the
     128-task array the task description calls out).
   - `full_s4` (`J4`) — `--dependency=afterok:$J3`.
   - `full_s5` (`J5`, returned as `J5_FULL`) — `--dependency=afterok:$J4`.
5. **`cwfm_s9`** (`J9`) — `--dependency=afterok:$J5_FULL:$J5_LIB`. Depends on
   **both** map chains' final stage-5 jobs, as required — this is the
   three-way comparison and it cannot start until both new maps have finished
   clustering.

`local TAG RESULTS CHUNKS THROTTLE DEP` inside `submit_map()` shadow the
outer scope cleanly; the function's status line (`echo "$TAG: s2=... s5=..."`)
is redirected to stderr (`>&2`) so that `J5_LIB=$(submit_map ...)` and
`J5_FULL=$(submit_map ...)` capture only the final job id from stdout, not the
status line. Verified this by reading the function body directly, not by
running it.

The chain is coherent: no dangling dependency reference, both map chains
branch off the split as siblings (not chained to each other), and the
comparison waits on both.

## Frozen parameters and `:census` — verbatim lines

Grep of both files confirms every frozen value from the brief is present
character-for-character:

```
probe_scaling.sh:28:        -c 0.8 --cov-mode 0 -e 10 --max-seqs 300 --alignment-type 1 -a --threads 16
submit_maps.sh:48:         --layer $LAYER --tm-floor 0.5 --epi-floor 0.6 --k 30 --mutual \
submit_maps.sh:27:        --map lib12k:build_lib12k:results_lib12k:census'") || exit 1
```

`stage3_search.sbatch` (unmodified, not part of this task) itself hardcodes
`-c 0.8 --cov-mode 0 -e 10 --max-seqs 300`, and `submit_maps.sh` passes
`ALN_TYPE=1` via `--export`, so `--alignment-type 1` reaches the actual stage-3
search calls too even though it isn't spelled out again in `submit_maps.sh`
itself.

## Read-only cluster verification performed

- `scontrol show partition normal.4h` — exists, `MaxTime=04:00:00`, matching
  every `--time=04:00:00` (and shorter) request in the scripts.
- `ls` confirmed `results_v6`, `results_sw_k30_f50_mut`, `logs/`, `slurm/`,
  `scripts/` all exist under `$FM`.
- Confirmed `slurm/stage1_extract_nested.sbatch`, `slurm/stage3_search.sbatch`,
  `slurm/stage_single.sbatch`, `scripts/stage1b_split.py`,
  `scripts/stage9_compare.py`, `scripts/stage0_sample.py` all exist on Euler.
- Confirmed `/cluster/work/wicky/lvalsecchi/binder_map_venv/bin/python` exists
  (symlink to `python3`).
- Read `scripts/stage1b_split.py`, `stage2_builddb.sh`,
  `stage3_split_queries.py`, `stage4_graph.py`, `stage5_cluster.py`,
  `stage9_compare.py`, `stage0_sample.py`, `stage_single.sbatch`,
  `stage1_extract_nested.sbatch`, `stage3_search.sbatch` in full to check
  every argument name (`--layer`, `--n-chunks`, `--results`, `--tm-floor`,
  `--epi-floor`, `--k`, `--mutual`, `--axis`) against what `submit_maps.sh`
  passes. All match.
- Ran the existing pytest suite on Euler
  (`PYTHONPATH=src:scripts python -m pytest -q`): exit code 0, progress
  characters were only `.` (pass) and `s` (skip) — no `F` or `E` — across
  ~173 collected tests. No Python source was touched by this task, so this
  was a confirmation, not a fix.
- Did **not** run `sbatch`, `probe_scaling.sh`, or `submit_maps.sh`.

## Things found in the brief worth flagging (not fixed, since I was told to
match the brief's script text and not modify stages 2/3/5)

1. **The `:census` omission** in the brief's `submit_maps.sh` code block —
   already covered above; fixed per the task's explicit instruction.
2. **The verification gates in Steps 4/4b/5/6 are manual, not wired into the
   job dependency graph.** `full_s3` (the expensive 128-task array) is
   `--dependency=afterok:$J2`, i.e. it launches automatically the moment
   `full_s2` (stage 2, `builddb`) succeeds — there's no automated pause for a
   human to run the coord-hash / id-format / TM-spot-check gates first. The
   brief is self-consistent about this: Step 4 says "before `full_s3`
   starts... If it fails, `scancel` the `full_s3`/`s4`/`s5` ids" — i.e. the
   intended safety mechanism is a human racing to `scancel` the just-submitted
   `full_s3` array job during the real wall-clock time `full_s2` (dedup +
   Foldseek `createdb` over the whole 500k set) takes to run, not a hard gate
   in the chain. This is not a bug I should silently fix (Stage 2/3 are
   explicitly off-limits, and restructuring the dependency graph to insert a
   manual-approval stage wasn't asked for), but it's worth the user knowing:
   **whoever runs `submit_maps.sh` needs to watch for `full_s2` finishing and
   be ready to `scancel` `full_s3` quickly if the gates in Steps 4/4b fail.**
3. **`probe_scaling.sh` can only be run meaningfully after `results_full`'s
   stage-2 DB exists** (it reads `results_full/design/db/uniq.lookup` by
   default, and needs at least 150,000 entries there for the `N=150000` probe
   to be a real measurement rather than a truncated `head`). Since
   `submit_maps.sh` auto-chains `full_s3` right after `full_s2`, running the
   probe "before committing a 128-task array" (as the task description frames
   it) means running it in the same narrow window as concern #2 above, and
   `scancel`-ing `full_s3` if the probe's projected CPU-h looks bad. This
   follows directly from the brief's own design; I did not alter it.
4. Step 3's "Expected: seven job ids printed" underclaims: the script actually
   submits 11 `sbatch` calls (1 + 1 + 4 + 4 + 1). This is a doc detail in the
   brief's narrative (Steps 3–7 are instructions for whoever runs the pipeline
   later, not code I authored), not something in the two scripts themselves,
   so I left it as-is and flag it here rather than editing Steps 3-7 prose.

None of the above required changing the scripts beyond the mandated
`:census` fix — they are operational/documentation notes for whoever executes
Steps 3–7 later.

## Files changed

- `/Users/leo/Documents/binder_map/probe_scaling.sh` (new, executable)
- `/Users/leo/Documents/binder_map/submit_maps.sh` (new, executable)

Both were synced to Euler via `bash sync.sh` (pushed to
`/cluster/work/wicky/lvalsecchi/rna_binder/analyses/foldseek_map/`) for the
syntax check; sync.sh is unmodified — it already listed both filenames.

## Self-review

- **Completeness**: both files exist, are executable, syntax-check clean
  locally and on Euler, and contain every frozen parameter the task demanded.
- **Correctness of the chain**: traced above, job by job — no dangling
  dependency, both map chains branch correctly off the shared split, the
  comparison waits on both final stage-5 jobs.
- **No overbuilding**: no extra scripts, no extra stages, no modification to
  stages 2/3/5 sbatch or python files, no test files added (shell scripts,
  as the task expected — "no new pytest coverage is expected"). The only
  content change versus the brief's literal text is the one-line `:census`
  fix plus two comment lines explaining it.
- **Baseline suite**: confirmed green on Euler (exit 0, all pass/skip, no
  fail/error) both before and after adding the two untracked shell files
  (they don't touch anything pytest collects).
- **No Slurm jobs submitted**: confirmed — every cluster interaction was
  `bash -n`, `ls`, `scontrol show partition`, or the existing `pytest` suite.
  `sbatch` was never invoked directly; it only appears as text inside the two
  driver scripts on disk.

## Concerns

- The manual-gate timing window described in point 2 above (racing to
  `scancel full_s3` before it starts consuming real GPU-cluster CPU-hours) is
  a genuine operational risk for whoever runs this pipeline for real. It's
  inherent to the brief's design, not something I introduced, but it's worth
  the user's attention before they run `submit_maps.sh`.
- Everything else is low-risk: the scripts are unmodified drivers matching an
  already-tested pipeline (stages 0–9 were built and tested in Tasks 1–7).
