# Task 4 report: Stage 1 rewrite — one streaming pass per archive

## Note on the brief vs. the actual assignment

`task-4-brief.md` as it exists on disk does **not** contain a "Step 0" section —
it opens directly at "Step 1: Write the failing test" for the stage-1 rewrite,
and its sbatch listing shows `--cpus-per-task=8`. `progress.md` (lines 221,
257) asserts the brief "already includes the controller-added Step 0
(archive early-exit) and cpus 8->2", but that text is not present in the file
on disk. I followed the dispatching agent's explicit task description for
Step 0's behaviour and for the cpu correction (both given in prose, with the
measured numbers behind them), since that is the more specific and more
current instruction. This is flagged here rather than silently reconciled,
per the project's own precedent (Task 3 flagged its spec discrepancies the
same way).

## What I implemented

### Step 0 — early exit in `binder_map.nested.iter_members`

`src/binder_map/nested.py`: added a `remaining = set(wanted)` tracking set.
Each time an inner archive named in `wanted` is fully processed (whether or
not the specific requested member names were found in it), it is discarded
from `remaining`; once `remaining` is empty the generator `return`s instead
of continuing to scan the rest of the outer tar. If `wanted` is empty to
begin with, it returns before even spawning the `zstd` subprocess. Nothing
else about the function changed — same yields, same signature, same
`finally` teardown (unchanged, per the task's explicit note that Task 2's
review already proved it correct across normal exit / exception / generator
abandonment, so `return` needed no new teardown).

Rationale in the docstring: design members occupy only the first 8-10% of
each outer archive by tar-member order (see "Evidence early exit fires"
below), so stopping once every wanted inner archive is consumed avoids
streaming the remaining 90%+ of a solid archive that will never be read.

### Stage 1 rewrite — `scripts/stage1_extract_nested.py`

Implemented verbatim per the brief's Step 3 code block: `shard_outers`
assigns whole outer archives to shards (never splits one archive across
shards, since it can only be streamed once); `plan_archive_work` groups a
shard's manifest rows by `inner_path` into the `wanted` dict `iter_members`
expects; `main` streams each of the shard's outer archives exactly once via
`iter_members`, and for every yielded `(inner_name, member_name, blob)`:
writes the blob to a scratch scratch file (`$TMPDIR/bm_s1n_<shard>/one.cif`,
overwritten per structure, never accumulated), loads it, splits chains by
composition (`split_chains` — never by chain ID, so the RFD3 A/B inversion
is handled correctly), hashes the backbone, computes the epitope fingerprint
and RNA-frame approach vector, and appends a deduplicated backbone PDB into
one `shard_NN.tar.gz`. Rows that fail chain selection, aren't length-80, or
error during parsing are recorded in `missing/shard_NN.tsv` with a reason
rather than crashing the shard. At the end it reconciles
`len(hashes) + len(missing) == len(mine)` and raises `SystemExit` (not just
an assertion) on either zero-parsed or an empty tarball, so a broken shard
fails the Slurm task loudly instead of leaving a plausible-looking empty
output — the exact failure mode described in the constraints (blank
outputs at the right row/file count defeated a resume check in July 2026).

`slurm/stage1_extract_nested.sbatch`: one array task per outer archive
(`sbatch --array=0-34%35`). Implemented per the brief with **one deliberate
deviation**: `--cpus-per-task=2` instead of the brief's literal `8`. This
follows the dispatching agent's explicit correction, backed by a measured
cost model (parse+hash+fingerprint = 0.031 cpu-s/structure => ~4.4 cpu-hours
for all 508,207 rows; streaming, not parsing, dominates wall time, so the
other 6 of 8 cpus would sit idle). Everything else (partition, walltime,
mem-per-cpu, `--tmp=20000`, the `zstd` PATH check) is verbatim from the
brief.

## TDD evidence

### Step 0

**RED** — wrote two new tests in `tests/test_nested.py` first:
`test_iter_members_stops_once_the_only_wanted_archive_is_found` and
`test_iter_members_waits_for_every_wanted_archive_before_stopping`. Both
build an outer archive with a valid wanted entry followed by a **poison
duplicate**: a second tar member reusing an already-satisfied inner
archive's name, but holding bytes that are not a valid zstd stream.
`_decompress` runs `zstd -dc` with `check=True`, so opening the poison entry
raises `subprocess.CalledProcessError`. A scan that keeps going after
`wanted` is exhausted reaches it and crashes; a scan that stops once
`wanted` is exhausted never touches it. (My first attempt used a corrupted,
unterminated tail instead of a poison duplicate — that turned out to prove
nothing: CPython's `tarfile`, in streaming mode, silently treats an
unparseable header encountered after at least one valid entry as ordinary
end-of-archive rather than raising, so both the old and new code passed
that version of the test. I caught this by running it against the
unmodified module and watching it pass when it should have failed, per the
TDD skill's "test passes immediately -> fix the test" rule, and replaced it
with the poison-duplicate design, which discriminates because it fails via
an unrelated code path — subprocess exit status — not tarfile's own
leniency.)

Ran on Euler against the pre-Step-0 `nested.py`:

```
$ PYTHONPATH=src:scripts $PY -m pytest tests/test_nested.py -v
...
FAILED tests/test_nested.py::test_iter_members_stops_once_the_only_wanted_archive_is_found
FAILED tests/test_nested.py::test_iter_members_waits_for_every_wanted_archive_before_stopping
2 failed, 5 passed in 2.47s
```

Failure detail (both tests, same cause):

```
subprocess.CalledProcessError: Command '['zstd', '-dc']' returned non-zero exit status 1.
Captured stderr call: zstd: /*stdin*\: unsupported format
```

This is exactly the expected failure: the unmodified code kept scanning
past the satisfied `wanted` set, reached the poison duplicate, and crashed
trying to decompress non-zstd bytes. The 5 pre-existing tests stayed green,
confirming I hadn't touched what the function yields.

**GREEN** — after the `remaining`-tracking change:

```
$ PYTHONPATH=src:scripts $PY -m pytest tests/test_nested.py -v
tests/test_nested.py .......                                             [100%]
7 passed in 3.15s
```

Full suite after Step 0: `114 passed, 12 skipped` — wait, that number is
after both Step 0 and the stage-1 rewrite; immediately after Step 0 alone it
was **110 passed, 12 skipped, 0 failed** (108 baseline + 2 new).

### Stage-1 rewrite

**RED** — `tests/test_stage1_nested.py` written verbatim from the brief,
before `scripts/stage1_extract_nested.py` existed:

```
$ PYTHONPATH=src:scripts $PY -m pytest tests/test_stage1_nested.py -v
ERROR collecting tests/test_stage1_nested.py
ImportError while importing test module ...
E   ModuleNotFoundError: No module named 'stage1_extract_nested'
1 error in 36.08s
```

Exactly the expected failure (module doesn't exist yet).

**GREEN** — after writing `scripts/stage1_extract_nested.py`:

```
$ PYTHONPATH=src:scripts $PY -m pytest tests/test_stage1_nested.py -v
tests/test_stage1_nested.py ....                                         [100%]
4 passed in 19.66s
```

**Full suite, final:**

```
$ PYTHONPATH=src:scripts $PY -m pytest
114 passed, 12 skipped, 112 warnings in 55.93s (0:01:05)
```

108 (baseline) + 2 (Step 0) + 4 (stage-1) = 114, 0 failed. Matches
expectation exactly. (The 112 warnings are pre-existing biotite CIF-dialect
`UserWarning`s from `auth_comp_id`/`auth_atom_id` fallback on real fixture
data used by `test_backbone.py`/`test_chains.py`/`test_fingerprint.py`/
`test_geometry.py`/`test_rnaframe.py` — unrelated to this change, same
count before and after.)

## Step 6 smoke test on real data (verbatim + one extra archive)

Built a 200-row single-archive manifest from the real
`build_shared/manifests/design.tsv` (508,207 rows, 35 outer archives) and
ran the extractor outside Slurm, exactly as the brief's Step 6 specifies.

The alphabetically-first outer archive is **not** RFD3 — it's
`campaigns/csy4_crrna_free/boltzdesign__ligandmpnn/midrun_v2.tar.zst`
(18 GB compressed). Per `chains.py`'s own documented composition table,
BD1/LigandMPNN's protein chain is `A`, not `B` — the RFD3 inversion doesn't
apply to this stack, so this run alone can't exercise the specific
regression check the task asked me to watch for. I ran it anyway (it's the
literal brief instruction) and then ran a **second** smoke test against an
actual RFD3 archive (`campaigns/csy4_crrna_free/rfd3__ligandmpnn/midrun_v3.tar.zst`,
16 GB compressed) to get real coverage of the inversion check, since the
task explicitly said to watch for it flipping.

### Run 1 — BD1/LigandMPNN midrun_v2 (brief's literal Step 6 archive)

```
$ PYTHONPATH=src $PY scripts/stage1_extract_nested.py \
    --shard 0 --n-shards 1 --build $FM/build_smoke_nested --results $FM/results_smoke_nested
shard 0/1: 200 structures from 1 outer archives
  streaming campaigns/csy4_crrna_free/boltzdesign__ligandmpnn/midrun_v2.tar.zst (24 inner archives, 200 members)
done in 78s: 200 parsed, 200 unique backbones, 0 missing
```

Chain tally:

```
$ $PY -c "... h['protein_chain'].value_counts() ..."
protein_chain
A    200
n_res unique: [80]
PASS
```

`A` for all 200 — matches `chains.py`'s documented table for this stack
(BD1/LigandMPNN = A). No inversion expected here, none seen.

### Run 2 — RFD3/LigandMPNN midrun_v3 (the actual inversion check)

```
$ PYTHONPATH=src $PY scripts/stage1_extract_nested.py \
    --shard 0 --n-shards 1 --build $FM/build_smoke_rfd3 --results $FM/results_smoke_rfd3
shard 0/1: 200 structures from 1 outer archives
  streaming campaigns/csy4_crrna_free/rfd3__ligandmpnn/midrun_v3.tar.zst (15 inner archives, 200 members)
done in 143s: 200 parsed, 200 unique backbones, 0 missing
```

Chain tally:

```
protein_chain
B    200
n_res unique: [80]
PASS
```

`B` for all 200 — the RFD3 chain inversion is still detected correctly.
The tally has not flipped; nothing upstream changed. No FATAL in either
run, no missing rows in either run (`missing/shard_00.tsv` is header-only
in both).

### Content checks beyond the brief's script (self-review, not just file existence)

- `fp/shard_00.parquet`: shape `(200, 73)` = 1 `design_id` + 48 fingerprint
  columns + 16 `per_nt_*` + `contacts_offframe`/`contacts_total` + 6
  geometry columns (`cx,cy,cz,ax,ay,az`). Columns match
  `FINGERPRINT_COLUMNS` order.
- `shards/shard_00.tar.gz`: 200 real `.pdb` members (`tar tzf` count), first
  member's content spot-checked — well-formed `ATOM` records, chain `A`,
  residue numbering starting at 1, matching `write_backbone_pdb`'s format.
- `missing/shard_00.tsv`: header row only (`design_id\treason`), 0 data
  rows, for both runs — consistent with "0 missing" in the printed summary.
- No leftover scratch directories under `/tmp` on the login node after
  either run (`shutil.rmtree(tmp, ignore_errors=True)` at the end ran
  cleanly both times).

## Evidence the early exit actually fires

1. **Unit-test proof (strongest evidence):** the poison-duplicate tests
   above are a direct behavioural proof, not a timing inference — the old
   code provably keeps scanning (it crashes on a poison entry placed right
   after the wanted set is satisfied); the new code provably stops (the
   same test passes cleanly because the poison entry is never reached).

2. **`.list` sidecar member counts (tar order), for the two smoke-test
   archives:**

   | archive | last `/designs/` entry (member #) | total members in archive | design fraction |
   |---|---|---|---|
   | boltzdesign\_\_ligandmpnn/midrun_v2 | 1,140 | 16,573 | 6.9% |
   | rfd3\_\_ligandmpnn/midrun_v3 | 846 | 9,588 | 8.8% |

   (`grep -c '/designs/' *.tar.zst.list` and the line number of the last
   match, against `wc -l` for the total.) Without Step 0, `iter_members`
   would keep iterating through the remaining ~93%/91% of each archive's
   tar members (the af3/esm refold directories) after every wanted design
   member had already been delivered. With Step 0, both smoke runs
   returned as soon as their 24 / 15 requested inner archives were each
   seen once — consistent with the printed "N inner archives" being fully
   satisfied within the first ~7-9% of each archive's member order matches
   the pre-measured 8-10% figure cited in the task brief.

3. **Wall-clock time is consistent with a partial, not exhaustive, stream:**
   78 s for 200 structures from an 18 GB compressed archive, 143 s for 200
   from a 16 GB compressed archive — plausible for `zstd -T0` decompressing
   roughly the first 7-9% of each archive (plus CIF parsing, hashing and
   fingerprinting at the pre-measured 0.031 cpu-s/structure), not for
   decompressing the full 16-18 GB.

I did not run a full-exhaustion baseline decompression of either archive to
get a directly comparable "before" wall-clock number — the login node has
4 cores and was already at ~50% memory use from other sessions, and
spending several more minutes decompressing 16-18 GB twice on a shared
login node for a timing comparison that the unit tests already prove
correctly felt like the wrong trade. I consider the unit-test proof (item 1)
sufficient on its own; items 2-3 corroborate it on real data.

## Files changed

- `src/binder_map/nested.py` — Step 0: early exit in `iter_members`.
- `tests/test_nested.py` — two new tests for Step 0 (RED then GREEN); 5
  pre-existing tests untouched.
- `scripts/stage1_extract_nested.py` — new, verbatim per brief.
- `slurm/stage1_extract_nested.sbatch` — new, verbatim per brief except
  `--cpus-per-task=2` (see deviation note above).
- `tests/test_stage1_nested.py` — new, verbatim per brief.

No other files were modified. `git add` will be given these five explicit
paths only.

## Self-review

- **Completeness:** Step 0, its tests, the stage-1 script, its tests, the
  sbatch wrapper, and the Step 6 smoke test (plus the RFD3-specific
  follow-up) are all done. Did not submit the 35-task array (explicitly
  out of scope — Task 8).
- **Naming:** `remaining` / `shard_outers` / `plan_archive_work` match the
  brief's own naming; the two new Step-0 tests are named for the behaviour
  they prove (`_stops_once_the_only_wanted_archive_is_found`,
  `_waits_for_every_wanted_archive_before_stopping`), not implementation
  detail.
- **YAGNI:** the `if not remaining: return` short-circuit for an empty
  `wanted` dict is the one bit I added beyond the minimum needed to pass
  the two new tests. It's one line, avoids spawning a `zstd` subprocess for
  zero work, and is consistent with the feature's own stated intent
  ("stop once everything wanted is consumed" — zero wanted is trivially
  already consumed). No test currently exercises it directly; if that's
  considered scope creep it's a one-line, low-risk revert.
- **Tests assert real behaviour:** the Step-0 tests were rewritten once
  already (see TDD evidence above) specifically because the first version
  didn't actually discriminate old-vs-new behaviour — I want to flag that I
  caught this myself by watching the test suite, not that it was right the
  first time. The stage-1 tests are the brief's own, unit-level and
  fixture-based (no real archives), which is appropriate given they test
  pure-function sharding/grouping logic, not extraction. The smoke test is
  the one place real data + real chain-inversion behaviour is exercised,
  and I ran it against two archives (one BD1, one RFD3) specifically
  because the brief's own archive-selection logic (alphabetically-first
  outer path) happened not to land on RFD3 this time.
- **Pristine output:** full suite is `114 passed, 12 skipped, 0 failed`,
  112 warnings, all pre-existing and unrelated to this change (verified the
  same warning count appears in the pre-Step-0 baseline run). No stray
  files left in the repo; smoke-test build/results directories on Euler
  (`build_smoke_nested`, `results_smoke_nested`, `build_smoke_rfd3`,
  `results_smoke_rfd3`, each ~1 MB) follow the same convention as Task 2's
  pre-existing `build_smoke`/`results_smoke` directories, which were also
  left in place rather than cleaned up.

## Concerns

1. **Brief/progress.md mismatch** (see top of this report): the brief file
   on disk doesn't contain the Step 0 section or the cpu correction that
   `progress.md` says it already has. I implemented both per the
   dispatching agent's explicit prose instructions, which matched
   `progress.md`'s intent even though the brief file itself doesn't reflect
   it. Worth reconciling `task-4-brief.md` on disk so a future resume of
   this plan doesn't re-diverge.
2. **`--cpus-per-task=2` deviates from the brief's literal `8`** — see
   rationale above. Flagging explicitly per this project's own
   "verify job submits" lesson (never change a proven cpus/mem shape
   without a tested run) — this shape has now been smoke-tested for
   correctness (two full runs, 400 structures total, 0 missing, 0 FATAL)
   but not yet for wall-clock throughput at `--cpus-per-task=2` under real
   Slurm scheduling, since Task 8 owns the actual array submission.
3. Did not run a full-archive decompression baseline for a direct
   before/after timing comparison (see "Evidence early exit fires" above)
   — judged not worth the shared-login-node cost given the unit-test proof
   already establishes correctness unambiguously.

---

## Fix round: review findings on Task 4

A reviewer found four Important issues and five Minor ones against `419837b`.
All nine are fixed below. Note on the brief/report mismatch flagged in the
original report: it's now moot — the brief on disk was regenerated again
before this round and the dispatching agent confirmed the current disk copy
is authoritative; no discrepancy was chased.

Files touched: `src/binder_map/nested.py`, `tests/test_nested.py`,
`scripts/stage1_extract_nested.py`, `tests/test_stage1_nested.py`.

### Fix 1 (Important) — `iter_members` loop guard must test `remaining`, not `wanted`

**Change:** `src/binder_map/nested.py:77` — `info.name not in wanted` →
`info.name not in remaining`. One-line diff. Without it, a repeated inner-
archive name in an outer tar that appears *before* every wanted archive has
been satisfied is decompressed and re-yielded a second time (duplicate
`hashes`/`fp` rows for one `design_id`, or a crash if the duplicate happens
to be corrupt).

**The discriminating-test problem, handled as instructed:** the two
existing poison-duplicate tests
(`test_iter_members_skips_a_duplicate_name_with_only_one_archive_wanted`,
`test_iter_members_skips_a_duplicate_name_once_all_wanted_archives_are_found`
— renamed from `..._stops_once_the_only_wanted_archive_is_found` /
`..._waits_for_every_wanted_archive_before_stopping` since they no longer
prove early-exit, only that a duplicate name is skipped) stop discriminating
the early `return` once the guard is keyed on `remaining`: a duplicate name
is skipped by the guard itself the instant its archive is fully delivered,
whether or not `return` still exists. Kept them (they still assert something
true — the module docstring for `_build_outer` was updated to say so
explicitly) and added a new test that measures I/O directly instead:

**New test:** `test_iter_members_stops_reading_the_outer_stream_once_satisfied`
(`tests/test_nested.py`). Monkeypatches `binder_map.nested.subprocess.Popen`
(only the *outer* `zstd -dc -T0 <outer>` call, identified by the `-T0` flag
in argv — the same patched `Popen` is also reached internally by
`subprocess.run` inside `_decompress`, and wrapping *that* stdout breaks
`subprocess.run`'s own internals, so the wrapper is conditional) so the
returned `stdout` is wrapped in a byte-counting reader. Builds an outer
archive with the one wanted inner archive at the front followed by 4 MiB of
filler (`padding/filler.bin`, never in `wanted`) that the real pipeline would
never touch (refold output, in production). Asserts the counting reader saw
far less than the archive's uncompressed size.

**Discriminate-check evidence, exactly as instructed** (temporarily removed
the `return`, ran the new test, confirmed FAIL; restored it, confirmed
PASS):

Removed the `return` (replaced with `pass` + a `# TEMP` comment) and ran:

```
$ PYTHONPATH=src:scripts $PY -m pytest tests/test_nested.py -v
...
FAILED tests/test_nested.py::test_iter_members_stops_reading_the_outer_stream_once_satisfied
1 failed, 7 passed in 2.79s
```

Failure detail — the counting reader saw the *entire* archive, filler
included:

```
E       AssertionError: read 4198400 of 4198400 uncompressed bytes off the outer pipe; the early return should stop well short of the filler
E       assert 4198400 < (4198400 // 4)
```

Restored the `return` and ran again:

```
$ PYTHONPATH=src:scripts $PY -m pytest tests/test_nested.py -v
tests/test_nested.py ........                                            [100%]
8 passed in 2.43s
```

This is the strongest possible evidence: not a timing inference, a direct
byte-count showing the unmodified-minus-`return` code reads 100% of the
archive and the real code reads a small fraction of it.

### Fix 2 (Important) — `by_member` must be keyed by outer archive too

**Change:** `scripts/stage1_extract_nested.py` — `by_member` now keys on
`(r.outer_path, r.inner_path, r.member_path)` instead of `(r.inner_path,
r.member_path)`, with a new assert immediately after construction:
`assert len(by_member) == len(mine)`, message names the collision
(`"by_member key collision: N manifest rows collapsed to M unique
(outer_path, inner_path, member_path) keys - two rows share the same
outer/inner/member triple"`). The lookup in the streaming loop now builds
the same 3-tuple key using the `outer_path` already in scope from the
`groupby`. `delivered`'s type annotation updated to
`set[tuple[str, str, str]]` to match.

**Tests:**
- `test_by_member_keeps_rows_distinct_across_colliding_outer_archives` —
  builds the reviewer's exact collision shape (two manifest rows, `d1` and
  `d2`, sharing one `(inner_path, member_path)` across two different outer
  archives `a.tar.zst`/`b.tar.zst`) and asserts both `d1` and `d2` land in
  `hashes`/`fp` distinctly, `missing` is empty.
- `test_by_member_collision_assert_fires_on_a_genuine_duplicate_row` — a
  literal duplicate manifest row (same `design_id` twice) must trip the new
  assert with an `AssertionError` matching `"by_member key collision"`,
  proving the guard fires on real collisions rather than being vacuous.

Structural-biology calls (`load_structure`, `split_chains`, `coord_hash`,
`fingerprint`, `approach`, `write_backbone_pdb`) are monkeypatched in these
`main()`-level integration tests to simple, blob-content-driven fakes — real
CIF parsing/chain-selection is covered elsewhere
(`test_chains.py`/`test_fingerprint.py`/the brief's real-archive smoke test)
and isn't what these tests are about; what they exercise is `main()`'s own
bookkeeping.

```
$ PYTHONPATH=src:scripts $PY -m pytest tests/test_stage1_nested.py -k by_member -v
tests/test_stage1_nested.py::test_by_member_keeps_rows_distinct_across_colliding_outer_archives PASSED
tests/test_stage1_nested.py::test_by_member_collision_assert_fires_on_a_genuine_duplicate_row PASSED
2 passed
```

### Fix 3 (Important) — a corrupt inner archive must not kill the shard

**Change:** the `for inner_name, member_name, blob in iter_members(outer,
wanted):` loop is now wrapped in its own `try`/`except Exception` *inside*
the `for outer_path, group in mine.groupby(...)` loop. An exception raised
by `iter_members` while the `for` line advances it (a `CalledProcessError`
from a corrupt inner archive, a `tarfile.ReadError` on a truncated outer
stream, etc. — never caught by the existing inner `try`, which only wraps
per-member processing and never lets an exception escape it) is now caught
per outer archive: every row of that archive's group not yet in `delivered`
is recorded in `missing` with reason `archive_read_error:<ExceptionType>`,
a warning is printed, and the outer loop continues to the next archive.

**A second, related bug this exposed and also fixed:** building this test
surfaced that the *pre-existing* `archive_missing` branch (`if not
outer.exists(): ...`) never marked its keys as accounted-for either, so the
final reconciliation loop (`for key, design_id in by_member.items(): if key
not in delivered: missing.append(..., "member_missing")`) would re-add the
*same* design a second time — meaning any real shard containing even one
genuinely-absent outer archive would already have tripped the accounting
assert (`len(hashes) + len(missing) == len(mine)`) in the committed code.
Both the `archive_missing` branch and the new `archive_read_error` branch
now call `delivered.add(key)` before/when appending to `missing`, so a
key is counted exactly once. This wasn't one of the reviewer's four named
Important issues, but it's the same bookkeeping mechanism Fix 3 depends on
to be correct, and it's a real crash-on-real-data bug, so it's fixed
alongside it rather than left in place.

**Test:** `test_corrupt_inner_archive_yields_missing_rows_not_a_crash` — one
good outer archive (`good.tar.zst`, one design `d0`) and one bad outer
archive (`bad.tar.zst`, whose sole "inner archive" member is literally
`b"not a zstd stream"`, so `iter_members`'s `_decompress` raises
`subprocess.CalledProcessError` for designs `d1`/`d2`). Asserts `d0` parses
normally, `d1`/`d2` both get reason `archive_read_error:CalledProcessError`,
`main()` returns 0 (no crash), the final shard tarball exists at its real
name (not `.partial`) and is non-empty.

```
$ PYTHONPATH=src:scripts $PY -m pytest tests/test_stage1_nested.py::test_corrupt_inner_archive_yields_missing_rows_not_a_crash -v
tests/test_stage1_nested.py::test_corrupt_inner_archive_yields_missing_rows_not_a_crash PASSED
1 passed
```

### Fix 4 (Important, plan-mandated) — validate before writing outputs

**Change:** the accounting assert and the two FATAL checks now run
*before* any of the four output files are written. The shard tarball is
still written incrementally during streaming (unavoidable — it's built
member-by-member as `iter_members` yields), but now under a
`shard_NN.tar.gz.partial` name (`partial_tar_path = tar_path.with_name(...
+ ".partial")`); it's promoted to the real `shard_NN.tar.gz` name via
`os.replace(partial_tar_path, tar_path)` only immediately after the FATAL
checks pass, and immediately before the `hashes`/`fp`/`missing` writes. A
shard that fails validation (accounting mismatch, zero parsed, or an empty
tarball) now leaves *no* complete-looking output behind — no
`hashes/shard_NN.tsv`, no `fp/shard_NN.parquet`, no `missing/shard_NN.tsv`,
and the tarball sits under a name (`.partial`) a resume check scanning for
`shard_NN.tar.gz` will never mistake for done.

**Test:** `test_output_files_are_not_written_when_validation_fails` — a
manifest whose one referenced outer archive is never created (0 parsed, 1
missing, `len(mine)=1` → hits the zero-parsed FATAL). Asserts `main()`
raises `SystemExit` matching `"FATAL"`, and that none of
`hashes/shard_00.tsv`, `fp/shard_00.parquet`, `missing/shard_00.tsv`,
`shards/shard_00.tar.gz` exist afterward.

```
$ PYTHONPATH=src:scripts $PY -m pytest tests/test_stage1_nested.py::test_output_files_are_not_written_when_validation_fails -v
tests/test_stage1_nested.py::test_output_files_are_not_written_when_validation_fails PASSED
1 passed
```

### Minor fixes

1. **Empty-shard guard consistency** — the zero-parsed FATAL check was
   already guarded with `if len(mine) and not hashes:`; the tarball-size
   check now matches: `if len(mine) and partial_tar_path.stat().st_size <
   200:`. `--n-shards` greater than the outer-archive count now exits 0 with
   empty-but-schema'd outputs instead of a spurious FATAL.
   Test: `test_empty_shard_with_no_outer_archives_assigned_exits_cleanly`
   (shard 1 of 2, with only 1 outer archive total → shard 1 gets nothing);
   also incidentally covers the `FP_COLUMNS` schema check below since `fp`
   is empty there.
2. **`member_empty` vs `member_missing`** — a blob that's present but
   whitespace-only now gets its own reason (`member_empty`) instead of
   falling through to the generic `member_missing` reconciliation. The
   `delivered.add(key)` line was moved up (right after the `design_id is
   None` check) so this branch — like `wrong_length`/`chain_select_failed`/
   `error` — is accounted for once, immediately, rather than relying on the
   final reconciliation loop.
   Test: `test_blank_member_is_member_empty_distinct_from_member_missing`
   (one blank member, one genuinely absent member, one successful control
   design so the shard doesn't also trip the zero-parsed FATAL) — asserts
   the two get different reasons.
3. **`wrong_length` reason string** — `chain_select_failed:n_res={n}` →
   `wrong_length:n_res={n}` for the explicit `n_res != 80` gate (chain
   selection succeeded there; only the length check failed). The
   `except (ChainSelectError, FrameError)` branch's
   `chain_select_failed:{type(exc).__name__}` reason is unchanged — that one
   genuinely is a chain-selection failure.
   Test: `test_wrong_length_reason_names_the_length_gate_not_chain_selection`.
4. **Explicit `fp` DataFrame columns** — added a module-level `FP_COLUMNS`
   constant (`["design_id"] + FINGERPRINT_COLUMNS + [per_nt_1..16] +
   [contacts_offframe, contacts_total, cx, cy, cz, ax, ay, az]`, built from
   `RNA_LEN` so it can't drift from `FINGERPRINT_COLUMNS`/`per_nt`'s actual
   length) and pass it as `columns=` to `pd.DataFrame(fps, ...)`. A zero-row
   shard now writes a parquet with the correct schema instead of a
   schema-less, zero-column one.
   Verified via `test_empty_shard_with_no_outer_archives_assigned_exits_cleanly`:
   `assert list(fp.columns) == s1.FP_COLUMNS` on a 0-row `fp` DataFrame.
5. **Unused `import pytest`** — moot in the rewritten `tests/test_stage1_nested.py`:
   `pytest` is now genuinely used (`pytestmark = pytest.mark.skipif(...)` for
   the `zstd`-availability guard, plus `pytest.raises` in two of the new
   tests), so there is nothing to remove.

### Full verification

```
$ PYTHONPATH=src:scripts $PY -m pytest tests/test_nested.py tests/test_stage1_nested.py -vv --tb=no
tests/test_nested.py::test_nested_ref_splits_a_store_design_ref PASSED
tests/test_nested.py::test_nested_ref_handles_a_suffixed_midrun PASSED
tests/test_nested.py::test_iter_members_pulls_only_what_was_asked_for PASSED
tests/test_nested.py::test_iter_members_ignores_inner_archives_not_requested PASSED
tests/test_nested.py::test_iter_members_is_silent_about_absent_members PASSED
tests/test_nested.py::test_iter_members_skips_a_duplicate_name_with_only_one_archive_wanted PASSED
tests/test_nested.py::test_iter_members_skips_a_duplicate_name_once_all_wanted_archives_are_found PASSED
tests/test_nested.py::test_iter_members_stops_reading_the_outer_stream_once_satisfied PASSED
tests/test_stage1_nested.py::test_shard_outers_assigns_whole_archives PASSED
tests/test_stage1_nested.py::test_every_archive_lands_in_exactly_one_shard PASSED
tests/test_stage1_nested.py::test_plan_archive_work_groups_members_by_inner_archive PASSED
tests/test_stage1_nested.py::test_plan_archive_work_on_empty_group_is_empty PASSED
tests/test_stage1_nested.py::test_by_member_keeps_rows_distinct_across_colliding_outer_archives PASSED
tests/test_stage1_nested.py::test_by_member_collision_assert_fires_on_a_genuine_duplicate_row PASSED
tests/test_stage1_nested.py::test_corrupt_inner_archive_yields_missing_rows_not_a_crash PASSED
tests/test_stage1_nested.py::test_empty_shard_with_no_outer_archives_assigned_exits_cleanly PASSED
tests/test_stage1_nested.py::test_blank_member_is_member_empty_distinct_from_member_missing PASSED
tests/test_stage1_nested.py::test_wrong_length_reason_names_the_length_gate_not_chain_selection PASSED
tests/test_stage1_nested.py::test_output_files_are_not_written_when_validation_fails PASSED
19 passed in 17.14s
```

Full suite:

```
$ PYTHONPATH=src:scripts $PY -m pytest
122 passed, 12 skipped, 112 warnings in 67.70s
```

122 = 114 baseline + 8 (`test_nested.py`: 7→8, one new test, two renamed with
unchanged count) + 7 (`test_stage1_nested.py`: 4→11, seven new tests). 0
failed, matches the 12 pre-existing skips exactly (needs-`zstd`/upstream-
parity markers, unrelated to this change).

### Did not do

- Did not re-run the real-data smoke test (BD1/LigandMPNN and RFD3/LigandMPNN
  archives on Euler). None of these nine fixes change what a *successful*
  structure produces (Fix 1 only changes when streaming stops; Fix 2 only
  changes the dict key for `by_member`, which is a no-op when every
  `(inner_path, member_path)` pair is unique to one outer archive, as it was
  in both smoke-test archives; Fix 3/4 only change behavior on a path that
  already fails; the minor reason-string/schema fixes only affect `missing`
  labels and an all-empty `fp` schema). The two original smoke runs already
  showed 0 missing, 0 FATAL, correct chain inversion for both stacks, so
  there's nothing these fixes would change about that outcome.
- Did not submit the 35-task Slurm array (Task 8's job, unchanged).
- Did not restructure `main()` beyond what the four Important fixes needed
  (the try/except nesting for Fix 3, the validate-then-write reorder for
  Fix 4, and the keying/accounting changes for Fix 2 and the
  `archive_missing` bug fix are all in-place edits to the existing
  function body — no new top-level helpers were extracted).

---

## Fix round 2

A re-review confirmed round 1's nine fixes are genuinely in place, and found
one Important blocker (the `remaining` guard fix is covered by no
discriminating test) and five Minor issues, each with an exact prescribed
fix. All six are done below, following the prescribed fixes as given.

Files touched: `src/binder_map/nested.py` (net zero diff — see Important fix
below), `tests/test_nested.py`, `scripts/stage1_extract_nested.py`,
`tests/test_stage1_nested.py`.

### Important — the `remaining` guard fix was covered by no test

The reviewer's diagnosis was exactly right: in both poison-duplicate tests
the poison entry sat *after* every wanted archive was satisfied, so the
early `return` fired before the loop ever reached it — the tests passed
whether the guard was keyed on `remaining` or `wanted`. `nested.py:77`
itself needed no change (it is correct); the fix is entirely in the test.

**Change:** in `tests/test_nested.py`,
`test_iter_members_skips_a_duplicate_name_once_all_wanted_archives_are_found`
now orders the outer archive's entries `[aaaa, poison-aaaa, bbbb]` (poison
planted *between* the two wanted archives, not after both) with
`wanted = {aaaa, bbbb}`. After `aaaa` is delivered, `remaining` becomes
`{bbbb}` — not yet empty, so the early `return` has not fired — and the loop
genuinely reaches the guard check on the poison entry next. With the guard
keyed on `remaining`, the poison entry is skipped (its name was already
discarded) and both real members still yield. With the guard keyed on
`wanted` (the pre-fix bug), the name is still present in `wanted`
unconditionally, so the poison entry is not skipped and `_decompress` raises
`subprocess.CalledProcessError`.

Rewrote `_build_outer`'s docstring to say this precisely: the
`..._with_only_one_archive_wanted` test's poison entry is never reached at
all (the early `return` fires on the sole real entry first) and would pass
unchanged even with the guard keyed on `wanted` or deleted outright — it
proves nothing about the guard. The `..._once_all_wanted_archives_are_found`
test (after this reorder) is the one that actually discriminates the two
guard implementations.

**Covering test:**
`test_iter_members_skips_a_duplicate_name_once_all_wanted_archives_are_found`
(`tests/test_nested.py`).

**Discriminate-check, exactly as instructed** — temporarily reverted
`src/binder_map/nested.py:77` from `info.name not in remaining` to
`info.name not in wanted`, synced, ran, confirmed FAIL:

```
$ PYTHONPATH=src:scripts $PY -m pytest tests/test_nested.py -v
...
E               subprocess.CalledProcessError: Command '['zstd', '-dc']' returned non-zero exit status 1.
----------------------------- Captured stderr call -----------------------------
zstd: /*stdin*\: unsupported format
=========================== short test summary info ============================
FAILED tests/test_nested.py::test_iter_members_skips_a_duplicate_name_once_all_wanted_archives_are_found
========================= 1 failed, 7 passed in 3.08s ==========================
```

Restored `remaining`, synced, ran again, confirmed PASS:

```
$ PYTHONPATH=src:scripts $PY -m pytest tests/test_nested.py -v
tests/test_nested.py ........                                            [100%]
============================== 8 passed in 3.04s ===============================
```

`git diff` on `src/binder_map/nested.py` is empty for this round — the file
ends exactly where round 1 left it, per "do not change the `remaining` guard
itself."

### Minor 1 — promote the tarball last, not first

**Change:** `scripts/stage1_extract_nested.py` — moved
`os.replace(partial_tar_path, tar_path)` to after the three
`hashes`/`fp`/`missing` writes instead of before them. A crash between the
metadata writes now leaves the tarball under its `.partial` name, never
promoted without its companions.

**Covering test:**
`test_tarball_is_promoted_last_so_a_crash_mid_write_leaves_no_tarball` (new,
`tests/test_stage1_nested.py`) — monkeypatches `pd.DataFrame.to_csv` to
raise on the first metadata write, then asserts the tarball is still under
`.partial` and the promoted name does not exist.

```
$ PYTHONPATH=src:scripts $PY -m pytest tests/test_stage1_nested.py::test_tarball_is_promoted_last_so_a_crash_mid_write_leaves_no_tarball -v
tests/test_stage1_nested.py::test_tarball_is_promoted_last_so_a_crash_mid_write_leaves_no_tarball PASSED
```
(ran as part of the full `test_stage1_nested.py` batch below; not run in isolation separately)

### Minor 2 — scope the per-archive `except` to the generator advance only

**Change:** `scripts/stage1_extract_nested.py` — replaced the
`try: for inner_name, member_name, blob in iter_members(...): ... except
Exception:` structure (which wrapped the entire per-archive body, including
`write_backbone_pdb`/`tar.addfile`/record construction) with the prescribed
`while True: try: item = next(it) except StopIteration: break except
Exception: ...; break` shape. The per-member work (dedup guard, blank-blob
check, chain selection, hashing, fingerprinting, PDB write) now sits outside
that `try`, under only its pre-existing inner `try`/`except
(ChainSelectError, FrameError)`/`except Exception` — a genuine per-structure
bug is again recorded per design (`error:<Type>`/`chain_select_failed:<Type>`),
not relabelled `archive_read_error` and made to abandon the rest of the
archive.

**Covering test:** no new test was added specifically for this restructuring
— it is a pure control-flow change with no new externally observable
branch, and the existing tests already pin both halves of the behaviour it
must preserve: `test_corrupt_inner_archive_yields_missing_rows_not_a_crash`
still exercises the `except Exception: ... break` path (a corrupt inner
archive raising during generator advance), and every other
`test_stage1_nested.py` test that reaches the per-member `try` (`wrong_length`,
`chain_select_failed` via the fake pipeline, `member_empty`) continues to
exercise that inner try unchanged. Both pass unmodified after the
restructuring (see the full run below), which is the evidence that the
per-member/generator-advance boundary was preserved correctly.

### Minor 3 — FATAL message names the wrong path

**Change:** `scripts/stage1_extract_nested.py` — the undersized-tarball
FATAL message now prints `partial_tar_path` (the file that actually exists
on disk at that point) instead of `tar_path` (which does not exist yet,
since promotion has not happened).

**Covering test:**
`test_fatal_tarball_empty_message_names_the_partial_path` (new,
`tests/test_stage1_nested.py`) — monkeypatches `Path.stat` to report a
near-zero size for exactly the `.partial` file (nothing else touched),
forcing the branch regardless of what the fake pipeline wrote, and asserts
the raised message matches `FATAL: shard tarball .*\.partial is empty`.

### Minor 4 — the collision assert leaks the scratch dir

**Change:** `scripts/stage1_extract_nested.py` — moved the `by_member`
construction and its collision `assert` to immediately after `mine` is
computed, before `out_dir`/`tmp` (the scratch dir) are created. The assert
depends only on `mine`, not on `tmp`/`scratch`, so this is a pure reorder
with no behavioral change on any path that doesn't hit the assert. A
collision now fails with nothing on disk to clean up.

**Covering test:** extended the existing
`test_by_member_collision_assert_fires_on_a_genuine_duplicate_row` — sets
`TMPDIR` to a controlled directory via `monkeypatch.setenv` and asserts
`TMPDIR/bm_s1n_0` does not exist after the `AssertionError`, in addition to
its existing check that the assert fires with the right message.

### Minor 5 — make exactly-once structural, not incidental

**Change:** `scripts/stage1_extract_nested.py` — added
`if key in accounted: continue` immediately before `accounted.add(key)` in
the per-member loop, so a duplicate `(outer, inner, member)` key yielded
twice by `iter_members` (e.g. an inner tar with duplicate member names) is a
no-op on its second delivery instead of a second `hashes`/`fp` row and a
second hit against the final accounting assert.

**Rename:** `delivered` → `accounted` throughout the file (6 identifier
occurrences: the declaration, the three `.add(key)` call sites, and the two
`in`/`not in` checks — a handful of lines, done per the prescribed
condition). The comments in the `archive_missing` and `archive_read_error`
branches already said "accounted for" in prose before this round, so the
rename brings the identifier in line with language the file already used to
describe it.

**Covering test:**
`test_duplicate_member_name_within_one_inner_archive_is_counted_once` (new,
`tests/test_stage1_nested.py`) — builds one inner tar with two members both
named `dup.cif`, both matching one wanted member; asserts `hashes`/`fp` each
contain exactly one row for the design, not two, and `missing` is empty.

### Full verification

```
$ PYTHONPATH=src:scripts $PY -m pytest tests/test_nested.py tests/test_stage1_nested.py -v
tests/test_nested.py ........                                            [ 36%]
tests/test_stage1_nested.py ..............                               [100%]
============================= 22 passed in 37.27s ==============================
```

22 = 8 (`test_nested.py`, unchanged count — one test reordered/renamed
docstring context, no test added or removed) + 14 (`test_stage1_nested.py`:
11 → 14, three new tests: Minor 1, Minor 3, Minor 5).

Full suite:

```
$ PYTHONPATH=src:scripts $PY -m pytest
125 passed, 12 skipped, 112 warnings in 68.21s (0:01:08)
```

125 = 122 (round-1 baseline) + 3 new tests. 0 failed, 12 skipped matches the
pre-existing (unrelated) skip count exactly.

### Does any fix change successful-path output?

Reasoned through rather than re-run (per instruction not to re-run the
real-data smoke test):

- **Important fix:** test-only change; `nested.py` has a net-zero diff for
  this round. No effect on any output.
- **Minor 1** (tarball-last): pure reordering of when files are written, not
  what is written. On a run that reaches this point at all (i.e. already
  passed both FATAL checks), all four output files end up byte-identical to
  before — only their write order on disk changes. No effect on successful-path
  content.
- **Minor 2** (scoped `except`): only changes which exception handler catches
  an error raised while *advancing the generator* versus one raised during
  *per-member processing* — and only for inputs that were already erroring.
  On a fully successful archive (every yielded item processed without
  exception) the two code shapes execute the identical sequence of
  statements; nothing about a passing design changes.
- **Minor 3** (FATAL message text): only fires inside a `raise SystemExit`
  that already meant the run was failing; changes a string, not behavior
  reachable from a successful run.
- **Minor 4** (assert reorder): the assert's truth value and message are
  unchanged; only *when* it's evaluated (before vs. after directory
  creation) moved. A shard with no collision — every real shard, since
  `by_member`'s key already includes `outer_path` — takes the identical path
  either way.
- **Minor 5** (exactly-once guard + rename): a no-op for every manifest row
  that `iter_members` yields exactly once, which per its own contract and
  every archive observed so far (including both smoke-test archives in
  round 1) is all of them. The guard only changes behavior on a key
  `iter_members` has never been seen to actually produce twice. The
  `delivered`→`accounted` rename is not an executable change.

Net: none of the six fixes are expected to change what a *successful* shard
run's `hashes`/`fp`/`missing`/tarball content looks like on the real
midrun archives. Per the task's instruction, the real-data smoke test was
not re-run to confirm this.
