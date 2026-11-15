# Task 11 + 12 report

Branch `foldseek-binder-map`. Two commits, in order:

- `0c3ccc6` — Add stage 0: per-layer manifests and curated explorer metadata (Task 11)
- `8781a76` — Add RNA-variant census as a hard gate before compute (Task 12)

## Task 11 — stage 0 manifests and explorer metadata

Files created:
- `/Users/leo/Documents/binder_map/scripts/stage0_select.py`
- `/Users/leo/Documents/binder_map/tests/test_stage0.py`

Both transcribed verbatim from `.superpowers/sdd/task-11-brief.md`. No bugs found
in the brief's code — checked the real `passing.csv` against every assumption
the code makes (see below) and all held.

### Step 1/2 — write test, confirm it fails

```
$ /Users/leo/.venvs/binder_map/bin/python -m pytest tests/test_stage0.py -v
...
ERROR collecting tests/test_stage0.py
ModuleNotFoundError: No module named 'stage0_select'
1 error in 0.94s
```

Matches the brief's expected failure exactly.

### Step 3 — implement `scripts/stage0_select.py`

Transcribed verbatim. Before transcribing, verified against the real CSV:

```
$ head -3 /Users/leo/Documents/binder_filter/exports/passing.csv
# binder_filter export  n=71866  store_n=6821957
# filter={"tiers": ...}
design_id,uuid8,length,stack,campaign,gen,invfold,...,design_cif,af3_cif,esm_cif
```

- All 24 `METADATA_COLUMNS` exist in the CSV header (checked by listing all 65
  columns and diffing).
- `design_cif`/`af3_cif`/`esm_cif` all contain exactly one `!` per row (0 rows
  missing a `!`, 0 rows with more than one) across all 71,866 rows, so the
  `str.split("!", n=1, expand=True)` in `build_manifest` never produces a
  `None` member_path.
- `stack` values: `RFD3+LiMPNN` 37,951, `BG+LiMPNN` 19,162, `BG+BIF` 9,881,
  `BD1+LiMPNN` 3,344, `HD+LiMPNN` 1,528 — matches the test's expected 5-stack
  set exactly.
- HD `design_id`s (`csy4free_HD_LiMPNN_len80_<uuid8>_recycleN_sampleM_packK`)
  have zero matches for `_seq\d+$` — confirms the brief's claim that HD ids
  carry no `_seqN` suffix.

### Step 4 — tests pass

```
$ /Users/leo/.venvs/binder_map/bin/python -m pytest tests/test_stage0.py -v
tests/test_stage0.py::test_reads_all_rows_past_the_two_comment_lines PASSED
tests/test_stage0.py::test_expected_stacks_present PASSED
tests/test_stage0.py::test_manifest_splits_archive_and_member PASSED
tests/test_stage0.py::test_manifest_is_sorted_by_archive PASSED
tests/test_stage0.py::test_metadata_columns_all_exist PASSED
tests/test_stage0.py::test_hd_designs_have_no_seq_suffix PASSED
6 passed in 2.28s
```

Full suite (no regressions):

```
$ /Users/leo/.venvs/binder_map/bin/python -m pytest
100 passed, 117 warnings in 45.63s
```

94 pre-existing + 6 new = 100. Warnings are the expected `auth_comp_id`/
`auth_atom_id` fallback `UserWarning`s from biotite plus one UMAP
`precomputed_knn` warning already present before this task; none suppressed.

### Step 5 — run stage 0 for real

```
$ /Users/leo/.venvs/binder_map/bin/python scripts/stage0_select.py --out build
read 71,866 rows from /Users/leo/Documents/binder_filter/exports/passing.csv
  design   71,866 rows,  8,931 archives -> build/manifests/design.tsv
  af3      71,866 rows,  8,931 archives -> build/manifests/af3.tsv
  esm      71,866 rows,  8,932 archives -> build/manifests/esm.tsv
wrote metadata for 71,866 designs

$ cat build/meta/counts.json
{
  "input_csv": "/Users/leo/Documents/binder_filter/exports/passing.csv",
  "input_sha256": "39ea0f79aef924386dcfd5ceda5085ff332218a71bf1d495ff0fcb2da323d873",
  "n_rows": 71866,
  "subsample_per_stack": null,
  "layer_counts": {"design": 71866, "af3": 71866, "esm": 71866},
  "stack_counts": {
    "RFD3+LiMPNN": 37951, "BG+LiMPNN": 19162, "BG+BIF": 9881,
    "BD1+LiMPNN": 3344, "HD+LiMPNN": 1528
  }
}
```

**Deviation from the brief's expectation:** the brief estimated "roughly
6,780 archives for design and similar for the refolds." The real count is
8,931 (design/af3) and 8,932 (esm) archives — about 32% higher than the
brief's guess. This is not a test assertion (no test checks archive count),
just a stated expectation, and it is a property of the real data (one archive
per `uuid8`, and there are more distinct uuids in the passing set than the
brief anticipated). No code change needed; noted here as the one place real
output differs from the brief's prose.

### Step 6 — commit

```
$ git add scripts/stage0_select.py tests/test_stage0.py
$ git commit -m "Add stage 0: per-layer manifests and curated explorer metadata"
[foldseek-binder-map 0c3ccc6] Add stage 0: per-layer manifests and curated explorer metadata
 2 files changed, 174 insertions(+)
```

## Task 12 — RNA-variant census (the frame gate)

Files created:
- `/Users/leo/Documents/binder_map/scripts/stage0b_rna_census.py`
- `/Users/leo/Documents/binder_map/slurm/stage_single.sbatch`

Both transcribed verbatim from `.superpowers/sdd/task-12-brief.md`.

### Bug review before running

Read the census script against `binder_map.chains` / `binder_map.rnaframe`
(both already implemented and tested in Tasks 1–10) line by line:

- `extract_member` shells out to `tar --zstd -xOf <RB>/<archive> <member>` —
  matches the stated constraint that this Python's `tarfile` has no zstd
  codec.
- `load_structure(tmp)` is given a bare `bm_census.cif` (no `.gz` suffix), and
  `tar -xOf` already yields decompressed bytes, so the gzip-detection branch
  in `load_structure` is correctly skipped.
- Exception handling only catches `(FileNotFoundError, ChainSelectError)`,
  which matches the documented exception surface of `chains.py`.
- `io` and `tarfile` are imported but unused (the extraction path went
  through the `tar` CLI instead of Python's `tarfile`, per the docstring's
  own explanation). This is a harmless leftover, not a functional bug, and
  the instruction was to transcribe verbatim absent a *real* bug — left as
  written.
- No other logic bugs found. Ran it for real rather than guessing further.

### Step 1 — `slurm/stage_single.sbatch`

Transcribed verbatim. Not exercised directly in this task (the brief's own
verification step invokes the census over `ssh` rather than through
`sbatch`), so it is untested here beyond being syntactically a valid sbatch
script; it exists as the shared template for later stages (2, 4, 5, 6, 6b, 7).

### Step 2 — `scripts/stage0b_rna_census.py`

Transcribed verbatim, per the review above.

### Step 3 — sync to Euler

```
$ ssh euler 'mkdir -p /cluster/work/wicky/lvalsecchi/rna_binder/analyses/foldseek_map/logs'
$ rsync -av --exclude '.git' --exclude 'tests/fixtures' --exclude 'build' \
    ./ euler:/cluster/work/wicky/lvalsecchi/rna_binder/analyses/foldseek_map/
$ rsync -av build/ euler:/cluster/work/wicky/lvalsecchi/rna_binder/analyses/foldseek_map/build/
sent 61184782 bytes  received 148 bytes  3706358 bytes/sec
```

Verified `src/` and `scripts/` (including the two new files) actually landed
on Euler by listing them remotely before running anything — both present.
Remote target directory did not exist beforehand (fresh checkout); created it
first.

### Step 4 — run the census

```
$ ssh euler 'cd /cluster/work/wicky/lvalsecchi/rna_binder/analyses/foldseek_map && \
    PYTHONPATH=src /cluster/work/wicky/lvalsecchi/binder_map_venv/bin/python \
    scripts/stage0b_rna_census.py --build build --per-group 3'
```

Output (warnings elided — all are the expected biotite
`auth_comp_id`/`auth_atom_id` fallback `UserWarning`s, ~100 of them, one pair
per structure read):

```
=== RNA variant census -> build/meta/rna_variants.tsv ===
rna_sequence        length  status
ACUGCCGUAUAGGCAG    16      ok        91
CACUGCCGUAUAGGCAGC  18      ok        14

all RNA variants align to the canonical frame
```

`EXIT_CODE=0` (captured explicitly after the ssh command). No `FRAME_MISMATCH`,
no `unreadable:` rows, no `FATAL`/`Traceback` anywhere in the log — checked
with an explicit grep, not just eyeballing the tail.

**Full `build/meta/rna_variants.tsv`** (105 rows = 3 layers x 5 stacks x 7
campaigns, `--per-group 3`), pulled back locally for the record:

| layer | stack | campaigns present | sequence | length | offset | status |
|---|---|---|---|---|---|---|
| design | BD1+LiMPNN | v2,v4,v5,v6,v7,v8,v9 | CACUGCCGUAUAGGCAGC | 18 | 1 | ok |
| design | BG+BIF | v2,v4,v5,v6,v7,v8,v9 | ACUGCCGUAUAGGCAG | 16 | 0 | ok |
| design | BG+LiMPNN | v2,v4,v5,v6,v7,v8,v9 | ACUGCCGUAUAGGCAG | 16 | 0 | ok |
| design | HD+LiMPNN | v3,v4,v5,v6,v7,v8,v9 | CACUGCCGUAUAGGCAGC | 18 | 1 | ok |
| design | RFD3+LiMPNN | v3,v4,v5,v6,v7,v8,v9 | ACUGCCGUAUAGGCAG | 16 | 0 | ok |
| af3 | all 5 stacks | as above | ACUGCCGUAUAGGCAG | 16 | 0 | ok |
| esm | all 5 stacks | as above | ACUGCCGUAUAGGCAG | 16 | 0 | ok |

Every row is `status=ok`. Two distinct sequences total, both align: the
canonical 16-mer at offset 0 (91 of 105 group/sequence rows — every af3 row,
every esm row, and the design-layer BG+BIF/BG+LiMPNN/RFD3+LiMPNN rows), and
the 18-nt `CACUGCCGUAUAGGCAGC` variant at offset 1 (14 rows — only
design-layer BD1+LiMPNN and HD+LiMPNN). This matches the brief's stated
expectation exactly, including which variant appears and at which offset.

Incidental finding, not a bug: the 18-nt closing-pair variant only shows up
in the **design** layer for BD1+LiMPNN/HD+LiMPNN; both refold layers (af3,
esm) show the canonical 16-mer for every stack including those two. Consistent
with the refold target being re-specified to the canonical 16-mer regardless
of what the design-layer cif carried. Not something this task needed to
change — the frame code already treats both as in-frame via `frame_offset`.

### Step 5 — commit

```
$ git add scripts/stage0b_rna_census.py slurm/stage_single.sbatch
$ git commit -m "Add RNA-variant census as a hard gate before compute"
[foldseek-binder-map 8781a76] Add RNA-variant census as a hard gate before compute
 2 files changed, 131 insertions(+)
```

## Self-review

- Both files' code bodies were checked character-for-character against the
  briefs before writing (diffed mentally line by line while transcribing);
  no discrepancies introduced.
- Task 11: only `scripts/stage0_select.py` and `tests/test_stage0.py` staged
  and committed — confirmed via `git status` immediately before `git add` and
  again after, in both commits. `build/` was never staged (it's pipeline
  output, not source).
- Task 12: only `scripts/stage0b_rna_census.py` and `slurm/stage_single.sbatch`
  staged and committed, same verification.
- Did not touch `RNA_CANONICAL`, `frame_offset`, or any matching logic in
  `src/binder_map/rnaframe.py` — the gate passed on real data with no need to
  loosen anything, so per the task instructions nothing there was touched.
- No unit tests exist for Task 12 by design (per the brief); verification is
  the actual Euler run above, including an explicit grep for
  `FRAME_MISMATCH|FATAL|Traceback` across the full census log rather than
  trusting the visible summary alone.
- One dead-code nit left as-is per "transcribe verbatim, fix only real bugs":
  `io` and `tarfile` imports in `stage0b_rna_census.py` are unused now that
  extraction goes through the `tar` CLI. Not a functional bug.
- The only deviation from either brief's stated expectations was the design
  numbers (~8,931 archives vs. the brief's "roughly 6,780" estimate), which
  is a fact about the real data, not a code issue — flagged above rather than
  silently accepted or "corrected" toward the brief's guess.
