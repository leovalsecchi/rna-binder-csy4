# Task 13 report: Stage 1 — extraction and fingerprinting

Commit: `5e93acc` on branch `foldseek-binder-map`.

## What was implemented

- `scripts/stage1_extract.py` — one pass per shard over its assigned whole
  archives: extracts backbone PDBs into a single `shard_%02d.tar.gz`,
  writes `coord_hash`/`protein_chain`/`n_res` to `hashes/shard_%02d.tsv`,
  the 48-dim fingerprint + `per_nt_1..16` + `contacts_offframe/total` +
  `cx/cy/cz/ax/ay/az` to `fp/shard_%02d.parquet`, and unreadable/rejected
  design_ids to `missing/shard_%02d.tsv`. Everything streams through
  `$TMPDIR`; no loose files land under `results/`.
- `slurm/stage1_extract.sbatch` — array wrapper, one task per shard,
  `normal.4h`, 4 cpus, 4G/cpu, used as
  `sbatch --job-name=bm_s1_esm --array=0-63%16 slurm/stage1_extract.sbatch esm`.

Followed the brief's code as the base, with two deviations (both bugs, not
style choices — details below).

## Deviations from the brief (two genuine bugs found and fixed)

### 1. `read_archive_members`: member misassignment on multi-member archives

The brief's version extracts all requested members with one `tar --zstd
-xOf archive -O member1 member2 ...` call, splits the stdout stream on
`\ndata_` (the mmCIF block boundary), and does `zip(members, blobs)` to
pair each member name back to its content. That assumes `tar -O` emits
matched members **in the order they were requested**. It does not — `tar`
streams the archive once and emits members in the archive's own on-disk
order, regardless of the order they were listed on the command line.
Verified directly:

```
$ tar -cf t.tar b.cif a.cif c.cif   # archive order: b, a, c
$ tar --zstd -xOf t.tar.zst a.cif c.cif b.cif   # requested: a, c, b
hello b
hello a
hello c
```

`build/manifests/design.tsv` gets its `member_path` order from
`passing.csv`'s row order (task 11's `build_manifest`, stable-sorted only
by `archive_path`), which has no relationship to the archive's internal
member order. 7,541 of the 8,931 archives (84%) hold more than one member,
so this was not a corner case — it would have quietly cross-assigned CIF
content between design_ids across most of the design layer. Reproduced
with realistic CIF-style content (`data_design<X>` headers) requested in
non-archive order:

```
designs/a.cif requested -> got: data_designB
designs/c.cif requested -> got: data_designA
designs/b.cif requested -> got: data_designC
```

— a full circular mislabeling with zero errors raised. The `if
proc.returncode != 0: return {}` guard was a second problem riding the
same function: a single genuinely-missing member in a requested batch (tar
exits non-zero, but still emits the members it did find) would discard the
*entire* archive's structures as `member_missing`, not just the missing
one.

**Fix**: extract requested members by name into real files under
`$TMPDIR/.../members/` (`tar --zstd -x -C workdir -f archive member...`),
then read each file by its exact `member_path`. This sidesteps the
ordering assumption entirely and makes partial misses precise (only the
files that didn't materialize are `member_missing`). Verified against the
same reordered synthetic archive — all three members came back correctly
matched by name, and a genuinely-missing member no longer nukes its
siblings. `shutil.rmtree` replaces the final `tmp.glob("*").unlink(...)`
cleanup (which would have raised `IsADirectoryError` on the new
`members/` subdirectory) with a single recursive removal of `$TMPDIR`'s
scratch dir.

### 2. Brief's own step-5 verification command: sort/cut ordering

The brief's step 5 pipes `sort -k1,1 -t"\t" build/manifests/design.tsv |
cut -f2,4` into `join`. `design.tsv`'s field 1 is `row_id`, not
`design_id` (which is field 2) — so this sorts by `row_id` and only
extracts `design_id`/`stack` *afterward*, handing `join` a stream that is
NOT sorted on the key it's told to join on (`-2 1`, now `design_id` post-
`cut`). `join` requires both streams to be pre-sorted on the actual join
key for its merge algorithm to work; violating that silently drops
non-adjacent matches instead of erroring. Running it as written produced
one spurious match out of 154 rows in the shard. Fixed by reordering to
`cut -f2,4 build/manifests/design.tsv | sort -k1,1 -t"\t"` (cut first,
then sort on the extracted `design_id`). This is a bug in the brief's
shell snippet, not in `stage1_extract.py`; the Python implementation was
unaffected and the real tally (below) came from the corrected pipeline.

Both bugs are documented in the commit message.

## Step 3: single-shard smoke test (esm layer)

```
$ ssh euler 'cd .../foldseek_map && \
  TMPDIR=/tmp PYTHONPATH=src .../binder_map_venv/bin/python \
  scripts/stage1_extract.py --layer esm --shard 0 --n-shards 512'
layer=esm shard=0/512: 166 structures from 18 archives
done in 22s: 166 parsed, 166 unique backbones, 0 missing
```

Matches expectation: a few hundred structures, 0 missing.

## Step 4: shard contents + inode footprint (esm)

```
=== shard members ===
csy4free_BD1_LiMPNN_len80_00bdadf5_design13_seq1.pdb
csy4free_BD1_LiMPNN_len80_6abda374_design54_seq6.pdb
csy4free_BD1_LiMPNN_len80_6abda374_design54_seq8.pdb
=== count === 166
=== inodes used by results/ === 9
=== protein chain tally (must be A for esm) ===
    166 A
      1 protein_chain      (header line self-match, not a data row)
=== fingerprint sanity ===
rows 166 cols 73
contacts_total: min 144 median 306.0
zero-contact rows: 0
offframe nonzero: 0
```

166 members named `<design_id>.pdb`, no dedup collisions in the esm layer
(refold coordinates rarely round-hash identical). 73 columns = `design_id`
+ 48 fingerprint + 16 `per_nt` + `contacts_offframe` + `contacts_total` +
6 geometry = 73. Median 306 contacts, all-A chain — as expected for a
refold layer.

## Step 5: design-layer RFD3 chain-inversion tally (the critical check)

```
$ ssh euler '... stage1_extract.py --layer design --shard 0 --n-shards 512'
layer=design shard=0/512: 154 structures from 18 archives
done in 11s: 154 parsed, 78 unique backbones, 0 missing
```

Join (corrected pipeline, see deviation #2 above) of `hashes/shard_00.tsv`
against `build/manifests/design.tsv`'s `design_id, stack`:

```
      4 A BD1+LiMPNN
     20 A BG+BIF
     45 A BG+LiMPNN
      6 A HD+LiMPNN
     79 B RFD3+LiMPNN
      1 protein_chain stack   (header self-match)
```

**4 + 20 + 45 + 6 + 79 = 154, matching the 154 parsed rows exactly.**
RFD3+LiMPNN is unanimously chain `B` (inverted); every other stack is
unanimously chain `A`. This is exactly the expected guard result — chain
selection by composition (not chain ID) is working correctly for the
inverted 53%-of-layer stack.

## Combined inode footprint (both smoke-test layers)

```
$ find results | wc -l
19
```

19 inodes total under `results/` for both the esm and design smoke-test
shards combined (8 output files + directories) — in the tens as required,
nowhere near the inode cap. `$TMPDIR` scratch (`/tmp/bm_esm_0`,
`/tmp/bm_design_0`) was confirmed cleaned up after each run (`ls
/tmp/bm_*` → no such file).

## Fingerprint sanity, design layer (extra, not required by the brief but checked)

```
rows 154 cols 73
contacts_total: min 16 median 57.0
zero-contact rows: 0
offframe nonzero: 5
```

Lower contact counts than esm (median 57 vs 306) is expected: three of the
five design stacks ship backbone-only CIFs with no sidechains. The 5
off-frame rows are consistent with the known BG+BIF length variant
(canonical 16-mer plus one flanking base pair, per `rnaframe.py`'s
docstring).

## Local pre-flight checks (before spending Euler time)

- Full existing suite still green: `pytest -q` → 100 passed, no changes
  under `src/`.
- Ran the extraction pipeline (chain split → hash → fingerprint →
  approach → backbone PDB) directly against all 6 local fixtures:
  `rfd3_design.cif.gz` and `rfd3_design_seq7.cif.gz` both report
  `chain=B`, identical `coord_hash` (`801f8372...`) despite being
  different seq variants of the same backbone — confirms the dedup
  rationale in the brief's comment. `af3_refold`, `bd1_design`,
  `bgbif_design`, `esm_refold` all report `chain=A`. All 6 report
  `n_res=80`. `bd1_design.cif.gz` legitimately has `contacts_total=0`
  (that specific fixture design just doesn't contact the RNA within
  cutoff) — a property of the fixture, not a code issue.
- Verified via direct unit-level test that `n_res == 80` is a safe literal
  for this dataset: `grep -c len80` on all three manifests (`design`,
  `af3`, `esm`) matches their total row counts exactly (71,866/71,866 for
  each) — no other length variant present.

## Self-review / concerns

- `pd.DataFrame(hashes)` / `pd.DataFrame(fps)` would produce a columnless
  empty frame if a shard parsed zero structures (e.g. an all-archive-
  missing shard). Not fixed — no shard in this dataset is expected to hit
  that, and the brief doesn't ask for it, but it's a latent rough edge
  worth knowing about if a shard ever comes back truly empty (its parquet
  write could error or lose the schema).
- The `from binder_map import RNA_LEN` import in the brief's script is
  unused in the body. Left as-is (harmless, matches the brief, no lint
  gate in this repo enforces it) rather than treating an unused import as
  a "genuine bug."
- I did not delete the smoke-test shard_00 outputs from Euler's
  `results/{esm,design}/` — they're exactly what the real array run will
  overwrite when shard 0 executes for real, and at 19 inodes total there's
  no cleanup need.
- Both fixes were verified against synthetic/local repros before touching
  Euler, then the real behavior was re-confirmed on Euler with actual
  archives (steps 3-5 above) — the fixes weren't just theoretical.
