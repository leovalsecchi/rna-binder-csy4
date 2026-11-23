# gpubench — per-stack GPU efficiency bench

`rna_binder/bench_gpueff/`. Built 2026-11-11. Label `gpubench_v1`.

Measures **GPU-h per design, per stack × role**, cleanly enough to rank the five
generation stacks against each other. Feed the result your campaign pass rates
to get GPU-h per passing design.

## Why not just use the campaign ledger

The 197,797 GPU-h `sacct` ledger (see `gpu_hours/`) cannot rank the stacks:

1. **32% of it is waste** — only 68% of campaign GPU-h ended `COMPLETED`; the
   rest was preempt, fail or timeout, unevenly across stacks.
2. **Different GPUs** — RFD3/BG/BD1 ran on Turing `eu-g3` (quadro_rtx_6000 /
   titan_rtx), HD on Ampere/Ada. That confounds exactly the comparison wanted.
3. **BoltzGen gen is unsplittable there** — one `boltzgen_multif.sh` job feeds
   both BIF and LiMPNN, so the ledger records it as `shared`.
4. **No clean denominator** — the ledger has no per-stack design counts.

The bench fixes all four by re-running campaign-identical jobs on one card model
and counting the designs each job actually produced.

## Design

**Generation** — whole campaign-identical jobs, because per-design gen cost
includes job-level amortization. 2 × RFD3 (640 designs each), 2 × BoltzGen
(600 BIF + 600 LiMPNN each, shared backbone), 2 × BD1 (~600), 4 × HD (150 each).

**Refolding** — a fixed **600 designs per stack**, one AF3 + one B2 + one ESM job
each. This is the part that has to be equal: AF3 and ESM load the model once and
run over `--input_dir` (`refold_af3_scratch.sh:262-290`), which is where their
10-13× speedup lives, so per-design cost depends on batch size. HD makes only 150
designs per uuid, so refolding per-uuid would give HD a 150-design batch and
everyone else 600 — inflating HD's refold cost for a reason unrelated to HD.
`bench_refold.py` therefore pools each stack's designs and uses the pool as the
refold `source_run`.

The pool is a **packed** source, not a symlink dir (corrected 2026-11-11 — see
*What the first run got wrong* below): a finished gen run has no `designs/` at
all, so the pool is built by extracting the chosen members out of each run's
`designs.tar.zst` and repacking them as one `pools/<stack>/designs.tar.zst`,
which the refold worker already knows how to read
(`refold_{af3,b2,esm}_scratch.sh` all source `archive_pipeline.sh` and unpack a
packed `source_run` to node-local scratch). Inode cost: 2 files per stack
instead of 600 symlinks, and nothing is written to Lustre per design.

**Why 600 and not a single common N:** two stacks are quantized. HD comes in
multiples of 150 (`num_recycles` 15 × `num_seqs` 10, and `num_designs` must stay
1 or de-novo trajectories collide on `random_init_NNN` tags). RFD3 comes in
multiples of 80 (`n_batches` × `diffusion_batch_size` 8 backbones × 10 seqs, and
the diffusion batch size is the GPU batch — changing it changes what's measured).
The smallest N both hit exactly is LCM(150,80) = 1200. Decoupling gen from refold
gets the fairness where it matters at half the cost.

## What the first run got wrong (2026-11-11)

Both of these were found *after* the gen jobs had already run, and both are now
fixed in the scripts. They are recorded because each one silently produces a
plausible-looking wrong number.

**1. A finished gen run has no `designs/` directory.** The finalize tail packs it
into `designs.tar.zst` + `.tar.zst.list` and deletes the tree — the inode
discipline the whole campaign runs under. Every bench script was written against
an unpacked `designs/`, so all three read zero designs: the pool builder would
have refused to submit, and the collector would have rejected every job for
"0 designs". All reads now go through `bench_lib.py`. The tarball doubles as the
completeness signal, since only the finalize tail publishes it — so a cancelled
or still-running job is excluded automatically rather than contributing a
partial denominator.

A trap inside that fix: the member filter must key on **extension alone**. Each
stack names its designs differently — RFD3/BG/BD1 use `_design<D>_seq<S>`, HD
uses `_recycle<N>_sample<M>_pack1` — so a name-shaped filter drops a whole stack.
An early `"_design" in name` test zeroed out HD completely while looking like it
worked for everyone else.

**2. HD gen was submitted without the DeepSpeed JIT fix, and one job hung for
good.** `launch.py` does not emit the per-task `TORCH_EXTENSIONS_DIR` block; it
is a post-generation patch that every new HD batch has to re-apply
(`ops/submit/patch_topup_resume_jit.py`). Without it, all four HD jobs raced for
the single lock in `containers/runtime_cache/torch_extensions/evoformer_attn`:

| uuid | start → DeepSpeed init | outcome |
|---|--:|---|
| e8724b45 | 6 min | COMPLETED 2.66 h |
| 9412f761 | 9 min | COMPLETED 3.09 h |
| d6f3d66d | 10 min | COMPLETED |
| 828c688a | **119 min** | stalled at recycle 1, killed after 4.8 h |

`bench_gen.sh` now applies the fix itself for any sbatch containing
`HalluDesign_run.py`, and aborts rather than submitting unpatched.

The clock corruption here is the dangerous part: `sacct` cannot see it. The hung
job would have reported `COMPLETED` on 1 GPU with a full design count, and its
4.8 h would have been averaged straight into HD's number. So `collect_gpuh.py`
now measures start → DeepSpeed init from the HD log and **rejects** any gen job
over 20 minutes, alongside the existing requeue/preempt/GPU-count rejections.
Clean jobs sit at 6–10 min; a lock-loser sits at two hours. The 6–10 min that
clean jobs *do* spend is a kernel rebuild the campaign's pre-seeded jobs never
pay, so bench HD gen carries ~4% of overhead the campaign did not — far too
small to touch a 50× gap, but it is there.

**3. The rewritten BG+BIF CIFs were unreadable by the refolder** (found 2026-11-12,
after the first refold round). All three BG+BIF refolds died in stage 1 on
`KeyError: '_atom_site.B_iso_or_equiv'`; the four stacks whose design files were
untouched refolded normally, which localised it immediately.

Two traps stacked. First, `st.b_factor = arr` is a **silent no-op in biotite** —
`__setattr__` only routes to an annotation that already exists, so a new name
becomes a plain Python attribute and the writer never emits the column.
`set_annotation("b_factor", ...)` is required. Second, the worker parses designs
with **BioPython**, not biotite (`refold_af3_scratch.sh:118-120`, inside
`boltzgen.sif` — not the host python, and not `alphafold3.sif`, none of whose
three interpreters has `Bio`). `MMCIFParser` hard-requires `B_iso_or_equiv` and
`occupancy`, which biotite is happy to read back without. So a biotite round-trip
looks clean and still produces a file the pipeline rejects.

`bench_refold.py` now runs `probe_biopython.py` over 10 files of **every** pool,
in `boltzgen.sif`, before packing, and aborts rather than submitting. The fixed
rewrite returns the identical 80-aa sequence as the original re-prediction —
which is the proof that the rewrite moves coordinates only, leaving refold cost
untouched.

## Deliberate deviations from the campaign

Model params, seeds policy, containers, cpus and mem are the **production configs
untouched** — only the `slurm:` block is overridden, exactly as the campaign
drivers did. Five changes, each required:

| change | why |
|---|---|
| `--gpus=rtx_4090:1`, no `--nodelist` | All Euler 4090s are `eu-g6` (80 nodes, all in `gpuhe.24h`), so the card model *is* the node pin. One architecture removes GPU generation as a confounder. Also dodges `eu-g4-032`, which is decommissioned but still baked into the refold configs — submitting with it fails `Invalid node name`. |
| non-preemptible partitions only, never `.bulk` | A preempted job's `ElapsedRaw` is a *wrong* clock, not a slow one. |
| `--requeue` / `--open-mode` stripped everywhere | All four gen templates and all three refold builders emit them. The AF3 worker skips already-done designs, so a requeued attempt measures only the remainder. |
| HD: `--ntasks=4` → `--cpus-per-task=4` | `launch.py:209` emits `--ntasks` with no `--cpus-per-task` for HD only. The campaign fixed this in the array bucketer, so this **reproduces** campaign conditions. Left alone, HD's single process gets 1 core instead of 4. |
| `refold_*.sh` → `refold_*_scratch.sh` | Every campaign driver does this rewrite before submit (`mv8_refold_bulk.py:227`), so the `_scratch` variants are what actually ran. Node-local scratch also keeps the bench off the Lustre inode quota. |

**Not** changed, deliberately: B2's `--ntasks-per-node=4`. The campaign never
fixed that one, so leaving it is fidelity. It is identical across stacks and
cannot affect the ranking.

Consequence to expect: **bench numbers will be faster than the ledger** for
RFD3/BG/BD1, because they ran on Turing. That is the point.

## The BG+BIF true-backbone fix

`boltzgen_multif.sh:209-246` harvests BG+BIF's "design" from BoltzGen's
`refold_cif/` — a Boltz-2 **re-prediction** of the designed sequence, not the
backbone the design step produced. The fold is right (0.82 Å median) but the pose
is not: median 4.4 Å RNA-aligned, 25% of designs >10 Å from where BoltzGen put the
binder. So every campaign BG+BIF `ligand_RMSD` / `sc_DockQ` was measured against a
different docking than the design step produced.

`bif_true_designs.py` rewrites each bench BG+BIF design as the real `spec_<D>.cif`
protein backbone (N,CA,C,O) + the RNA, with residue names transplanted
positionally from the BIF sequence. The spec cannot be used directly — it carries
BoltzGen's own co-designed sequence, ~50% identical to the BIF sequence, which
trips `lrmsd_rna_aligned.py`'s `MIN_REGISTER=0.9` gate and makes DockQ's internal
alignment a coin flip. Hence the transplant, reusing the already-validated
`analyses/bif_true_bb_rescore/rescore_lib.py` rather than a second implementation.

Backbone-only is correct, not a shortcut: BG+LiMPNN design files are themselves
backbone-only (LigandMPNN writes no sidechains without `--pack_side_chains`), so
this makes the two inverse folders **symmetric**, which is the point of a stack
comparison. `ligand_RMSD` is backbone-only anyway; only `fnat` moves on stripping.

Generation is untouched: `--steps inverse_folding folding` stays exactly as
production, so the BoltzGen job's GPU cost is unchanged and still comparable to
the campaign. The refolders read only the **sequence**, so their cost is unchanged
too. Only the stored coordinates — the scoring reference — change. The
re-prediction is preserved in `designs_refold_cif/` so both references stay
computable and the delta is quantifiable.

The rewrite is applied to the **staged** copy during the pool build, not to the
gen run on Lustre: a finished run is packed, so rewriting in place would mean
unpacking and repacking a production artefact to no purpose. Production output is
never mutated, and the re-prediction stays exactly where the gen job left it.

### The self-check
`test_transplant.py` rebuilds all 20 BG+LiMPNN designs in
`analyses/bif_backbone_probe/bundle/` from their specs and requires exact
agreement — a BG+LiMPNN design file *is* the spec backbone, so any drift means
the transplant, the spec→design index mapping or the chain handling is wrong.

Result 2026-11-11: **20/20 reproduced, worst coord deviation 0.00050 Å**, matching
the 2026-09-30 probe — and re-confirmed at the same 0.00050 Å against this
bench's own packed BG output, which also proves the tar reading and the
spec→design index mapping. The same run measured the stored BIF designs at a raw CA
RMSD of 7.54–25.57 Å (median 18.87) from their true backbones — the bug, in one
number. `bif_true_designs.py --validate` gates `--apply`, so a broken transplant
can never reach the design files.

## Running it

```bash
cd /cluster/work/wicky/lvalsecchi/rna_binder/bench_gpueff
python3 test_transplant.py           # must print PASS before anything else
bash bench_gen.sh                    # dry run; prints every sbatch header
bash bench_gen.sh --submit           # 10 gen jobs -> gen_manifest.tsv
```

Phase 2 is chained with `--dependency=afterany` on the gen jobs and fires by
itself (`afterany`, not `afterok`, so a failed gen job gets *reported* rather than
hanging the chain). It runs the true-backbone rewrite, then the equal-N pools and
refold submit. To do it by hand instead:

```bash
python3 bif_true_designs.py --label gpubench_v1 --validate   # gate
python3 bench_refold.py --label gpubench_v1 --submit         # pools + 15 jobs
python3 collect_gpuh.py
```

Re-running a single stack's gen (e.g. to replace a lost job) uses the same
submitter — `REPS`, `HD_UUIDS` and `HD_SEED_BASE` are env-overridable, and the
seed must be fresh:

```bash
REPS=0 HD_UUIDS=1 HD_SEED_BASE=900005 bash bench_gen.sh --submit
```

## Reading the output

`collect_gpuh.py` → `gpuh_per_design.tsv` plus a printed table.

It **rejects** rather than averages: any job not `COMPLETED`, any job with more
than one `sacct -D` attempt (requeued), any job not on exactly 1 GPU, any job with
0 designs on disk. Rejections are listed, never silently dropped. `-D` matters —
without it requeued attempts vanish and the surviving row under-reports.

BoltzGen gen is reported as its own `bg_shared` row, with the measured stage
split plus 50/50 and full-charge-each. The split is **two-way, not three-way**:
`boltzgen_multif.sh` prints its stage banners with no timestamps and production
is left untouched, so the only datable boundary is the end of backbone design —
the `_backbones/BG_len80_<uuid8>/spec_*` files, which survive unpacked on Lustre.
The finer boltzif-vs-ligandmpnn split is not recoverable: `bif_native/` is never
published, and both inverse folders' trees are packed into tarballs written at
the very end of the job, so their mtimes are identical. Measured: **17% of a BG
gen job is genuinely shared** backbone design, 83% is inverse folding for the two
stacks together.

## Measured gen results (2026-11-11, n=2 clean jobs per stack)

| stack | GPU-h/design (gen) | ×best |
|---|--:|--:|
| RFD3+LiMPNN | 0.000353 | 1.0 |
| BG+BIF / BG+LiMPNN | 0.001081 each (50/50) · 0.002162 (full charge each) | 3.1 / 6.1 |
| BD1+LiMPNN | 0.006217 | 17.6 |
| HD+LiMPNN | 0.019176 | **54.3** |

Faster than the campaign ledger across the board, as predicted (4090 vs Turing),
but the *spread* is wider than the ledger suggested: HD/RFD3 is 54× here against
41× from campaign medians. Refold numbers follow once phase 3 runs.

## Expected, from campaign medians

The original sanity targets, kept for comparison. Gen medians are Turing for
RFD3/BG/BD1, so the bench should come in at or below these; refold medians are
already Ampere/Ada.

| stack | gen h/job | designs/job | gen h/design | refold h/design | total | ×best |
|---|--:|--:|--:|--:|--:|--:|
| RFD3+LiMPNN | 0.55 | 640 | 0.00086 | 0.0058 | 0.0067 | 1.0 |
| BG+BIF | 2.40 shared | 600 | 0.0020 | 0.0058 | 0.0078 | 1.2 |
| BG+LiMPNN | shared | 600 | 0.0020 | 0.0058 | 0.0078 | 1.2 |
| BD1+LiMPNN | 6.40 | ~600 | 0.0107 | 0.0058 | 0.0165 | 2.5 |
| HD+LiMPNN | 5.30 | 150 | 0.0353 | 0.0058 | 0.0411 | 6.1 |

Refolding is stack-independent (same backends, same 80-mer, same 18-nt RNA), so
the ranking lives almost entirely in generation — HD ≈ 41× RFD3 there.

Bench cost ≈ 58 GPU-h (gen ~40, refold ~18) = 0.03% of the campaign.

## Gotchas found while building this

- **`launch.py` has no standalone `refold_esm` branch.** `config_kind: refold_af3`
  (line 1356) and `refold_b2` (1426) exist; an ESM config falls through to the
  generator path and dies on `KeyError: 'generator'`. All three
  `build_refold_*_sbatch` functions do exist, so `bench_refold.py` imports them
  and drives all three uniformly.
- **Euler routes by walltime, not by the partition you ask for.** RFD3 at
  `--time=04:00:00` lands in `gpuhe.4h` even with `--partition=gpuhe.24h`. Both
  tiers are the same eu-g6 4090s, so the measurement is unaffected.
- **`eu-g6-064`** is on the campaign refold driver's permanent exclude list and is
  inside the pinned pool — excluded here too.
- **BoltzGen multif creates two dirs sharing one uuid** (boltzif + ligandmpnn) and
  only one carries `run.sbatch`; resolve by the file, not the first match.
- **spec padding is inconsistent** (`spec_7` vs `spec_07`) — always match on the
  parsed int.

## Files

| file | role |
|---|---|
| `bench_lib.py` | reads/writes the packed `designs.tar.zst`; the finished-run gate |
| `probe_biopython.py` | parses designs with the worker's own parser; the pool gate |
| `bench_gen.sh` | gen submits, header edits, HD JIT fix, submit gate |
| `bench_refold.py` | equal-N pools + all 3 refold backends |
| `bif_true_designs.py` | BG+BIF true-backbone rewrite (`--validate` self-check; the rewrite itself runs inside the pool build) |
| `test_transplant.py` | the self-check, against 20 real campaign pairs |
| `collect_gpuh.py` | GPU-h/design table, with rejections |
| `bench_phase2.sbatch` | the dependency-chained phase 2 |
| `gen_manifest.tsv`, `refold_manifest.tsv` | job provenance the collector reads |
