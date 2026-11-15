# csy4 6kseq_4h — concurrency bug fixes, BIF Route B, and recovery (2026-06-11)

Single-session as-built. Three independent problems found + fixed, plus a full
recovery of the `csy4_crrna_free/6kseq_4h` campaign. All code edits keep `.bak_*`
backups next to the originals.

---

## 1. Finalize-tail canonicalize concurrency abort (gen jobs)  — FIX DEPLOYED

**Symptom:** 20/60 `rfd3__ligandmpnn/6kseq_4h` gen jobs `FAILED 2:0` at ~12–17 min
(well under the 4 h cap) after writing all 64×10 designs, leaving 704 loose files
(`640 .pdb + 64 .fa`) in the shared `designs/` and no `<uuid>.tar.zst` slice. All
designs intact; only the finalize tail died.

**Cause:** the post-generation tail (`canonicalize_designs.py` → purge → `pack_slice`)
runs under `set -euo pipefail` with the `canonicalize` `apptainer exec` **unguarded**.
Under ~60-way concurrent launch the container/FS transiently fails; the non-zero exit
aborts the job before the `|| true`-guarded pack. Re-running canonicalize succeeds
(transient). `canonicalize_designs.py` also did `sorted(designs.iterdir())` over the
*shared* dir while 59 siblings churned it — the race surface.

**Fix** (`*.bak_pre_canonretry`):
- All 4 `*_multif.sh` + `postprocess_halludesign.sh`: the canonicalize `apptainer exec`
  is wrapped in an `until`-retry (5 tries, jittered backoff `sleep (RANDOM%15)+5*try`;
  `exit 1` only after all fail → designs stay loose+recoverable). set-e-safe (validated).
- `canonicalize_designs.py`: scoped mode now `sorted(designs.glob(f"*{match}*"))`
  instead of `iterdir()`+filter → a scoped job never lists siblings' churn.

---

## 2. boltzgen mpnn_native shared-path collision (launch.py arg order)  — FIX DEPLOYED

**Symptom:** 14/30 `boltzgen__ligandmpnn/6kseq_4h` UUIDs crashed *mid*-MPNN
(`exit 1`, partial 0–480 designs), 11 with
`FileNotFoundError … designs/mpnn_native_design31//seqs/spec_31.fa`.

**Cause:** `boltzgen_multif.sh` puts each backbone's LiMPNN workdir at
`$MPNN_RUN_DIR/mpnn_native_design<idx>` and runs `purge "$MPNN_RUN_DIR"/mpnn_native_design*`.
That is only safe if `MPNN_RUN_DIR` is per-UUID. `launch.py`'s
`build_boltzgen_multif_sbatch` passed the **shared** designs dir as positional `$2/$3`
(`BIF_RUN_DIR`/`MPNN_RUN_DIR`) instead of the per-UUID `bif_run_dir`/`mpnn_run_dir`
→ all 30 jobs shared `designs/mpnn_native_design<idx>/`, and one job's purge deleted a
sibling's in-flight `seqs/`. (rfd3 was already correct.)

**Fix** (`launch.py.bak_pre_bgargfix`, lines ~614–615): pass `{bif_run_dir}` and
`{mpnn_run_dir}` (per-UUID) as `$2/$3`; `$6/$7` stay the shared designs dirs. The
script's own cleanup comments already expected per-UUID `BIF_RUN_DIR/MPNN_RUN_DIR`.

> NB: these 14 designs are genuinely incomplete (mid-MPNN), NOT recoverable from disk —
> they get re-run. The completeness gate (below) refused to pack them.

---

## 3. BIF Route B — 10 designs/backbone parity with LiMPNN  — LIVE (validating)

**Problem:** BoltzIF produced **1 design/backbone** (64 for 6kseq) vs LiMPNN's
**10/backbone** (640) — systematic 10× imbalance in every production run (3kseq 32:320,
6kseq 64:640). The two inverse-folding stacks were never comparable.

**Cause:** BoltzGen *already* runs inverse-folding 10×/backbone
(`bif_cli.inverse_fold_num_sequences=10` → `multiplicity:10`; the IF DataLoader runs
20 items for 2 backbones). The compute is spent — but `DesignWriter`
(`/app/src/boltzgen/task/predict/writer.py`) collapses them at write time:
`total_files = multiplicity * n_samples` reads `multiplicity=1` (bad
`getattr(trainer.datamodule.cfg, "multiplicity", 1)` lookup), so it uses the
non-unique `{stem}.cif` name → 10 replicas overwrite → 1 survives. The other 9 are
computed then discarded.

**Fix (Route B, deployed via bind-mount — NO container rebuild):**
1. `…/rna_binder/boltzgen_patch/writer.py` — 2-line patch: take the unique
   `{stem}_{global_idx}` branch when `"data_sample_idx" in batch` (the datamodule sets
   it per replica), widen `num_digits`. → emits `spec_<bb>_<seq>` (10/bb).
2. `boltzgen_multif.sh` (`*.bak_pre_routeb`): BIF `apptainer exec` adds
   `--bind …/boltzgen_patch/writer.py:/app/src/boltzgen/task/predict/writer.py`; the
   harvest copy-loop maps `spec_<bb>_<seq>` → `_design<bb>_seq<seq>` (back-compat with
   `spec_<bb>`). Every boltzgen multif run now auto-produces 10/bb.

**Cost (measured, rtx_4090):**
- IF compute already paid → +0.
- Folding: ~5.6–8.3 s/structure × 740 ≈ **~100 min** → **~+1.7 h GPU/UUID**.
  Whole job ~40 min → ~2.5–3 h (fits gpuhe.4h thinly; 24 h safer).
- **VRAM unchanged ~4.5 GB peak** (folding is sequential; the 20 GB request is ~4×
  over-provisioned). Route B is pure wall-clock, zero extra VRAM.
- Refold load also ×10 on the BIF side (740/UUID) — decide refold scope deliberately.

**Validation:** standalone IF → 10/bb confirmed (`spec_00_00..spec_00_09`). End-to-end
pipeline = job **3039490** (`_routeb_val`, num_designs=4, expect BIF≈40/LiMPNN≈40) —
queued at session end; confirm before the campaign rerun.

---

## Recovery + campaign state

- Recovered `6kseq_4h` to: **rfd3 60/60, boltzif 30/30, bg_limpnn 16/16** packed
  (`canonicalize + pack_slice`, completeness-gated so no short slice). Log:
  `…/rna_binder/_recovery_2026-06-11.log`.
- **Recovery must run as an sbatch job** (`normal.4h`, `es_wicky`) — the login-node
  process reaper kills the CPU-heavy canonicalize+pack loop after ~15 min. `pgrep` is
  node-local, so a "dead?" check on another login node lies (two runs raced once;
  harmless — `pack_slice` is flock'd + gated).
- **Old boltzgen `6kseq_4h` archived** (reversible) to
  `…/csy4_crrna_free/_archive/pre_routeb_2026-06-11/{boltzgen__boltzif,boltzgen__ligandmpnn}__6kseq_4h`,
  and the live dirs cleared → user reruns Route B via the GUI into a clean `6kseq_4h`
  once 3039490 is green. rfd3 slices untouched.

## Files
- Patch: `…/rna_binder/boltzgen_patch/writer.py` (+ `.orig`)
- Backups: `…/code/scripts/stages/{*_multif.sh.bak_pre_canonretry, boltzgen_multif.sh.bak_pre_routeb}`,
  `canonicalize_designs.py.bak_pre_canonretry`, `…/scripts/launch.py.bak_pre_bgargfix`
- Recovery log: `…/rna_binder/_recovery_2026-06-11.log`

---

## Route B — UPDATE (validated, with a 2nd fix the writer patch alone missed)

The writer patch + bind alone did NOT work in the pipeline: `--only_inverse_fold`
collapses to 1/backbone (its IF DataLoader runs N=multiplicity *total* items, not
per-backbone — per the boltzgen README it expects "a fully specified structure").
The standalone test that "worked" had used `--steps inverse_folding folding`, which
fans IF out per design-backbone × multiplicity = 10/bb.

**2nd fix** (`boltzgen_multif.sh`, backup `.bak_pre_stepsfix`): switched the BIF call
`--only_inverse_fold` -> `--steps inverse_folding folding`, and harvest from
`intermediate_designs_inverse_folded/refold_cif` (where `--steps` routes outputs).

**VALIDATED 2026-06-11** (v4 job 3053806, num_designs=2 -> 2 backbones): BIF=20,
LiMPNN=20, both 10/bb, named `_design<bb>_seq<seq>`, BIF==LiMPNN. Green light for rerun.

**Runtime (measured):** BIF `folding` AND `design_folding` both process the 10x set
(~6-8 s/struct, high node variance) -> num_designs=64 ~= 3.4 h. USE gpuhe.24h, NOT 4h.
VRAM unchanged ~4.5 GB. (gpu_capacity_calc.py boltzgen sec_per_unit 63 -> 190;
boltzdesign 370 -> 108 + "use 24h"; rfd3 dbs=104 flagged extrapolated/untested.)
