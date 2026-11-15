# midrun_v6 — AS-RUN conditions & lessons (the working recipe)

**Written 2026-07-08.** Authoritative record of what **midrun_v6** actually ran as. It supersedes
`docs/MIDRUN_V4_AS_RUN_CONDITIONS_2026-07-05.md` (v4/v5) — the array **mechanism**, node-local
scratch, and bucketer pattern are unchanged and still correct; the **only real deltas** are the
routing (v6 moved the fast generators to the high-cap `gpu.4h` tier), the measured runtimes that
justified it, and the v6 file names. Read the v4 doc for the full derivation of each fix; this doc
records the v6 specifics + what changed.

- Repo root `RB=/cluster/work/wicky/lvalsecchi/rna_binder`; `code/` is **not** git-tracked.
- Target: `csy4_crrna_free` (16-nt csy4 crRNA epitope).

---

## 0. TL;DR — the recipe (v6)

1. **Generate** per-uuid `run.sbatch` with `launch.py --num-runs N` (uses `--cpus-per-task`), then
   **cancel the broken native array** it submits, and **re-bucket** into self-contained Slurm
   arrays — one per routing signature — with **`ops/submit/bucket_midrun_v6.py`**. **Never `--ntasks`.**
2. **Route the fast Turing stacks to `gpu.4h`** — NEW in v6: RFD3 (~24 min) and BoltzGen (~2.7 h max)
   fit inside 4 h, and the `gpu.4h` QOS has ~2× the concurrency budget of the 24 h tier (≈96 GPU /
   1770 G). **BD1 (6.4 h median) can't fit 4 h → stays `gpu.24h` + `gpu.120h`.** **HD** (sm_80+ only)
   → `gpuhe`/`gpupr` `.24h` (unchanged from v5), spread onto `.120h` tiers for concurrency.
3. **Refold on the `.bulk` tiers with node-local scratch workers** (`refold_*_scratch.sh`,
   `--tmp=30G`): AF3/ESM → `gpuhe.bulk`, B2 → `gpu.bulk`, all at `--time=06:00:00`. Node-local
   scratch keeps project inodes flat (~6 inodes/leaf) — *the* inode fix.
4. **Two self-chaining loops** auto-fill refolds for every newly-packed uuid; **HD's per-uuid
   packs are consolidated** into the shared pool by `hd_combine.sh`.
5. **Exclude the black-hole node `eu-g6-064`** everywhere.

---

## 1. Design allocation (equal ~90k/stack, 5 stacks / 4 generators)

| Stack | Config (`code/configs/…`) | Designs/uuid | Target uuids | GPU floor | **Measured runtime/uuid (v6)** |
|---|---|---|---|---|---|
| RFD3 + LigandMPNN | `rfd3__multif/default.yaml` | 640 (64 bb × 10 seq) | **141** | Turing sm_75 OK | **~24 min** → `gpu.4h` |
| BoltzGen + BoltzIF | `boltzgen__multif/default.yaml` | ~593 | **152** (shared) | Turing sm_75 OK | **~2.3 h med / 2.7 h max** → `gpu.4h` |
| BoltzGen + LigandMPNN | *(same BG job, Route-B)* | ~595 | *(same 152)* | — | — |
| BoltzDesign (BD1) + LigandMPNN | `boltzdesign__multif/default.yaml` | ~600 (60 bb × 10) | **150** | Turing sm_75 OK | **~6.4 h med / 11 h max** → `gpu.24h`+`.120h` |
| HalluDesign (de-novo) + LigandMPNN | `halludesign/random_init.yaml` | ~150 (15 rec × 10) | **600** | **Ampere/Ada sm_80+ ONLY** | ~3–7 h |

**BoltzGen's ONE job feeds BOTH the boltzif and ligandmpnn stacks** (Route-B multif). Total
gen-uuids ≈ **1,043**. Per-stack config invariants (cys-free, seed handling, HD worktree +
3 symlinks, BD1 filter DISABLED, static `csy4_rna_only.pdb`) are **identical to v4** — see the v4
doc §1; nothing changed for v6 except HD seeds are **70001–70600**.

---

## 2. Array mechanism — the `--ntasks` trap (critical, UNCHANGED)

Identical to v4. `launch.py`'s native `--num-runs>1` array emits `#SBATCH --ntasks=4` →
`SLURM_NTASKS=4` → PyTorch-Lightning stacks (RFD3/BoltzGen/BD1) hard-abort. **Fix:** generate per-uuid
`run.sbatch`, cancel the broken native array, re-bucket with **`ops/submit/bucket_midrun_v6.py`**:
- groups per-uuid `run.sbatch` by routing signature `(partition, time, nodelist, gres, mem)` — so
  **routing is whatever the submit script baked into each `run.sbatch`**, not hardcoded in the bucketer.
- `#SBATCH --cpus-per-task=4` — **NO `--ntasks`** → `SLURM_NTASKS=1`.
- explicit `#SBATCH --qos=es_wicky/{base}/{tier}` derived from the partition, `--chdir=code`.
- inline `SCRIPTS=(…)` indexed by `$SLURM_ARRAY_TASK_ID`, dispatched with the guard
  `env -u SLURM_NTASKS -u SLURM_NPROCS -u SLURM_NTASKS_PER_NODE -u SLURM_TASKS_PER_NODE bash "$S"`.
- `--exclude=eu-g6-064` (the one black-hole node flagged this campaign).

HD gen is a **`--num-runs 1` loop with a unique per-uuid seed** (`submit_mv6_gen_hd.sh`), then bucketed the same way.

---

## 3. Routing (final, as-run) — **the main change from v4**

Turing stacks → **public `gpu.*`** (18× eu-g3 24G Turing, off the depleted gpuhe fairshare); HD →
**gpuhe** Ampere/Ada + **gpupr** A100. To land on the public pool: **`gpu:1` (drop `gpumem:24g`)** +
pin the eu-g3 nodelist + explicit QOS. **v6 change: the fast stacks moved to `gpu.4h`** — its QOS
caps ~96 GPU / 1770 G, ≈2× the 24 h tier, so RFD3/BG fold ~twice as fast there.

| Stack | Partition(s) | Nodelist | gres | Walltime |
|---|---|---|---|---|
| RFD3 (141) | **`gpu.4h`** | eu-g3-[062-077,083-084] (all 18) | gpu:1 | 04:00:00 |
| BoltzGen (152) | **`gpu.4h`** | eu-g3-[062-077,083-084] (all 18) | gpu:1 | 04:00:00 |
| BD1 (75+75) | `gpu.24h` / `gpu.120h` | eu-g3-[062-070] / eu-g3-[071-077,083-084] | gpu:1 | 23:30:00 / 48:00:00 |
| HD (400) | `gpuhe.24h` | eu-g4-[001-031],eu-g6-[001-080] | gpumem:24g | 23:30:00 |
| HD (200) | `gpupr.24h` | eu-a65-[02-07],eu-ts-[01-02] | gpumem:40g | 23:30:00 |

QOS = `es_wicky/{base}/{tier}` (e.g. `es_wicky/gpu/4`, `es_wicky/gpuhe/24`). HD is redistributed
across the **`.24h` and `.120h`** tiers of gpuhe/gpupr post-submit to multiply concurrency (each
partition-tier is a separate QOS memory budget — observed v6 spread: gpuhe.24h/120h + gpupr.24h/120h).
**HD is still hard fairshare-capped** (gpuhe ≈ 0.197) → ~40–80 concurrent regardless of config →
**HD is always the long pole (~1–2 days)**. No config/mem/pool lever beats it (an HD gpuhe QOS
mem-cap raise is **not available** — cannot be changed), so plan HD as the tail. **Never route BD1
to a100/gpupr** (reserved-forever ReqNodeNotAvail trap).

---

## 4. Refold system — node-local scratch + `.bulk` (UNCHANGED mechanism)

Node-local scratch workers `code/scripts/stages/refold_{af3,b2,esm}_scratch.sh` (require `--tmp=NG`,
all heavy I/O under `$TMPDIR/rf_$SLURM_JOB_ID`, copy back only `designs.tar.zst*` — the `*` copies
the `.list` manifest too) → **~6 Lustre inodes/leaf instead of ~6000**. This is *the* inode fix;
never run the old Lustre-scratch workers at scale.

**Refold routing** (`ops/refold/mv6_refold_bulk.py` `PLACE`), safe because node-local:
- **AF3 → `gpuhe.bulk`**, **ESM → `gpuhe.bulk`**, **B2 → `gpu.bulk`**, all `--time=06:00:00`,
  QOS `es_wicky/{gpuhe,gpu}/bulk`. Walltime must be set **at submit** (`scontrol update TimeLimit`
  can only LOWER a running job) — the `.bulk` branch bakes partition+qos+time. 6 h (was 4 h → AF3
  TIMEOUT).
- Each refold sbatch is `refold_X_scratch.sh`, gets `--tmp=30G`, `--exclude=eu-g6-064`, `--array=…%40`.
- IPSAE (`ipsae_*`) + tip-contact sidecars emitted at refold.

**⚠ `.bulk` submit caveat (observed 2026-07-08).** A one-off **single (non-array) job** submitted to
`gpuhe.bulk` with partition + `--qos=es_wicky/gpuhe/bulk` + `--time=06:00:00` **all baked in** was
still **rerouted by the submit plugin to `gpuhe.24h`** (QOS became `.../gpuhe/24`, blocked on
`QOSMaxMemoryPerUser`). Forcing bulk required **post-submit** `scontrol update jobid=<J>
partition=gpuhe.bulk qos=es_wicky/gpuhe/bulk` (same pool → no Features/gres change needed). The
refold **arrays** land on bulk directly, so the reroute appears specific to ad-hoc single-job submits
— use the scontrol fallback for those.

**Two auto-fill loops** (self-chaining `normal.4h`, ~600 s tick): `ops/refold/mv6_refold_loop.slurm`
(driver `mv6_refold_bulk.py`, all 5 stacks; STOP via `_mv6_refold_loop.STOP`) and
`hd_refold_bulk_loop` (driver `hd_refold_bulk.py`, reads the **shared** HD pool). HD packs per-uuid
(`<stack>/midrun_v6_<uuid>/designs.tar.zst`) → must be consolidated into the shared pool with
`ops/refold/hd_combine.sh midrun_v6 --go` (schedule it, don't wait to discover HD is invisible).

---

## 5. Lessons (carried from v4, still binding) + v6 additions

Carried forward (all still true): node-local scratch = the inode fix; `.list` copyback
(`designs.tar.zst*`); HD per-uuid → `hd_combine.sh`; black-hole `--exclude`; `.bulk` safe only with
node-local + requeue; `%40` refold throttle; array-aware gap-fill (`gen_gap_fill.py`); fairshare is
the HD wall; clean CSV needs all-3-refolds present. See v4 doc §5 for the reasoning.

**v6 additions:**
- **`gpu.4h` is the throughput lever for sub-4h stacks.** RFD3/BG on `gpu.4h` get ~2× the 24h-tier
  concurrency. Measure the true per-uuid runtime first (RFD3 24 m, BG 2.7 h max) and leave margin —
  BD1 at 6.4 h median genuinely can't use it.
- **`.bulk` reroute on single-job submits** (see §4 caveat) — bake it into the array path; scontrol
  fallback for one-offs.

---

## 6. Current state (2026-07-08, in progress)

Campaign underway across all 5 stacks: **674 packed design-slices** of the ~1,043 uuid target as of
2026-07-08 (RFD3/BG/BD1 on `gpu.4h`/`.24h`/`.120h`, HD generating under the fairshare cap). Not yet at
steady state — HD is the tail as always. Run `tally -r` + `gap_audit.py` for live counts.

---

## 7. File inventory (v6)

| Purpose | Path |
|---|---|
| Gen array bucketer | `ops/submit/bucket_midrun_v6.py` |
| Turing-stack gen submit | `ops/submit/submit_mv6_gen_turing.sh` (RFD3/BG → gpu.4h; BD1 → 24h/120h) |
| HD seeded gen submit | `ops/submit/submit_mv6_gen_hd.sh` (+ `ops/submit/mv6_hd_seed_manifest.txt`, seeds 70001–70600) |
| Array-aware gen gap-fill | `ops/submit/gen_gap_fill.py` (point `LABEL=midrun_v6`) |
| Gap audit (gen+refold vs queue) | `ops/submit/gap_audit.py` |
| Refold driver + placement | `ops/refold/mv6_refold_bulk.py` (af3/esm→gpuhe.bulk, b2→gpu.bulk, 6h) |
| Refold auto-loop | `ops/refold/mv6_refold_loop.slurm` |
| HD refold driver / loop | `ops/refold/hd_refold_bulk.py` |
| HD per-uuid → shared consolidator | `ops/refold/hd_combine.sh midrun_v6 --go` |
| Node-local scratch refold workers | `code/scripts/stages/refold_{af3,b2,esm}_scratch.sh` |
| CSV analysis pipeline | `code/scripts/run_analysis.sh` |
| Live tally | `code/scripts/tally.sh` (`tally -r` to rescan) |

**Infra loops that must be running:** `mv6_refold_loop`, `hd_refold_bulk_loop`, `inode_watcher`,
`hold_governor`, `gpu_watcher`.

---

## 8. Next-campaign checklist (bake it in from the start)

1. Set `LABEL=midrun_v7` in `bucket_midrun_v6.py`, `submit_mv7_gen_*.sh`, `gen_gap_fill.py`,
   `gap_audit.py`, `mv6_refold_bulk.py`, `hd_refold_bulk.py` (or clone to `*_v7.*`). Confirm
   analysis/tally auto-detect the new run.
2. **Pre-flight (unchanged configs):** `csy4_rna_only.pdb` present; HD worktree + 3 symlinks present;
   all 4 configs cys-free & walltime-correct; **fresh HD seeds** (don't reuse 70001–70600).
3. **Gen:** measure runtimes → route sub-4h stacks to `gpu.4h`, ≥4h stacks to `.24h`/`.120h`.
   `launch.py --num-runs N` per lane → cancel broken native arrays → `bucket_midrun_v7.py` → submit.
   **1-task/stack pre-flight test first.**
4. **Refold workers already node-local** — verify `refold_*_scratch.sh` + `--tmp`/`--exclude`/`%40`/
   `cp …tar.zst*`. Start `mv7_refold_loop` + `hd_refold_bulk_loop`.
5. **HD combine** on a timer from day 1.
6. **Black-hole excludes** on every array.
7. Monitor with `gap_audit.py` + `tally -r`; keep the inode watcher on.
8. `run_analysis.sh midrun_v7 <5 run dirs>` for the clean CSV **once BD1 refolds finish**.
