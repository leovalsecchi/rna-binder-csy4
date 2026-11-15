# midrun_v4 — AS-RUN conditions & lessons (the working recipe)

**Written 2026-07-05, after v4 reached steady state.** This is the *authoritative* record of
what midrun_v4 **actually ran as** — the launch conditions **plus** every operational fix we
discovered mid-run that made it work. It is the blueprint for **midrun_v5** (identical setup).

- Companion (as-planned, full config detail): `docs/MIDRUN_V4_LAUNCH_CONDITIONS_2026-07-03.md`.
- Repo root `RB=/cluster/work/wicky/lvalsecchi/rna_binder`; `code/` is **not** git-tracked.
- Target: `csy4_crrna_free` (16-nt csy4 crRNA epitope).

---

## 0. TL;DR — the recipe that worked

1. **Generate** per-uuid `run.sbatch` with `launch.py` (which uses `--cpus-per-task`), then
   **re-bucket** into self-contained Slurm arrays — one per (stack × lane). **Never `--ntasks`.**
2. **Route Turing-OK stacks** (RFD3 / BoltzGen / BD1) to the **public `gpu.*` eu-g3 pool**
   (off the depleted gpuhe fairshare); **HD** (sm_80+ only) to **gpuhe** (eu-g4/eu-g6 Ampere/Ada)
   + **gpupr** (A100).
3. **Refold on the `.bulk` tiers with node-local scratch workers** (`refold_*_scratch.sh`,
   `--tmp=30G`). This is *the* fix — it kept project inodes flat (~6 inodes/leaf instead of
   ~6000) so the campaign never wedged the Lustre inode quota.
4. **Two self-chaining loops** auto-fill refolds for every newly-packed uuid; **HD's per-uuid
   packs are consolidated** into the shared pool by `hd_combine.sh` so tally + the refolder see them.
5. **Exclude the known black-hole nodes** everywhere.

Outcome at steady state: **742 uuids, 375,220 designs, 1,053,950 refolds**, project inodes
**5.92M / 8.25M hard** (stable, never breached). 4 of 5 stacks done/near-done; HD the long pole
(fairshare-throttled).

---

## 1. Design allocation (equal ~90k/stack, 5 stacks / 4 generators)

| Stack | Config (`code/configs/…`) | Designs/uuid | Target uuids | GPU floor | Runtime/uuid |
|---|---|---|---|---|---|
| RFD3 + LigandMPNN | `rfd3__multif/default.yaml` | 640 (64 bb × 10 seq) | **141** | Turing sm_75 OK | ~50 min |
| BoltzGen + BoltzIF | `boltzgen__multif/default.yaml` | ~593 | **150** (shared) | Turing sm_75 OK | ~3.4 h |
| BoltzGen + LigandMPNN | *(same BG job, Route-B)* | ~595 | *(same 150)* | — | — |
| BoltzDesign (BD1) + LigandMPNN | `boltzdesign__multif/default.yaml` | ~600 (60 bb × 10) | **150** | Turing sm_75 OK | ~7 h → **24h tier** |
| HalluDesign (de-novo) + LigandMPNN | `halludesign/random_init.yaml` | ~150 (15 rec × 10) | **600** | **Ampere/Ada sm_80+ ONLY** | ~3–7 h |

**BoltzGen's ONE job feeds BOTH the boltzif and ligandmpnn stacks** (Route-B multif) → 150 BG
uuids ≈ 90k boltzif + 90k limpnn. Total gen-uuids ≈ **1,043**.

### Key per-stack config invariants (all cys-free; verified, unedited)
- **RFD3**: `n_batches=8 × diffusion_batch_size=8` → 64 bb, len 80; `mpnn_cli.omit_AA=C`.
  `chains_to_design=A` is **dead** — `rfd3_multif.sh` auto-detects the real protein chain at
  runtime (the seq-collapse fix; validated 0 cys / 13,440 designs in v3).
- **BoltzGen**: `protocol=protein-anything`, `num_designs=60`, `bif_cli.inverse_fold_avoid=C`,
  `mpnn_cli.omit_AA=C`. **Requires static input** `targets/csy4_crrna_free/csy4_rna_only.pdb`
  (521 atoms / chain B / res 4-21) — `cp`'d by the stage script; regen from the .cif if missing.
- **BD1**: `design_samples=60`, len 80, **`i_ptm_cutoff=0.0` + `complex_plddt_cutoff=0.0`
  (filter DISABLED — must stay 0 or ~0 designs survive)**, `omit_AA=C`. **Never route to
  a100/gpupr** (reserved-forever ReqNodeNotAvail trap).
- **HD**: `random_init.yaml` (**NOT** `default.yaml` = old helix mode). Runs from git worktree
  `/cluster/work/wicky/lvalsecchi/halludesign_upgrade` (branch `upgrade-1b08312`, commit
  `42ccb00`, de-novo). `random_init_chain_spec="A:80"`, `num_designs=1` (**must stay 1**),
  `rna=CACUGCCGUAUAGGCAGC` (18-nt), `num_seqs=10`, `num_recycles=15`. **Seed overridden per
  uuid** (50001–50600) — `--num-runs N` reuses one seed → identical designs. Needs its 3 runtime
  symlinks (`LigandMPNN/model_params`, `Protenix/release_data`, `cutlass`).

---

## 2. Array mechanism — the `--ntasks` trap (critical)

`launch.py`'s native `--num-runs>1` array path emits `#SBATCH --ntasks=4` (L925 / HD builder
L209). That exports `SLURM_NTASKS=4`, which makes **PyTorch-Lightning stacks (RFD3, BoltzGen,
BD1) hard-abort** (`RuntimeError: --ntasks=4 … not supported`). HD tolerates it (not Lightning-DDP).

**The fix (proven pattern):**
1. `launch.py --num-runs N` to *generate* per-uuid `run.sbatch` (uses `--cpus-per-task`, safe) —
   then **cancel** the broken native array it submits.
2. Re-bucket into one self-contained array per (stack × lane), via **`ops/submit/bucket_midrun_v4.py`**:
   - `#SBATCH --cpus-per-task=4` — **NO `--ntasks`** → SLURM_NTASKS=1
   - explicit `#SBATCH --qos=es_wicky/{base}/{tier}`, `--chdir=<code>`
   - inline `SCRIPTS=(…)` indexed by `$SLURM_ARRAY_TASK_ID`
   - dispatch with the guard:
     `env -u SLURM_NTASKS -u SLURM_NPROCS -u SLURM_NTASKS_PER_NODE -u SLURM_TASKS_PER_NODE bash "$S"`
   - `--exclude=eu-g6-064` on gpuhe/gpupr lanes (black-hole)

HD gen itself is run as a **`--num-runs 1` loop with a unique per-uuid seed** (flat run dirs
`<stack>/midrun_v4_<uuid>/`), then those run.sbatch are bucketed the same way.

---

## 3. Routing (final, as-run)

Turing stacks → **public `gpu.*`** (18× eu-g3 24G Turing, **off the depleted gpuhe fairshare**);
HD → **gpuhe** Ampere/Ada + **gpupr** A100. To land on the public pool: **drop `gpumem:24g` →
`gpu:1`** + pin the eu-g3 nodelist + explicit QOS. Spreading over `.24h`/`.120h` (and gpuhe
`.4h`/`.120h` for redistribution) multiplies concurrency — **each partition-tier is a separate
QOS memory budget**, so more tiers = more parallel slots.

| Stack | Partition(s) | Nodelist | gres | QOS |
|---|---|---|---|---|
| RFD3 (141) | `gpu.24h` | eu-g3-[062-077,083-084] | gpu:1 | es_wicky/gpu/24 |
| BoltzGen (76+76) | `gpu.24h` / `gpu.120h` | eu-g3-[062-070] / eu-g3-[071-077,083-084] | gpu:1 | gpu/24, gpu/120 |
| BD1 (75+75) | `gpu.24h` / `gpu.120h` | (same eu-g3 split) | gpu:1 | gpu/24, gpu/120 |
| HD (400) | `gpuhe.24h` | eu-g4-[001-031],eu-g6-[001-080] | gpumem:24g | es_wicky/gpuhe/24 |
| HD (200) | `gpupr.24h` | eu-a65-[02-07],eu-ts-[01-02] | gpumem:40g | es_wicky/gpupr/24 |

`gpu.120h` is forced by a >24h walltime (48:00:00). **HD is fairshare-capped** (gpuhe ≈ 0.197)
→ ~40–80 concurrent regardless of config → **HD is always the long pole (~1–2 days)**. Only real
speed-ups: RawShares admin bump, or lean harder on gpupr (0.42 fairshare).

---

## 4. Refold system — node-local scratch (THE lesson)

**Problem in earlier runs:** each refold leaf wrote ~6000 tiny files to Lustre; at campaign scale
this blew the **project inode quota** (`lfs quota -p 510551 /cluster/work/wicky`, soft 7.5M /
hard 8.25M) and every pack started failing.

**The fix — node-local scratch workers** (`code/scripts/stages/refold_{af3,b2,esm}_scratch.sh`):
- Require `--tmp=NG` → SLURM provisions per-job node-local NVMe at `$TMPDIR`
  (`/scratch/tmp.<jobid>.lvalsecchi`, 452 GB, auto-deleted on job end, **doesn't touch Lustre
  quota**).
- `WORK="$TMPDIR/rf_${SLURM_JOB_ID}"`; **all** heavy I/O (INPUTS / OUTPUTS / CANONICAL_DIR,
  Triton/CUEQ caches) lives under `$WORK`.
- `apptainer exec … --bind "$TMPDIR"` so the container sees node-local scratch.
- **Copy back only the packed tarball**: `cp -f "$WORK"/designs.tar.zst* "$REFOLD_DIR/"` — the
  `*` is load-bearing: it copies **both** `designs.tar.zst` **and** its `.list` manifest.
  → **~6 Lustre inodes/leaf instead of ~6000.**

**Final refold routing** (`ops/refold/mv4_refold_bulk.py` `PLACE`), safe *because* node-local:
- **All three refold backends submit *directly* to their `.bulk` partition at `--time=06:00:00`**:
  `af3`/`esm → gpuhe.bulk`, `b2 → gpu.bulk`. The base sbatch is gpu.4h/4h (rejects >4h); direct
  submit works because both `gpu.bulk` and `gpuhe.bulk` are MaxTime=15d + AllowQos=ALL and the
  `es_wicky/{gpu,gpuhe}/bulk` QOS have no MaxWall. **`scontrol update TimeLimit` cannot RAISE** a
  running job’s limit (users may only lower it), so the walltime must be set *at submit* — hence
  the direct-to-bulk branch (`elif part.endswith(".bulk") and wall:` bakes partition+qos+time).
  Observed a rfaf3 4h TIMEOUT before the bump. (AF3 was on public `gpu.*` at launch; moved to
  `gpuhe.bulk` once node-local made bulk inode-safe and preemption harmless.)
- Every refold sbatch is rewritten `refold_X.sh → refold_X_scratch.sh`, gets `--tmp=30G`,
  `--exclude=eu-lo-g2-028,eu-lo-g3-059,eu-g6-064`, and **`--array=…%40`** (throttle — see §5).
- IPSAE (`ipsae_*`) + tip-contact sidecars emitted at refold.

**Two auto-fill loops** (self-chaining `normal.4h`, ~600 s tick):
- `ops/refold/mv4_refold_loop.slurm` (job driver `mv4_refold_bulk.py`) — all 5 stacks; STOP via
  `_mv4_refold_loop.STOP`.
- `hd_refold_bulk_loop` (driver `hd_refold_bulk.py`) — HD-specific; reads the **shared** HD pool.
- Each tick submits AF3+ESM+B2 for every packed uuid missing that backend. Newly-packed uuids are
  picked up on the next tick (normal fresh-pack lag, not a gap).

---

## 5. Lessons learned (mid-run fixes, in priority order)

1. **Node-local scratch = the inode fix** (§4). Without it the campaign wedges the project quota.
   Bake it in from day 1 for v5; never run the old Lustre-scratch refold workers at scale.
2. **`.list` sidecar must be copied back** (`designs.tar.zst*`). tally counts refolds from the
   `.list`; if you copy only the tarball, refolds show **0** in tally though the data is fine.
3. **HD per-uuid layout → must be combined.** HD packs `<stack>/midrun_v4_<uuid>/designs.tar.zst`;
   everything else uses shared `<stack>/midrun_v4/designs/<uuid>.tar.zst`. Run
   **`ops/refold/hd_combine.sh midrun_v4 --go`** (periodically) to consolidate — otherwise HD is
   invisible to both tally and the refolder. (Consider cron-ing this for v5.)
4. **Black-hole nodes** (`eu-lo-g2-028, eu-lo-g3-059, eu-g6-064`) wedge GPUs → fast-fail
   `CUDA unknown error`. `--exclude` them on **every** array (gen + refold).
5. **`.bulk` is safe *only* with node-local scratch + requeue.** Preemption then just requeues
   with no loose Lustre files. AF3 on bulk without node-local = ~12% success (don't).
6. **`%40` array throttle** on refolds — stops a whole batch hitting CPU-bound Stage-3
   (canonicalize/pack) at once (thundering-herd contention that stalled B2).
7. **Gen gap-fill must be array-aware** (`ops/submit/gen_gap_fill.py`): only resubmit uuids with a
   `run.sbatch`, **no packed slice, AND not covered by a live array** (parse each live array's
   `SCRIPTS=()` + direct `run.sbatch` from `squeue -o %o`). Otherwise you double-submit. In v4 the
   "hundreds missing" were almost all *in-flight in live arrays* (fairshare-pending) — only ~18
   genuinely fell out.
8. **Fairshare is the HD wall**, not config. No walltime/mem/pool lever beats it. Plan HD as the
   ~1–2 day tail.
9. **Clean CSV needs all-3-refolds present** — run the analysis pipeline (§7) *after* refolds
   fill, or BD1 (slowest AF3) is under-represented in the clean set.

---

## 6. Outcome & health (steady state, 2026-07-05)

| stack | uuids | designs | AF3 | ESM | B2 |
|---|---|---|---|---|---|
| rfd3/limpnn | 141 | 90,240 | ~99% | ~99% | 100% |
| boltzgen/boltzif | 152 | 91,200 | ~99% | 100% | 100% |
| boltzgen/limpnn | 152 | 91,200 | ~96% | 100% | 100% |
| boltzdesign/limpnn | ~134 | ~80k | filling | ~85% | ~85% |
| halludesign/limpnn | 167→600 | growing | 100% of packed | 100% | 100% |
| **TOTAL** | **742** | **375,220** | **326,580** | **363,552** | **363,818** |

Refolds total **1,053,950**. Project inodes **5.92M / 8.25M** (never breached — node-local held).
Gen genuine gaps: **0** (all done or in a live array). Refold genuine gaps: **0** stuck (only
fresh-pack lag). HD still generating under the fairshare cap = the remaining wall-clock.

---

## 7. File inventory (reuse verbatim for v5, swap the label)

| Purpose | Path |
|---|---|
| Gen array bucketer | `ops/submit/bucket_midrun_v4.py` |
| Turing-stack gen submit | `ops/submit/submit_mv4_gen_turing.sh` |
| HD seeded gen submit | `ops/submit/submit_mv4_gen_hd.sh` (+ `mv4_hd_seed_manifest.txt`) |
| Array-aware gen gap-fill | `ops/submit/gen_gap_fill.py` |
| Gap audit (gen+refold vs queue) | `ops/submit/gap_audit.py` |
| Refold driver + placement | `ops/refold/mv4_refold_bulk.py` |
| Refold auto-loop | `ops/refold/mv4_refold_loop.slurm` |
| HD refold driver / loop | `ops/refold/hd_refold_bulk.py` |
| HD per-uuid → shared consolidator | `ops/refold/hd_combine.sh` |
| Node-local scratch refold workers | `code/scripts/stages/refold_{af3,b2,esm}_scratch.sh` |
| CSV analysis pipeline | `code/scripts/run_analysis.sh` (collect → clean → sample metrics) |
| Live tally | `code/scripts/tally.sh` (`tally -r` to rescan) |

**Infra loops that must be running:** `mv4_refold_loop`, `hd_refold_bulk_loop`, `inode_watcher`,
`hold_governor`, `gpu_watcher`.

---

## 8. midrun_v5 checklist (bake the lessons in from the start)

Everything below is identical to v4 except the label — the point of v5 is **zero mid-run
firefighting** because every fix is already in place before launch.

1. Set `LABEL=midrun_v5` in `bucket_midrun_v4.py`, `gen_gap_fill.py`, `gap_audit.py`,
   `mv4_refold_bulk.py`, `hd_refold_bulk.py` (or clone to `*_v5.*`). Confirm the analysis/tally
   auto-detect the new run.
2. **Pre-flight (unchanged configs):** confirm `csy4_rna_only.pdb` present; HD worktree +
   3 symlinks present; all 4 configs still cys-free & walltime-correct.
3. **Gen:** `launch.py --num-runs N` per Turing lane + HD seeded loop → cancel the broken native
   arrays → `bucket_midrun_v4.py` → submit. **1-task/stack pre-flight test first.**
4. **Refold workers already node-local** — verify `refold_*_scratch.sh` are the ones referenced
   and `--tmp`/`--exclude`/`%40`/`cp …tar.zst*` are all present (§4). Start `mv4_refold_loop` +
   `hd_refold_bulk_loop`.
5. **HD combine:** schedule `hd_combine.sh midrun_v5 --go` on a timer (don't wait to discover HD
   is invisible).
6. **Black-hole excludes** on every array from the start.
7. Monitor with `gap_audit.py` (gen+refold vs queue) and `tally -r`; keep the inode watcher on.
8. Run `run_analysis.sh midrun_v5 <5 run dirs>` for the clean CSV **once BD1 refolds finish**.
