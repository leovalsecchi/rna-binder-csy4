# midrun_v4 — exact launch conditions (2026-07-03)

Campaign label: **`midrun_v4`**. Target: **450,000 designs**, **5 stacks**, equal **90k/stack**.
Launched as **Slurm arrays** via the proven external bucketer (`--cpus-per-task`, never
`--ntasks`). Target: `csy4_crrna_free` (16-nt csy4 crRNA epitope). Repo root
`RB=/cluster/work/wicky/lvalsecchi/rna_binder`; `code/` is NOT git-tracked.

All facts below were verified live from the configs + stage scripts + memory before launch
(5-agent audit). No pipeline code was changed; no config files were edited.

---

## 1. Design allocation (equal 90k/stack)

| Stack | Config | Designs/uuid | uuids | Yield |
|---|---|---|---|---|
| RFD3 (RFdiffusion3 + LigandMPNN) | `rfd3__multif` | 640 (64 bb × 10 seq) | **141** | ~90k |
| HalluDesign (de-novo + LigandMPNN) | `halludesign/random_init.yaml` | 150 (15 recycles × 10 seq) | **600** | ~90k |
| BoltzGen + BoltzIF | `boltzgen__multif` | ~593/uuid | **152** (shared) | ~90k |
| BoltzGen + LigandMPNN | `boltzgen__multif` (same job, Route-B) | ~595/uuid | (same 152) | ~90k |
| BoltzDesign (BD1) + LigandMPNN | `boltzdesign__multif` | ~600 (60 bb × 10 seq) | **150** | ~90k |

**Total gen-uuids: 1,043** (441 Turing-stack array tasks + 600 HD + 1 test each = the 4 test
uuids from pre-flight run to completion as part of the campaign; see §6). BoltzGen's ONE job
feeds BOTH boltzif + ligandmpnn stacks (Route-B multif).

---

## 2. Exact per-stack config (the yaml files, unedited)

### RFD3 — `code/configs/rfd3__multif/default.yaml`
- `generator_params`: `n_batches=8`, `diffusion_batch_size=8` → **64 backbones**; length **80** (launch.py default). Diffusion seed = random per uuid (distinct backbones).
- `generator_params.mpnn_cli`: `omit_AA=C` (cys-free), `chains_to_design=A` **← dead/ignored; `rfd3_multif.sh` auto-detects the real protein chain at runtime** (the seq-collapse fix), `batch_size=10`, `number_of_batches=1` → 10 seq/bb, `temperature=0.1`, `seed=42`.
- container: `foundry.sif` (+ `LigandMPNN.sif`, weights `ligandmpnn_v_32_010_25.pt`).
- slurm defaults (overridden at launch): `rtx_4090:1`, `gpuhe.4h`, `04:00:00`, cpus=4, mem_per_cpu=5G. **Uses `--cpus-per-task` (SLURM_NTASKS=1).**
- Cys-free: **VALIDATED** — midrun_v3 = 0 cys / 13,440 designs. GPU floor: **Turing sm_75 OK**. Runtime ~50 min/uuid.

### HalluDesign — `code/configs/halludesign/random_init.yaml`  (⚠ NOT `default.yaml`, which is old helix mode)
- `generator_params.workdir`: `/cluster/work/wicky/lvalsecchi/halludesign_upgrade` (git worktree, branch `upgrade-1b08312`, commit `42ccb00` = de-novo build). Runs de-novo straight out of the worktree (NOT promoted to the live `~/halludesign` = 25f1045). Its 3 runtime symlinks verified present pre-launch: `LigandMPNN/model_params`, `Protenix/release_data`, `cutlass`.
- `generator_params.params`: `HalluDesign_model=cross_model`, **`random_init_chain_spec="A:80"`** (de-novo, no input PDB), `num_designs=1` (MUST stay 1), `seed` **overridden per-uuid** (see §6), `rna=CACUGCCGUAUAGGCAGC` (**18-nt full csy4 crRNA**), `num_seqs=10`, `num_recycles=15`, `ref_time_steps=60`, `mpnn=ligand_mpnn`, `mpnn_num_packs=1`, `fix_chain_index=B`.
- container: `af3_protenix_esm.sif` (`apptainer run --nv`). Uses ESM2 LM internally (NOT ESMFold2).
- slurm defaults (overridden): `a100_80gb:1`, `ntasks=4`, `mem_per_cpu=7G`, `gpuhe.24h`, `23:30:00`. Per-uuid run.sbatch carries `--ntasks=4` but HD is NOT Lightning-DDP (tolerates it); once bucketed the array's `--cpus-per-task` governs.
- Cys-free: **YES** (0% in midrun_v2, LigandMPNN default omit). GPU floor: **Ampere/Ada sm_80+ ONLY — NO Turing** (silent 19KB junk on Turing). Runtime ~3–7 h/uuid.

### BoltzGen — `code/configs/boltzgen__multif/default.yaml`  (Route-B → 2 stacks)
- `generator_params.bg_cli`: `protocol=protein-anything`, **`num_designs=60`**, `diffusion_batch_size=1`. Spec: protein chain D len 80 + RNA chain B from `csy4_rna_only.pdb`.
- `generator_params.bif_cli`: **`inverse_fold_avoid=C`** (cys-free BIF), `inverse_fold_num_sequences=10`.
- `generator_params.mpnn_cli`: **`omit_AA=C`** (cys-free LiMPNN), `chains_to_design=A`, `batch_size=10`, `number_of_batches=1`, `temperature=0.1`, `seed=42`.
- container: `boltzgen.sif` (+ `LigandMPNN.sif`). **Requires static input** `targets/csy4_crrna_free/csy4_rna_only.pdb` (521 atoms / chain B / res 4-21) — verified present pre-launch (restored 2026-07-03 from the .cif).
- slurm defaults (overridden): `rtx_4090:1`, `gpuhe.24h`, `23:30:00`, cpus=4, mem_per_cpu=16G. **Uses `--cpus-per-task`.**
- Cys-free: **VALIDATED** this session (bif_fill ran on this exact config). GPU floor: **Turing sm_75 OK**. Runtime ~3.4 h/uuid at num_designs=60.

### BoltzDesign (BD1) — `code/configs/boltzdesign__multif/default.yaml`
- `generator_params`: `target_name=csy4_rna`, `pdb_path=.../csy4_rna_only.pdb`, `target_type=rna`, `pdb_target_ids=B`, **`design_samples=60`**, `length_min=length_max=80`, `inter_chain_cutoff=10`, **`i_ptm_cutoff=0.0` + `complex_plddt_cutoff=0.0` (filter DISABLED — must stay 0, else ~0 designs)**, `run_ligandmpnn/alphafold/rosetta=False`, `boltz_checkpoint=boltz1_conf.ckpt`, `ccd_path=ccd.pkl`.
- `generator_params.mpnn_cli`: `chains_to_design=A`, `batch_size=10`, `number_of_batches=1`, `temperature=0.1`, `seed=42`, **`omit_AA=C`**.
- container: `boltz_design1.sif` (+ `LigandMPNN.sif`).
- slurm defaults (overridden): `rtx_4090:1`, `gpuhe.24h`, `23:30:00`, cpus=4, mem_per_cpu=5G, tmp=25G. **Uses `--cpus-per-task`.**
- Cys-free: **YES** (0% midrun_v2). GPU floor: **Turing sm_75 OK** (also only stack that runs sm_120). Runtime ~7 h/uuid → 24h tier mandatory. **Never route to a100/gpupr** (reserved-forever trap).

---

## 3. Routing (throughput-max: Turing stacks → off-fairshare public pool; HD → sm_80+ pools)

Overrides applied per lane via `--override slurm.{gpus,gres,partition,time,nodelist}` (+ HD seed).

| Stack | Lane | Partition | Nodelist | gres | Walltime | QOS |
|---|---|---|---|---|---|---|
| RFD3 (141) | public | `gpu.24h` | eu-g3-[062-077,083-084] | gpu:1 | 23:30:00 | es_wicky/gpu/24 |
| BoltzGen (76+76) | public A/B | `gpu.24h` / `gpu.120h` | eu-g3-[062-070] / eu-g3-[071-077,083-084] | gpu:1 | 23:30 / 48:00 | gpu/24, gpu/120 |
| BD1 (75+75) | public A/B | `gpu.24h` / `gpu.120h` | eu-g3-[062-070] / eu-g3-[071-077,083-084] | gpu:1 | 23:30 / 48:00 | gpu/24, gpu/120 |
| HD (400) | owner | `gpuhe.24h` | eu-g4-[001-031],eu-g6-[001-080] | gpumem:24g | 23:30:00 | es_wicky/gpuhe/24 |
| HD (200) | overflow | `gpupr.24h` | eu-a65-[02-07],eu-ts-[01-02] | gpumem:40g | 23:30:00 | es_wicky/gpupr/24 |

Rationale: public `gpu.*` (18× eu-g3 24G Turing) is OFF the depleted gpuhe fairshare → Turing
stacks (RFD3/BG/BD) go there; HD (sm_80+ only) gets gpuhe Ampere/Ada + gpupr A100 (fresh 0.42
fairshare). To land on the public pool: **drop gpumem:24g → `gpu:1`** + pin eu-g3 + explicit QOS.
`gpu.120h` forced by a >24h walltime (48h). `--exclude=eu-g6-064` (black-hole node) on gpuhe lanes.

---

## 4. Array mechanism (the critical detail)

**Never used `launch.py`'s native `--num-runs>1` array path** — it emits `#SBATCH --ntasks=4`
(L925), which exports `SLURM_NTASKS=4` and makes PyTorch-Lightning stacks (RFD3, BoltzGen, BD1)
hard-abort. Instead, the proven pattern (== midrun_v3 `convert_to_arrays.sh` / this session's
`rebucket_bif_fill.py`):

1. Generate per-uuid `run.sbatch` via `launch.py` (which uses `--cpus-per-task`, safe).
2. Re-bucket into one self-contained array per (stack × lane) with:
   - `#SBATCH --cpus-per-task=4` (**NO `--ntasks`** → SLURM_NTASKS=1)
   - explicit `#SBATCH --qos=es_wicky/{base}/{tier}`
   - `#SBATCH --chdir=<code>`, `--exclude=eu-g6-064` on gpuhe/gpupr
   - inline `SCRIPTS=(...)` indexed by `$SLURM_ARRAY_TASK_ID`
   - dispatch: `env -u SLURM_NTASKS -u SLURM_NPROCS -u SLURM_NTASKS_PER_NODE -u SLURM_TASKS_PER_NODE bash "$S"`

Bucketer: `ops/submit/bucket_midrun_v4.py`. Pre-flight test (1 task/stack) confirmed all 5 stacks
run clean via this mechanism (zero `--ntasks`/Lightning errors).

---

## 5. Generation method

- **RFD3 / BoltzGen / BD1** (`ops/submit/submit_mv4_gen_turing.sh`): `launch.py --num-runs N` per
  lane (creates N nested run dirs `<stack>/midrun_v4/run_<uuid>/`), then **cancelled** the broken
  native array launch.py submits. Diffusion self-randomizes → distinct backbones per uuid.
  Cancelled broken arrays: 5649276, 5649303, 5649316, 5649317, 5649326.
- **HalluDesign** (`ops/submit/submit_mv4_gen_hd.sh`): `--num-runs 1` LOOP with a **unique seed per
  uuid** (`--override generator_params.params.seed=$SEED`, seeds **50001–50600**), flat run dirs
  `<stack>/midrun_v4_<uuid>/`, single job cancelled per iteration. (Cannot use `--num-runs N` — it
  reuses one seed → identical designs.) Seed manifest: `ops/submit/mv4_hd_seed_manifest.txt`.

---

## 6. Submitted arrays (gen)

Pre-flight test uuids (running to completion as first of campaign): rfd3 `4b1043ef`, boltzgen
`dc7fe23b`, bd1 `fcdc1a0d`, hd `46db14f3` (test arrays 5647814/815/816 + 5647505).

Turing arrays (submitted 2026-07-03):

| Array | Stack | Partition | Tasks |
|---|---|---|---|
| **5649974** | rfd3 | gpu.24h | 141 |
| **5649977** | boltzgen | gpu.24h | 76 |
| **5649976** | boltzgen | gpu.120h | 76 |
| **5649979** | bd1 | gpu.24h | 75 |
| **5649978** | bd1 | gpu.120h | 75 |

HD arrays (submitted after its 600-uuid seeded gen completed):

| Array | Stack | Partition | Tasks |
|---|---|---|---|
| **5651620** | hd | gpuhe.24h | 400 |
| **5651619** | hd | gpupr.24h | 200 |

Refold loop: **5650147** (`mv4_refold_loop`).

**Total: 9 gen arrays, 1,043 uuids** (+ 4 pre-flight test uuids running to completion).

---

## 7. Refold trigger

- `ops/refold/mv4_refold_bulk.py` (clone of bif_fill; `DEFAULT_RUN=midrun_v4`, `--stack` arg).
- `ops/refold/mv4_refold_loop.slurm` — self-chaining `normal.4h` loop, every 600 s submits refolds
  for any newly-packed uuid across ALL 5 stacks. Placement: **AF3 → public gpu.24h/.120h (eu-g3
  24G); ESM → gpuhe.bulk; B2 → gpu.bulk.** STOP via `_mv4_refold_loop.STOP` sentinel.
- Started: job **5650147**. Idles until designs pack (no auto-refold at gen, per policy).

---

## 8. Pre-flight verifications (all ✅)
- BoltzGen/BD1 target `csy4_rna_only.pdb` present (521 atoms, chain B, res 4-21).
- HD worktree `halludesign_upgrade` + its 3 runtime symlinks present.
- All 4 configs cys-free and walltime-correct as-is (no edits).
- Array mechanism validated by 1-task/stack test — all 5 stacks generate, zero `--ntasks` errors.
- `aggregate.py` carries the gen/invfold seq-key fix + AF3 `_gap` regex fix (use for v4 aggregation).

## 9. Caveats / watch-items
- **Public pool contention**: bif_fill AF3 refolds are camping the same eu-g3 public nodes →
  midrun_v4 Turing arrays queue behind them; drains as bif_fill winds down.
- **gpupr A100 trap**: HD gpupr.24h lane may stick on `ReqNodeNotAvail` if A100s are reserved →
  reroute those uuids to gpuhe.24h if they don't schedule.
- **HD is the wall-clock long pole** (600 uuids, fairshare-capped ~40–80 concurrent) → ~2–3 days.
- **Black-hole node** eu-g6-064 excluded on gpuhe/gpupr lanes.
- HD `COMPLETED ≠ success` — gap-audit every HD uuid for a 150-line slice; resubmit gaps w/ a new seed.
- **Inode safety**: `ops/watchers/inode_watcher` (job 5635976) blanket-holds ALL pending jobs at
  **7.4M** files, releases at **7.2M** (raised from 7.0M this launch — baseline was ~7.05M and rising,
  so the old 7.0M release line would have dead-locked a hold). Current ~7.05M / soft 7.5M / hard 8.25M.
  Refold loops excluded from holds. Shut off: `echo 0 > ops/watchers/inode_watcher/enabled`.

## 10. Key paths
- Gen: `ops/submit/submit_mv4_gen_turing.sh`, `ops/submit/submit_mv4_gen_hd.sh`
- Bucketer: `ops/submit/bucket_midrun_v4.py`; arrays under `campaigns/csy4_crrna_free/_mv4_arrays/`
- Refold: `ops/refold/mv4_refold_bulk.py`, `ops/refold/mv4_refold_loop.slurm`
- HD seeds: `ops/submit/mv4_hd_seed_manifest.txt`
