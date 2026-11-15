# GPU compatibility — AF3, Boltz-2 & ESMFold2 refolding on Euler

Which Euler GPU types can run the rna_binder refold workers (AlphaFold3 and
Boltz-2). Tested 2026-06-10 by refolding **one** design from `3kseq` (rfd3),
MSA-off, via the production `refold_af3.sh` / `refold_b2.sh` workers.

> Status: **12/14 confirmed.** Only `a100_80gb` is still queued (its gpupr nodes
> are under a reservation — `ReqNodeNotAvail`); it is the same Ampere sm_80 as the
> `a100`, which passed both, so both are expected to pass.

## Result matrix

| GPU type (`--gpus`) | VRAM | Arch (CC) | Partition | AF3 | Boltz-2 |
|---|---|---|---|---|---|
| `rtx_2080_ti`      | 11 GB | Turing (sm_75) | gpu.4h   | ❌ **fail** | ✅ ok (2:26) |
| `nvidia_titan_rtx` | 24 GB | Turing (sm_75) | gpu.4h   | ❌ **fail** | ✅ ok (2:44) |
| `rtx_6000`         | 24 GB | Turing (sm_75) | gpu.4h   | ❌ **fail** | ✅ ok (2:00) |
| `rtx_3090`         | 24 GB | Ampere (sm_86) | gpuhe.4h | ✅ ok (2:34) | ✅ ok (3:21) |
| `rtx_4090`         | 24 GB | Ada (sm_89)    | gpuhe.4h | ✅ ok (2:10) | ✅ ok (2:44) |
| `a100`             | 40 GB | Ampere (sm_80) | gpupr.4h | ✅ ok (3:19) | ✅ ok (3:11) |
| `a100_80gb`        | 80 GB | Ampere (sm_80) | gpupr.4h | ⏳ queued¹ | ⏳ queued¹ |
| `rtx_pro_6000`     | 96 GB | Blackwell (sm_120) | cuda13pr.* | ⛔ skipped | ⛔ skipped |

✅ = COMPLETED · ❌ = FAILED · ⏳ = not finished yet · ⛔ = not tested (see below).
¹ nodes reserved at test time; expected ✅/✅ (identical sm_80 to `a100`).

**Bottom line:** use **Ampere or newer** for AF3 (rtx_3090, rtx_4090, a100,
a100_80gb). Boltz-2 runs on **any** of these GPUs, Turing included. Avoid the three
Turing cards (rtx_2080_ti, nvidia_titan_rtx, rtx_6000) for AF3, and avoid
rtx_pro_6000 (Blackwell) for both.

## Key finding

**AF3 requires Ampere or newer (compute capability ≥ sm_80).** On all three Turing
cards (`rtx_2080_ti`, `nvidia_titan_rtx`, `rtx_6000`, all sm_75) AF3 fails
immediately — within ~1 min — with the identical error:

```
ValueError: implementation='triton' for FlashAttention is unsupported on this
GPU generation. Please use implementation='xla' instead.
```

AF3's default Triton FlashAttention kernel is Ampere+. **This is now resolved via the
`flash_attention_implementation=xla` path** — `refold_af3.sh` auto-picks `xla` when
`nvidia-smi` reports compute_cap < 80, so AF3 runs on the Turing 24 GB cards
(`quadro_rtx_6000`, `nvidia_titan_rtx`). The required CC-7.x flag
`XLA_FLAGS=--xla_disable_hlo_passes=custom-kernel-fusion-rewriter` was already set.

**Validated 2026-06-11** (10-design, 3-arm refold — jobs 2999125 / 2999126 / 2999128):
titan_rtx + quadro_rtx_6000 confidences match the 3090 within AF3's own sampling noise
(ranking_score Δ≤0.04, ptm/iptm ≤0.02, no pass/fail flips) **and ran faster** than the
3090 (11–13 min vs 24 min for 10 designs). **This is now the DEFAULT for AF3
refolds** — both the `refold__af3` config nodelist and the GUI default (`turing_24g`)
span all usable ≥24 GB cards (3090/4090 + titan/quadro); pick `ampere_24g` to
restrict to 3090/4090.
The 11 GB `rtx_2080_ti` stays excluded by the 24 GB VRAM floor (not by architecture).
The `❌ FAIL` marks in the table above are *out-of-the-box* (default triton) behaviour.

**Boltz-2 runs on every GPU tested**, including the 11 GB Turing `rtx_2080_ti` —
Boltz uses different attention kernels and the design is small (80-aa protein +
16-nt RNA), so VRAM is not the limiter at this size.

`rtx_pro_6000` (Blackwell, sm_120) was **deliberately not tested**: PyTorch in the
current containers has no sm_120 support, so it cannot run AF3 or Boltz-2 as-is.

## ESMFold2 (`refold_esm`, added 2026-06-17)

Local ESMFold2 (`esmc.sif`, `biohub/ESMFold2` + ESMC-6B backbone) via the
`refold_esm.sh` worker. **Needs Ampere+ (sm_80+)** — its fused pair-bias attention
is a Triton kernel with the *same* sm_80 floor as AF3's, so the Turing 24 GB cards
fail (no working xla fallback for ESM). Config + GUI pin `ampere_24g` (3090/4090);
a100 also works. VRAM ~21.5 GB on a 24 GB card for csy4-size complexes (16 nt RNA
+ ~80–150 aa binder); much larger complexes can OOM.

**Validated 2026-06-17** on `rtx_3090` (job 3743276, 3 real csy4 designs): one-time
startup ~200 s (25 GB node-local weight stage + ESMC-6B load + Triton compile),
then **~4 s/design** steady-state → 600 designs ≈ 43 min, fits `gpuhe.4h`. Only the
3090 has been run; 4090/a100 (also sm_80+) expected to pass. **Key gotcha:** the
worker bulk-copies the 25 GB weights cache to node-local `$TMPDIR` first — mmap'ing
them off Lustre crawls at ~1 MB/s on busy nodes (see `ESMFOLD2_SPEEDFIX.md`).

## Method

- **Input:** one design `csy4free_RFD3_LiMPNN_len80_3ffd6899_design0_seq0.cif`
  extracted from `campaigns/csy4_crrna_free/rfd3__ligandmpnn/3kseq/designs/3ffd6899.tar.zst`
  into a 1-design source `campaigns/csy4_crrna_free/_gputest/src/designs/`.
- **Workers:** the production `code/scripts/stages/refold_af3.sh` and
  `refold_b2.sh`, MSA-off, `num_diffusion_samples=1`. Each prints `nvidia-smi -L`
  so the actual card is logged.
- **Jobs:** 7 GPU types × {AF3, B2} = 14 SLURM jobs, `--time=01:00:00`,
  4 CPU × 8 GB. GPU pinned with `--gpus=<type>:1`; partition per type (below).
- **Partition mapping** (from `sinfo`; a type only schedules in a partition that
  hosts it):
  - `gpu.4h` → rtx_2080_ti, nvidia_titan_rtx, rtx_6000
  - `gpuhe.4h` → rtx_3090, rtx_4090
  - `gpupr.4h` → a100, a100_80gb (private; accessible, no `--account` needed)
- **Locations** (on Euler): jobs under
  `campaigns/csy4_crrna_free/_gputest/refold_{af3,b2}/<key>/` with `run.sbatch`,
  `slurm.log`, `slurm.err`; job→GPU manifest in `_gputest/jobs.tsv`.

## Re-check / re-run

```bash
# states + reasons (run on Euler, or via the SSH wrapper)
sacct -j $(cut -f1 .../_gputest/jobs.tsv | paste -sd, -) -X -o JobName,State,ExitCode,Elapsed
# the failure cause lives in each <dir>/slurm.err
```
SLURM gres/partition strings were pre-validated with `sbatch --test-only`; both
short (`rtx_4090`) and long (`nvidia_geforce_rtx_4090`) `--gpus` names are accepted.

## GPU selection strategy (implemented 2026-06-10)

Adopted from the supervisor's (mlampis) BindCraft scripts: **select by GPU memory
floor**, not by pinning a single type — plus, for AF3, an **Ampere+ node-range**
to dodge the Turing cards that fail it.

His mechanisms:
- sbatch header `--gpus=1 --gres=gpumem:24g` → any ≥24 GB GPU, no `--partition`/type
  (`gpumem` is a `no_consume` constraint gres).
- optional node-range pins via `-w`: `eu-g4-[001-032],eu-g6-[001-080]` (3090+4090)
  or `eu-g6-[001-080]` (4090 only).
- `USE_BULK=1 → -p bulk --requeue`; packs N workers onto one GPU via
  `XLA_PYTHON_CLIENT_MEM_FRACTION=0.95/N`.

Node → GPU map (from `sinfo`):

| node range | GPU | VRAM | AF3 |
|---|---|---|---|
| `eu-lo-g2-*` | rtx_2080_ti | 11 GB | ❌ Turing |
| `eu-g3-*` | titan_rtx / quadro_rtx_6000 | 24 GB | ❌ Turing |
| `eu-g4-*` | rtx_3090 | 24 GB | ✅ Ampere |
| `eu-g6-*` | rtx_4090 | 24 GB | ✅ Ada |
| `eu-a65-*` | a100 / a100_80gb | 40/80 GB | ✅ Ampere |
| `eu-g7-*` | rtx_pro_6000 | 96 GB | ⛔ Blackwell |

**Implemented in the rna_binder refold pipeline** (cluster repo is not git;
backups `*.bak_pre_gpusel`):
- `code/scripts/launch.py` — both refold sbatch builders now emit an optional
  `#SBATCH --nodelist=` when `slurm.nodelist` is set (mirrors the existing `gres`).
- `configs/refold__af3/default.yaml` — `gpus: 1`, `gres: gpumem:24g`,
  `nodelist: eu-g4-[001-032],eu-g6-[001-080]`. AF3 needs sm_80+, so this pins the
  Ampere/Ada 24 GB pool and excludes the Turing `eu-g3` cards. **112 eligible
  nodes** (32 × 3090 + 80 × 4090) vs 80 when pinned to `rtx_4090:1`.
- `configs/refold__b2/default.yaml` — `gpus: 1`, `gres: gpumem:24g`, no nodelist
  (Boltz-2 runs on every generation → widest pool / shortest queue).

Verified with `sbatch --test-only` (nodelist acts as a candidate pool for a 1-GPU
job) and `launch.py --dry-run` (emits the expected `--gpus`/`--gres`/`--nodelist`).
To pin AF3 even tighter (e.g. 4090 only) set `nodelist: eu-g6-[001-080]`.

**a100 / large complexes (>24 GB).** a100 (40 GB) and a100_80gb (80 GB) live in
the `gpupr.*` partition, NOT gpuhe — so a nodelist on a gpuhe job can't reach
them. Select them by switching the partition and setting a matching VRAM floor:
`partition: gpupr.<4h|24h|120h>` + `gpus: 1` + `gres: gpumem:40g` (or `80g`); no
nodelist needed (gpupr holds only a100s). Verified via `--test-only`
(`gpupr.4h`/`24h` honored, not rewritten; `gpumem:40g`→any a100, `80g`→80 GB SKU).
Note: gpupr has no `bulk` tier, so a100 can't run preemptible.

## GUI (Streamlit) — refold GPU dropdown

The launcher's refold GPU selectors are **strategy presets**, not pinned types
(`stacks.REFOLD_GPU_STRATEGIES` → expanded by `app.expand_refold_gpu`):

| strategy | gpus | gres | nodelist / partition |
|---|---|---|---|
| `ampere_24g` (AF3 default) | 1 | gpumem:24g | nodelist eu-g4+eu-g6 |
| `rtx4090_24g` | 1 | gpumem:24g | nodelist eu-g6 |
| `any_24g` (B2 default) | 1 | gpumem:24g | — |
| `a100_40g` | 1 | gpumem:40g | partition gpupr.\<tier\> |
| `a100_80g` | 1 | gpumem:80g | partition gpupr.\<tier\> |

AF3 only offers sm_80+ strategies (no `any_24g`). The a100 tier follows the
walltime dropdown (4h/24h/120h → gpupr.4h/24h/120h). So GUI auto-refolds and
CLI/manual refolds now emit identical GPU selection. Generation-stack GPU
dropdowns still pin a type (unchanged).
