# Generator GPU compatibility

Which generation stacks run on which GPUs. Sibling of `GPU_REFOLD_COMPAT.md`
(which covers AF3/Boltz-2 refolding). Tested 2026-06-11.

## TL;DR

Unlike AF3 (which needed `--flash_attention_implementation=xla` to run on Turing),
**the generators run on Turing sm_75 OUT OF THE BOX — no flag needed.** RFD3,
BoltzGen and BoltzDesign all completed a smoke generation on a Titan RTX with zero
Triton/sm_75/CUDA errors. So they can join the 24 GB titan/quadro pool with only a
GPU-selection change (no code/flag change to the generators themselves).

## Matrix

| Stack | Container | Turing 24 G (titan/quadro, sm_75) | Blackwell sm_120 | Ampere/Ada/a100 |
|---|---|---|---|---|
| **RFD3** (`rfd3`) | `foundry.sif` (torch 2.7.1+cu126) | ✅ runs, no flag | ❌ build lacks sm_120 | ✅ (current pin: rtx_4090) |
| **BoltzGen** | `boltzgen.sif` | ✅ runs, no flag | ❌ | ✅ (current pin: rtx_3090) |
| **BoltzDesign** (BD1) | `boltz_design1.sif` (torch 2.8.0+cu128) | ✅ runs, no flag | ⚠️ build *has* sm_120 but untested | ✅ (current pin: rtx_4090) |
| **HalluDesign** | `af3_protenix_esm.sif` (torch 2.6.0+cu126) | ✅ fits 24 G — prod config ran clean 3h+ no OOM (pilot 3034847); a100_80gb pin was over-spec | ❌ build lacks sm_120 | ✅ (current pin) |
| LigandMPNN (inverse folder) | `LigandMPNN.sif` (torch 2.2.1+cu121) | ✅ | ❌ | ✅ |

`✅` = verified-or-expected runnable; `⚠️` = caveat; `❌` = not without a container rebuild.

## How it was tested

**Phase 0 — static arch probe (one short GPU job).** `torch.cuda.get_arch_list()`
inside each container, run on a real GPU node (get_arch_list returns `[]` on a
login node — it gates on `torch.cuda.is_available()`). Results:

| Container | torch / CUDA | compiled arch_list | Triton |
|---|---|---|---|
| foundry (RFD3) | 2.7.1 / cu126 | sm_50,60,70,**75**,80,86,90 | 3.3.1 |
| boltz_design1 | 2.8.0 / cu128 | sm_70,**75**,80,86,90,100,**120** | 3.4.0 |
| af3_protenix_esm | 2.6.0 / cu126 | sm_50,60,70,**75**,80,86,90 | 3.1.0 |
| LigandMPNN | 2.2.1 / cu121 | sm_50,60,70,**75**,80,86,90 | 2.2.0 |

Every build includes **sm_75** (Turing). Only `boltz_design1` (cu128) includes
**sm_120** (Blackwell) — the rest cap at sm_90, so Blackwell = rebuild-only. All
carry Triton, so the open risk was a Triton flash-attn kernel that's sm_80+ only
(the AF3 failure mode) — Phase 1 settled that.

**Phase 1 — smoke generation on Titan RTX (eu-g3).** Each generator launched at
smoke scale, GPU pinned to Turing (`--gpus=nvidia_titan_rtx:1`, which routes to an
eu-g3 node):

| Stack | Job | Node | State | Elapsed | Notes |
|---|---|---|---|---|---|
| RFD3 (`rfd3__multif`, 2 backbones) | 3025194 | eu-g3-070 | COMPLETED | 4:15 | 2 backbones + 20 MPNN designs, packed |
| BoltzGen (`num_designs=1`) | 3025197 | eu-g3-070 | COMPLETED | 8:10 | design + inverse_folding OK |
| BoltzDesign (`design_samples=2`) | 3025200 | eu-g3-076 | FAILED* | 36:20 | *ran full model; failure = smoke-scale quality filter (iptm 0.29 < 0.5 cutoff), NOT a GPU error |

No Triton / sm_75 / CUDA / OOM errors in any log.

## GPU pools / node ranges (Euler)

- rtx_3090 (Ampere) = `eu-g4-[001-032]`
- rtx_4090 (Ada) = `eu-g6-[001-080]`
- **Turing 24 G** (sm_75): quadro_rtx_6000 = `eu-g3-[062-066],eu-g3-[083-084],eu-lo-g3-[049-061]`; titan_rtx = `eu-g3-[067-079]`
- a100 / a100_80gb = `eu-a65-*` (in the gpupr partition, NOT gpuhe)
- excluded everywhere: 11 GB 2080 Ti (`eu-lo-g2`, < 24 G), Blackwell rtx_pro_6000 (`eu-g7`, sm_120 — breaks every build except boltz_design1)

## Enablement

The GUI generation GPU dropdown gains **"24 GB · all usable (3090/4090/titan/quadro)"**
(key `turing_24g`). `app.pair_gpumem` expands it to `gpus=1 + gpumem:20g +
nodelist=ALL_24G_NODES` (= eu-g4 + eu-g6 + eu-g3 + eu-lo-g3), so SLURM grabs
whichever 24 GB card is free first (most convenient / shortest queue). It is the
**default for RFD3 / BoltzGen / BoltzDesign**. HalluDesign's GUI default still pins
`a100_80gb`, **but it is NOT VRAM-bound** — its prod config ran clean on a 24 GB rtx_4090
for 3h+ with no OOM (pilot 3034847), so the a100 pin was over-spec; HD is a candidate to
flip to `turing_24g` too (use a 24 h walltime — on the 4090 it ran slower than on a100).

Pinned by **nodelist**, not GPU type — on Euler a typed `--gpus=<turing>:1` routes
loosely (a quadro request landed on a titan node). The launch.py generation builders
now emit `--nodelist` (additive; no-op unless set; deployed to the Euler copy).
No generator flag change is needed (the generators run on sm_75 natively).
