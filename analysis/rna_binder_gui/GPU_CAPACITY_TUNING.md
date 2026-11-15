# RNA-binder pipeline — GPU runtime, VRAM & capacity tuning

**Scope:** csy4_crrna_free target (80-aa binder + 16-nt RNA), measured from the
**3kseq** campaign on Euler (2026-06-09/10) plus the 200/600seq runs. All numbers
come from Slurm accounting — peak VRAM and GPU util live in the `.batch` step's
`TRESUsageInMax` (`gres/gpumem`, `gres/gpuutil`); see `euler-sacct-vram-gpuutil`.

> Pull any job's real VRAM/util:
> `sacct -j <id> --format=JobID,State,ElapsedRaw,TRESUsageInMax%120 -p`
> (gpumem/gpuutil appear on the `<id>.batch` line, not the main line.)

---

## 1. The two levers

GPU usage has two independent dimensions, each driven by a different kind of knob:

| Lever | Knob type | Queue-dependent? |
|---|---|---|
| **Fill VRAM** | *depth / parallelism* (diffusion_batch_size, samples, recycles) | No — it's about the card |
| **Fill walltime** | *count* (n_batches, num_designs, design_samples, num_seqs) | Yes — **4 h → 24 h → 120 h = 1× → 6× → 30×** |

Target ~90 % of the walltime to leave room for model-load + the finalize tail.
The *count* knobs only change runtime, **never** VRAM. Confirmed: rfd3 with
`n_batches` 7 / 8 / **60** all peak at **4306 MB**.

---

## 2. Generation — measured (3kseq)

| Stack | GPU | avg runtime | peak VRAM | util | notes |
|---|---|---|---|---|---|
| **rfd3** | RTX 4090 | 12.9 min | **4.3 GB** | 70 % | 64 backbones (n_batches 8 × dbs 8); 3/5 "failed" in the finalize tail (GPU work fine) |
| **boltzgen** | RTX 4090 | 23.5 min | **4.5 GB** | 84 % | one job feeds **both** boltzif (32 designs) + ligandmpnn (320) |
| **boltzdesign** | A100 80 GB | 83.5 min | **6.7 GB** | 87 % | **wastes the A100 (8 %)** — see right-sizing |
| **halludesign** | A100 80 GB | 79.5 min | **79.5 GB** | 68 % | nearly fills the 80 GB card |

## 3. Refolds — measured (3kseq, in-flight)

| Refold | GPU | runtime (completed) | peak VRAM | util |
|---|---|---|---|---|
| **AF3** | RTX 4090 | ~50 min avg (scales with #designs) | **~23 GB** (fills 4090) | 90 % |
| **B2** | RTX 4090 | ~28 min avg | **~2.7 GB** (wastes 4090) | 39 % |

**rfd3 AF3 timeout root cause (not a bug):** rfd3's refolds were submitted with
`--partition=gpuhe.4h` / `--time=04:00:00`; every other stack used `bulk` /
`24:00:00`. With ~600–640 designs/slice either way (HD 600, rfd3 640), 4 h isn't
enough → rfd3 TIMEOUT, HD completes. **Design count is a red herring; the queue
walltime is the cause.** Root: `configs/refold__af3/default.yaml` defaults to
`gpuhe.4h`/4 h (its `.bak_pre_1sample` was `gpuhe.24h`/12 h); the GUI overrides
`slurm.time` to 24 h, so only raw/CLI launches hit this. **Also:** AF3 actually
uses ~23 GB but requests `gpumem:20g` — bump the request to `24g`.

> Euler nuance: `job_submit.lua` derives the `gpuhe.{4h,24h,120h}` tier from
> `--time` and **overrides** any explicit `--partition`. `--time` is the real
> lever; `bulk` is the special preemptible 24 h alias (`-p bulk`, needs
> `--requeue`, which the refold template already has).

---

## 4. What drives VRAM (evidence)

**Moves VRAM (per-example size & depth):**
- **Complex length / token count** — dominant lever; fixed here (80 aa + 16 nt).
- **RFD3 `diffusion_batch_size`** — **measured (probes §8):** 2/8/32/64/96 →
  3.5 / 4.3 / 8.2 / 12.8 / 18.3 GB; **dbs=128 → 23.2 GB = OOM on a 4090.** Fit
  ≈ **+157 MB/unit** on a ~3.1 GB base; **4090 ceiling ≈ dbs 96–104.** Raising
  dbs also lifts GPU util **70 % → 99 %**.
- **HalluDesign `num_recycles` / `num_packs` / `num_seqs`** — HD is the hog (~79.5 GB).
- **Refold `num_diffusion_samples` / `diffusion_samples`, recycles, `use_rna_msa`**
  — pinned to 1 sample / MSA-off here (why msaon ≈ msaoff at ~23 GB: RNA MSA is empty).

**Does NOT move VRAM (throughput / scheduling only):**
- **RFD3 `n_batches`** — 7 / 8 / 60 → identical 4306 MB.
- **LiMPNN `batch_size`/`number_of_batches`, BoltzGen `num_designs`, BoltzDesign `design_samples`**.
- **Number of designs to refold** — processed serially → flat VRAM (this is why rfd3 *timed out* rather than OOM'd).
- **`gres:gpumem` request + GPU-type selection** — a *reservation*, not usage:
  boltzdesign used 6.8 GB whether it requested 20 / 40 / **80 g**.

**Per-stack baselines (len80+16nt, depth pinned):** rfd3 4.3 · boltzgen 4.5 ·
boltzdesign 6.8 · halludesign ~79.5 · AF3 refold ~23 · B2 refold ~2.7 GB.

---

## 5. Measured per-unit costs (for sizing count knobs)

| Stack | per-unit time | unit → designs |
|---|---|---|
| RFD3 | **~17 s / backbone** | 1 backbone → `mpnn batch_size` designs (10) |
| BoltzGen | **~63 s / `num_designs`** (marginal; 44 s single-point) | 1 → `bif + mpnn_bs` designs |
| BoltzDesign | **~370 s / `design_sample`** (marginal; 313 s single-point) | 1 → `mpnn_bs` designs (10) |
| HalluDesign | **~8 s / design** | designs = `num_seqs × num_recycles × num_packs` |
| AF3 refold | **~88 s / design** (measured: 163/640 in 4 h; boltzif 32 in ~2733 s) | serial → **≤144 designs fits 4 h** |
| B2 refold | **~5.5 s / design + ~150 s overhead** | serial |

Anchors: rfd3 64 bb/772 s & 480 bb/8097 s; bg 32/1412 s & (600seq) /3168 s;
bd1 16/5011 s & (600seq) /21288 s; hd 600/4773 s; rfd3 b2 640/3540 s.

---

## 6. Max-out recommendations ("cost no concern")

90 %-walltime targets: 4 h = 12960 s · 24 h = 77760 s · 120 h = 388800 s.

### RFD3 (RTX 4090)
- **Fill VRAM (validated):** `diffusion_batch_size` 8 → **96–104** (dbs=96 measured
  18.3 GB @ 99 % util; **128 OOMs**). This also takes util 70 % → 99 %. Also bump
  LiMPNN `batch_size` 10 → ~50.
- **Fill walltime:** `n_batches = ⌈backbones / diffusion_batch_size⌉` (dbs=104)

| Queue | backbones | designs (×10) | n_batches @ dbs104 |
|---|---|---|---|
| 4 h | ~760 | ~8,300 | ~8 |
| 24 h | ~4,570 | ~45,800 | ~44 |
| 120 h | ~22,870 | ~228,800 | ~220 |

(At dbs≈96 the per-backbone time roughly halved vs dbs=8, so these n_batches are
conservative — re-measure one multi-batch run at dbs=96 to push further.)

### BoltzGen (RTX 4090) — stays VRAM-light
- **VRAM:** ~3–5 GB and **cannot be filled** by exposed knobs — `bg_cli.diffusion_batch_size`
  is a *util* knob (measured 54 % → 81 % at 2→8), not a VRAM filler. Set it to ~8
  for compute; accept the card is memory-underused (or use a smaller GPU).
- **Fill walltime:** `bg_cli.num_designs` ≈ **205 / 1,235 / 6,170** (4 h / 24 h / 120 h).

### BoltzDesign — **move to RTX 4090** (no depth knob, 6.8 GB fixed)
- **Fill walltime:** `design_samples` ≈ **35 / 210 / 1,050** (4 h / 24 h / 120 h).

### HalluDesign (A100 80 GB) — VRAM already maxed
- **Keep A100 80 GB; do NOT raise `num_recycles`** (VRAM driver, would OOM).
- **Fill walltime:** raise `num_seqs` so total designs ≈ **1,600 / 9,700 / 48,600**.

### GPU right-sizing
- **BoltzDesign 4090 not A100** (uses 6.8 GB; compute-bound).
- **B2 refold** wastes the 4090 (2.7 GB) — bump `diffusion_samples` to use it.
- **AF3 refold / HalluDesign** already fill their cards.

---

## 7. AF3 / B2 refold strategy

**Q: just run them on bulk?** Bulk (24 h, preemptible, `-p bulk`) is the safe
zero-code default and is already what bd1/bg/hd refolds use — it never times out.
**But refolds are embarrassingly parallel across designs**, so given cost is no
concern the better play is **shard + fan-out**:

- **AF3** ≈ **88 s/design** (measured: rfd3 reached [163/640] in 4 h before TIMEOUT;
  32-design boltzif slice finished in ~2733 s) → a single 4 h job holds ~**144 designs**.
  Split a 640-design slice into **⌈640/144⌉ = 5 shards** on the **fast 4 h queue**
  (non-preemptible, higher priority) running in parallel → ~5× faster wall-clock
  *and* no timeout *and* no preemption restarts. Beats one ~16 h bulk job.
- **B2** ≈ 5.5 s/design → even 4 h holds ~2,300 designs; **never needs sharding**.
  Just use 4 h. Optionally bump `diffusion_samples` 1 → 5+ to use the idle VRAM.
- **Quality option:** with cost no concern, AF3 `num_diffusion_samples` 1 → 5
  (paper default, best-of-5) uses the otherwise-idle headroom — but 5× the time,
  so shard more.

**Decision rule:**
| Case | Queue |
|---|---|
| B2, any slice | 4 h (single job) |
| AF3, ≤ ~128 designs/slice | 4 h (single job — no shard/bulk needed) |
| AF3, 600–640 designs | **5 shards × 4 h** (preferred) *or* 1× 24 h/bulk (simple) |
| AF3, 120 h single job | ~4,400 designs |

Sharding needs a small `--num-shards N` on the refold launcher (split the per-uuid
design list into N chunks, submit one refold job per chunk). Not yet implemented.

---

## 8. VRAM probe jobs (submitted + completed 2026-06-10) — RESULTS

Tiny jobs (`n_batches=1` / few designs, 1 h, label `vramprobe`).

| jobid | test | result | verdict |
|---|---|---|---|
| 2885269 | rfd3 dbs=32 | 8.2 GB, util 97 % | ✅ |
| 2885271 | rfd3 dbs=64 | 12.8 GB, util 99 % | ✅ |
| 2885273 | rfd3 dbs=96 | **18.3 GB**, util 99 % | ✅ safe ceiling |
| 2885274 | rfd3 dbs=128 | 23.2 GB | ❌ **OOM / FAILED** |
| 2885275/76/77 | boltzgen bg_dbs=2/4/8 | 2.8 / 2.8 / 2.9 GB, util 54/57/81 % | VRAM flat; util knob only |
| 2885279 | boltzdesign on RTX 4090 | 6.8 GB, util 70 %, done in 8 min | ✅ right-sizing confirmed |

**Conclusions:** rfd3 dbs slope ≈ +157 MB/unit, base ~3.1 GB, **4090 ceiling
dbs 96–104** (128 OOMs); dbs also saturates compute (util → 99 %). BoltzGen VRAM
can't be filled via exposed knobs. BoltzDesign runs fine on a 4090.
**TODO: delete the `vramprobe_*` run dirs.**

---

## 9. Caveats
1. Extrapolated VRAM targets (rfd3 dbs, bg dbs) come from ≤2 points — confirm via §8 before a 120 h run.
2. No-depth-knob stacks (BoltzDesign, B2) can't fill a big card — right-size the GPU instead.
3. **Finalize-tail risk:** a 120 h job that dies in the unguarded canonicalize→pack tail loses 5 days. For rfd3/HD, many medium jobs can beat one giant one — the risk is wall-clock, not cost.
4. BoltzGen/BoltzDesign count columns are ceilings (marginal rates a touch higher than single-point); monitor the first job of any new config.

---

*Companion tool:* `gpu_capacity_calc.py` computes these count/depth knobs for any
stack × queue. See `--help`.
