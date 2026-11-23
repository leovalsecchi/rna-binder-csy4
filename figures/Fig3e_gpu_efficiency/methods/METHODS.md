# Methods — GPU efficiency of the five design stacks

Source material for the methods section. Written 2026-11-12, immediately after
the run, so nothing here is reconstructed from memory. Everything it refers to
is copied into `sources/` beside this file; paths in the text are relative to
that folder.

The analysis answers one question: **how much GPU time does each design stack
cost per design, and per design that survives the campaign filter?** It produced
panel Fig3e.

---

## 1. Why a dedicated benchmark rather than the campaign's own accounting

The campaign's SLURM accounting (197,797 GPU-h over `sacct`) cannot rank the
stacks against each other, for four reasons:

1. **A third of it is waste.** Only 68% of campaign GPU-h ended `COMPLETED`; the
   rest was preemption, failure or timeout, distributed unevenly across stacks.
2. **The stacks ran on different GPUs.** RFD3, BoltzGen and BoltzDesign ran
   largely on Turing cards (`quadro_rtx_6000`, `titan_rtx`), HalluDesign on
   Ampere and Ada. GPU generation is confounded with stack identity, which is
   precisely the comparison wanted.
3. **BoltzGen generation is one job serving two stacks** and the ledger records
   it as shared, with no split.
4. **There is no per-stack design denominator** in the accounting.

A dedicated benchmark (`gpubench_v1`) re-ran campaign-identical jobs on a single
GPU model and counted the designs each job actually produced.

## 2. Hardware and scheduling

All benchmark jobs were pinned to NVIDIA RTX 4090 cards by requesting
`--gpus=rtx_4090:1`. On this cluster every RTX 4090 sits in one homogeneous
partition group (80 nodes, `eu-g6`), so the card model is a complete node pin;
no node list was needed. One node on the campaign's permanent exclude list was
excluded here too.

Only non-preemptible partitions were used. A preempted job's elapsed time is not
a slow clock but a **wrong** one, and `sacct` cannot distinguish the two after
the fact.

Every job ran with exactly one GPU, verified from `AllocTRES` at analysis time
(`bench/sacct_all_bench_jobs.psv`).

## 3. Generation

Each stack was run as whole, unmodified production jobs, because per-design
generation cost legitimately includes job-level amortisation (model load,
container start, backbone batch size). Model parameters, seeds policy,
containers, CPU and memory requests are the **production configuration
untouched**; only the SLURM block was overridden, exactly as the campaign's own
submission drivers did. The resolved configuration of one job per stack is in
`configs/*.resolved_config.yaml`.

Replicates and design counts per job:

| stack | jobs | designs per job |
|---|--:|--:|
| RFD3 + LigandMPNN | 3 | 640 |
| BoltzGen (shared) | 3 | 600 BoltzIF + 600 LigandMPNN |
| BoltzDesign + LigandMPNN | 2 (3rd running) | 600 |
| HalluDesign + LigandMPNN | 4 | 150 |

HalluDesign produces 150 designs per job (15 recycles x 10 sequences; its
`num_designs` must remain 1, since de-novo trajectories otherwise collide on a
shared tag), so four jobs were run to reach a comparable total.

## 4. Refolding, and why it is measured on an equal pool

Each design was re-predicted with three backends: AlphaFold3, Boltz-2 and
ESMFold2, which is the campaign's own three-backend consensus procedure.

AlphaFold3 and ESMFold2 load their model **once** and iterate over an input
directory, which is where their order-of-magnitude speed-up over per-design
invocation comes from. Per-design refold cost therefore depends on how many
designs share one model load. Refolding per generation job would have given
HalluDesign a 150-design batch and every other stack ~600, inflating
HalluDesign's apparent refold cost for a reason that has nothing to do with
HalluDesign.

Refolding was therefore measured on a **fixed pool of exactly 600 designs per
stack**, round-robined across that stack's benchmark jobs so no single job
dominates, and submitted as one AlphaFold3 + one Boltz-2 + one ESMFold2 job per
stack — fifteen jobs over five identically-sized batches.

That this was the right control is visible in the result: refold cost came out
**stack-independent at 0.00364–0.00370 GPU-h per design, a 1.6% spread** across
five stacks. Any larger spread would have indicated the batching, not the stack.

## 5. BoltzGen: attributing one allocation to two stacks

BoltzGen's generation job produces backbones once and then inverse-folds them
twice, with BoltzIF and with LigandMPNN, emitting two separate design sets from
a single allocation. Splitting that job 50/50 is wrong in two independent ways:
it credits each stack with a backbone saving that only exists if you want both
products, and it assumes the two inverse folders cost the same.

The stage boundaries are not recoverable from a normal run — the worker prints
its stage banners without timestamps, and both design trees are archived at job
end with identical modification times. One additional generation job was
therefore run with its standard output passed through a timestamper, changing
nothing about the computation (that job's total, 4,658 s, sits inside the range
of the two uninstrumented jobs, 4,672 s and 4,800 s). Measured boundaries are in
`logs/boltzgen_stage_timing_20640689.txt`:

| stage | seconds | share |
|---|--:|--:|
| 1 — backbone design (**shared**) | 765 | 16.7% |
| 2 — BoltzIF inverse folding + folding | 3,336 | 72.6% |
| 3 — LigandMPNN over the 60 backbones | 492 | 10.7% |

BoltzIF costs **6.8x** LigandMPNN, because it runs a structure prediction per
design (`--steps inverse_folding folding`) while LigandMPNN only assigns
sequence to a fixed backbone.

Each BoltzGen stack is therefore charged **the full shared backbone stage plus
its own inverse folding** — what that stack would cost if run alone, which is
what a comparison between stacks is asking. The two shares sum to more than the
job, which is correct: the backbone stage is not divisible.

## 6. Which jobs were counted

GPU-hours are elapsed wall time x allocated GPUs, from `sacct -D -X`. The `-D`
flag matters: without it, requeued attempts are hidden and the surviving row
under-reports.

A job was **rejected rather than averaged in** if it was not `COMPLETED`, had
more than one attempt (i.e. was requeued), was not on exactly one GPU, or
produced no designs. Rejections are listed, never silently dropped
(`logs/phase3b_*.log`).

One further rejection rule was added after a failure mode that `sacct` cannot
see. HalluDesign jobs contend for a single shared JIT-compilation lock for a
DeepSpeed attention kernel; a job that loses that race can sit for hours at 0%
GPU utilisation and still finish, reporting `COMPLETED` on one GPU with a full
design count and a wall clock that is simply wrong. The gap between job start
and framework initialisation separates the two cleanly — 6, 9 and 10 minutes for
healthy jobs against 119 minutes for one that stalled — so any generation job
exceeding 20 minutes to initialisation is rejected. The stalled job was
discarded and re-run.

## 7. Pass rate

Pass rate is the fraction of a stack's **own** generated designs that satisfy the
campaign filter. The filter admitted 55,706 of 6,992,515 designs overall (0.80%).

The filter itself is in `gate/campaign_gate_55706.json`, extracted verbatim from
the header of the campaign's exported passing set (see
`gate/campaign_gate_provenance.txt`). In words:

- ipTM >= 0.6, taking the **maximum** of the AlphaFold3 and ESMFold2 values
- self-consistency RMSD <= 3 A in **either** backend
- self-consistency DockQ >= 0.7 in **either** backend
- longest single-residue run <= 5
- hairpin-foot clash <= 5 binder heavy atoms in **both** backends

⚠️ When parsing that JSON: the ipTM criterion is a **top-level** `iptm` key
(`{"on": true, "value": "max", "mode": "threshold", "threshold": 0.6}`),
combined with the tiers by the top-level `"combine": "AND"`. It is *not* inside
`tiers`, so walking only the tiers finds seven criteria and silently misses the
one that removes the most designs. `iptm_max` is also not a stored column — it
is the element-wise maximum of the two backends' ipTM.

Per-stack counts are in `gate/cascade_by_stack.csv`, variant `old_bif` — the
variant that sums to exactly 55,706. That variant scores BoltzGen + BoltzIF on
its **stored** values, which are what the campaign actually gated and ranked on.
A separate rescore of those designs against their true design-step backbones
exists and gives a slightly different count (9,515 rather than 9,042); it
answers "what are these backbones really like", not "what did the campaign do",
and must not be mixed into campaign-history numbers.

| stack | generated | passing | pass rate |
|---|--:|--:|--:|
| RFD3 + LigandMPNN | 1,524,796 | 27,881 | 1.829% |
| BoltzGen + LigandMPNN | 1,519,958 | 14,752 | 0.971% |
| BoltzGen + BoltzIF | 1,480,140 | 9,042 | 0.611% |
| BoltzDesign + LigandMPNN | 1,444,614 | 3,229 | 0.224% |
| HalluDesign + LigandMPNN | 1,023,007 | 802 | 0.078% |
| **total** | **6,992,515** | **55,706** | **0.797%** |

Note the two measurements come from different places on purpose: GPU-hours from
the benchmark (so the hardware is controlled), pass rates from the full campaign
(so the denominators are large). Pass rate is a property of the designs, not of
the GPU they were made on, so the pairing is sound.

## 8. The derived quantity

    GPU-h per passing design = (GPU-h per design) / (pass rate)

## 9. Results

| stack | gen | refold | total GPU-h/design | designs/GPU-h | pass rate | **GPU-h per passing design** |
|---|--:|--:|--:|--:|--:|--:|
| RFD3 + LigandMPNN | 0.000339 | 0.003657 | 0.003996 | 250 | 1.829% | **0.22** |
| BoltzGen + LigandMPNN | 0.000592 | 0.003696 | 0.004288 | 233 | 0.971% | **0.44** |
| BoltzGen + BoltzIF | 0.001930 | 0.003666 | 0.005596 | 179 | 0.611% | **0.92** |
| BoltzDesign + LigandMPNN | 0.006217 | 0.003638 | 0.009855 | 101 | 0.224% | **4.41** |
| HalluDesign + LigandMPNN | 0.019176 | 0.003647 | 0.022823 | 44 | 0.078% | **29.1** |

Per refold backend, averaged over the five pools: ESMFold2 0.00067, AlphaFold3
0.00147, Boltz-2 0.00151 GPU-h per design. ESMFold2 is ~2.2x cheaper than
either of the others.

Three things worth stating in the text:

- **Refolding dominates for the cheap stacks.** It is 91% of RFD3's total cost
  and 77% of BoltzGen + LigandMPNN's. Reporting generation cost alone would
  overstate the penalty for the expensive stacks by roughly tenfold.
- **The spread compresses, then re-expands.** Generation spans 57x between RFD3
  and HalluDesign; total cost per design spans only 5.7x once the common refold
  cost is added; cost per *passing* design spans **133x** once pass rate is
  applied. Pass rate, not throughput, is what separates the stacks.
- **HalluDesign generation is high-variance.** Its four accepted jobs gave
  0.0166, 0.0177, 0.0206 and 0.0584 GPU-h per design on identical hardware; the
  median is reported. The outlier is CPU contention on a shared node, a real
  cluster condition rather than a measurement fault.

## 10. Deliberate deviations from campaign conditions

All are recorded in `bench/README.md` with the reasoning. Summarised:

| deviation | reason |
|---|---|
| single GPU model (RTX 4090) | removes GPU generation as a confounder; makes benchmark numbers faster than the campaign ledger for the three stacks that ran on Turing |
| non-preemptible partitions only | a preempted elapsed time is a wrong clock |
| requeue disabled | a requeued attempt measures only the remaining work |
| node-local scratch worker variants | what the campaign drivers actually submitted |
| HalluDesign CPU request corrected | the launcher emits a task count without a per-task CPU count for this stack only; the campaign's submission layer fixed it, so this reproduces campaign conditions rather than departing from them |

One thing deliberately **not** corrected: Boltz-2's task-count request, which the
campaign never fixed either. It is identical across stacks and so cannot affect
the ranking.

## 11. Caveats

- HalluDesign benchmark jobs each paid a one-off JIT kernel rebuild (6–10 min on
  a ~2.7 h job, so ~4% overhead) that the campaign's pre-seeded jobs did not.
  Far too small to affect a 50x ranking, but it is in the numbers.
- BoltzDesign has two accepted replicates at the time of writing; a third was
  running and will move the median by well under one percent.
- Pass rates are from the campaign, which ran across several hardware
  generations and software revisions. They are design properties, but they are
  not from the same jobs the GPU-hours came from.

---

## File manifest

| path | what it is |
|---|---|
| `bench/README.md` | the benchmark's own protocol document, including the three measurement traps found while building it |
| `bench/bench_gen.sh` | generation submitter: header edits, the HalluDesign JIT fix, the submit gate |
| `bench/bench_refold.py` | builds the equal-N pools and submits all three backends |
| `bench/bench_lib.py` | reads the archived design sets; the finished-run gate |
| `bench/bif_true_designs.py`, `bench/test_transplant.py` | the BoltzIF true-backbone correction and its self-check |
| `bench/probe_biopython.py` | parses designs with the refold worker's own parser; the pool gate |
| `bench/collect_gpuh.py` | the accounting: rejections, per-stack medians, the BoltzGen stage split |
| `bench/gen_manifest.tsv`, `bench/refold_manifest.tsv` | every job this analysis submitted, with its run directory |
| `bench/sacct_all_bench_jobs.psv` | raw `sacct -D -X` for all of them |
| `bench/gpuh_per_design.tsv` | the accepted per-job numbers behind the results table |
| `configs/*.resolved_config.yaml` | the exact production configuration each stack resolved, containers included |
| `configs/refold_*.default.yaml` | the three refold backend configurations |
| `gate/campaign_gate_55706.json` | THE campaign filter, as machine-readable JSON |
| `gate/campaign_gate_provenance.txt` | where that JSON was extracted from |
| `gate/cascade_by_stack.csv` | per-stack funnel; variant `old_bif` is the one summing to 55,706 |
| `logs/boltzgen_stage_timing_20640689.txt` | the timestamped stage banners behind the BoltzGen split |
| `logs/phase2_*.log`, `logs/phase3*.log` | the run's own output, including every rejection |
