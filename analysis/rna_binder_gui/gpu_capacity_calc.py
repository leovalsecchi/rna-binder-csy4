#!/usr/bin/env python3
"""GPU capacity calculator for the RNA-binder pipeline.

Given a stack and a queue tier (4h / 24h / 120h), compute:
  - the *count* knob value that fills the walltime (≈90% by default), and
  - the *depth* knob value that fills the GPU's VRAM (where one exists),
plus expected designs, runtime and peak VRAM.

All rates/VRAM models are measured from the 3kseq + 200/600seq runs on Euler
(see GPU_CAPACITY_TUNING.md). They are anchored to the csy4_crrna_free target
(80-aa binder + 16-nt RNA); change --sec-per-unit if your target differs.

2026-06-11 update (6kseq_4h sacct, 65 boltzgen / 131 rfd3 / 55 bd / 51 hd / 582
refold jobs): added measured peak host-RAM per stack (the `QOSMaxMemoryPerUser`
limiter) + right-sized `mem_per_cpu` requests; VRAM models confirmed. boltzgen
runtime updated for **Route B** (BIF switched to `--steps inverse_folding` + a patched
DesignWriter so it keeps 10/backbone like LiMPNN — `--only_inverse_fold` collapsed to
1/bb — folding now does 10x the structures, ~+100 min/UUID; VRAM unchanged).

Examples:
  python3 gpu_capacity_calc.py --all                 # full matrix, all stacks/queues
  python3 gpu_capacity_calc.py --stack rfd3 --queue 24h
  python3 gpu_capacity_calc.py --refold af3 --designs 640   # queue/shard advice
  python3 gpu_capacity_calc.py --mem                 # measured RAM + right-sized --mem
"""
import argparse

# 90%-of-walltime is the default fill target (model-load + finalize-tail margin).
QUEUES = {"4h": 14400, "24h": 86400, "120h": 432000}

# Usable VRAM target per GPU (MB) — capacity minus ~2 GB headroom.
GPU_USABLE = {
    "rtx_4090": 22000, "ampere_24g": 22000, "any_24g": 22000,
    "a100": 38000, "a100_80gb": 78000,
}
GPU_CAP = {"rtx_4090": 24564, "ampere_24g": 24564, "any_24g": 24564,
           "a100": 40960, "a100_80gb": 81920}

# --- stack models -----------------------------------------------------------
# sec_per_unit : measured seconds per 'count unit'
# count_path   : GUI/launch override path for the walltime knob
# unit         : name of the count unit
# designs_per_unit : designs produced per count unit (default mpnn batch_size=10)
# depth        : (path, base_MB, per_unit_MB) for the VRAM knob, or None
# vram_fixed   : measured peak VRAM (MB) when there is no depth knob
# gpu          : recommended GPU
# ram_mb       : measured peak host RAM (MB) — the QOSMaxMemoryPerUser limiter.
#                Right-size mem_per_cpu to ~2x this (see --mem); 2026-06-11 sacct.
STACKS = {
    "rfd3": dict(
        # VRAM model fit to MEASURED probes (dbs 32/64/96 = 8.2/12.8/18.3 GB);
        # dbs=128 OOMs a 4090, so cap at 104. Raising dbs also lifts util 70%->99%.
        sec_per_unit=17.0, count_path="generator_params.n_batches", unit="backbone",
        designs_per_unit=10, gpu="rtx_4090", vram_fixed=4306, ram_mb=10100,
        depth=("generator_params.diffusion_batch_size", 3138, 157.5), depth_default=8,
        depth_cap=104,
        note="EXTRAPOLATED: production only ever ran dbs=8 (~5 GB VRAM); VRAM probes maxed at "
             "dbs=96 (18.3 GB), so dbs=104 is unvalidated and may OOM a 4090 — verify before "
             "trusting the high-count rows. RAM 10.1 GB peak -> config cut 32->20 GB 2026-06-11"),
    "boltzgen": dict(
        # Route B VALIDATED 2026-06-11 (v4 job 3053806: 2 bb -> 20 BIF + 20 LiMPNN = 10/bb,
        # matched, _design<bb>_seq<seq>). BIF keeps 10 seqs/bb via --steps inverse_folding +
        # patched DesignWriter; BOTH folding AND design_folding then process the 10x-bigger
        # set (~6-8 s/struct, high node variance) -> the two fold passes dominate runtime.
        # sec_per_unit 63 -> 190 (num_designs=64 ~= 3.4 h, was ~40 min).
        sec_per_unit=190.0, count_path="generator_params.bg_cli.num_designs",
        unit="num_designs", designs_per_unit=10, gpu="rtx_4090", vram_fixed=4510, ram_mb=39000,
        depth=None,
        note="USE 24h NOT 4h under Route B: num_designs=64 ~= 3.4 h (4h holds only ~61). VRAM "
             "unchanged ~4.5 GB (folding sequential). RAM 39 GB -> KEPT 64 GB. Now 10/bb on "
             "BOTH boltzif + ligandmpnn (matched)."),
    "boltzdesign": dict(
        # Runtime RE-FIT 2026-06-11: production design_samples=200 -> ~6 h (a bd1 run hit
        # 5h54m, cancelled on the 4h cap; 72/104 bd1 jobs cancelled). ~108 s/sample in the
        # production regime (+~250 s fixed overhead seen on small 2-5 sample smoke runs).
        sec_per_unit=108.0, count_path="generator_params.design_samples",
        unit="design_sample", designs_per_unit=10, gpu="rtx_4090", vram_fixed=6848, ram_mb=9500,
        depth=None, note="USE 24h NOT 4h — production design_samples=200 needs ~6 h (measured). "
                         "VRAM 6.8 GB, RAM 9.5 GB -> config cut 40->20 GB 2026-06-11"),
    "halludesign": dict(
        sec_per_unit=8.0, count_path="generator_params.params.num_seqs",
        unit="design", designs_per_unit=1, gpu="a100_80gb", vram_fixed=79500, ram_mb=27500,
        depth=None,
        note="VRAM maxed (a100/80GB); RAM 27.5 GB -> leave 32 GB. raise num_seqs (NOT "
             "num_recycles); designs=num_seqs*num_recycles*num_packs"),
}

# refolds: serial over designs -> overhead + per-design.  ram_mb = measured peak host RAM.
REFOLDS = {
    "af3": dict(overhead=300, sec_per_design=88.0, vram=23000, gpu="ampere_24g", ram_mb=9300),
    "b2": dict(overhead=150, sec_per_design=5.5, vram=2700, gpu="any_24g", ram_mb=10400),
    # ESMFold2: ~200s one-time startup (25GB node-local weight stage + ESMC-6B
    # load + Triton compile), ~4s/design steady-state (measured 2026-06-17, 3090).
    "esm": dict(overhead=200, sec_per_design=4.0, vram=22000, gpu="ampere_24g", ram_mb=12000),
}


def fmt_hms(s):
    s = int(s); return f"{s//3600}h{(s%3600)//60:02d}m"


def depth_fill(stack, gpu=None):
    s = STACKS[stack]; gpu = gpu or s["gpu"]
    if not s.get("depth"):
        probe = s.get("depth_probe")
        return None, (f"no depth knob — fixed ~{s['vram_fixed']/1024:.1f} GB"
                      + (f"; probe `{probe}`" if probe else ""))
    path, base, per = s["depth"]
    usable = GPU_USABLE.get(gpu, 22000)
    val = int((usable - base) / per)
    cap = s.get("depth_cap")
    capped = cap is not None and val > cap
    if capped:
        val = cap
    pred = base + per * val
    tag = f"capped at {cap} (128 OOMs)" if capped else "measured fit"
    return val, f"{path} ≈ {val}  (predicts {pred/1024:.1f} GB on {gpu}; {tag})"


def report_stack(stack, queue, margin):
    s = STACKS[stack]; wall = QUEUES[queue] * margin
    count = int(wall / s["sec_per_unit"])
    designs = count * s["designs_per_unit"]
    dval, dmsg = depth_fill(stack)
    if stack == "rfd3":
        dbs = dval or s["depth_default"]
        nb = -(-count // dbs)  # ceil
        countstr = f"n_batches={nb} (×dbs {dbs} = {nb*dbs} backbones)"
        designs = nb * dbs * s["designs_per_unit"]
    elif stack == "halludesign":
        nr, npk = 20, 1
        ns = -(-designs // (nr * npk))
        countstr = f"num_seqs={ns} (×{nr} recycles ×{npk} packs)"
        designs = ns * nr * npk
    else:
        countstr = f"{s['count_path'].split('.')[-1]}={count}"
    return dict(stack=stack, queue=queue, count=countstr, designs=designs,
                runtime=fmt_hms(designs / s["designs_per_unit"] * s["sec_per_unit"]),
                vram=dmsg, gpu=s["gpu"])


def print_stack_rows(rows):
    print(f"{'STACK':12s} {'QUEUE':5s} {'COUNT KNOB (fill walltime)':40s} {'~designs':>9s}  VRAM-fill knob")
    print("-" * 110)
    for r in rows:
        print(f"{r['stack']:12s} {r['queue']:5s} {r['count']:40s} {r['designs']:>9,d}  {r['vram']}")
    print(f"\nGPU per stack: " + ", ".join(f"{k}={v['gpu']}" for k, v in STACKS.items()))
    shown = list(dict.fromkeys(r["stack"] for r in rows))
    notes = [(s, STACKS[s].get("note")) for s in shown if STACKS[s].get("note")]
    if notes:
        print("\nNOTES — read before sizing a real run:")
        for s, n in notes:
            print(f"  [{s}] {n}")


def report_refold(kind, designs, margin):
    r = REFOLDS[kind]
    print(f"\n=== {kind.upper()} refold — {designs} designs "
          f"(~{r['sec_per_design']} s/design + {r['overhead']} s overhead, "
          f"~{r['vram']/1024:.1f} GB on {r['gpu']}) ===")
    runtime = r["overhead"] + designs * r["sec_per_design"]
    print(f"single job runtime ≈ {fmt_hms(runtime)}")
    for q in ("4h", "24h", "120h"):
        cap = int((QUEUES[q] * margin - r["overhead"]) / r["sec_per_design"])
        fits = "OK (single job)" if runtime <= QUEUES[q] * margin else \
               f"NO — shard into {-(-designs // max(1,cap))} jobs of ≤{cap}"
        print(f"  {q:5s}: holds ~{cap:>6,d} designs/job  ->  {fits}")
    if kind == "af3":
        print("  recommendation: shard to 4h (non-preemptible, parallel) > one 24h/bulk job.")
    else:
        print("  recommendation: 4h single job is plenty; bump diffusion_samples to use idle VRAM.")


def report_mem():
    print(f"{'STACK':14s}{'peak RAM used':>15s}{'~2x request':>14s}   GPU")
    print("-" * 60)
    rows = [(k, v.get("ram_mb"), v.get("gpu", "")) for k, v in STACKS.items()]
    rows += [(f"refold_{k}", v.get("ram_mb"), v.get("gpu", "")) for k, v in REFOLDS.items()]
    for name, ram, gpu in rows:
        if not ram:
            continue
        rec = max(16, round(ram * 2 / 1024 / 4) * 4)  # ~2x peak, nearest 4 GB, floor 16
        print(f"{name:14s}{ram/1024:>13.1f} G{rec:>12d} G   {gpu}")
    print("\nRAM is the QOSMaxMemoryPerUser limiter — request ~2x measured peak, not the old\n"
          "32-64 GB defaults (halving them ~doubles concurrent jobs under the cap). Applied\n"
          "2026-06-11: rfd3 & boltzdesign 32/40 -> 20 GB; boltzgen KEPT 64 GB (39 GB peak +\n"
          "Route B headroom); halludesign stays 32 GB (27.5 GB peak, works).")


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--stack", choices=list(STACKS))
    p.add_argument("--queue", choices=list(QUEUES))
    p.add_argument("--all", action="store_true", help="full matrix: all stacks × all queues")
    p.add_argument("--refold", choices=list(REFOLDS), help="refold queue/shard advice")
    p.add_argument("--designs", type=int, default=640, help="design count for --refold")
    p.add_argument("--margin", type=float, default=0.9, help="walltime fill fraction (default 0.9)")
    p.add_argument("--mem", action="store_true",
                   help="measured peak host-RAM per stack + right-sized --mem (QOS limiter)")
    a = p.parse_args()

    if a.mem:
        report_mem(); return
    if a.refold:
        report_refold(a.refold, a.designs, a.margin); return
    if a.all or (not a.stack):
        rows = [report_stack(s, q, a.margin) for s in STACKS for q in QUEUES]
        print_stack_rows(rows)
    else:
        queues = [a.queue] if a.queue else list(QUEUES)
        print_stack_rows([report_stack(a.stack, q, a.margin) for q in queues])
    if not a.stack and not a.refold:
        print("\n(refold advice: --refold af3 --designs N)")


if __name__ == "__main__":
    main()
