#!/usr/bin/env python3
"""Fill midrun_v5 GEN gaps: resubmit uuids that have a run.sbatch, NO packed designs
slice, and are NOT covered by a live gen array (so no double-runs). Bucketed into
--cpus-per-task arrays per (stack x lane), respecting arch floor + walltime, %40 throttle.
Write-only unless --submit."""
import re, os, sys, glob, subprocess, collections

RB = "/cluster/work/wicky/lvalsecchi/rna_binder"
C = RB + "/campaigns/csy4_crrna_free"; CODE = RB + "/code"
OUT = C + "/_mv5_gapfill_arrays"; LOGDIR = OUT + "/logs"; os.makedirs(LOGDIR, exist_ok=True)
EXCLUDE = "eu-lo-g2-028,eu-lo-g3-059,eu-g6-064"; LABEL = "midrun_v5"
SUBMIT = "--submit" in sys.argv
only = None
for a in sys.argv:
    if a.startswith("--stacks="):
        only = set(a.split("=", 1)[1].split(","))

# stack -> (stack_key, layout)  layout: "multif" (shared designs/) or "hd" (per-uuid)
STACKS = collections.OrderedDict([
    ("bg",  ("boltzgen__boltzif", "multif")),
    ("bd",  ("boltzdesign__ligandmpnn", "multif")),
    ("hd",  ("halludesign__ligandmpnn", "hd")),
])
AMPERE = "eu-g4-[001-031],eu-g6-[001-080]"
G3A = "eu-g3-[062-070]"; G3B = "eu-g3-[071-077,083-084]"
LANES = {
    "bg": [("gpu.24h","es_wicky/gpu/24","23:30:00","gpu:1",G3A),
           ("gpu.120h","es_wicky/gpu/120","48:00:00","gpu:1",G3B),
           ("gpuhe.120h","es_wicky/gpuhe/120","48:00:00","gpumem:24g",AMPERE)],
    "bd": [("gpu.24h","es_wicky/gpu/24","23:30:00","gpu:1",G3A),
           ("gpu.120h","es_wicky/gpu/120","48:00:00","gpu:1",G3B),
           ("gpuhe.24h","es_wicky/gpuhe/24","23:30:00","gpumem:24g",AMPERE)],
    "hd": [("gpuhe.24h","es_wicky/gpuhe/24","23:30:00","gpumem:24g",AMPERE),
           ("gpuhe.120h","es_wicky/gpuhe/120","48:00:00","gpumem:24g",AMPERE),
           ("gpupr.24h","es_wicky/gpupr/24","23:30:00","gpumem:40g","eu-a65-[02-07],eu-ts-[01-02]")],
}
TEST_UUIDS = {"4b1043ef", "dc7fe23b", "fcdc1a0d", "46db14f3"}


def sh(*a):
    return subprocess.run(a, capture_output=True, text=True).stdout


def inflight_runsbatch():
    """run.sbatch paths covered by a live gen job (array SCRIPTS or a direct run.sbatch)."""
    live = set()
    for o in sh("squeue", "-u", "lvalsecchi", "-h", "-o", "%o").splitlines():
        o = o.strip()
        if not o.endswith(".sbatch"):
            continue
        if o.endswith("/run.sbatch"):
            live.add(o); continue
        try:
            txt = open(o).read()
            m = re.search(r"SCRIPTS=\((.*?)\n\)", txt, re.DOTALL)
            if m:
                for ln in m.group(1).splitlines():
                    ln = ln.strip()
                    if ln.endswith("run.sbatch"):
                        live.add(ln)
        except OSError:
            pass
    return live


def hdr(p, k):
    pat = re.compile(r"^#SBATCH --%s=(.*)$" % re.escape(k))
    for ln in open(p):
        m = pat.match(ln.strip())
        if m:
            return m.group(1)
    return None


def uuid_of(p):
    m = re.search(r"(?:run_|%s_)([0-9a-f]{8})" % LABEL, p)
    return m.group(1) if m else None


LIVE = inflight_runsbatch()
print("live gen run.sbatch (in-flight): %d" % len(LIVE))
made = []
for stack, (sk, layout) in STACKS.items():
    if only and stack not in only:
        continue
    if layout == "multif":
        runs = glob.glob(f"{C}/{sk}/{LABEL}/run_*/run.sbatch")
        packed = lambda u: os.path.exists(f"{C}/{sk}/{LABEL}/designs/{u}.tar.zst")
    else:  # hd per-uuid
        runs = glob.glob(f"{C}/{sk}/{LABEL}_*/run.sbatch")
        packed = lambda u: os.path.exists(f"{os.path.dirname(_)}/designs.tar.zst")
    gaps = []
    for _ in sorted(runs):
        u = uuid_of(_)
        if not u or u in TEST_UUIDS:
            continue
        if _ in LIVE:
            continue                     # in-flight -> skip
        if layout == "multif":
            if os.path.exists(f"{C}/{sk}/{LABEL}/designs/{u}.tar.zst"):
                continue                 # done
        else:
            if os.path.exists(f"{os.path.dirname(_)}/designs.tar.zst"):
                continue                 # done
        gaps.append(_)
    lanes = LANES[stack]
    buckets = [[] for _ in lanes]
    for i, p in enumerate(gaps):
        buckets[i % len(lanes)].append(p)
    print(f"\n{stack}: gap={len(gaps)} -> {len([b for b in buckets if b])} lanes")
    mem = hdr(gaps[0], "mem-per-cpu") if gaps else "8G"
    for (part, qos, wall, gres, nodelist), scr in zip(lanes, buckets):
        name = part.replace('.','_')
        if not scr:
            continue
        n = len(scr)
        print(f"    {part:11s} n={n:<4d} qos={qos}")
        sb = os.path.join(OUT, f"arr_{stack}_{name}.sbatch")
        L = ["#!/bin/bash", f"#SBATCH --job-name=mv5_{stack}", f"#SBATCH --partition={part}",
             f"#SBATCH --qos={qos}", f"#SBATCH --time={wall}", "#SBATCH --cpus-per-task=4",
             f"#SBATCH --mem-per-cpu={mem or '8G'}", "#SBATCH --gpus=1",
             f"#SBATCH --gres={gres}", f"#SBATCH --nodelist={nodelist}",
             f"#SBATCH --exclude={EXCLUDE}", "#SBATCH --requeue", "#SBATCH --open-mode=append",
             f"#SBATCH --array=0-{n-1}%40", f"#SBATCH --chdir={CODE}",
             f"#SBATCH --output={LOGDIR}/{stack}_{name}_%A_%a.log",
             f"#SBATCH --error={LOGDIR}/{stack}_{name}_%A_%a.err", "", "set -euo pipefail", "SCRIPTS=("]
        L += [f"  {s}" for s in scr]
        L += [")", 'S="${SCRIPTS[$SLURM_ARRAY_TASK_ID]}"', 'RUNDIR=$(dirname "$S")',
              f'echo "[gap] {stack} $SLURM_ARRAY_TASK_ID $(hostname) -> $S"',
              'env -u SLURM_NTASKS -u SLURM_NPROCS -u SLURM_NTASKS_PER_NODE -u SLURM_TASKS_PER_NODE \\',
              '    bash "$S" > "$RUNDIR/slurm.log" 2> "$RUNDIR/slurm.err"']
        open(sb, "w").write("\n".join(L) + "\n")
        made.append(sb)

print(f"\n=== {len(made)} gap-fill arrays ===")
if SUBMIT:
    for sb in made:
        r = subprocess.run(["sbatch", "--parsable", sb], capture_output=True, text=True)
        print("submitted %-24s -> %s%s" % (os.path.basename(sb), r.stdout.strip(), r.stderr.strip()))
else:
    print("(dry-run; add --submit)")
