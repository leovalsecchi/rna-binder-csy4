#!/usr/bin/env python3
"""Redistribute midrun_v4 PENDING gen uuids across more partitions for throughput.

Each partition tier is a separate QOS budget; spreading pending work across gpu.24h/.120h
+ gpuhe.24h/.120h/.4h + gpupr.24h lets far more run concurrently. Respects per-stack
"right conditions":
  * BoltzGen (Turing-OK, ~3.4h): gpu.24h/.120h + gpuhe.120h + gpuhe.4h (eu-g6 4090 only, fits 4h)
  * BoltzDesign (Turing-OK, ~7h -> >=24h): gpu.24h/.120h + gpuhe.24h/.120h (NOT 4h)
  * HalluDesign (Ampere/Ada ONLY, ~3-7h -> >=24h): gpuhe.24h/.120h + gpupr.24h (NOT Turing, NOT 4h)

PENDING = run dir has no slurm.log (never dispatched). Running/done uuids (have slurm.log)
are left untouched. Emits one --cpus-per-task array per (stack x lane), env -u guard.
Write-only unless --submit.
"""
import re, os, sys, glob, subprocess, collections

RB = "/cluster/work/wicky/lvalsecchi/rna_binder"
C = RB + "/campaigns/csy4_crrna_free"
CODE = RB + "/code"
OUT = C + "/_mv4_redist_arrays"; LOGDIR = OUT + "/logs"; os.makedirs(LOGDIR, exist_ok=True)
EXCLUDE = "eu-g6-064"; LABEL = "midrun_v4"
SUBMIT = "--submit" in sys.argv
TEST_UUIDS = {"4b1043ef", "dc7fe23b", "fcdc1a0d", "46db14f3"}

STACKS = {"bg": "boltzgen__boltzif", "bd": "boltzdesign__ligandmpnn", "hd": "halludesign__ligandmpnn"}
AMPERE = "eu-g4-[001-031],eu-g6-[001-080]"   # 3090 + 4090 (HD-eligible)
G3A = "eu-g3-[062-070]"; G3B = "eu-g3-[071-077,083-084]"   # public 24G Turing bins
# (name, partition, qos, walltime, gres, nodelist)
LANES = {
    "bg": [
        ("gpu24",   "gpu.24h",   "es_wicky/gpu/24",   "23:30:00", "gpu:1",      G3A),
        ("gpu120",  "gpu.120h",  "es_wicky/gpu/120",  "48:00:00", "gpu:1",      G3B),
        ("ghe120",  "gpuhe.120h","es_wicky/gpuhe/120","48:00:00", "gpumem:24g", AMPERE),
        ("ghe4",    "gpuhe.4h",  "es_wicky/gpuhe/4",  "03:50:00", "gpumem:24g", "eu-g6-[001-080]"),
    ],
    "bd": [
        ("gpu24",   "gpu.24h",   "es_wicky/gpu/24",   "23:30:00", "gpu:1",      G3A),
        ("gpu120",  "gpu.120h",  "es_wicky/gpu/120",  "48:00:00", "gpu:1",      G3B),
        ("ghe24",   "gpuhe.24h", "es_wicky/gpuhe/24", "23:30:00", "gpumem:24g", AMPERE),
        ("ghe120",  "gpuhe.120h","es_wicky/gpuhe/120","48:00:00", "gpumem:24g", AMPERE),
    ],
    "hd": [
        ("ghe24",   "gpuhe.24h", "es_wicky/gpuhe/24", "23:30:00", "gpumem:24g", AMPERE),
        ("ghe120",  "gpuhe.120h","es_wicky/gpuhe/120","48:00:00", "gpumem:24g", AMPERE),
        ("gpupr24", "gpupr.24h", "es_wicky/gpupr/24", "23:30:00", "gpumem:40g", "eu-a65-[02-07],eu-ts-[01-02]"),
    ],
}


def hdr(path, key):
    pat = re.compile(r"^#SBATCH --%s=(.*)$" % re.escape(key))
    with open(path) as fh:
        for line in fh:
            m = pat.match(line.strip())
            if m:
                return m.group(1)
    return None


def uuid_of(p):
    m = re.search(r"(?:run_|%s_)([0-9a-f]{8})" % LABEL, p)
    return m.group(1) if m else None


made = []
print("=== redistribution plan (PENDING = no slurm.log) ===")
for stack, sk in STACKS.items():
    paths = sorted(set(glob.glob(f"{C}/{sk}/{LABEL}/run_*/run.sbatch") +
                       glob.glob(f"{C}/{sk}/{LABEL}_*/run.sbatch")))
    pending = [p for p in paths
               if not os.path.exists(os.path.join(os.path.dirname(p), "slurm.log"))
               and uuid_of(p) not in TEST_UUIDS]
    lanes = LANES[stack]
    buckets = [[] for _ in lanes]
    for i, p in enumerate(pending):
        buckets[i % len(lanes)].append(p)
    print(f"\n{stack}: {len(pending)} pending -> {len(lanes)} lanes")
    mem = hdr(pending[0], "mem-per-cpu") if pending else "8G"
    for (name, part, qos, wall, gres, nodelist), scripts in zip(lanes, buckets):
        if not scripts:
            continue
        n = len(scripts)
        print(f"    {part:11s} qos={qos:18s} n={n:<3d} gres={gres}")
        sb = os.path.join(OUT, f"arr_{stack}_{name}.sbatch")
        L = [
            "#!/bin/bash",
            f"#SBATCH --job-name=mv4_{stack}",
            f"#SBATCH --partition={part}",
            f"#SBATCH --qos={qos}",
            f"#SBATCH --time={wall}",
            "#SBATCH --cpus-per-task=4",
            f"#SBATCH --mem-per-cpu={mem or '8G'}",
            "#SBATCH --gpus=1",
            f"#SBATCH --gres={gres}",
            f"#SBATCH --nodelist={nodelist}",
        ]
        if part.startswith("gpuhe") or part.startswith("gpupr"):
            L.append(f"#SBATCH --exclude={EXCLUDE}")
        L += [
            "#SBATCH --requeue", "#SBATCH --open-mode=append",
            f"#SBATCH --array=0-{n-1}",
            f"#SBATCH --chdir={CODE}",
            f"#SBATCH --output={LOGDIR}/{stack}_{name}_%A_%a.log",
            f"#SBATCH --error={LOGDIR}/{stack}_{name}_%A_%a.err",
            "", "set -euo pipefail", "SCRIPTS=(",
        ]
        L += [f"  {s}" for s in scripts]
        L += [
            ")",
            'S="${SCRIPTS[$SLURM_ARRAY_TASK_ID]}"',
            'RUNDIR=$(dirname "$S")',
            f'echo "[mv4-redist] stack={stack} lane={name} task=$SLURM_ARRAY_TASK_ID host=$(hostname) -> $S"',
            'env -u SLURM_NTASKS -u SLURM_NPROCS -u SLURM_NTASKS_PER_NODE -u SLURM_TASKS_PER_NODE \\',
            '    bash "$S" > "$RUNDIR/slurm.log" 2> "$RUNDIR/slurm.err"',
        ]
        with open(sb, "w") as fh:
            fh.write("\n".join(L) + "\n")
        made.append(sb)

print(f"\n=== {len(made)} arrays built ===")
if SUBMIT:
    for sb in made:
        r = subprocess.run(["sbatch", "--parsable", sb], capture_output=True, text=True)
        print("submitted %-24s -> %s%s" % (os.path.basename(sb), r.stdout.strip(), r.stderr.strip()))
else:
    print("(dry-run; add --submit)")
