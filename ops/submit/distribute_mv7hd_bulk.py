#!/usr/bin/env python3
"""Redistribute the mv7_hd HD tail off the starved gpuhe.bulk: move 80 -> gpuhe.24h and
80 -> gpuhe.120h (higher-priority tiers that schedule), rest stay bulk. Rebuild all as clean
compact --cpus-per-task arrays (bucketer syntax). Writes 3 arr files (no cancel/submit)."""
import os

AR = "/cluster/work/wicky/lvalsecchi/rna_binder/campaigns/csy4_crrna_free/_mv7_arrays"
RB = "/cluster/work/wicky/lvalsecchi/rna_binder"
LOGDIR = AR + "/logs"; os.makedirs(LOGDIR, exist_ok=True)

runs = []
for src in ["arr_mv7hd_final.sbatch", "arr_mv7hd_final_gpuhe_24h.sbatch"]:
    p = os.path.join(AR, src)
    if not os.path.exists(p):
        continue
    for l in open(p):
        s = l.strip()
        if s.startswith("/cluster") and s.endswith("run.sbatch"):
            runs.append(s)
runs = sorted(set(runs))
N = len(runs)
print("pooled mv7_hd run.sbatch:", N)

MOVE = 80
n120 = MOVE
n24 = MOVE + 19            # 80 moved from bulk + the 19 already on 24h
plan = [
    ("gpuhe.120h", "es_wicky/gpuhe/120", "47:30:00", runs[0:n120]),
    ("gpuhe.24h",  "es_wicky/gpuhe/24",  "23:30:00", runs[n120:n120 + n24]),
    ("gpuhe.bulk", "es_wicky/gpuhe/bulk", "23:30:00", runs[n120 + n24:]),
]

def write(part, qos, wall, rlist):
    sig = part.replace(".", "_")
    sb = "%s/arr_mv7hd_dist_%s.sbatch" % (AR, sig)
    n = len(rlist)
    L = ["#!/bin/bash", "#SBATCH --job-name=mv7_hd", "#SBATCH --partition=%s" % part,
         "#SBATCH --qos=%s" % qos, "#SBATCH --time=%s" % wall, "#SBATCH --cpus-per-task=4",
         "#SBATCH --mem-per-cpu=7G", "#SBATCH --gpus=1", "#SBATCH --gres=gpumem:24g",
         "#SBATCH --nodelist=eu-g4-[001-031],eu-g6-[001-080]", "#SBATCH --exclude=eu-g6-064",
         "#SBATCH --requeue", "#SBATCH --open-mode=append", "#SBATCH --array=0-%d" % (n - 1),
         "#SBATCH --chdir=%s/code" % RB,
         "#SBATCH --output=%s/mv7hd_dist_%s_%%A_%%a.log" % (LOGDIR, sig),
         "#SBATCH --error=%s/mv7hd_dist_%s_%%A_%%a.err" % (LOGDIR, sig),
         "", "set -euo pipefail", "SCRIPTS=("]
    L += ["  " + r for r in rlist]
    L += [")", 'S="${SCRIPTS[$SLURM_ARRAY_TASK_ID]}"', 'RUNDIR=$(dirname "$S")',
          'echo "[mv7hd_dist] task=$SLURM_ARRAY_TASK_ID host=$(hostname) -> $S"',
          'env -u SLURM_NTASKS -u SLURM_NPROCS -u SLURM_NTASKS_PER_NODE -u SLURM_TASKS_PER_NODE \\',
          '    bash "$S" > "$RUNDIR/slurm.log" 2> "$RUNDIR/slurm.err"']
    open(sb, "w").write("\n".join(L) + "\n")
    print("wrote %-32s tasks=%-4d part=%s" % (os.path.basename(sb), n, part))

for part, qos, wall, rlist in plan:
    write(part, qos, wall, rlist)
print("distribution: 120h=%d  24h=%d  bulk=%d  (total %d)" % (n120, n24, N - n120 - n24, N))
