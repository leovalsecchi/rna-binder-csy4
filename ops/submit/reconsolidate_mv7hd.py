#!/usr/bin/env python3
"""Consolidate fragmented mv7_hd singletons on a given tier into ONE compact array.
Usage: reconsolidate_mv7hd.py <partition> <time>   e.g. gpuhe.24h 23:30:00
Writes the array sbatch + prints the cancel IDs to stdout (caller cancels + submits)."""
import subprocess, re, os, sys
from collections import defaultdict

PART = sys.argv[1]           # e.g. gpuhe.24h  (submit tier; for .bulk submit routes to .24h then base-id move)
WALL = sys.argv[2] if len(sys.argv) > 2 else "23:30:00"
base, tier = PART.split(".")[0], PART.split(".")[1].rstrip("h")
QOS = "es_wicky/%s/%s" % (base, tier)

RB = "/cluster/work/wicky/lvalsecchi/rna_binder"
OUTDIR = RB + "/campaigns/csy4_crrna_free/_mv7_arrays"
LOGDIR = OUTDIR + "/logs"; os.makedirs(LOGDIR, exist_ok=True)

out = subprocess.run(["squeue","-u","lvalsecchi","-h","-r","-o","%i|%j|%P|%T"],
                     capture_output=True, text=True).stdout
tasks = [l.split("|")[0] for l in out.splitlines()
         if len(l.split("|")) >= 4 and l.split("|")[1] == "mv7_hd"
         and l.split("|")[2] == PART and l.split("|")[3] == "PENDING"]
byarr = defaultdict(list)
for t in tasks:
    a, tid = t.rsplit("_", 1); byarr[a].append(int(tid))

t2r = {}
for a, tids in byarr.items():
    cmd = subprocess.run(["scontrol","show","job","%s_%d" % (a, tids[0])],
                         capture_output=True, text=True).stdout
    m = re.search(r"Command=(\S+)", cmd)
    if not m: continue
    scripts = [ln.strip() for ln in open(m.group(1))
               if ln.strip().startswith("/cluster") and ln.strip().endswith("run.sbatch")]
    for tid in tids:
        if 0 <= tid < len(scripts): t2r["%s_%d" % (a, tid)] = scripts[tid]

runs = sorted({r for r in t2r.values() if os.path.exists(r)})
sig = PART.replace(".", "_")
sb = OUTDIR + "/arr_mv7hd_recon_%s.sbatch" % sig
n = len(runs)
if n == 0:
    print("NOTHING to consolidate on %s" % PART); sys.exit(0)
L = ["#!/bin/bash", "#SBATCH --job-name=mv7_hd", "#SBATCH --partition=%s" % PART,
     "#SBATCH --qos=%s" % QOS, "#SBATCH --time=%s" % WALL, "#SBATCH --cpus-per-task=4",
     "#SBATCH --mem-per-cpu=7G", "#SBATCH --gpus=1", "#SBATCH --gres=gpumem:24g",
     "#SBATCH --nodelist=eu-g4-[001-031],eu-g6-[001-080]", "#SBATCH --exclude=eu-g6-064",
     "#SBATCH --requeue", "#SBATCH --open-mode=append", "#SBATCH --array=0-%d" % (n-1),
     "#SBATCH --chdir=%s/code" % RB,
     "#SBATCH --output=%s/mv7hd_recon_%s_%%A_%%a.log" % (LOGDIR, sig),
     "#SBATCH --error=%s/mv7hd_recon_%s_%%A_%%a.err" % (LOGDIR, sig),
     "", "set -euo pipefail", "SCRIPTS=("]
L += ["  %s" % r for r in runs]
L += [")", 'S="${SCRIPTS[$SLURM_ARRAY_TASK_ID]}"', 'RUNDIR=$(dirname "$S")',
      'echo "[mv7hd_recon] task=$SLURM_ARRAY_TASK_ID host=$(hostname) -> $S"',
      'env -u SLURM_NTASKS -u SLURM_NPROCS -u SLURM_NTASKS_PER_NODE -u SLURM_TASKS_PER_NODE \\',
      '    bash "$S" > "$RUNDIR/slurm.log" 2> "$RUNDIR/slurm.err"']
open(sb, "w").write("\n".join(L) + "\n")
print("ARRAY_FILE=%s" % sb)
print("ARRAY_SIZE=%d" % n)
print("CANCEL_IDS=" + " ".join(tasks))
