#!/usr/bin/env python3
"""Re-consolidate the fragmented mv7_hd singletons on gpuhe.bulk into ONE compact array.
Collect each pending task's run.sbatch (via its source array's SCRIPTS list), write a
compact --cpus-per-task array (Lightning-safe), and emit the cancel list. Does NOT cancel
or submit — caller does that explicitly."""
import subprocess, re, os
from collections import defaultdict

RB = "/cluster/work/wicky/lvalsecchi/rna_binder"
OUTDIR = RB + "/campaigns/csy4_crrna_free/_mv7_arrays"
LOGDIR = OUTDIR + "/logs"
os.makedirs(LOGDIR, exist_ok=True)

out = subprocess.run(["squeue","-u","lvalsecchi","-h","-r","-o","%i|%j|%P|%T"],
                     capture_output=True, text=True).stdout
tasks = []
for l in out.splitlines():
    p = l.split("|")
    if len(p) >= 4 and p[1] == "mv7_hd" and p[2] == "gpuhe.bulk" and p[3] == "PENDING":
        tasks.append(p[0])

byarr = defaultdict(list)
for t in tasks:
    arr, tid = t.rsplit("_", 1)
    byarr[arr].append(int(tid))

# map task-index -> run.sbatch using each source array's inline SCRIPTS list
task_to_run = {}
for arr, tids in byarr.items():
    cmd = subprocess.run(["scontrol","show","job","%s_%d" % (arr, tids[0])],
                         capture_output=True, text=True).stdout
    m = re.search(r"Command=(\S+)", cmd)
    if not m:
        print("!! no Command for array", arr); continue
    arrfile = m.group(1)
    scripts = [ln.strip() for ln in open(arrfile)
               if ln.strip().startswith("/cluster") and ln.strip().endswith("run.sbatch")]
    for tid in tids:
        if 0 <= tid < len(scripts):
            task_to_run["%s_%d" % (arr, tid)] = scripts[tid]
        else:
            print("!! task %s_%d out of range (%d scripts)" % (arr, tid, len(scripts)))

runs = sorted(set(task_to_run.values()))
missing = [r for r in runs if not os.path.exists(r)]
runs = [r for r in runs if os.path.exists(r)]
print("pending tasks=%d  source arrays=%s  unique run.sbatch=%d  (dropped missing=%d)"
      % (len(tasks), list(byarr.keys()), len(runs), len(missing)))

n = len(runs)
sb = OUTDIR + "/arr_mv7hd_reconsolidated.sbatch"
L = [
  "#!/bin/bash",
  "#SBATCH --job-name=mv7_hd",
  "#SBATCH --partition=gpuhe.bulk",
  "#SBATCH --qos=es_wicky/gpuhe/bulk",
  "#SBATCH --time=23:30:00",
  "#SBATCH --cpus-per-task=4",
  "#SBATCH --mem-per-cpu=7G",
  "#SBATCH --gpus=1",
  "#SBATCH --gres=gpumem:24g",
  "#SBATCH --nodelist=eu-g4-[001-031],eu-g6-[001-080]",
  "#SBATCH --exclude=eu-g6-064",
  "#SBATCH --requeue",
  "#SBATCH --open-mode=append",
  "#SBATCH --array=0-%d" % (n - 1),
  "#SBATCH --chdir=%s/code" % RB,
  "#SBATCH --output=%s/mv7hd_recon_%%A_%%a.log" % LOGDIR,
  "#SBATCH --error=%s/mv7hd_recon_%%A_%%a.err" % LOGDIR,
  "",
  "set -euo pipefail",
  "SCRIPTS=(",
]
L += ["  %s" % r for r in runs]
L += [
  ")",
  'S="${SCRIPTS[$SLURM_ARRAY_TASK_ID]}"',
  'RUNDIR=$(dirname "$S")',
  'echo "[mv7hd_recon] task=$SLURM_ARRAY_TASK_ID host=$(hostname) -> $S"',
  'env -u SLURM_NTASKS -u SLURM_NPROCS -u SLURM_NTASKS_PER_NODE -u SLURM_TASKS_PER_NODE \\',
  '    bash "$S" > "$RUNDIR/slurm.log" 2> "$RUNDIR/slurm.err"',
]
open(sb, "w").write("\n".join(L) + "\n")
open("/tmp/mv7bulk_cancel.txt", "w").write("\n".join(tasks) + "\n")
print("wrote %s  (array 0-%d)" % (sb, n - 1))
print("cancel list -> /tmp/mv7bulk_cancel.txt (%d tasks)" % len(tasks))
