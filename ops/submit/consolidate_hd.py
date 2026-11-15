#!/usr/bin/env python3
"""Consolidate fragmented HD gen (mv7_hd + mv8_hd) pending records back into ONE clean compact
array per (name,tier). Split-off task records (JobID != ArrayJobID, from per-task scontrol ops)
get folded back in. Skips already-packed uuids. Writes arr files + cancel list + submit plan
to a SHARED-FS dir (node-safe). Does NOT cancel/submit — the driver shell does that."""
import subprocess, re, os, glob
from collections import defaultdict

RB = "/cluster/work/wicky/lvalsecchi/rna_binder"
CAMP = RB + "/campaigns/csy4_crrna_free"
HD = CAMP + "/halludesign__ligandmpnn"
OUT = CAMP + "/_consolidate_arrays"; os.makedirs(OUT + "/logs", exist_ok=True)
NAMES = {"mv7_hd", "mv8_hd"}

# already-packed uuids -> skip (don't re-run finished work)
packed = set()
for t in glob.glob(HD + "/midrun_v*/designs/*.tar.zst"):
    packed.add(os.path.basename(t).split(".")[0])
for t in glob.glob(HD + "/midrun_v*_*/designs.tar.zst"):
    m = re.search(r"midrun_v[0-9]+_([0-9a-f]{8})", t)
    if m: packed.add(m.group(1))

def expand_k(k):
    k = k.split("%")[0]; idx = []
    for part in k.split(","):
        if "-" in part:
            a, b = part.split("-")[:2]
            try: idx.extend(range(int(a), int(b) + 1))
            except: pass
        elif part.isdigit(): idx.append(int(part))
    return idx

_arr = {}
def scripts_of(f):
    if f not in _arr:
        try: _arr[f] = [l.strip() for l in open(f) if l.strip().startswith("/cluster") and l.strip().endswith("run.sbatch")]
        except: _arr[f] = []
    return _arr[f]

out = subprocess.run(["squeue","-u","lvalsecchi","-h","-r","-t","PENDING","-o","%j|%P|%K|%o|%i"],
                     capture_output=True, text=True).stdout
groups = defaultdict(list); cancel_ids = set(); seen = set(); skipped_packed = 0
for line in out.splitlines():
    p = line.split("|")
    if len(p) < 5: continue
    name, part, kk, arrfile, jobid = p[0], p[1], p[2], p[3], p[4]
    if name not in NAMES: continue
    cancel_ids.add(jobid)
    sc = scripts_of(arrfile)
    for k in expand_k(kk):
        if 0 <= k < len(sc):
            rs = sc[k]
            m = re.search(r"midrun_v[0-9]+_([0-9a-f]{8})", rs)
            if m and m.group(1) in packed: skipped_packed += 1; continue
            if rs not in seen and os.path.exists(rs):
                seen.add(rs); groups[(name, part)].append(rs)

def qos_of(part):
    b, t = part.split(".")[0], part.split(".")[1].rstrip("h")
    return "es_wicky/%s/%s" % (b, t)

def write_array(name, part, runs):
    n = len(runs); sig = (name + "_" + part).replace(".", "_")
    b, t = part.split(".")[0], part.split(".")[1].rstrip("h")
    wall = "47:30:00" if t == "120" else "23:30:00"
    gres = "gpumem:24g" if b == "gpuhe" else "gpumem:40g"
    nodes = "eu-g4-[001-031],eu-g6-[001-080]" if b == "gpuhe" else "eu-a65-[02-03,05-07],eu-ts-[01-02]"
    sb = OUT + "/arr_%s.sbatch" % sig
    L = ["#!/bin/bash", "#SBATCH --job-name=%s" % name, "#SBATCH --partition=%s" % part,
         "#SBATCH --qos=%s" % qos_of(part), "#SBATCH --time=%s" % wall, "#SBATCH --cpus-per-task=4",
         "#SBATCH --mem-per-cpu=7G", "#SBATCH --gpus=1", "#SBATCH --gres=%s" % gres,
         "#SBATCH --nodelist=%s" % nodes, "#SBATCH --exclude=eu-g6-064", "#SBATCH --requeue",
         "#SBATCH --open-mode=append", "#SBATCH --array=0-%d" % (n - 1), "#SBATCH --chdir=%s/code" % RB,
         "#SBATCH --output=%s/logs/%s_%%A_%%a.log" % (OUT, sig),
         "#SBATCH --error=%s/logs/%s_%%A_%%a.err" % (OUT, sig),
         "", "set -euo pipefail", "SCRIPTS=("]
    L += ["  %s" % r for r in runs]
    L += [")", 'S="${SCRIPTS[$SLURM_ARRAY_TASK_ID]}"', 'RUNDIR=$(dirname "$S")',
          'echo "[consol] task=$SLURM_ARRAY_TASK_ID host=$(hostname) -> $S"',
          'env -u SLURM_NTASKS -u SLURM_NPROCS -u SLURM_NTASKS_PER_NODE -u SLURM_TASKS_PER_NODE \\',
          '    bash "$S" > "$RUNDIR/slurm.log" 2> "$RUNDIR/slurm.err"']
    open(sb, "w").write("\n".join(L) + "\n")
    return sb, n

print("pending HD task-RECORDS to cancel: %d   (skipped already-packed: %d)" % (len(cancel_ids), skipped_packed))
plan = []
for (name, part), runs in sorted(groups.items()):
    if not runs: continue
    sb, n = write_array(name, part, runs)
    plan.append("%s|%s|%s|%d" % (name, part, sb, n))
    print("  %-8s %-12s -> %4d uuids   %s" % (name, part, n, os.path.basename(sb)))
open(OUT + "/cancel_ids.txt", "w").write(" ".join(sorted(cancel_ids)) + "\n")
open(OUT + "/submit_plan.txt", "w").write("\n".join(plan) + "\n")
print("wrote %s/{cancel_ids.txt,submit_plan.txt} + %d arr files" % (OUT, len(plan)))
