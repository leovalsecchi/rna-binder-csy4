#!/usr/bin/env python3
"""Re-bucket the existing bif_fill per-uuid run.sbatch into corrected Slurm arrays.

The launch.py array wrapper set `#SBATCH --ntasks=4`, which exports SLURM_NTASKS=4
into the job env; BoltzGen (PyTorch-Lightning) then refuses to start. This rebuilds
one array per lane using `--cpus-per-task=4` (SLURM_NTASKS=1) plus an `env -u
SLURM_NTASKS ...` guard on the inner dispatch. Reuses existing run dirs (no new uuids).

Write-only by default; pass --submit to sbatch them.
"""
import re, glob, os, sys, subprocess, collections

LABEL_DIR = "/cluster/work/wicky/lvalsecchi/rna_binder/campaigns/csy4_crrna_free/boltzgen__boltzif/bif_fill_2026-07-03"
REARRAY = os.path.join(LABEL_DIR, "_rearray")
LOGDIR = os.path.join(REARRAY, "array_logs")
os.makedirs(LOGDIR, exist_ok=True)

# per-lane cap -> 18+18+12+12+15 = 75 uuids
CAPS = {"gpu.24h": 18, "gpu.120h": 18, "gpuhe.24h": 12, "gpuhe.120h": 12, "gpuhe.4h": 15}
SUBMIT = "--submit" in sys.argv


def hdr(path, key):
    pat = re.compile(r"^#SBATCH --%s=(.*)$" % re.escape(key))
    with open(path) as fh:
        for line in fh:
            m = pat.match(line.strip())
            if m:
                return m.group(1)
    return None


runs = sorted(glob.glob(os.path.join(LABEL_DIR, "run_*", "run.sbatch")))
buckets = collections.OrderedDict((p, []) for p in CAPS)
for r in runs:
    part = hdr(r, "partition")
    if part in CAPS and len(buckets[part]) < CAPS[part]:
        buckets[part].append(r)

made = []
for part, scripts in buckets.items():
    if not scripts:
        print("!! no run.sbatch found for %s" % part)
        continue
    rep = scripts[0]
    time_ = hdr(rep, "time")
    nodelist = hdr(rep, "nodelist")
    gres = hdr(rep, "gres")
    mem = hdr(rep, "mem-per-cpu") or "16G"
    n = len(scripts)
    sig = part.replace(".", "_")
    sb = os.path.join(REARRAY, "array_%s.sbatch" % sig)
    L = [
        "#!/bin/bash",
        "#SBATCH --job-name=bg_bif_fill_%s" % sig,
        "#SBATCH --requeue",
        "#SBATCH --open-mode=append",
        "#SBATCH --partition=%s" % part,
        "#SBATCH --time=%s" % time_,
        "#SBATCH --cpus-per-task=4",         # NOT --ntasks: keeps SLURM_NTASKS=1
        "#SBATCH --mem-per-cpu=%s" % mem,
        "#SBATCH --gpus=1",
        "#SBATCH --nodelist=%s" % nodelist,
        "#SBATCH --gres=%s" % gres,
        "#SBATCH --array=0-%d" % (n - 1),
        "#SBATCH --output=%s/%%A_%%a.log" % LOGDIR,
        "#SBATCH --error=%s/%%A_%%a.err" % LOGDIR,
        "",
        "set -euo pipefail",
        "SCRIPTS=(",
    ]
    L += ["  %s" % s for s in scripts]
    L += [
        ")",
        'S="${SCRIPTS[$SLURM_ARRAY_TASK_ID]}"',
        'RUNDIR=$(dirname "$S")',
        'echo "[gen-array] task=$SLURM_ARRAY_TASK_ID host=$(hostname) part=%s -> $S"' % part,
        "# strip SLURM_NTASKS et al so BoltzGen/Lightning sees a single-process run",
        'env -u SLURM_NTASKS -u SLURM_NPROCS -u SLURM_NTASKS_PER_NODE -u SLURM_TASKS_PER_NODE \\',
        '    bash "$S" > "$RUNDIR/slurm.log" 2> "$RUNDIR/slurm.err"',
    ]
    with open(sb, "w") as fh:
        fh.write("\n".join(L) + "\n")
    print("wrote %-42s tasks=%-3d part=%-11s time=%-9s gres=%s" % (os.path.basename(sb), n, part, time_, gres))
    made.append(sb)

total = sum(len(v) for v in buckets.values())
print("=== %d arrays, %d total uuids ===" % (len(made), total))
if SUBMIT:
    for sb in made:
        out = subprocess.run(["sbatch", "--parsable", sb], capture_output=True, text=True)
        print("submitted %s -> %s%s" % (os.path.basename(sb), out.stdout.strip(), out.stderr.strip()))
else:
    print("(dry-run; re-run with --submit to launch)")
