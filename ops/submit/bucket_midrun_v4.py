#!/usr/bin/env python3
"""Bucket midrun_v4 per-uuid run.sbatch into corrected Slurm arrays (proven pattern).

Globs each stack's run dirs (nested `LABEL/run_*/` from --num-runs N, and flat
`LABEL_*/` from the --num-runs 1 HD loop), groups by routing signature, and emits one
self-contained inline-SCRIPTS array per (stack x lane). Headers use --cpus-per-task
(NEVER --ntasks -> SLURM_NTASKS=1, Lightning-safe), explicit --qos per tier, --chdir=code,
--exclude for black-holes, and an `env -u SLURM_NTASKS ...` dispatch guard.

  python3 bucket_midrun_v4.py [--stacks rfd3,boltzgen,bd1] [--submit]   (CAP=N env optional)
"""
import re, os, sys, glob, subprocess, collections

RB = "/cluster/work/wicky/lvalsecchi/rna_binder"
CODE = RB + "/code"
C = RB + "/campaigns/csy4_crrna_free"
OUT = C + "/_mv4_arrays"
LOGDIR = OUT + "/logs"
os.makedirs(LOGDIR, exist_ok=True)
EXCLUDE = "eu-g6-064"          # black-hole node flagged this session
LABEL = "midrun_v4"
CAP = int(os.environ.get("CAP", "0"))

STACKS = collections.OrderedDict([
    ("rfd3", "rfd3__ligandmpnn"),
    ("boltzgen", "boltzgen__boltzif"),
    ("bd1", "boltzdesign__ligandmpnn"),
    ("hd", "halludesign__ligandmpnn"),
])
# the 4 test-array uuids already RUNNING (must not be re-submitted)
EXCLUDE_UUIDS = {"4b1043ef", "dc7fe23b", "fcdc1a0d", "46db14f3"}

SUBMIT = "--submit" in sys.argv
only = None
for a in sys.argv:
    if a.startswith("--stacks="):
        only = set(a.split("=", 1)[1].split(","))


def hdr(path, key):
    pat = re.compile(r"^#SBATCH --%s=(.*)$" % re.escape(key))
    with open(path) as fh:
        for line in fh:
            m = pat.match(line.strip())
            if m:
                return m.group(1)
    return None


def uuid_of(rs):
    m = re.search(r"(?:run_|%s_)([0-9a-f]{8})" % LABEL, rs)
    return m.group(1) if m else None


made = []
for stack, sk in STACKS.items():
    if only and stack not in only:
        continue
    paths = sorted(set(glob.glob(f"{C}/{sk}/{LABEL}/run_*/run.sbatch") +
                       glob.glob(f"{C}/{sk}/{LABEL}_*/run.sbatch")))
    paths = [p for p in paths if uuid_of(p) not in EXCLUDE_UUIDS]
    if not paths:
        print("  %s: no run.sbatch found" % stack)
        continue
    groups = collections.OrderedDict()
    for rs in paths:
        key = (hdr(rs, "partition"), hdr(rs, "time"), hdr(rs, "nodelist"),
               hdr(rs, "gres"), hdr(rs, "mem-per-cpu") or "8G")
        groups.setdefault(key, []).append(rs)
    for (part, time_, nodelist, gres, mem), scripts in groups.items():
        if CAP:
            scripts = scripts[:CAP]
        n = len(scripts)
        base, tier = part.split(".")[0], part.split(".")[1].rstrip("h")
        qos = "es_wicky/%s/%s" % (base, tier)
        sig = ("%s_%s" % (stack, part)).replace(".", "_")
        sb = os.path.join(OUT, "arr_%s.sbatch" % sig)
        L = [
            "#!/bin/bash",
            "#SBATCH --job-name=mv4_%s" % stack,
            "#SBATCH --partition=%s" % part,
            "#SBATCH --qos=%s" % qos,
            "#SBATCH --time=%s" % time_,
            "#SBATCH --cpus-per-task=4",
            "#SBATCH --mem-per-cpu=%s" % mem,
            "#SBATCH --gpus=1",
            "#SBATCH --gres=%s" % gres,
            "#SBATCH --nodelist=%s" % nodelist,
        ]
        if base in ("gpuhe", "gpupr"):
            L.append("#SBATCH --exclude=%s" % EXCLUDE)
        L += [
            "#SBATCH --requeue",
            "#SBATCH --open-mode=append",
            "#SBATCH --array=0-%d" % (n - 1),
            "#SBATCH --chdir=%s" % CODE,
            "#SBATCH --output=%s/%s_%%A_%%a.log" % (LOGDIR, stack),
            "#SBATCH --error=%s/%s_%%A_%%a.err" % (LOGDIR, stack),
            "",
            "set -euo pipefail",
            "SCRIPTS=(",
        ]
        L += ["  %s" % s for s in scripts]
        L += [
            ")",
            'S="${SCRIPTS[$SLURM_ARRAY_TASK_ID]}"',
            'RUNDIR=$(dirname "$S")',
            'echo "[mv4] stack=%s task=$SLURM_ARRAY_TASK_ID host=$(hostname) -> $S"' % stack,
            'env -u SLURM_NTASKS -u SLURM_NPROCS -u SLURM_NTASKS_PER_NODE -u SLURM_TASKS_PER_NODE \\',
            '    bash "$S" > "$RUNDIR/slurm.log" 2> "$RUNDIR/slurm.err"',
        ]
        with open(sb, "w") as fh:
            fh.write("\n".join(L) + "\n")
        print("wrote %-28s stack=%-9s tasks=%-4d part=%-11s qos=%-18s gres=%s"
              % (os.path.basename(sb), stack, n, part, qos, gres))
        made.append(sb)

print("=== %d arrays built ===" % len(made))
if SUBMIT:
    for sb in made:
        r = subprocess.run(["sbatch", "--parsable", sb], capture_output=True, text=True)
        print("submitted %-28s -> %s%s" % (os.path.basename(sb), r.stdout.strip(), r.stderr.strip()))
else:
    print("(dry-run; add --submit to launch)")
