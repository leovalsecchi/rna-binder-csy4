#!/usr/bin/env python3
"""Bucket midrun_v9 per-uuid run.sbatch into corrected Slurm arrays (proven pattern).

Globs each stack's run dirs (nested `LABEL/run_*/` from --num-runs N, and flat
`LABEL_*/` from the --num-runs 1 HD loop), groups by routing signature, and emits one
self-contained inline-SCRIPTS array per (stack x lane). Headers use --cpus-per-task
(NEVER --ntasks -> SLURM_NTASKS=1, Lightning-safe), explicit --qos per tier, --chdir=code,
--exclude for black-holes, and an `env -u SLURM_NTASKS ...` dispatch guard.

  python3 bucket_midrun_v9.py [--stacks rfd3,boltzgen,bd1] [--submit]   (CAP=N env optional)
"""
import re, os, sys, glob, subprocess, collections

RB = "/cluster/work/wicky/lvalsecchi/rna_binder"
CODE = RB + "/code"
C = RB + "/campaigns/csy4_crrna_free"
OUT = C + "/_mv9_hd_topup_arrays"
LOGDIR = OUT + "/logs"
os.makedirs(LOGDIR, exist_ok=True)
EXCLUDE = "eu-g6-064"          # black-hole node flagged this session
LABEL = "midrun_v9"
CAP = int(os.environ.get("CAP", "0"))

STACKS = collections.OrderedDict([
    ("rfd3", "rfd3__ligandmpnn"),
    ("boltzgen", "boltzgen__boltzif"),
    ("bd1", "boltzdesign__ligandmpnn"),
    ("hd", "halludesign__ligandmpnn"),
])
# the 4 test-array uuids already RUNNING (must not be re-submitted)
# Load ALREADY-COMPLETE uuids so finished work is never resubmitted. The bucketer globs
# every run.sbatch and groups only by routing signature — without this it would re-run the
# 427 done HD (~7 h each) and 195 done BD1 (~6.4 h each) uuids, and BD1's multif self-cleans
# a uuid's work dir on restart, regenerating over completed slices.
_DONE_FILE = os.environ.get("MV9_DONE_FILE", RB + "/ops/submit/mv9_done_uuids.txt")
EXCLUDE_UUIDS = set()
if os.path.exists(_DONE_FILE):
    with open(_DONE_FILE) as _fh:
        EXCLUDE_UUIDS = {l.strip() for l in _fh
                         if l.strip() and not l.startswith("#")}
    print("excluding %d already-complete uuids (%s)" % (len(EXCLUDE_UUIDS), _DONE_FILE))
else:
    sys.exit("FATAL: %s missing — refusing to run without a completion filter" % _DONE_FILE)

_ALLOW_FILE = os.environ.get("MV9_TOPUP_MANIFEST",
                            RB + "/ops/submit/mv9_hd_topup_seed_manifest.txt")
if not os.path.exists(_ALLOW_FILE):
    sys.exit("FATAL: %s missing - refusing to bucket without an allowlist" % _ALLOW_FILE)
ALLOW = set()
with open(_ALLOW_FILE) as _fh:
    for _l in _fh:
        _p = _l.split()
        if len(_p) >= 2 and _p[1] != "NA":
            ALLOW.add(_p[1])
print("ALLOWLIST: %d topup uuids (%s)" % (len(ALLOW), _ALLOW_FILE))

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
    paths = [p for p in paths if uuid_of(p) in ALLOW]
    print("  %s: matched %d of %d topup uuids" % (stack, len(paths), len(ALLOW)))
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
            "#SBATCH --job-name=mv9_%s_topup" % stack,
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
            'echo "[mv7] stack=%s task=$SLURM_ARRAY_TASK_ID host=$(hostname) -> $S"' % stack,
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
