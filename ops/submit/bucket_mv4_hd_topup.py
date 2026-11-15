#!/usr/bin/env python3
"""Bucket ONLY the midrun_v4 HD topup uuids into Slurm arrays — same proven pattern as
bucket_midrun_v4.py (--cpus-per-task, never --ntasks; explicit --qos; env-guard dispatch;
black-hole exclude), but restricted to an ALLOWLIST read from the topup seed manifest so
the 600 already-completed HD gens are NEVER re-bucketed/resubmitted. Writes to a separate
_mv4_hd_topup_arrays/ dir.

  python3 bucket_mv4_hd_topup.py [--submit]
"""
import re, os, sys, glob, subprocess, collections

RB = "/cluster/work/wicky/lvalsecchi/rna_binder"
CODE = RB + "/code"
C = RB + "/campaigns/csy4_crrna_free"
SK = "halludesign__ligandmpnn"
OUT = C + "/_mv4_hd_topup_arrays"
LOGDIR = OUT + "/logs"
os.makedirs(LOGDIR, exist_ok=True)
EXCLUDE = "eu-g6-064"
LABEL = "midrun_v4"
MANIFEST = RB + "/ops/submit/mv4_hd_topup_seed_manifest.txt"
SUBMIT = "--submit" in sys.argv

# allowlist = uuids from the topup manifest (col 2), the ONLY dirs we touch
ALLOW = set()
with open(MANIFEST) as fh:
    for line in fh:
        parts = line.split("\t")
        if len(parts) >= 2 and re.fullmatch(r"[0-9a-f]{8}", parts[1].strip()):
            ALLOW.add(parts[1].strip())
print("allowlist uuids (%d): %s" % (len(ALLOW), ",".join(sorted(ALLOW))))


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


paths = sorted(set(glob.glob(f"{C}/{SK}/{LABEL}_*/run.sbatch")))
paths = [p for p in paths if uuid_of(p) in ALLOW]
print("matched %d run.sbatch (of allowlist %d)" % (len(paths), len(ALLOW)))
missing = ALLOW - {uuid_of(p) for p in paths}
if missing:
    print("  WARNING: no run.sbatch for uuids: %s" % ",".join(sorted(missing)))

groups = collections.OrderedDict()
for rs in paths:
    key = (hdr(rs, "partition"), hdr(rs, "time"), hdr(rs, "nodelist"),
           hdr(rs, "gres"), hdr(rs, "mem-per-cpu") or "8G")
    groups.setdefault(key, []).append(rs)

made = []
for (part, time_, nodelist, gres, mem), scripts in groups.items():
    n = len(scripts)
    base, tier = part.split(".")[0], part.split(".")[1].rstrip("h")
    qos = "es_wicky/%s/%s" % (base, tier)
    sig = ("hd_topup_%s" % part).replace(".", "_")
    sb = os.path.join(OUT, "arr_%s.sbatch" % sig)
    L = [
        "#!/bin/bash",
        "#SBATCH --job-name=mv4_hd_topup",
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
        "#SBATCH --output=%s/hd_topup_%%A_%%a.log" % LOGDIR,
        "#SBATCH --error=%s/hd_topup_%%A_%%a.err" % LOGDIR,
        "",
        "set -euo pipefail",
        "SCRIPTS=(",
    ]
    L += ["  %s" % s for s in scripts]
    L += [
        ")",
        'S="${SCRIPTS[$SLURM_ARRAY_TASK_ID]}"',
        'RUNDIR=$(dirname "$S")',
        'echo "[mv4_hd_topup] task=$SLURM_ARRAY_TASK_ID host=$(hostname) -> $S"',
        'env -u SLURM_NTASKS -u SLURM_NPROCS -u SLURM_NTASKS_PER_NODE -u SLURM_TASKS_PER_NODE \\',
        '    bash "$S" > "$RUNDIR/slurm.log" 2> "$RUNDIR/slurm.err"',
    ]
    with open(sb, "w") as fh:
        fh.write("\n".join(L) + "\n")
    print("wrote %-30s tasks=%-3d part=%-11s qos=%-18s gres=%s"
          % (os.path.basename(sb), n, part, qos, gres))
    made.append(sb)

print("=== %d topup arrays built ===" % len(made))
if SUBMIT:
    for sb in made:
        r = subprocess.run(["sbatch", "--parsable", sb], capture_output=True, text=True)
        print("submitted %-30s -> %s%s" % (os.path.basename(sb), r.stdout.strip(), r.stderr.strip()))
else:
    print("(dry-run; add --submit to launch)")
