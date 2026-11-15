#!/usr/bin/env python3
"""hd_consolidate — self-healing HD-gen array de-fragmenter (loop worker).

WHY. Per-task `scontrol requeue` of a stuck HD task (hd_watcher restarting an
evoformer-JIT / wedged-node hang) PEELS that task into its own JobID. Slurm's main
scheduler evaluates only the top ~250 pending jobs by priority, so a pile of peeled
singletons clogs that window and starves gpuhe/gpupr. `squeue -o %i` HIDES the peels
by re-folding them into `_[range]` displays, but the scheduler and job_summary count
them as separate records.

DETECT (fold-proof). records = distinct JobID via `squeue -r -o %A`; bases = distinct
ArrayJobID `%F`. peels = records - bases. (A compact array's tasks all share the base
JobID; every peeled task carries a distinct JobID but its parent's ArrayJobID.) `%i`
folding cannot fool this.

HEAL. When peels >= --threshold, re-pack ALL pending mv*_hd into ONE compact array per
(name,partition): collect pending tasks -> map task-idx to its per-uuid run.sbatch via
the live array's inline SCRIPTS=() -> skip already-packed uuids -> cancel by NAME+STATE
(server-side, never touches RUNNING) -> resubmit one compact array per group -> base-id
scontrol-move the .bulk groups to their bulk partition -> nice the non-primary campaign
below the primary. All scontrol ops are on the BASE array id (never a task element), so
the pack itself never peels.

Prints the loop-framework summary token: `needing-work=0` when clean, `submitted=<N>`
(and `moved_to_bulk=<M>`) when it re-packs. DRY-RUN unless --apply.
"""
import subprocess, re, os, glob, time, argparse
from collections import defaultdict

RB   = "/cluster/work/wicky/lvalsecchi/rna_binder"
CAMP = RB + "/campaigns/csy4_crrna_free"
HD   = CAMP + "/halludesign__ligandmpnn"
OUT  = CAMP + "/_consolidate_arrays"
USER = "lvalsecchi"
NAME_RE = re.compile(r"^mv\d+_hd$")


def sh(*a):
    return subprocess.run(a, capture_output=True, text=True).stdout


# ---- fold-proof fragmentation detector -------------------------------------
def measure():
    recs = set()
    for x in sh("squeue", "-u", USER, "-h", "-t", "PENDING", "-r", "-o", "%A|%j").splitlines():
        p = x.split("|")
        if len(p) == 2 and NAME_RE.match(p[1]):
            recs.add(p[0])
    bases = set()
    for x in sh("squeue", "-u", USER, "-h", "-t", "PENDING", "-o", "%F|%j").splitlines():
        p = x.split("|")
        if len(p) == 2 and NAME_RE.match(p[1]):
            bases.add(p[0])
    return len(recs), len(bases), len(recs) - len(bases)


# ---- collect + write compact arrays (proven in the 2026-07-23 manual pack) ---
def expand_k(k):
    k = k.split("%")[0]; idx = []
    for part in k.split(","):
        if "-" in part:
            a, b = part.split("-")[:2]
            try: idx.extend(range(int(a), int(b) + 1))
            except ValueError: pass
        elif part.isdigit():
            idx.append(int(part))
    return idx


_arr = {}
def scripts_of(f):
    if f not in _arr:
        try:
            _arr[f] = [l.strip() for l in open(f)
                       if l.strip().startswith("/cluster") and l.strip().endswith("run.sbatch")]
        except OSError:
            _arr[f] = []
    return _arr[f]


def packed_set():
    s = set()
    for t in glob.glob(HD + "/midrun_v*/designs/*.tar.zst"):
        s.add(os.path.basename(t).split(".")[0])
    for t in glob.glob(HD + "/midrun_v*_*/designs.tar.zst"):
        m = re.search(r"midrun_v[0-9]+_([0-9a-f]{8})", t)
        if m: s.add(m.group(1))
    return s


def qos_of(part):
    b, t = part.split(".")[0], part.split(".")[1].rstrip("h")
    return "es_wicky/%s/%s" % (b, t)


def collect():
    """Read live pending mv*_hd tasks -> {(name,part): [run.sbatch,...]} (packed uuids skipped).
    Caches each live array's inline SCRIPTS=() BEFORE any cancel, so overwriting the arr
    file afterward is race-free."""
    packed = packed_set()
    groups = defaultdict(list); seen = set(); skipped = 0
    for line in sh("squeue", "-u", USER, "-h", "-r", "-t", "PENDING",
                   "-o", "%j|%P|%K|%o").splitlines():
        p = line.split("|")
        if len(p) < 4: continue
        name, part, kk, arrfile = p[0], p[1], p[2], p[3]
        if not NAME_RE.match(name): continue
        sc = scripts_of(arrfile)
        for k in expand_k(kk):
            if 0 <= k < len(sc):
                rs = sc[k]
                m = re.search(r"midrun_v[0-9]+_([0-9a-f]{8})", rs)
                if m and m.group(1) in packed:
                    skipped += 1; continue
                if rs not in seen and os.path.exists(rs):
                    seen.add(rs); groups[(name, part)].append(rs)
    return groups, skipped


def write_array(packdir, name, part, runs):
    n = len(runs); sig = (name + "_" + part).replace(".", "_")
    b, t = part.split(".")[0], part.split(".")[1].rstrip("h")
    wall = "47:30:00" if t == "120" else "23:30:00"
    gres = "gpumem:24g" if b == "gpuhe" else "gpumem:40g"
    nodes = "eu-g4-[001-031],eu-g6-[001-080]" if b == "gpuhe" else "eu-a65-[02-03,05-07],eu-ts-[01-02]"
    # UNIQUE per-pack path: never overwrite an arr file a still-live older array's tasks
    # reference (a re-pended/preempted old task would otherwise read the new SCRIPTS -> wrong uuid).
    sb = packdir + "/arr_%s.sbatch" % sig
    L = ["#!/bin/bash", "#SBATCH --job-name=%s" % name, "#SBATCH --partition=%s" % part,
         "#SBATCH --qos=%s" % qos_of(part), "#SBATCH --time=%s" % wall, "#SBATCH --cpus-per-task=4",
         "#SBATCH --mem-per-cpu=7G", "#SBATCH --gpus=1", "#SBATCH --gres=%s" % gres,
         "#SBATCH --nodelist=%s" % nodes, "#SBATCH --exclude=eu-g6-064", "#SBATCH --requeue",
         "#SBATCH --open-mode=append", "#SBATCH --array=0-%d" % (n - 1), "#SBATCH --chdir=%s/code" % RB,
         "#SBATCH --output=%s/logs/%s_%%A_%%a.log" % (packdir, sig),
         "#SBATCH --error=%s/logs/%s_%%A_%%a.err" % (packdir, sig),
         "", "set -euo pipefail", "SCRIPTS=("]
    L += ["  %s" % r for r in runs]
    L += [")", 'S="${SCRIPTS[$SLURM_ARRAY_TASK_ID]}"', 'RUNDIR=$(dirname "$S")',
          'echo "[consol] task=$SLURM_ARRAY_TASK_ID host=$(hostname) -> $S"',
          'env -u SLURM_NTASKS -u SLURM_NPROCS -u SLURM_NTASKS_PER_NODE -u SLURM_TASKS_PER_NODE \\',
          '    bash "$S" > "$RUNDIR/slurm.log" 2> "$RUNDIR/slurm.err"']
    open(sb, "w").write("\n".join(L) + "\n")
    return sb, n


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--threshold", type=int, default=20,
                    help="re-pack only when peeled records reach this many (bounds churn)")
    ap.add_argument("--primary", default="mv7_hd",
                    help="campaign kept at full priority; all other mv*_hd niced below it")
    ap.add_argument("--nice", type=int, default=1500)
    ap.add_argument("--apply", action="store_true")
    a = ap.parse_args()

    recs, bases, peels = measure()
    print("[hd_consolidate] %s records=%d bases=%d peels=%d thresh=%d apply=%s"
          % (time.strftime("%F %T"), recs, bases, peels, a.threshold, a.apply))
    if peels < a.threshold:
        print("  peels<thresh -> needing-work=0")
        return

    groups, skipped = collect()
    if not groups:
        print("  fragmented but 0 unpacked uuids to re-pack -> needing-work=0")
        return
    ntot = sum(len(v) for v in groups.values())
    print("  fragmented: re-pack %d uuids -> %d compact arrays (skipped already-packed=%d)"
          % (ntot, len(groups), skipped))
    for (name, part), runs in sorted(groups.items()):
        print("    %-8s %-12s %4d uuids" % (name, part, len(runs)))
    if not a.apply:
        print("  DRY-RUN -> needing-work=%d (would re-pack)" % ntot)
        return

    # cancel by NAME+STATE (server-side; never touches RUNNING), then drain
    names = sorted({n for (n, _) in groups})
    for nm in names:
        sh("scancel", "-u", USER, "-t", "PENDING", "--name", nm)
    for _ in range(25):
        time.sleep(6)
        r, _b, _p = measure()
        if r == 0:
            break
        for nm in names:
            sh("scancel", "-u", USER, "-t", "PENDING", "--name", nm)

    # write compact arrays to a UNIQUE per-pack dir so a still-live older array's arr file
    # is never overwritten (re-pended/preempted old tasks keep reading their correct SCRIPTS)
    packdir = OUT + "/packs/" + time.strftime("%Y%m%d_%H%M%S")
    os.makedirs(packdir + "/logs", exist_ok=True)
    nsub = moved = 0
    for (name, part), runs in sorted(groups.items()):
        sb, n = write_array(packdir, name, part, runs)
        jid = sh("sbatch", "--parsable", sb).strip()
        if not jid.isdigit():
            print("    submit FAILED %-8s %-12s (%s)" % (name, part, jid[:60]))
            continue
        nsub += n
        tags = []
        if part.endswith(".bulk"):
            base = part.split(".")[0]
            sh("scontrol", "update", "jobid=%s" % jid, "partition=%s" % part,
               "qos=es_wicky/%s/bulk" % base)
            moved += 1; tags.append("bulk")
        if name != a.primary:
            sh("scontrol", "update", "jobid=%s" % jid, "nice=%d" % a.nice)
            tags.append("nice%d" % a.nice)
        print("    submitted %-8s %-12s (%4d) -> %s%s"
              % (name, part, n, jid, (" [" + ",".join(tags) + "]") if tags else ""))
    print("  submitted=%d moved_to_bulk=%d arrays=%d" % (nsub, moved, len(groups)))


if __name__ == "__main__":
    main()
