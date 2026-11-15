#!/usr/bin/env python3
"""array_packer — re-pack PENDING jobs that peeled out of arrays back into compact arrays.

Slurm budgets (Euler): default_queue_depth=250, bf_max_job_assoc=50, bf_max_job_user_part=50.
Peeled singletons (from preempt-requeue or per-task scontrol) each consume one slot, so a
fragmented array starves everything behind it -> free GPUs while jobs sit at Reason=Priority.

Per base array: read its array.sbatch (%o) + manifest.txt, map PENDING task indices -> leaf
run.sbatch, drop already-done leaves, then regroup by (jobname,partition,qos,timelimit) and
submit ONE compact array per signature. All ops on BASE ids so the pack never peels.
"""
import argparse, re, subprocess
from collections import defaultdict
from datetime import datetime
from pathlib import Path

def sh(*a): return subprocess.run(a, capture_output=True, text=True)

ap = argparse.ArgumentParser()
ap.add_argument("--user", default="lvalsecchi")
ap.add_argument("--name", default="", help="only pack jobs whose name matches this substring")
ap.add_argument("--min-records", type=int, default=20, help="only pack a base array with >= this many pending RECORDS")
ap.add_argument("--apply", action="store_true")
a = ap.parse_args()

# Enumerate from `scontrol show job -o`, NOT squeue. squeue re-groups records that share a
# display prefix back into ONE compact line (e.g. "9073455_[0-1,3-7,22-23,...]"), so it
# reported 37 lines where the scheduler actually held 483 records -- a 13x under-count that
# kept --min-records from ever tripping while gpuhe.bulk sat at 417 peeled records.
# scontrol prints exactly ONE line per REAL record, which is what the scheduler budgets
# (default_queue_depth=250, bf_max_job_assoc / bf_max_job_user_part=50) are spent on.
# Do NOT "fix" this with `squeue -r` instead: that flattens to TASKS, making a healthy
# compact array (1 record, 642 tasks) look identical to 642 peels, and the packer would
# start cancelling and resubmitting perfectly good arrays.
WANT = ("JobId", "ArrayTaskId", "JobName", "Partition", "QOS", "TimeLimit", "Command",
        "ReqNodeList")
base = defaultdict(list)
for ln in sh("scontrol", "show", "job", "-o").stdout.splitlines():
    if f"UserId={a.user}(" not in ln or "JobState=PENDING" not in ln: continue
    # A HELD job is PENDING with Reason=JobHeldUser. Repacking it would cancel it and
    # resubmit an UN-HELD array, silently defeating `scontrol hold`. Held = deliberate.
    if "Reason=JobHeldUser" in ln: continue
    f = {}
    for tok in ln.split():
        k, sep, v = tok.partition("=")
        if sep and k in WANT: f.setdefault(k, v)
    jid, idx, name = f.get("JobId", ""), f.get("ArrayTaskId", ""), f.get("JobName", "")
    part, qos = f.get("Partition", ""), f.get("QOS", "")
    tl, cmd = f.get("TimeLimit", ""), f.get("Command", "")
    nodelist = f.get("ReqNodeList", "")
    if nodelist in ("(null)", "None"): nodelist = ""
    if a.name and a.name not in name: continue
    # ANY .sbatch/.sh array script, not just the manifest-shaped "array.sbatch" — the
    # bucketer's gen arrays and the apo array are equally packable (see header).
    if not (cmd.endswith(".sbatch") or cmd.endswith(".sh")): continue
    # KEY on the array.sbatch path: a peeled task gets a NEW ArrayJobId, so keying on that
    # would split one fragmented array into N groups of 1. The command path is stable.
    base[(name, part, qos, tl, cmd, nodelist)].append((jid, idx))

def expand(idx):
    out = []
    for tok in idx.split(","):
        tok = tok.split("%")[0]
        if "-" in tok:
            s, e = tok.split("-")[:2]
            if s.isdigit() and e.isdigit(): out += list(range(int(s), int(e)+1))
        elif tok.isdigit(): out.append(int(tok))
    return out

def _inline_scripts(path):
    """Leaves from an inline `SCRIPTS=( ... )` block (bucketer-style arrays)."""
    try: txt = Path(path).read_text()
    except OSError: return None
    m = re.search(r"^SCRIPTS=\(\s*$(.*?)^\)\s*$", txt, re.M | re.S)
    if not m: return None
    out = [l.strip() for l in m.group(1).splitlines() if l.strip()]
    return out or None


def compact(idxs):
    """[0,1,2,5,7,8] -> '0-2,5,7-8' (keeps the --array spec short)."""
    out, i = [], 0
    while i < len(idxs):
        j = i
        while j + 1 < len(idxs) and idxs[j + 1] == idxs[j] + 1: j += 1
        out.append(str(idxs[i]) if j == i else f"{idxs[i]}-{idxs[j]}")
        i = j + 1
    return ",".join(out)


groups = defaultdict(list); victims = []; index_jobs = []
for (name, part, qos, tl, cmd, nodelist), tasks in base.items():
    recs = len(tasks)
    idxs = sorted({i for _j, k in tasks for i in expand(k)})
    if recs < a.min_records: continue
    rec_ids = [j for j, _k in tasks]
    man = Path(cmd).parent / "manifest.txt"
    try: leaves = man.read_text().splitlines(); mode = "manifest"
    except OSError:
        leaves = _inline_scripts(cmd); mode = "inline"
        if not leaves:
            # index mode: no leaf list exists (e.g. apo drives off chunk_${i}.tsv). Resubmit
            # the SAME script with just the pending indices — they have not run, so nothing
            # is lost and no done-detection is needed.
            spec = compact(idxs)
            print(f"  {name} {part}: records={recs} tasks={len(idxs)} -> INDEX repack  [{Path(cmd).name}]")
            index_jobs.append((name, part, qos, tl, cmd, nodelist, spec, len(idxs)))
            victims += rec_ids
            continue
    keep = []
    for i in idxs:
        if 0 <= i < len(leaves):
            lp = leaves[i].strip()
            if lp and not (Path(lp).parent / "designs.tar.zst").exists(): keep.append(lp)
    print(f"  {name} {part}: records={recs} tasks={len(idxs)} -> repack {len(keep)} ({mode})  [{Path(cmd).parent.name}]")
    if keep:
        groups[(name, part, qos, tl, cmd, nodelist)] += keep
        victims += rec_ids

if not groups and not index_jobs:
    print("nothing to pack"); raise SystemExit(0)
print(f"\nWould cancel {len(victims)} pending records -> submit "
      f"{len(groups) + len(index_jobs)} compact arrays")
if not a.apply:
    print("DRY-RUN — add --apply"); raise SystemExit(0)

for c in range(0, len(victims), 300): sh("scancel", *victims[c:c+300])
import time; time.sleep(5)
ts = datetime.now().strftime("%Y%m%d_%H%M%S")
for (name, part, qos, tl, cmd, nodelist, spec, n) in index_jobs:
    args = ["sbatch", "--parsable", f"--array={spec}", f"--partition={part}",
            f"--qos={qos}", f"--time={tl}"]
    if nodelist: args.append(f"--nodelist={nodelist}")
    args.append(cmd)
    r = sh(*args)
    jid = r.stdout.strip().split(";")[0]
    if not jid.isdigit():
        print(f"  {name} {part}: SUBMIT_ERR {r.stderr.strip()[:120]}"); continue
    mv = sh("scontrol", "update", f"jobid={jid}", f"partition={part}", f"qos={qos}")
    ok = "moved" if mv.returncode == 0 else f"MOVE_FAIL {mv.stderr.strip()[:60]}"
    print(f"  {name} {part}: {n} tasks (index) -> jobid {jid} ({ok})")

for (name, part, qos, tl, cmd, nodelist), leaves in groups.items():
    seen=set(); leaves=[l for l in leaves if not (l in seen or seen.add(l))]
    old = Path(cmd).parent
    src_sbatch = old / "array.sbatch"
    if not src_sbatch.exists(): src_sbatch = Path(cmd)   # inline-SCRIPTS arrays
    new = old.parent / f"packed_{re.sub(r'[^A-Za-z0-9]','_',part)}_{ts}"
    (new/"array_logs").mkdir(parents=True, exist_ok=True)
    (new/"manifest.txt").write_text("".join(l+"\n" for l in leaves))
    hdr=[]
    for ln in src_sbatch.read_text().splitlines():
        if not ln.startswith("#SBATCH"): continue
        if ln.startswith("#SBATCH --array"): ln=f"#SBATCH --array=0-{len(leaves)-1}%200"
        elif ln.startswith("#SBATCH --output"): ln=f"#SBATCH --output={new}/array_logs/%A_%a.log"
        elif ln.startswith("#SBATCH --error"): ln=f"#SBATCH --error={new}/array_logs/%A_%a.err"
        hdr.append(ln)
    body=('set -euo pipefail\nSCRIPTS=(\n'+"\n".join(f"  {l}" for l in leaves)+'\n)\n'
          'S="${SCRIPTS[$SLURM_ARRAY_TASK_ID]}"\nRUNDIR=$(dirname "$S")\n'
          'echo "[packed] task=$SLURM_ARRAY_TASK_ID host=$(hostname) -> $S"\n'
          'bash "$S" > "$RUNDIR/slurm.log" 2> "$RUNDIR/slurm.err"\n')
    sb = new/"array.sbatch"; sb.write_text("#!/bin/bash\n"+"\n".join(hdr)+"\n\n"+body); sb.chmod(0o755)
    r = sh("sbatch","--parsable",str(sb))
    jid = r.stdout.strip().split(";")[0]
    if not jid.isdigit():
        print(f"  {name} {part}: SUBMIT_ERR {r.stderr.strip()[:100]}"); continue
    # Euler remaps *.bulk -> *.24h at submit; move the BASE id back to the intended tier
    # (base id only, never a task element, so the pack never peels).
    mv = sh("scontrol","update",f"jobid={jid}",f"partition={part}",f"qos={qos}")
    ok = "moved" if mv.returncode==0 else f"MOVE_FAIL {mv.stderr.strip()[:60]}"
    print(f"  {name} {part}: {len(leaves)} tasks -> jobid {jid} ({ok})")
