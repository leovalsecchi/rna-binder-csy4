#!/usr/bin/env python3
"""hd_watcher — detect STUCK HalluDesign gen jobs (RUNNING but making no progress,
e.g. evoformer-JIT / wedged-node / MPS hangs) and requeue them.

The campaign array wrapper runs each per-uuid run.sbatch as
    bash "$S" > "$RUNDIR/slurm.log" 2> "$RUNDIR/slurm.err"
so ALL live progress lands in RUNDIR/slurm.log (startup logs, per-cycle prints)
and RUNDIR/designs/ (a recycle every ~7-11 min). Progress signal = newest mtime
of {RUNDIR/slurm.log, RUNDIR/designs}. No progress for >STALE_MIN while RUNNING
=> stuck. GRACE_MIN skips fresh jobs; MAX_REQUEUE caps retries per task.

RUNDIR is resolved authoritatively from the job's real StdOut (scontrol) whose
first line is '[..] task=N host=.. -> <RUNDIR>/run.sbatch'. (squeue %o gives the
array WRAPPER script, not the per-task log — do not use it.)

Only touches campaign HD gen jobs (names 'mvN_hd' / 'hd_<uuid8>'); hd_ns1_test
and refolds are never matched. DRY-RUN by default; --apply to requeue.

FRAGMENTATION / base-id discipline. `scontrol requeue <arrayjob>_<task>` is per-task and
PEELS that task into its own JobID (Slurm has no "restart one array element in place").
Peeled singletons clog the top-250 scheduler window. This is INTENTIONAL and healed
downstream by the hd_consolidate loop (ops/loops/hd_consolidate.py), which re-folds peels
back into compact arrays each tick. Requeue (not scancel) is deliberate: a hung uuid must
restart even if the consolidate loop is down. Do NOT add per-task hold/update here — those
also peel; operate on the BASE array id.
"""
import subprocess, os, re, time, json, argparse

NAME_RE = re.compile(r"^(mv\d+_hd|hd_[0-9a-f]{8})$")
FIRST_RE = re.compile(r"->\s+(\S+)/run\.sbatch")
STDOUT_RE = re.compile(r"StdOut=(\S+)")

def sh(cmd):
    return subprocess.run(cmd, capture_output=True, text=True).stdout

def parse_elapsed(s):
    d = 0
    if "-" in s:
        d, s = s.split("-", 1); d = int(d)
    p = [int(x) for x in s.split(":")]
    while len(p) < 3: p = [0] + p
    return d*86400 + p[0]*3600 + p[1]*60 + p[2]

def mtime(p):
    try: return os.path.getmtime(p)
    except OSError: return 0.0

def resolve_rundir(jid):
    """(stdout_log, run_dir) from the job's real StdOut first line."""
    m = STDOUT_RE.search(sh(["scontrol", "show", "job", "-o", jid]))
    if not m: return None, None
    so = m.group(1)
    try:
        with open(so) as f:
            first = f.readline()
    except OSError:
        return so, None
    mm = FIRST_RE.search(first)
    return so, (mm.group(1) if mm else None)

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--stale-min", type=float, default=45.0)
    ap.add_argument("--grace-min", type=float, default=20.0)
    ap.add_argument("--max-requeue", type=int, default=3)
    ap.add_argument("--max-per-run", type=int, default=15,
                    help="abort without acting if more than this many look stuck "
                         "(a mass hit almost always means a detection bug, not reality)")
    ap.add_argument("--user", default="lvalsecchi")
    ap.add_argument("--state", default=
        "/cluster/work/wicky/lvalsecchi/rna_binder/ops/watchers/hd_watcher/state.json")
    ap.add_argument("--apply", action="store_true")
    a = ap.parse_args()
    now = time.time()
    try: state = json.load(open(a.state))
    except Exception: state = {}

    rows = sh(["squeue", "-u", a.user, "-h", "-t", "RUNNING",
               "-o", "%i|%j|%M|%N"]).splitlines()
    n_hd = 0; stuck = []; skip_nord = 0
    for line in rows:
        f = line.split("|")
        if len(f) < 4: continue
        jid, name, el, node = f[0], f[1], f[2], f[3]
        if not NAME_RE.match(name): continue
        n_hd += 1
        elapsed = parse_elapsed(el)
        if elapsed < a.grace_min * 60: continue
        so, rd = resolve_rundir(jid)
        if not rd:
            skip_nord += 1; continue
        prog = max(mtime(os.path.join(rd, "slurm.log")),
                   mtime(os.path.join(rd, "designs")))
        if prog == 0:                       # no run-dir output at all yet -> conservative skip
            skip_nord += 1; continue
        idle = (now - prog) / 60
        if idle > a.stale_min:
            stuck.append((jid, name, node, elapsed/3600.0, idle, rd))

    print(f"[hd_watcher] {time.strftime('%F %T')} RUNNING_HD={n_hd} "
          f"stuck={len(stuck)} skip_norundir={skip_nord} "
          f"(stale>{a.stale_min:g}m grace>{a.grace_min:g}m) apply={a.apply}")
    if a.apply and len(stuck) > a.max_per_run:
        print(f"[hd_watcher] ABORT: {len(stuck)} stuck > --max-per-run={a.max_per_run}; "
              f"refusing to act (likely detection bug / FS mtime issue). Investigate.")
        return
    acted = 0
    for jid, name, node, eh, idle, rd in sorted(stuck, key=lambda x: -x[4]):
        n = state.get(jid, {}).get("requeues", 0)
        tag = "REQUEUE" if n < a.max_requeue else "SKIP(maxed)"
        uuid = os.path.basename(rd)
        print(f"  STUCK {jid:>14} {name:8} {node:11} elapsed={eh:4.1f}h "
              f"idle={idle:5.0f}m rq={n} -> {'' if a.apply else 'DRY:'}{tag}  {uuid}")
        if a.apply and n < a.max_requeue:
            r = sh(["scontrol", "requeue", jid]).strip()
            if r: print(f"      scontrol: {r}")
            st = state.setdefault(jid, {})
            st.update(requeues=n+1, last=now, node=node, name=name)
            acted += 1
    if a.apply:
        os.makedirs(os.path.dirname(a.state), exist_ok=True)
        json.dump(state, open(a.state, "w"))
        print(f"[hd_watcher] requeued {acted}; state saved.")

main()
