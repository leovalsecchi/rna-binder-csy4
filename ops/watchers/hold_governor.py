#!/usr/bin/env python3
"""
hold_governor.py -- per-partition admission-control governor for the SLURM queue.

WHY
  Euler's main scheduler evaluates only the top `default_queue_depth` (=250)
  pending jobs by priority each cycle; backfill budget is per-user-PER-PARTITION
  (`bf_max_job_user_part`=50). A big OWN backlog of generation jobs piled on ONE
  partition both floods the global window AND blows past that partition's backfill
  budget, starving higher-value jobs. Proven 2026-06-29: holding ~543 pending gen
  let 15 gpu-tier AF3 refolds start the instant the hold finished.

WHAT (mechanism #2 -- per-partition admission control + spread)
  Keep ~target/N un-held gen on EACH HD-gen-capable lane (default 200 -> 50/lane).
    - per lane: too many active -> HOLD newest excess; room -> RELEASE oldest held
  Spreading the active set evenly across lanes means (a) no single lane hoards the
  window, (b) each lane fills its own QOS cap in parallel, and (c) each lane stays
  within its 50-job backfill budget. Refolds live on gpu.* lanes (NOT in the spread
  set) so they are never touched and get their own per-partition backfill budget.

SAFE
  - Dry-run by default; only acts with --apply.
  - Only touches THIS user's gen jobs (GEN_RE) on the SPREAD lanes. HD gen is
    Ampere/Ada/a100-only & ~3-7h, so the lane set excludes gpuhe.4h (timeout) and
    gpu.* (Turing sm_75). Refolds/loops/other partitions untouched.
  - Idempotent. <= --max-ops scontrol calls per pass (Euler rate-limits ~0.5-1s/
    call). STOP sentinel honored. Logs to _hold_governor.log.

ARRAY-SAFETY / base-id discipline. This governor is deliberately array-safe: pending()
keeps only pure-integer JobIDs (jid.isdigit()), so it NEVER touches array elements, and
GEN_RE matches only standalone hd_<uuid8> jobs -- not the mvN_hd campaign ARRAYS. It
therefore cannot peel array tasks (unlike a per-task scontrol op). If ever extended to
manage arrays, hold/release by BASE array id only; peeled tasks are healed by
ops/loops/hd_consolidate.py.
"""
import argparse
import re
import subprocess
from pathlib import Path

USER = "lvalsecchi"
RB = Path("/cluster/work/wicky/lvalsecchi/rna_binder/ops/watchers")
STOP = RB / "_hold_governor.STOP"
LOG = RB / "_hold_governor.log"

GEN_RE = re.compile(r"^hd_[0-9a-f]{8}$")  # managed (holdable) jobs
# HD-gen-capable lanes only: Ampere/Ada/a100, >=24h walltime.
SPREAD_PARTITIONS = ["gpuhe.24h", "gpuhe.120h", "gpupr.24h", "gpupr.120h"]


def sh(*a):
    return subprocess.run(a, capture_output=True, text=True)


def pending():
    out = sh("squeue", "-u", USER, "-t", "PD", "-h", "-o", "%i|%j|%P|%r").stdout
    rows = []
    for ln in out.splitlines():
        p = ln.split("|")
        if len(p) != 4:
            continue
        jid, name, part, reason = (x.strip() for x in p)
        if jid.isdigit():
            rows.append((int(jid), name, part, reason.startswith("JobHeld")))
    return rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--target", type=int, default=200,
                    help="total un-held gen jobs to keep, split evenly across lanes")
    ap.add_argument("--max-ops", type=int, default=50,
                    help="max scontrol hold/release calls per pass (RPC-throttle guard)")
    ap.add_argument("--apply", action="store_true", help="act for real (default: dry-run)")
    args = ap.parse_args()

    if STOP.exists():
        print("[hold_governor] STOP sentinel present -> no-op")
        return

    rows = pending()
    n = len(SPREAD_PARTITIONS)
    base, rem = divmod(args.target, n)
    per = {p: base + (1 if i < rem else 0) for i, p in enumerate(SPREAD_PARTITIONS)}

    stamp = sh("date", "+%Y-%m-%dT%H:%M:%S").stdout.strip()
    ops_left = args.max_ops
    plan, loglines = [], []

    for p in SPREAD_PARTITIONS:
        active = sorted(j for (j, nm, pp, h) in rows if pp == p and GEN_RE.match(nm) and not h)
        held = sorted(j for (j, nm, pp, h) in rows if pp == p and GEN_RE.match(nm) and h)
        tgt = per[p]
        act, ids = "noop", []
        if len(active) > tgt:
            act = "hold"
            ids = sorted(active, reverse=True)[:len(active) - tgt][:max(0, ops_left)]
        elif len(active) < tgt and held:
            act = "release"
            ids = sorted(held)[:min(tgt - len(active), len(held))][:max(0, ops_left)]
        plan.append((p, len(active), tgt, len(held), act, len(ids)))
        if ids and ops_left > 0:
            ops_left -= len(ids)
            if args.apply:
                res = sh("scontrol", act, ",".join(map(str, ids)))
                loglines.append(f"[{stamp}] {p} {act} {len(ids)} ({ids[0]}..{ids[-1]}) rc={res.returncode}")
                if res.returncode != 0:
                    print("  err", p, res.stderr.strip()[:150])

    print(f"[{stamp}] target={args.target} per-lane~{base}" + ("" if args.apply else "  [DRY-RUN]"))
    for (p, a_n, tgt, h_n, act, nid) in plan:
        print(f"  {p:12s} active={a_n:3d} tgt={tgt:3d} held={h_n:3d} -> {act} {nid}")

    if args.apply and loglines:
        with open(LOG, "a") as fh:
            fh.write("\n".join(loglines) + "\n")


if __name__ == "__main__":
    main()
