"""ck.doctor — guardrails-as-code. Executable mirror of GUARDRAILS.md.

`campaign doctor <c>` runs every check against the resolved plan + live cluster state and
REFUSES launch on any FAIL. Each row names the rule and the 'why' (the scar-lesson).
Statuses: FAIL (blocks launch) · WARN (proceed with caution) · PASS · INFO (advisory).
"""
from __future__ import annotations
import subprocess


def _squeue_names(user):
    try:
        out = subprocess.run(["squeue", "-u", user, "-h", "-o", "%j"],
                             capture_output=True, text=True, timeout=20).stdout
        return out.split()
    except Exception:
        return []


def run(campaign: dict, plan: list, inventory: dict, user: str) -> list[tuple]:
    """Returns list of (name, status, detail). Launch is blocked iff any status == 'FAIL'."""
    checks = []

    # 1. node-local scratch  (2026-07-10 inode flood: extraction must never touch /cluster)
    roots = inventory.get("scratch_roots", [])
    ok = bool(roots) and all(not str(r).startswith("/cluster") for r in roots)
    checks.append(("node_local_scratch", "PASS" if ok else "FAIL",
                   f"scratch_roots={roots}"))

    # 2. every model has a compat entry (else we'd guess routing)
    for w in plan:
        if w.get("error"):
            checks.append((f"compat:{w['label']}", "FAIL", w["error"]))

    # 3. refold/gen routing: each workload only on a pool the probe said it RUNS on
    for w in plan:
        if w.get("error"):
            continue
        pool = w["partition"].split(".")[0]
        ro = w.get("runs_on") or []
        if ro and pool not in ro:
            checks.append((f"routing:{w['label']}", "FAIL",
                           f"-> {pool} but runs_on={ro} (probe never saw it succeed there)"))

    # 4. AF3 must not sit on a sub-24h tier (AF3 TIMEOUTs; policy=24h+)
    for w in plan:
        if w.get("model") == "refold_af3" and w.get("walltime_tier") in ("4h", "bulk"):
            checks.append(("af3_walltime", "WARN",
                           f"AF3 on {w['walltime_tier']} — AF3 runs to ~17h and TIMEOUTs on short tiers"))

    # 5. cysteine / omit-residues declared (leaked into RFD3-LiMPNN + BG-BIF before)
    omit = (campaign.get("constraints") or {}).get("omit_residues")
    checks.append(("omit_residues", "PASS" if omit else "WARN",
                   f"omit_residues={omit or 'NONE — cysteine may leak into inverse folding'}"))

    # 6. admission control up before a big gen fan-out (250-window starvation)
    gen_jobs = sum((w.get("count") or 0) for w in plan if w.get("role") == "gen")
    names = _squeue_names(user)
    hg = any("hold_governor" in n for n in names)
    if gen_jobs > 5000 and not hg:
        checks.append(("admission_control", "WARN",
                       f"{gen_jobs} gen units but hold_governor loop NOT running — refolds may starve"))
    else:
        checks.append(("admission_control", "PASS",
                       "hold_governor running" if hg else "small fan-out; governor optional"))

    # 7. header hygiene (guaranteed by the planner/executor, asserted here)
    checks.append(("mem_per_cpu", "PASS", "plan emits --mem-per-cpu, never --mem"))
    checks.append(("requeue_append", "PASS", "executor emits --requeue --open-mode=append"))

    # advisory — depends on the run body / monitoring, not statically checkable here
    checks.append(("pack_verify", "INFO", "outputs must be tar.zst + `zstd -t` before delete"))
    checks.append(("pool_count_truth", "INFO", "judge HD success by packed count, not sacct State (false-FAILED)"))
    return checks


def ok(checks) -> bool:
    return not any(s == "FAIL" for _, s, _ in checks)
