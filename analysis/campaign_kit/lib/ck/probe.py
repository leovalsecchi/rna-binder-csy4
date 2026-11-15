"""ck.probe — the compatibility probe.

Two ways to measure a model:
  * HISTORY-HARVEST (implemented): read peak VRAM / GPU-util / runtime / which-pool back
    out of sacct for the model's ALREADY-COMPLETED jobs. No new submissions. This is how
    an existing model earns its compat.yaml entry, and it self-corrects legacy myths by
    measuring (AF3 really runs on Turing; HD peak ~23G not 79.5G).
  * CALIBRATION-SUBMIT (stubbed below): for a NEW model with no history — submit a tiny
    job per candidate pool and harvest that. Same measurement engine; needs a launch.py hook.

The one primitive nothing in the repo had: reading gres/gpumem + gres/gpuutil out of
sacct's TRESUsageInMax (documented in 3 places, coded in none).
"""
from __future__ import annotations
import fnmatch
import re
import statistics as st
import subprocess

STD_GRES = [8, 11, 16, 24, 40, 48, 80]            # gres=gpumem:<g> ladder to round a measured peak up to
TIER_SECONDS = {"4h": 14400, "24h": 86400, "120h": 432000}


def _to_gib(v):
    """'23.4G' / '24000M' / '25165824K' / bare-bytes -> GiB float (None if unparseable)."""
    v = str(v).strip()
    mult = {"K": 1 / 1024 / 1024, "M": 1 / 1024, "G": 1.0, "T": 1024.0}
    if v and v[-1] in mult:
        try:
            return round(float(v[:-1]) * mult[v[-1]], 2)
        except ValueError:
            return None
    try:
        return round(float(v) / 1024 / 1024 / 1024, 2)
    except ValueError:
        return None


def _parse_tres(tres: str) -> dict:
    """Pull gres/gpumem (peak VRAM), gres/gpuutil, mem (host RAM) out of a TRESUsageInMax value."""
    out = {"gpumem_g": None, "gpuutil": None, "mem_g": None}
    for kv in tres.split(","):
        if "=" not in kv:
            continue
        k, v = kv.split("=", 1)
        if k == "gres/gpumem":
            out["gpumem_g"] = _to_gib(v)
        elif k == "gres/gpuutil":
            v = v.rstrip("%")
            out["gpuutil"] = int(v) if v.isdigit() else None
        elif k == "mem":
            out["mem_g"] = _to_gib(v)
    return out


def read_peak_vram_util(jobid: str) -> dict:
    """Single-job peak VRAM/util (used by the calibration-submit path). Reads the .batch step."""
    out = subprocess.run(
        ["sacct", "-j", jobid, "-p", "--noheader",
         "--format=JobID,State,ElapsedRaw,TRESUsageInMax%256"],
        capture_output=True, text=True).stdout
    res = {"jobid": jobid, "state": None, "elapsed_s": None,
           "gpumem_g": None, "gpuutil": None, "mem_g": None}
    for line in out.splitlines():
        f = line.split("|")
        if len(f) < 4:
            continue
        jid, state, elapsed, tres = f[0], f[1], f[2], f[3]
        if "." not in jid:
            res["state"] = (state.split() or [state])[0]
            res["elapsed_s"] = int(elapsed) if elapsed.strip().isdigit() else None
        elif jid.endswith(".batch") and tres:
            res.update(_parse_tres(tres))
    return res


def harvest(user: str, days: int = 14) -> list[dict]:
    """Every finished job in the window, joining the allocation row (name/state/elapsed/node)
    with its .batch step's TRESUsageInMax. One sacct call. Returns list of job dicts."""
    out = subprocess.run(
        ["sacct", "-u", user, "-S", f"now-{days}days", "-E", "now", "-P", "--noheader",
         "-o", "JobID,JobName,State,ElapsedRaw,NodeList,TRESUsageInMax%256"],
        capture_output=True, text=True).stdout
    alloc, tres = {}, {}
    for line in out.splitlines():
        f = line.split("|")
        if len(f) < 6:
            continue
        jid, name, state, elapsed, node, t = f[:6]
        base = jid.split(".")[0]
        if "." not in jid:
            alloc[base] = {"jobid": base, "name": name,
                           "state": (state.split() or [state])[0] if state else state,
                           "elapsed_s": int(elapsed) if elapsed.strip().isdigit() else None,
                           "node": node}
        elif jid.endswith(".batch") and t:
            tres[base] = _parse_tres(t)
    rows = []
    for base, a in alloc.items():
        a.update(tres.get(base, {}))
        rows.append(a)
    return rows


def pool_of_node(node: str, pools: list[dict]):
    """Which pool a NodeList ran on, by matching the first node token to each pool's globs."""
    if not node or node.lower().startswith("none"):
        return None
    n = re.split(r"[,\[]", node)[0]
    for p in pools:
        for g in p.get("node_globs", []):
            if fnmatch.fnmatch(n, g):
                return p["name"]
    return None


def _p95(vals):
    vals = sorted(v for v in vals if v is not None)
    if not vals:
        return None
    return vals[min(len(vals) - 1, int(0.95 * len(vals)))]


def summarize_by_pool(rows: list[dict], pools: list[dict]) -> dict:
    """Per-pool {done, fail, vram_p95, med_runtime_s} — for the validation printout."""
    by = {}
    for r in rows:
        pl = pool_of_node(r.get("node"), pools)
        if not pl:
            continue
        b = by.setdefault(pl, {"done": 0, "fail": 0, "_vram": [], "_rt": []})
        if r.get("state") == "COMPLETED":
            b["done"] += 1
            if r.get("gpumem_g"):
                b["_vram"].append(r["gpumem_g"])
            if r.get("elapsed_s"):
                b["_rt"].append(r["elapsed_s"])
        else:
            b["fail"] += 1
    for b in by.values():
        b["vram_p95"] = _p95(b.pop("_vram"))
        rt = b.pop("_rt")
        b["med_runtime_s"] = int(st.median(rt)) if rt else None
    return by


GRAB_FRAC = 0.9   # if a job's peak VRAM >= 90% of the card, the model grabs the whole card
                  # (JAX/Torch preallocation) — that reading is a card-size artifact, not a true need.


def _vram_floor(by: dict, pool_vram: dict):
    """True VRAM requirement, corrected for whole-card-grabbing models.

    A model that grabs the card reports VRAM ~= card size (e.g. HalluDesign: 22G on a 24G
    card but 78G on an 80G A100). The real floor is 'the smallest card it demonstrably fit in'.
      - If it ran on any pool WITHOUT grabbing -> use the max non-grab peak (its true need).
      - If it grabbed on every pool -> it fit in the smallest of those cards -> that card size.
    Returns (floor_g, memory_elastic_bool, vram_by_pool).
    """
    by_pool, non_grab, grab_cards = {}, [], []
    for pl, b in by.items():
        v, card = b.get("vram_p95"), pool_vram.get(pl)
        if b["done"] == 0 or v is None:
            continue
        grabbing = bool(card and v >= GRAB_FRAC * card)
        by_pool[pl] = {"vram_p95_g": v, "card_g": card, "grabbing": grabbing}
        (grab_cards.append(card) if grabbing else non_grab.append(v))
    if non_grab:
        return round(max(non_grab), 2), bool(grab_cards), by_pool
    if grab_cards:
        return float(min(grab_cards)), True, by_pool           # fit in the smallest grabbed card
    return None, False, by_pool


def _req_gpumem(floor, elastic):
    """gres=gpumem:<g> to request. If the model is elastic it demonstrably FIT in `floor`
    (a card size) -> request exactly that; otherwise round the true need up a rung + headroom."""
    if floor is None:
        return None
    if elastic and floor in STD_GRES:
        return int(floor)
    return next((g for g in STD_GRES if g >= floor + 1), 80)


def derive(model: dict, rows: list[dict], pools: list[dict], caps_by_pool: dict):
    """rows = jobs already filtered to THIS model. -> a compat.yaml entry (or None if no data)."""
    done = [r for r in rows if r.get("state") == "COMPLETED"]
    if not done:
        return None
    by = summarize_by_pool(rows, pools)
    runs_on = sorted(pl for pl, b in by.items() if b["done"] > 0)
    pool_vram = {p["name"]: p.get("vram_g") for p in pools}
    util = [r["gpuutil"] for r in done if r.get("gpuutil") is not None]
    memg = [r["mem_g"] for r in done if r.get("mem_g")]
    elapsed = [r["elapsed_s"] for r in done if r.get("elapsed_s")]
    units = (model.get("io", {}) or {}).get("units_per_job", 1) or 1

    floor, elastic, vram_by_pool = _vram_floor(by, pool_vram)
    req_gpumem = _req_gpumem(floor, elastic)
    med_job = int(st.median(elapsed)) if elapsed else None
    p95_job = _p95(elapsed)                                   # tier off the TAIL, not the median (AF3 runs to ~17h)
    tier = "bulk"
    if p95_job:
        need = p95_job * 1.2
        for t in ("4h", "24h", "120h"):
            if TIER_SECONDS[t] * 0.95 >= need:
                tier = t
                break
    return {
        "runs_on": runs_on,
        "n_jobs_measured": len(done),
        "vram_floor_g": floor,
        "req_gpumem_g": req_gpumem,
        "memory_elastic": elastic,                            # grabs the whole card (VRAM != true need)
        "vram_by_pool": vram_by_pool,
        "gpuutil_pct": int(st.median(util)) if util else None,
        "host_mem_peak_g": _p95(memg),
        "median_job_s": med_job,
        "p95_job_s": p95_job,
        "units_per_job": units,
        "sec_per_unit": round(med_job / units, 1) if med_job else None,   # approx (refolds may batch)
        "walltime_tier": tier,
        "qos_concurrency": {pl: caps_by_pool.get(pl, {}).get(tier) for pl in runs_on},
    }


# ---- calibration-submit path: probe a model with NO history from cold ----
import time as _time                                   # normal python (Euler login), not a workflow

_TERMINAL = {"COMPLETED", "FAILED", "TIMEOUT", "CANCELLED", "OUT_OF_MEMORY",
             "NODE_FAIL", "BOOT_FAIL", "PREEMPTED", "DEADLINE"}


def _arch_int(s):
    try:
        return int(str(s).replace("sm_", ""))
    except (TypeError, ValueError):
        return 0


def candidate_pools(model, inventory):
    """Pools whose arch >= model.compat_hints.min_arch (all pools if the hint is absent —
    an unknown model is tried everywhere and fails-and-learns)."""
    need = _arch_int((model.get("compat_hints") or {}).get("min_arch"))
    return [p for p in inventory["pools"] if _arch_int(p.get("arch")) >= need]


def poll_until_done(jobid, timeout_s=3600, interval_s=30, sleep=None):
    """Block until a job reaches a terminal state (or timeout), then return its measurement."""
    sleep = sleep or _time.sleep
    waited = 0
    while waited <= timeout_s:
        r = read_peak_vram_util(jobid)
        st = (r.get("state") or "").split()[0] if r.get("state") else ""
        if st in _TERMINAL:
            return r
        sleep(interval_s)
        waited += interval_s
    return {"jobid": jobid, "state": "POLL_TIMEOUT"}


def derive_from_calibration(model, measurements, inventory, caps_by_pool, calib_units):
    """measurements: [{pool,state,gpumem_g,gpuutil,elapsed_s}]. -> a compat entry (or None)."""
    pool_vram = {p["name"]: p.get("vram_g") for p in inventory["pools"]}
    by = {}
    for m in measurements:
        b = by.setdefault(m["pool"], {"done": 0, "fail": 0, "vram_p95": None, "med_runtime_s": None})
        if m.get("state") == "COMPLETED":
            b["done"] = 1
            b["vram_p95"] = m.get("gpumem_g")
            b["med_runtime_s"] = m.get("elapsed_s")
        else:
            b["fail"] = 1
    done = [m for m in measurements if m.get("state") == "COMPLETED"]
    if not done:
        return None
    runs_on = sorted(pl for pl, b in by.items() if b["done"] > 0)
    floor, elastic, vram_by_pool = _vram_floor(by, pool_vram)
    req = _req_gpumem(floor, elastic)
    units_pj = (model.get("io", {}) or {}).get("units_per_job", 1) or 1
    per_unit = [m["elapsed_s"] / max(1, calib_units) for m in done if m.get("elapsed_s")]
    spu = st.median(per_unit) if per_unit else None
    job_s = spu * units_pj if spu else None
    tier = "bulk"
    if job_s:
        need = job_s * 1.2
        for t in ("4h", "24h", "120h"):
            if TIER_SECONDS[t] * 0.95 >= need:
                tier = t
                break
    return {
        "runs_on": runs_on,
        "n_jobs_measured": len(done),
        "vram_floor_g": floor,
        "req_gpumem_g": req,
        "memory_elastic": elastic,
        "vram_by_pool": vram_by_pool,
        "gpuutil_pct": int(st.median([m["gpuutil"] for m in done if m.get("gpuutil") is not None]))
        if any(m.get("gpuutil") is not None for m in done) else None,
        "sec_per_unit": round(spu, 1) if spu else None,
        "units_per_job": units_pj,
        "walltime_tier": tier,
        "qos_concurrency": {pl: caps_by_pool.get(pl, {}).get(tier) for pl in runs_on},
        "source": f"calibration ({calib_units}u/pool)",
    }
