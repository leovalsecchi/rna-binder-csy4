"""ck.planner — campaign spec + compat.yaml -> submission plan.

Generalizes the two hand-assembled pieces:
  * ops/submit/submit_mv<N>_gen_*.sh  (the retyped per-stack routing table)
  * ops/submit/bucket_midrun_v<N>.py  (arrays by routing signature)
Every routing decision comes from the model's MEASURED compat entry — no hand-typed tuples.
"""
from __future__ import annotations
import math

# Pool preference: cheapest / highest-cap first. The public Turing pool `gpu` is preferred
# whenever a model's runs_on includes it (RFD3/BoltzGen/BoltzDesign/B2); Ampere `gpuhe` or
# A100 `gpupr` only when the model needs them (HalluDesign is sm_80-only, so its runs_on
# excludes `gpu` and it lands on gpuhe/gpupr automatically).
POOL_PREF = ["gpu", "gpuhe", "gpupr"]


def choose_pool(entry: dict, prefer: str | None = None) -> str | None:
    ro = entry.get("runs_on") or []
    if prefer and prefer in ro:
        return prefer
    for p in POOL_PREF:
        if p in ro:
            return p
    return ro[0] if ro else None


def mem_per_cpu_g(entry: dict, cpus: int = 4) -> str:
    hm = entry.get("host_mem_peak_g")
    if not hm:
        return "8G"
    return f"{max(2, math.ceil((hm * 1.2) / cpus))}G"


def plan_workload(label, role, model, count, compat, pools_by_name, prefer=None, cpus=4):
    e = compat.get(model)
    if not e:
        return {"label": label, "role": role, "model": model,
                "error": f"no compat entry — run: campaign probe --model {model} --write"}
    pool = choose_pool(e, prefer)
    tier = e.get("walltime_tier") or "24h"
    tiern = "bulk" if tier == "bulk" else tier.rstrip("h")
    qos_base = pools_by_name.get(pool, {}).get("qos_base", f"es_wicky/{pool}")
    upj = e.get("units_per_job") or 1
    return {
        "label": label, "role": role, "model": model,
        "partition": f"{pool}.{tier}",
        "walltime_tier": tier,
        "qos": f"{qos_base}/{tiern}",
        "gres": f"gpumem:{e['req_gpumem_g']}g" if e.get("req_gpumem_g") else "",
        "mem_per_cpu": mem_per_cpu_g(e, cpus),
        "cpus_per_task": cpus,
        "count": count,
        "units_per_job": upj,
        "n_jobs": math.ceil(count / upj) if count else None,
        "runs_on": e.get("runs_on"),
        "memory_elastic": e.get("memory_elastic"),
    }


def plan_campaign(campaign: dict, compat: dict, inventory: dict) -> list[dict]:
    pools_by_name = {p["name"]: p for p in inventory["pools"]}
    plan = []
    for st in campaign.get("stacks", []):
        gen, inv = st.get("gen"), st.get("inverse_fold")
        label = f"gen:{gen}/{inv}" if inv else f"gen:{gen}"
        plan.append(plan_workload(label, "gen", gen, st.get("count"), compat, pools_by_name))
    total = sum((s.get("count") or 0) for s in campaign.get("stacks", []))
    for b in (campaign.get("refold") or {}).get("backends", []):
        plan.append(plan_workload(f"refold:{b}", "refold", f"refold_{b}", total, compat, pools_by_name))
    return plan
