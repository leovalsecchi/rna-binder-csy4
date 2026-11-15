#!/usr/bin/env python3
"""Re-bake bulk-biased tiers into midrun_v8 HD run.sbatch (#SBATCH --partition/--time only;
nodelist/gres stay pool-fixed; bucketer derives --qos from partition). Then bucket → 6 COMPACT
arrays (bulk lanes route to .24h at submit; move whole-array to bulk by base-id afterward).
Bulk-biased: gpuhe bulk1152/24h480/120h288, gpupr bulk288/24h120/120h72."""
import glob, re, sys
from collections import Counter

HD = "/cluster/work/wicky/lvalsecchi/rna_binder/campaigns/csy4_crrna_free/halludesign__ligandmpnn"
files = sorted(glob.glob("%s/midrun_v8_*/run.sbatch" % HD))
gpuhe, gpupr = [], []
for f in files:
    t = open(f).read()
    if "eu-g4-" in t or "eu-g6-" in t: gpuhe.append(f)
    elif "eu-a65" in t or "eu-ts" in t: gpupr.append(f)
    else: print("UNCLASSIFIED:", f)
print("found: gpuhe=%d gpupr=%d total=%d" % (len(gpuhe), len(gpupr), len(files)))

# tier plans (part, time, count) — must sum to pool size
PLAN = {
    "gpuhe": [("gpuhe.bulk", "23:30:00", 1152), ("gpuhe.24h", "23:30:00", 480), ("gpuhe.120h", "47:30:00", 288)],
    "gpupr": [("gpupr.bulk", "23:30:00", 288),  ("gpupr.24h", "23:30:00", 120), ("gpupr.120h", "47:30:00", 72)],
}
def check(name, lst, plan):
    s = sum(n for _,_,n in plan)
    if s != len(lst):
        print("!! %s plan sums to %d but pool has %d — ADJUSTING last tier" % (name, s, len(lst)))
        p,tm,n = plan[-1]; plan[-1] = (p, tm, n + (len(lst)-s))
    return plan
PLAN["gpuhe"] = check("gpuhe", gpuhe, PLAN["gpuhe"])
PLAN["gpupr"] = check("gpupr", gpupr, PLAN["gpupr"])

assign = {}
for pool, lst in (("gpuhe", gpuhe), ("gpupr", gpupr)):
    i = 0
    for part, tm, n in PLAN[pool]:
        for f in lst[i:i+n]: assign[f] = (part, tm)
        i += n

n_ok = 0
for f, (part, tm) in assign.items():
    t = open(f).read()
    t2 = re.sub(r"^#SBATCH --partition=.*$", "#SBATCH --partition=%s" % part, t, count=1, flags=re.M)
    t2 = re.sub(r"^#SBATCH --time=.*$", "#SBATCH --time=%s" % tm, t2, count=1, flags=re.M)
    if ("--partition=%s" % part) not in t2 or ("--time=%s" % tm) not in t2:
        print("FAILED rewrite:", f); continue
    open(f, "w").write(t2); n_ok += 1
print("rebaked=%d" % n_ok)
print("tier counts:", dict(Counter(v[0] for v in assign.values())))
sys.exit(0 if n_ok == len(files) else 1)
