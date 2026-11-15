#!/usr/bin/env python3
"""verify_filter_counts.py — independent pass-count check for the binder-filter tool.

Deliberately shares NO code with store.py / filters.py: it streams the master CSV row by row
with stdlib csv and plain Python floats, and re-implements the filter semantics from the spec.
If this and the tool disagree, one of them has a wrong column mapping or a wrong NaN rule —
which is exactly the class of bug a test built on the same code path cannot catch.

    python3 verify_filter_counts.py --csv <master> --spec '<json>'

The spec is the same JSON the tool POSTs to /api/filter. Rank-mode ipTM gates need the whole
column in memory to rank, so they are computed in a second pass over the collected values.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
import sys

csv.field_size_limit(10 ** 7)

# tool metric -> master CSV column. Written out again on purpose: a typo here shows up as a
# disagreement, which is the point.
COL = {
    ("rmsd", "af3"): "dockq_LRMSD_af3", ("rmsd", "esm"): "dockq_LRMSD_esm",
    ("dockq", "af3"): "dockq_self_af3", ("dockq", "esm"): "dockq_self_esm",
    ("iptm", "af3"): "af3_msaoff_iptm", ("iptm", "esm"): "esm_iptm",
    ("hb", "af3"): "n_h_bonds_interface_af3", ("hb", "esm"): "n_h_bonds_interface_esm",
    ("pis", "af3"): "pi_stack_score_af3", ("pis", "esm"): "pi_stack_score_esm",
}


def num(s):
    if s is None or s == "":
        return None
    try:
        v = float(s)
    except ValueError:
        return None
    return None if math.isnan(v) else v


def cmp_ok(v, op, thr):
    if v is None:
        return False                       # NaN fails closed
    return {"<=": v <= thr, "<": v < thr, ">=": v >= thr, ">": v > thr}[op]


def tier_ok(row, tier):
    atoms = [a for a in tier.get("atoms", []) if a.get("on", True)]
    if not atoms:
        return True                        # empty tier does not filter
    res = []
    for a in atoms:
        key = (a["metric"], a.get("backend"))
        if key not in COL:
            print(f"  (skipping unsupported atom {key} — extend COL to check it)",
                  file=sys.stderr)
            continue
        res.append(cmp_ok(num(row.get(COL[key])), a.get("cmp", "<="), float(a["value"])))
    if not res:
        return True
    return all(res) if tier.get("op", "AND").upper() == "AND" else any(res)


def iptm_value(row, which):
    a, e = num(row.get(COL[("iptm", "af3")])), num(row.get(COL[("iptm", "esm")]))
    if which == "af3":
        return a
    if which == "esm":
        return e
    if a is None or e is None:
        return None
    return min(a, e) if which == "min" else (a + e) / 2.0


ap = argparse.ArgumentParser()
ap.add_argument("--csv", required=True)
ap.add_argument("--spec", required=True)
ap.add_argument("--limit", type=int, default=0)
a = ap.parse_args()
spec = json.loads(a.spec)
gate = spec.get("iptm") or {}
tiers = spec.get("tiers", {})
combine = spec.get("combine", "AND").upper()

n = n_base = 0
base_iptm = []          # ipTM of base-passing rows, for a rank gate
all_iptm = []           # every row's ipTM, for a global rank gate
with open(a.csv, newline="") as f:
    for row in csv.DictReader(f):
        n += 1
        t1 = tier_ok(row, tiers.get("t1", {}))
        t3 = tier_ok(row, tiers.get("t3", {}))
        ok = (t1 and t3) if combine == "AND" else (t1 or t3)
        v = iptm_value(row, gate.get("value", "min")) if gate.get("on") else None
        if gate.get("on") and gate.get("mode", "threshold") == "rank":
            if gate.get("scope", "filtered") == "global" and v is not None:
                all_iptm.append(v)
        if ok:
            n_base += 1
            if gate.get("on"):
                if gate.get("mode", "threshold") == "threshold":
                    if not cmp_ok(v, ">=", float(gate.get("threshold", 0.8))):
                        n_base -= 1
                else:
                    base_iptm.append(v)
        if a.limit and n >= a.limit:
            break

if gate.get("on") and gate.get("mode", "threshold") == "rank":
    vals = [x for x in base_iptm if x is not None]
    scope = gate.get("scope", "filtered")
    pool = len(vals) if scope == "filtered" else len(all_iptm)
    if gate.get("rank_kind", "n") == "pct":
        k = int(math.ceil(pool * float(gate.get("rank_pct", 20.0)) / 100.0))
    else:
        k = int(gate.get("rank_n", 10000))
    if scope == "filtered":
        n_pass = min(k, len(vals))
    else:
        cut = sorted(all_iptm, reverse=True)[:k]
        thr = cut[-1] if cut else float("inf")
        n_pass = sum(1 for x in vals if x >= thr)
        print(f"  (global rank cut at ipTM >= {thr:.6g}; ties may make this a "
              f"slight over-count vs the tool's exact top-k)", file=sys.stderr)
else:
    n_pass = n_base

print(json.dumps({"n_rows": n, "n_pass": n_pass}))
