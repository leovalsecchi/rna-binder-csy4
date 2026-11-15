#!/usr/bin/env python3
"""add_topup_store_cols.py — fill the 60 columns build_filter_store.py does not emit into the
v8-HD-topup mini store, so concat_stores.py produces a store with the SAME 85 columns as the
live one.

This matters more than it looks: concat_stores NaN-fills any column missing from one side, and
the filter treats NaN as FAIL. Skipping this step would leave all 170,558 new HD designs
failing every dsasa / foot / footn / footnb filter — they would be in the store yet invisible
to the gate, which is worse than not adding them at all.

Sources (each verified against the live store before use):
  master v8hdtopup_master.csv : rog, fbase, dsasa*, packstat, sc_value, delta_unsat, foot
  foot_clash/out2/            : footn_<tag>  <- foot_tp_f<tag>   (probe atoms, denom 129)
                                footnb_<tag> <- foot_tp_b<tag>   (binder atoms, denom ~600)
  pairrmsd_chunks_v8hdtopup/  : pose_disagree <- rmsd_pose_af3_esm
                                fold_disagree <- rmsd_prot_af3_esm
b2 never enters: only *_af3 / *_esm are read, and the pairrmsd b2 pair columns are ignored.
"""
import csv, glob, json, os, sys
import numpy as np

csv.field_size_limit(10 ** 7)
CC = "/cluster/work/wicky/lvalsecchi/rna_binder/analyses/cross_campaign/"
STORE = CC + "filter_store_v8hdtopup"
MASTER = CC + "v8hdtopup_master.csv"
FOOT_LO, FOOT_HI = 2138, 2206
TAGS = ("005","010","015","020","025","030","034","040","050","075","100")
BE = ("af3", "esm")

MASTER_MAP = {}
for b in BE:
    MASTER_MAP[f"rog_{b}"]                 = f"radius_of_gyration_ca_{b}"
    MASTER_MAP[f"fbase_{b}"]               = f"frac_base_contacts_{b}"
    MASTER_MAP[f"dsasa_{b}"]               = f"dsasa_total_{b}"
    MASTER_MAP[f"dsasa_polar_{b}"]         = f"dsasa_polar_{b}"
    MASTER_MAP[f"dsasa_hphobic_{b}"]       = f"dsasa_hphobic_{b}"
    MASTER_MAP[f"packstat_{b}"]            = f"packstat_{b}"
    MASTER_MAP[f"sc_value_{b}"]            = f"sc_value_{b}"
    MASTER_MAP[f"delta_unsat_hbonds_{b}"]  = f"delta_unsat_hbonds_{b}"
    MASTER_MAP[f"foot_{b}"]                = f"foot_{b}"

ids = np.load(os.path.join(STORE, "design_id.npy"), allow_pickle=False)
n = ids.shape[0]
pos = {(b.decode() if isinstance(b, bytes) else str(b)): i for i, b in enumerate(ids)}
assert len(pos) == n, "duplicate design_id in mini store"
print(f"mini store: {n:,} rows", flush=True)

def f32(v):
    if v is None or v == "" or v == "nan":
        return np.float32(np.nan)
    try:
        return np.float32(float(v))
    except ValueError:
        return np.float32(np.nan)

out = {k: np.full(n, np.nan, dtype=np.float32) for k in MASTER_MAP}
need = set(MASTER_MAP.values())
with open(MASTER, newline="") as f:
    r = csv.DictReader(f)
    miss = need - set(r.fieldnames)
    if miss:
        sys.exit(f"ABORT: master lacks {sorted(miss)}")
    seen = 0
    for row in r:
        i = pos.get(row["design_id"])
        if i is None:
            continue
        seen += 1
        for k, src in MASTER_MAP.items():
            out[k][i] = f32(row[src])
print(f"master columns filled from {seen:,} rows", flush=True)

# ---- foot counts: per-backend out2 sidecars -----------------------------------------
for b in BE:
    for t in TAGS:
        out[f"footn_{t}_{b}"] = np.full(n, np.nan, dtype=np.float32)
        out[f"footnb_{t}_{b}"] = np.full(n, np.nan, dtype=np.float32)
    got = 0
    for idx in range(FOOT_LO, FOOT_HI + 1):
        p = f"{CC}foot_clash/out2/foot_{b}_{idx:05d}.csv"
        if not os.path.exists(p):
            continue
        with open(p, newline="") as f:
            for row in csv.DictReader(f):
                i = pos.get(row["design_id"])
                if i is None:
                    continue
                got += 1
                for t in TAGS:
                    out[f"footn_{t}_{b}"][i] = f32(row.get(f"foot_tp_f{t}"))
                    out[f"footnb_{t}_{b}"][i] = f32(row.get(f"foot_tp_b{t}"))
    print(f"foot counts {b}: {got:,} rows matched", flush=True)

# ---- pairrmsd ------------------------------------------------------------------------
out["pose_disagree"] = np.full(n, np.nan, dtype=np.float32)
out["fold_disagree"] = np.full(n, np.nan, dtype=np.float32)
got = 0
for p in sorted(glob.glob(CC + "pairrmsd_chunks_v8hdtopup/pairrmsd_*.csv")):
    with open(p, newline="") as f:
        for row in csv.DictReader(f):
            i = pos.get(row["design_id"])
            if i is None:
                continue
            got += 1
            out["pose_disagree"][i] = f32(row.get("rmsd_pose_af3_esm"))
            out["fold_disagree"][i] = f32(row.get("rmsd_prot_af3_esm"))
print(f"pairrmsd: {got:,} rows matched", flush=True)

for k, v in sorted(out.items()):
    np.save(os.path.join(STORE, k + ".npy"), v)
print(f"\nwrote {len(out)} columns", flush=True)
print(f"{'column':28} {'non-null':>10}  {'pct':>7}")
for k in sorted(out):
    nn = int(np.isfinite(out[k]).sum())
    print(f"  {k:26} {nn:10,}  {100*nn/n:6.1f}%")

m = json.load(open(os.path.join(STORE, "meta.json")))
m["float_cols"] = {**m.get("float_cols", {}), **MASTER_MAP}
m["placeholder_cols"] = []
json.dump(m, open(os.path.join(STORE, "meta.json"), "w"), indent=1)
print("\nmeta.json updated; placeholder_cols cleared")
