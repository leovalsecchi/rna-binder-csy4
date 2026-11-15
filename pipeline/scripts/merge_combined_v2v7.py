#!/usr/bin/env python3
"""Merge additional versions' clean rows into the v2-v5 combined_clean master.

The master (combined_clean.csv) has 4 derived columns the raw collect->clean output lacks
(campaign, consensus_iptm, mean_iptm, has_consensus). We derive those for each added
version, align to the master's column set, concat, and dedup by design_id (master/earlier
frames win on any collision -> freezes v2-v5).

  python3 merge_combined_v2v7.py <master.csv> <out.csv> <ver:clean.csv> [<ver:clean.csv> ...]
    e.g.  merge_combined_v2v7.py combined_clean.csv combined_clean_v2v7.csv \
              v6:analyses/midrun_v6/collected_clean.csv v7:analyses/midrun_v7/collected_clean.csv
"""
import sys
import pandas as pd

master_path, out_path = sys.argv[1], sys.argv[2]
adds = sys.argv[3:]
if not adds:
    sys.exit("need at least one <ver:clean.csv> to add")

master = pd.read_csv(master_path, low_memory=False)
allcols = list(master.columns)
IPTM = ["af3_msaoff_iptm", "b2_iptm", "esm_iptm"]
HAS = ["has_af3_msaoff", "has_b2", "has_esm"]
print(f"master {master_path}: {len(master)} rows, {len(allcols)} cols")

frames = [master]
for spec in adds:
    ver, path = spec.split(":", 1)
    d = pd.read_csv(path, low_memory=False)
    ipt = d[IPTM].apply(pd.to_numeric, errors="coerce")
    d["consensus_iptm"] = ipt.min(axis=1)
    d["mean_iptm"] = ipt.mean(axis=1)
    d["has_consensus"] = d[HAS].astype(bool).all(axis=1)
    d["campaign"] = ver
    missing = [c for c in allcols if c not in d.columns]
    if missing:
        sys.exit(f"{ver}: still missing master cols after derive: {missing}")
    frames.append(d[allcols])
    print(f"  + {ver} ({path}): {len(d)} rows")

merged = pd.concat(frames, ignore_index=True)
before = len(merged)
merged = merged.drop_duplicates("design_id", keep="first")
print(f"concat {before} -> dedup by design_id {len(merged)}  ({before - len(merged)} dup design_ids dropped)")
merged.to_csv(out_path, index=False)
print(f"wrote {out_path}: {len(merged)} rows x {len(merged.columns)} cols")
print("per-version:")
print(merged["campaign"].value_counts().sort_index().to_string())
