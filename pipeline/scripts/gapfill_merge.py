#!/usr/bin/env python3
"""Gap-fill by REPLACING whole (campaign, stack) slices with freshly re-aggregated
clean sets, so late-completing (Type-A gap) designs are captured while sequence-dedup
integrity is preserved (fresh clean is internally seq-deduped).

For each <ver:clean.csv>:
  - derive the 4 master-only cols (campaign, consensus_iptm, mean_iptm, has_consensus)
  - find the stack value(s) present in the fresh clean
  - DROP master rows where campaign==ver AND stack in those values  (report the count;
    it should match the gap report's `in_table` for that stack)
  - append the fresh, aligned rows

Then design_id-dedup (keep first) and write. Nothing is clobbered: out != master.

  python3 gapfill_merge.py <master.csv> <out.csv> <ver:clean.csv> [<ver:clean.csv> ...]
"""
import sys
import pandas as pd

master_path, out_path = sys.argv[1], sys.argv[2]
adds = sys.argv[3:]
if not adds:
    sys.exit("need at least one <ver:clean.csv>")

# expected in_table (drop) counts from the gap report, for the safety check
EXPECT_DROP = {
    ("v2", "rfd3__ligandmpnn"): 0,
    ("v3", "halludesign__ligandmpnn"): 183590,
    ("v4", "halludesign__ligandmpnn"): 88840,
    ("v5", "halludesign__ligandmpnn"): 66850,
}

master = pd.read_csv(master_path, low_memory=False)
allcols = list(master.columns)
IPTM = ["af3_msaoff_iptm", "b2_iptm", "esm_iptm"]
HAS = ["has_af3_msaoff", "has_b2", "has_esm"]
print(f"master {master_path}: {len(master)} rows, {len(allcols)} cols")
print(f"master per-campaign:\n{master['campaign'].value_counts().sort_index().to_string()}\n")

fresh_frames = []
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
        sys.exit(f"{ver} ({path}): missing master cols after derive: {missing}")
    stacks = sorted(d["stack"].dropna().unique().tolist())
    if len(stacks) != 1:
        print(f"  !! WARN {ver} ({path}): fresh clean has {len(stacks)} stack values: {stacks}")
    # DROP the master slice(s) for this (campaign, stack)
    mask = (master["campaign"] == ver) & (master["stack"].isin(stacks))
    dropped = int(mask.sum())
    key = (ver, stacks[0] if stacks else "?")
    exp = EXPECT_DROP.get(key)
    flag = "" if exp is None else ("  OK" if dropped == exp else f"  <<< MISMATCH expected {exp}")
    print(f"  {ver}/{stacks}: fresh_clean={len(d)}  drop_from_master={dropped}{flag}")
    master = master[~mask]
    fresh_frames.append(d[allcols])

merged = pd.concat([master] + fresh_frames, ignore_index=True)
before = len(merged)
merged = merged.drop_duplicates("design_id", keep="first")
print(f"\nconcat {before} -> dedup by design_id {len(merged)}  ({before - len(merged)} dup design_ids dropped)")
merged.to_csv(out_path, index=False)
print(f"wrote {out_path}: {len(merged)} rows x {len(merged.columns)} cols")
print(f"\nper-campaign:\n{merged['campaign'].value_counts().sort_index().to_string()}")
