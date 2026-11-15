#!/usr/bin/env python3
"""Investigate why v2 rfd3__ligandmpnn is 100% absent from combined_clean.
Distinguish real gap vs stem/design_id format mismatch vs intentional exclusion."""
import glob, os, re
import pandas as pd

RB = "/cluster/work/wicky/lvalsecchi/rna_binder"
STACK = f"{RB}/campaigns/csy4_crrna_free/rfd3__ligandmpnn/midrun_v2"
MASTER = f"{RB}/analyses/cross_campaign/combined_clean.csv"

def design_stems(list_path):
    s = []
    for ln in open(list_path):
        b = os.path.basename(ln.strip())
        if any(x in b for x in ("_ESM_", "_AF3_", "_B2_")):
            continue
        for e in (".cif.gz", ".cif", ".pdb"):
            if b.endswith(e):
                s.append(b[:-len(e)]); break
    return s

lists = glob.glob(f"{STACK}/designs/*.tar.zst.list")
print(f"v2 rfd3 stack dir exists: {os.path.isdir(STACK)}")
print(f"n design .list files: {len(lists)}")
for be in ("af3", "b2", "esm"):
    rl = glob.glob(f"{STACK}/refold_{be}/msaoff_*/designs.tar.zst.list")
    print(f"  refold_{be} .list files: {len(rl)}")

stems = []
for dl in lists[:5]:
    stems += design_stems(dl)
print(f"\nsample v2 rfd3 design stems ({len(stems)} from first 5 lists):")
for s in stems[:8]:
    print("   ", s)

m = pd.read_csv(MASTER, usecols=["design_id", "campaign"], low_memory=False)
ids = set(m["design_id"])
print(f"\ncombined_clean: {len(m)} rows")
print("campaign value_counts:")
print(m["campaign"].value_counts().sort_index().to_string())

# Exact-match check
hits = [s for s in stems if s in ids]
print(f"\nexact stem-in-table hits (first {len(stems)}): {len(hits)}")

# Do ANY v2 design_ids look like rfd3? Inspect v2 id shapes
v2ids = m.loc[m.campaign == "v2", "design_id"].astype(str)
print(f"\nv2 rows in table: {len(v2ids)}")
print("v2 design_id samples:")
for x in v2ids.head(8):
    print("   ", x)
# token/prefix overlap: compare uuid-like prefixes
def uuid8(x):
    mm = re.search(r"[0-9a-f]{8}", x)
    return mm.group(0) if mm else None
stem_uuids = {uuid8(s) for s in stems if uuid8(s)}
v2_uuids = {uuid8(x) for x in v2ids if uuid8(x)}
print(f"\nstem 8-hex prefixes (sample stems): {len(stem_uuids)}")
print(f"v2-table 8-hex prefixes: {len(v2_uuids)}")
print(f"overlap of those stem-uuids with ANY v2-table uuid: {len(stem_uuids & v2_uuids)}")
print(f"overlap with ANY table uuid (all campaigns): "
      f"{len(stem_uuids & {uuid8(x) for x in m['design_id'].astype(str) if uuid8(x)})}")
