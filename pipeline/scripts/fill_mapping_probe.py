#!/usr/bin/env python3
"""Probe the info needed to gap-fill precisely: id-prefix per (campaign,stack),
the sequence column name, and the target stacks' dirs/refold coverage."""
import glob, os, re
import pandas as pd

RB = "/cluster/work/wicky/lvalsecchi/rna_binder"
CAMP = f"{RB}/campaigns/csy4_crrna_free"
OUT = f"{RB}/analyses/cross_campaign/combined_clean_v2v7.csv"

# 1) columns of the current deliverable
hdr = pd.read_csv(OUT, nrows=3, low_memory=False)
print("combined_clean_v2v7 columns (%d):" % len(hdr.columns))
print("  " + ", ".join(hdr.columns))
seqcol = [c for c in hdr.columns if "seq" in c.lower()]
print("seq-like cols:", seqcol)

# 2) id-prefix map: reload campaign+design_id, derive prefix = up to _len token
d = pd.read_csv(OUT, usecols=["design_id", "campaign"], low_memory=False)
def pref(x):
    m = re.match(r"(csy4free_[A-Za-z0-9]+_[A-Za-z0-9]+)_len", str(x))
    return m.group(1) if m else str(x)[:24]
d["pref"] = d["design_id"].map(pref)
print("\n(campaign, id-prefix) -> row count:")
g = d.groupby(["campaign", "pref"]).size().sort_index()
print(g.to_string())

# 3) target stacks: dirs + refold coverage
targets = [("midrun_v2", "rfd3__ligandmpnn"),
           ("midrun_v3", "halludesign__ligandmpnn"),
           ("midrun_v4", "halludesign__ligandmpnn"),
           ("midrun_v5", "halludesign__ligandmpnn")]
print("\ntarget stacks:")
for ver, stack in targets:
    sp = f"{CAMP}/{stack}/{ver}"
    nd = len(glob.glob(f"{sp}/designs/*.tar.zst.list"))
    cov = {be: len(glob.glob(f"{sp}/refold_{be}/msaoff_*/designs.tar.zst.list"))
           for be in ("af3", "b2", "esm")}
    # sample a stem to see the id-prefix
    stem = None
    for dl in glob.glob(f"{sp}/designs/*.tar.zst.list")[:1]:
        for ln in open(dl):
            b = os.path.basename(ln.strip())
            if any(x in b for x in ("_ESM_", "_AF3_", "_B2_")):
                continue
            for e in (".cif.gz", ".cif", ".pdb"):
                if b.endswith(e):
                    stem = b[:-len(e)]; break
            if stem:
                break
    print(f"  {ver}/{stack}: exists={os.path.isdir(sp)} designUUIDs={nd} refold={cov} sample_prefix={pref(stem) if stem else None}")
