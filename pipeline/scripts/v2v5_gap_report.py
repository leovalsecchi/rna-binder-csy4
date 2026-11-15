#!/usr/bin/env python3
"""Report v2-v5 Type-A gaps: designs that are refold-complete NOW (present in all 3
backends' .list sidecars) but NOT in combined_clean. Reads .list sidecars only (no CIFs).
Report only -- fills nothing. Per (version, stack): refold_complete vs in_table vs gap.

  python3 v2v5_gap_report.py
"""
import glob
import os
import re

import pandas as pd

RB = "/cluster/work/wicky/lvalsecchi/rna_binder"
CAMP = f"{RB}/campaigns/csy4_crrna_free"
MASTER = f"{RB}/analyses/cross_campaign/combined_clean.csv"
VERSIONS = ["midrun_v2", "midrun_v3", "midrun_v4", "midrun_v5"]
VER_TAG = {"midrun_v2": "v2", "midrun_v3": "v3", "midrun_v4": "v4", "midrun_v5": "v5"}
TAG = {"af3": "AF3", "b2": "B2", "esm": "ESM"}
BACKENDS = ["af3", "b2", "esm"]


def design_stems(list_path):
    s = set()
    for ln in open(list_path):
        b = os.path.basename(ln.strip())
        if "_ESM_" in b or "_AF3_" in b or "_B2_" in b:
            continue
        for e in (".cif.gz", ".cif", ".pdb"):
            if b.endswith(e):
                s.add(b[:-len(e)])
                break
    return s


def refold_stems(list_path, tag):
    s = set()
    if not os.path.exists(list_path):
        return s
    marker = f"_{tag}_"
    for ln in open(list_path):
        b = os.path.basename(ln.strip())
        i = b.find(marker)
        if i > 0:
            s.add(b[:i])
    return s


m = pd.read_csv(MASTER, usecols=["design_id", "campaign"], low_memory=False)
in_table_ids = set(m["design_id"])
print(f"combined_clean: {len(m)} rows, {len(in_table_ids)} unique design_ids")
print(f"{'ver':4} {'stack':26} {'refold_complete':>15} {'in_table':>9} {'gap':>8}")

tot_gap = 0
rows = []
for V in VERSIONS:
    vt = VER_TAG[V]
    for stackpath in sorted(glob.glob(f"{CAMP}/*/{V}")):
        stack = stackpath.split("/")[-2]
        complete = set()
        for dl in glob.glob(f"{stackpath}/designs/*.tar.zst.list"):
            mm = re.search(r"([0-9a-f]{8})\.tar", dl)
            if not mm:
                continue
            u = mm.group(1)
            rc = design_stems(dl)
            for be in BACKENDS:
                rc &= refold_stems(f"{stackpath}/refold_{be}/msaoff_{u}/designs.tar.zst.list", TAG[be])
            complete |= rc
        gap = complete - in_table_ids
        in_tab = len(complete & in_table_ids)
        tot_gap += len(gap)
        rows.append((vt, stack, len(complete), in_tab, len(gap)))
        print(f"{vt:4} {stack:26} {len(complete):>15} {in_tab:>9} {len(gap):>8}")

print(f"\n=== TOTAL Type-A gap candidates (refold-complete, absent from combined_clean): {tot_gap} ===")
print("(candidates -- some may be legit cysteine/sequence-dup drops; real gaps <= this)")
