#!/usr/bin/env python3
"""Design-level refold completion for midrun_v8, per stack.
Non-HD: <stack>/midrun_v8/{designs/*.list, refold_<be>/msaoff_<uuid>/designs.tar.zst.list}
HD:     <stack>/midrun_v8_<uuid>/designs.tar.zst.list  (+ refold_<be> if present)
Reports total designs vs per-backend refolded vs all-3-complete."""
import glob, os, re

RB = "/cluster/work/wicky/lvalsecchi/rna_binder"
C = f"{RB}/campaigns/csy4_crrna_free"
TAG = {"af3": "AF3", "b2": "B2", "esm": "ESM"}
BES = ["af3", "b2", "esm"]

def design_stems(list_path):
    s = set()
    for ln in open(list_path):
        b = os.path.basename(ln.strip())
        if any(x in b for x in ("_ESM_", "_AF3_", "_B2_")):
            continue
        for e in (".cif.gz", ".cif", ".pdb"):
            if b.endswith(e):
                s.add(b[:-len(e)]); break
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

def report(name, designs, comp):
    tot = len(designs)
    a, b, e = (len(comp[x]) for x in BES)
    all3 = len(comp["af3"] & comp["b2"] & comp["esm"])
    pct = (100 * all3 / tot) if tot else 0.0
    print(f"  {name:24} designs={tot:>7}  af3={a:>7} b2={b:>7} esm={e:>7}  all3={all3:>7}  ({pct:5.1f}%)")
    return tot, all3

print("=== midrun_v8 refold completion (design-level, from .list) ===\n")
grand_t = grand_c = 0

# --- non-HD stacks: standard midrun_v8 layout ---
for st in ["rfd3__ligandmpnn", "boltzgen__boltzif", "boltzgen__ligandmpnn", "boltzdesign__ligandmpnn"]:
    base = f"{C}/{st}/midrun_v8"
    designs = set()
    comp = {x: set() for x in BES}
    for dl in glob.glob(f"{base}/designs/*.tar.zst.list"):
        mm = re.search(r"([0-9a-f]{8})\.tar", dl)
        if not mm:
            continue
        u = mm.group(1)
        ds = design_stems(dl)
        designs |= ds
        for be in BES:
            comp[be] |= ds & refold_stems(f"{base}/refold_{be}/msaoff_{u}/designs.tar.zst.list", TAG[be])
    t, c = report(st, designs, comp)
    grand_t += t; grand_c += c

# --- HD: per-uuid dirs midrun_v8_<uuid> ---
hd_dirs = sorted(glob.glob(f"{C}/halludesign__ligandmpnn/midrun_v8_*"))
hd_dirs = [d for d in hd_dirs if os.path.isdir(d)]
designs = set()
comp = {x: set() for x in BES}
gen_done_uuids = 0
for d in hd_dirs:
    u = d.split("midrun_v8_")[-1]
    dl = f"{d}/designs.tar.zst.list"
    if os.path.exists(dl):
        ds = design_stems(dl)
        if ds:
            gen_done_uuids += 1
        designs |= ds
        for be in BES:
            comp[be] |= ds & refold_stems(f"{d}/refold_{be}/msaoff_{u}/designs.tar.zst.list", TAG[be])
print(f"  (HD per-uuid dirs: {len(hd_dirs)} total, {gen_done_uuids} with generated designs)")
t, c = report("halludesign__ligandmpnn", designs, comp)
grand_t += t; grand_c += c

print(f"\n=== v8 TOTAL: {grand_c}/{grand_t} designs all-3-refold-complete "
      f"({100*grand_c/grand_t if grand_t else 0:.1f}%) ===")
