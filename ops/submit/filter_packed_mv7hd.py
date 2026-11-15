#!/usr/bin/env python3
"""Rebuild the two mv7_hd consolidated arrays keeping ONLY uuids that are NOT already packed
(shared-pool OR run-dir tar). Writes new arr files (does not cancel/submit)."""
import glob, os, re
HD = "/cluster/work/wicky/lvalsecchi/rna_binder/campaigns/csy4_crrna_free/halludesign__ligandmpnn"
AR = "/cluster/work/wicky/lvalsecchi/rna_binder/campaigns/csy4_crrna_free/_mv7_arrays"

packed = set()
for t in glob.glob(HD + "/midrun_v7/designs/*.tar.zst"):
    packed.add(os.path.basename(t)[:-len(".tar.zst")])
for t in glob.glob(HD + "/midrun_v7_*/designs.tar.zst"):
    m = re.search(r"midrun_v7_([0-9a-f]{8})", t)
    if m:
        packed.add(m.group(1))
print("packed uuids:", len(packed))

def uuid_of(path):
    m = re.search(r"midrun_v7_([0-9a-f]{8})", path)
    return m.group(1) if m else None

for src, dst in [("arr_mv7hd_reconsolidated.sbatch", "arr_mv7hd_final.sbatch"),
                 ("arr_mv7hd_recon_gpuhe_24h.sbatch", "arr_mv7hd_final_gpuhe_24h.sbatch")]:
    sp = os.path.join(AR, src)
    if not os.path.exists(sp):
        print("MISSING", src); continue
    lines = open(sp).read().split("\n")
    header, scripts, tail, mode = [], [], [], 0
    for l in lines:
        if mode == 0:
            header.append(l)
            if l.strip() == "SCRIPTS=(":
                mode = 1
        elif mode == 1:
            if l.strip() == ")":
                tail.append(l); mode = 2
            elif l.strip().startswith("/cluster"):
                scripts.append(l.strip())
            else:
                tail.append(l)
        else:
            tail.append(l)
    kept = [s for s in scripts if uuid_of(s) and uuid_of(s) not in packed]
    dropped = len(scripts) - len(kept)
    if not kept:
        print("%s -> ALL %d already packed, nothing to submit" % (src, len(scripts))); continue
    header = [re.sub(r"--array=0-\d+", "--array=0-%d" % (len(kept) - 1), h) for h in header]
    out = "\n".join(header + ["  " + s for s in kept] + tail).rstrip("\n") + "\n"
    open(os.path.join(AR, dst), "w").write(out)
    print("%s -> %s : %d -> %d  (dropped %d already-packed)" % (src, dst, len(scripts), len(kept), dropped))
