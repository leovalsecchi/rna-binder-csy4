#!/usr/bin/env python3
"""Pull the two BoltzGen call-outs' ESMFold2 complexes out of the repacked tars.

The store's path templates are stale: the 2026-09-11 repack collapsed every
midrun directory into one solid tar.zst, so a refold now sits DOUBLE-NESTED

    <midrun>.tar.zst ! <midrun>/refold_esm/msaoff_<uuid>/designs.tar.zst
                     ! designs/<design_id>_ESM_sample0.cif

and a solid zstd stream has no random access. Both members sit ~90% of the way
in (16.1 GB and 7.6 GB), so each outer archive is streamed ONCE, sequentially,
and abandoned the moment its member is in hand. The inner archive is
decompressed in memory and dropped -- expanding these to scratch is what filled
/local/scratch in July 2026.
"""
import io
import os
import sys
import tarfile

import zstandard as zstd

B = ("/cluster/work/wicky/lvalsecchi/rna_binder/campaigns/csy4_crrna_free/"
     "boltzgen__ligandmpnn")
OUT = "/cluster/work/wicky/lvalsecchi/fig2b_esm"
JOBS = [
    ("bg_a", B + "/midrun_v2.tar.zst",
     "midrun_v2/refold_esm/msaoff_aebd763e/designs.tar.zst",
     "csy4free_BG_LiMPNN_len80_aebd763e_design11_seq0_ESM_sample0.cif"),
    ("bg_b", B + "/midrun_v4.tar.zst",
     "midrun_v4/refold_esm/msaoff_e1e86ec3/designs.tar.zst",
     "csy4free_BG_LiMPNN_len80_e1e86ec3_design9_seq0_ESM_sample0.cif"),
]

os.makedirs(OUT, exist_ok=True)
dctx = zstd.ZstdDecompressor()

for panel, outer, inner_name, cif in JOBS:
    print("\n=== %s  %s" % (panel, os.path.basename(outer)), flush=True)
    dst = os.path.join(OUT, cif)
    if os.path.exists(dst):
        print("  already have it, skipping", flush=True)
        continue
    blob = None
    seen = 0
    with open(outer, "rb") as fh, dctx.stream_reader(fh) as r:
        # r| : sequential streaming tar, the only mode a zstd stream supports
        with tarfile.open(fileobj=r, mode="r|") as tf:
            for m in tf:
                seen += 1
                if seen % 2000 == 0:
                    print("  ...%d members" % seen, flush=True)
                if m.name == inner_name:
                    print("  found inner archive at member %d (%.1f MB)"
                          % (seen, m.size / 1048576.0), flush=True)
                    blob = tf.extractfile(m).read()
                    break
    if blob is None:
        print("  INNER ARCHIVE NOT FOUND after %d members" % seen, flush=True)
        continue
    with dctx.stream_reader(io.BytesIO(blob)) as r2:
        data = io.BytesIO(r2.read())
    del blob
    with tarfile.open(fileobj=data, mode="r:") as tf2:
        names = tf2.getnames()
        hit = [n for n in names if n.endswith(cif)]
        if not hit:
            print("  CIF NOT in inner archive. %d entries, e.g. %s"
                  % (len(names), names[:3]), flush=True)
            continue
        payload = tf2.extractfile(hit[0]).read()
    with open(dst, "wb") as fh:
        fh.write(payload)
    print("  -> %s  (%d bytes)" % (dst, len(payload)), flush=True)

print("\n=== %s ===" % OUT, flush=True)
for f in sorted(os.listdir(OUT)):
    print("  ", f, os.path.getsize(os.path.join(OUT, f)), flush=True)
