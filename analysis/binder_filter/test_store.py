"""test_store.py — plain-assert tests for the columnar store loader.

Builds a tiny synthetic store in a temp dir and checks the parts that would fail silently in
production: path reconstruction (a corrupted ref only surfaces much later as a failed ChimeraX
open), the all-NaN placeholder detection, and subsample stability.

Run:  ~/miniconda3/bin/python test_store.py
"""
import json
import os
import shutil
import tempfile

import numpy as np

from store import Store

checks = 0


def ck(cond, msg):
    global checks
    checks += 1
    assert cond, msg


# a packed layout and a LOOSE one — both occur in the real campaign, which is why the store
# keys on a whole template triple rather than assuming a common prefix
TPL_PACKED = [
    "campaigns/x/bg__bif/midrun_v2/refold_af3/msaoff_<U>/designs.tar.zst!designs/<ID>_AF3_seed1_sample0.cif",
    "campaigns/x/bg__bif/midrun_v2/refold_esm/msaoff_<U>/designs.tar.zst!designs/<ID>_ESM_sample0.cif",
    "campaigns/x/bg__bif/midrun_v2/designs/<U>.tar.zst!<ID>.cif",
]
TPL_LOOSE = [
    "campaigns/x/hd__limpnn/midrun_v3/refold_af3/msaoff_<U>/designs/<ID>_AF3_seed1_sample0.cif",
    "campaigns/x/hd__limpnn/midrun_v3/refold_esm/msaoff_<U>/designs/<ID>_ESM_sample0.cif",
    "campaigns/x/hd__limpnn/midrun_v3/designs/<U>/<ID>.cif",
]

IDS = [b"csy4free_BG_BIF_len80_1b05319b_design52_seq3",
       b"csy4free_BG_BIF_len80_1b05319b_design52_seq4",
       b"csy4free_HD_LiMPNN_len80_aabbccdd_design1_seq0",
       b"csy4free_HD_LiMPNN_len80_aabbccdd_design2_seq0"]
UUS = [b"1b05319b", b"1b05319b", b"aabbccdd", b"aabbccdd"]
TPL = [0, 0, 1, 1]

d = tempfile.mkdtemp(prefix="fstore_test_")
try:
    np.save(f"{d}/design_id.npy", np.array(IDS, dtype="S64"))
    np.save(f"{d}/uuid8.npy", np.array(UUS, dtype="S16"))
    np.save(f"{d}/tpl_idx.npy", np.array(TPL, dtype=np.int16))
    np.save(f"{d}/stack.npy", np.array([0, 0, 1, 1], dtype=np.int16))
    np.save(f"{d}/campaign.npy", np.array([0, 0, 1, 1], dtype=np.int16))
    np.save(f"{d}/length.npy", np.array([80, 80, 80, 80], dtype=np.int16))
    np.save(f"{d}/iptm_af3.npy", np.array([.9, .8, .7, .6], dtype=np.float32))
    np.save(f"{d}/iptm_esm.npy", np.array([.5, .85, .65, .55], dtype=np.float32))
    np.save(f"{d}/rmsd_af3.npy", np.array([1., 2., 3., np.nan], dtype=np.float32))
    np.save(f"{d}/rog_af3.npy", np.full(4, np.nan, dtype=np.float32))
    with open(f"{d}/meta.json", "w") as f:
        json.dump({"n": 4, "built_utc": "2026-08-07T00:00:00Z",
                   "path_cols": ["af3_msaoff_cif_path", "esm_cif_path", "designed_cif_path"],
                   "templates": [TPL_PACKED, TPL_LOOSE],
                   "vocabs": {"stack": ["BG+BIF", "HD+LiMPNN"], "campaign": ["v2", "v3"]}}, f)

    s = Store(d)
    ck(s.n == 4, "n from meta")

    # ---- path reconstruction, both layouts
    p = s.paths(0)
    ck(p["af3_msaoff_cif_path"] ==
       "campaigns/x/bg__bif/midrun_v2/refold_af3/msaoff_1b05319b/designs.tar.zst!"
       "designs/csy4free_BG_BIF_len80_1b05319b_design52_seq3_AF3_seed1_sample0.cif",
       f"packed af3 path: {p['af3_msaoff_cif_path']}")
    ck(p["designed_cif_path"] ==
       "campaigns/x/bg__bif/midrun_v2/designs/1b05319b.tar.zst!"
       "csy4free_BG_BIF_len80_1b05319b_design52_seq3.cif",
       f"packed design path: {p['designed_cif_path']}")
    p = s.paths(2)
    ck("tar.zst" not in p["af3_msaoff_cif_path"], "loose layout has no archive")
    ck(p["af3_msaoff_cif_path"].endswith(
        "designs/csy4free_HD_LiMPNN_len80_aabbccdd_design1_seq0_AF3_seed1_sample0.cif"),
       f"loose af3 path: {p['af3_msaoff_cif_path']}")

    # the uuid8 is a SUBSTRING of the design_id — substituting in the wrong order corrupts
    # the path, so pin that the id survives intact
    for i in range(4):
        did = s.col("design_id")[i].decode()
        for v in s.paths(i).values():
            ck(did in v, f"row {i}: design_id mangled in {v}")

    # ---- derived consensus is computed here, never read from a b2-tainted column
    ck(np.allclose(s.col("iptm_min"), [.5, .8, .65, .55]), "iptm_min = min(af3, esm)")
    ck(np.allclose(s.col("iptm_mean"), [.7, .825, .675, .575]), "iptm_mean")
    ck(not any("b2" in c for c in s.numeric_cols()), "no b2 column may exist in the store")

    # ---- placeholder detection
    ck(s.has_data("iptm_af3"), "real column has data")
    ck(not s.has_data("rog_af3"), "all-NaN placeholder reports no data")
    ck(s.has_data("rmsd_af3"), "a column with SOME NaN still has data")

    # ---- labels
    ck(s.label("stack", 1) == "HD+LiMPNN", "vocab lookup")
    ck(s.label("stack", 99) == "", "out-of-range code is empty, not an exception")

    # ---- subsample is stable across calls (thresholds must not shuffle the point cloud)
    a = s.subsample(3, seed=7)
    b = s.subsample(3, seed=7)
    ck(np.array_equal(a, b), "same seed gives the same subsample")
    ck(a.size == 3 and np.all(np.diff(a) > 0), "subsample is sorted and unique")
    ck(s.subsample(999).size == 4, "subsample larger than the store is clamped to n")

    st = s.subsample(2, stratify="stack")
    ck(set(s.col("stack")[st].tolist()) == {0, 1}, "stratified sample covers every stack")

    ck(s.design_ids([0, 2]) == [IDS[0].decode(), IDS[2].decode()], "design_ids decodes")
    print(f"all {checks} checks passed")
finally:
    shutil.rmtree(d, ignore_errors=True)
