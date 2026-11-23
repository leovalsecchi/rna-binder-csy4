#!/usr/bin/env python3
"""Re-pick the two BoltzGen call-outs of Fig 2b so AF3 reproduces the design.

WHY. The figure renders and measures every call-out on the AF3 complex, but the
two BoltzGen panels were picked on the BETTER of the two refolds, and for them
the better one is ESMFold2: AF3 sits 29.0 A (bg_a, phosphate-reading) and 8.5 A
(bg_b, hairpin-head) from the design, DockQ 0.031 and 0.450. So those cells draw
a pose the design does not have. Fixed at the selection step, not per panel:
every candidate here must be reproduced by BOTH backends.

Mirrors code/Fig2b_map_callouts/pick_examples.py for the pool, the fold-map join
and the mode scores - same tables, same floors - and changes only the filter and
the scope (BG+LiMPNN, in the 12k).

Layout windows: the composite splits the map at x = 4.40 and orders slots by y,
which is what keeps 0 leader crossings. To replace bg_a in slot L2 and bg_b in
slot R1 without re-deriving the whole layout:
    bg_a  x < 4.40  and  2.60 < y < 7.08     (between L3 and L1)
    bg_b  x > 4.40  and  y > 5.45            (above R2)
Candidates outside those windows are still printed, flagged, because a better
binder may be worth re-doing the slots for.
"""
import json
import numpy as np
import pandas as pd

FM = "/cluster/work/wicky/lvalsecchi/rna_binder/analyses/foldseek_map"
STORE = ("/cluster/work/wicky/lvalsecchi/rna_binder/analyses/cross_campaign/"
         "filter_store_v2v9_hdtop")
D = f"{FM}/results_bbmap/design"

LOOP = [8, 9, 10]
STEM = [2, 3, 4, 5, 6, 7, 11, 12, 13, 14, 15, 16]
MOI = ("phosphate", "sugar", "base")
SHOWN_CLUSTERS = {7, 8, 5, 3, 4, 14}        # the six call-outs that stay
# AF3 complexes already on the Mac, so a pick among these needs no fetch
LOCAL = {"csy4free_BG_LiMPNN_len80_1caaae08_design23_seq0",
         "csy4free_BG_LiMPNN_len80_5de20158_design12_seq0",
         "csy4free_BG_LiMPNN_len80_81878aac_design15_seq0",
         "csy4free_BG_LiMPNN_len80_97eb1c8c_design25_seq0",
         "csy4free_BG_LiMPNN_len80_aebd763e_design11_seq0"}


def col(i, m):
    return "nt%02d_%s" % (i, m)


fp = pd.read_parquet("%s/fp/canonical.parquet" % D)
em = pd.read_parquet("%s/embed/fold/embed.parquet" % D)
nod = pd.read_csv("%s/dedup/nodes.tsv" % D, sep="\t")
can = pd.read_csv("%s/dedup/canonical.tsv" % D, sep="\t")
cl = pd.read_parquet("%s/clusters/fold/clusters.parquet" % D)[["node_id", "cluster_r10"]]

t = (nod.merge(can, on="struct_name").merge(em, on="node_id")
        .merge(cl, on="node_id").merge(fp, on="node_id"))


def hashes(tag):
    p = "%s/results_bbmap_%s/design/dedup/canonical.tsv" % (FM, tag)
    return set(pd.read_csv(p, sep="\t")["coord_hash"])


hp, hl = hashes("pass"), hashes("lib")
t["stage"] = np.where(t["coord_hash"].isin(hl), 2,
             np.where(t["coord_hash"].isin(hp), 1, 0)).astype("int8")

e = pd.read_parquet("%s/graph/fold/edges_cluster.parquet" % D, columns=["i", "j"])
deg = pd.concat([e["i"], e["j"]]).value_counts()
t["deg0"] = ~t["node_id"].isin(deg.index)

t = t[t["stage"] >= 1].copy()               # the gate, as pick_examples.py used
print("pool: %d gated backbones (stage>=1)" % len(t))

for m in MOI:
    t["c_%s" % m] = t[[col(i, m) for i in range(1, 17)]].sum(axis=1)
t["c_total"] = t[["c_%s" % m for m in MOI]].sum(axis=1)
t = t[t["c_total"] >= 50].copy()
per_nt = {i: t[[col(i, m) for m in MOI]].sum(axis=1) for i in range(1, 17)}
t["c_loop"] = sum(per_nt[i] for i in LOOP)
t["c_stem"] = sum(per_nt[i] for i in STEM)
t["c_stem_base"] = t[[col(i, "base") for i in STEM]].sum(axis=1)
t["f_backbone"] = (t["c_phosphate"] + t["c_sugar"]) / t["c_total"]
t["f_base"] = t["c_base"] / t["c_total"]
t["f_loop"] = t["c_loop"] / t["c_total"]
t["f_stem_base"] = t["c_stem_base"] / t["c_total"]

want = ["design_id", "uuid8", "tpl_idx", "campaign", "stack", "seq_run",
        "length", "iptm_af3", "iptm_esm", "dockq_af3", "dockq_esm",
        "rmsd_af3", "rmsd_esm", "hb_af3", "hb_esm", "pis_af3", "pis_esm",
        "foot_af3", "foot_esm"]
vocab = json.load(open("%s/meta.json" % STORE))["vocabs"]
ids = np.load("%s/design_id.npy" % STORE, mmap_mode="r")
assert ids.dtype.kind == "S", ids.dtype
target = pd.Index(t["design_id"].str.encode("ascii").values)
row = np.full(len(target), -1, dtype=np.int64)
CH = 1_000_000
for s in range(0, len(ids), CH):
    chunk = np.asarray(ids[s:s + CH])
    loc = target.get_indexer(chunk)
    hit = loc >= 0
    row[loc[hit]] = np.flatnonzero(hit) + s
found = row >= 0
print("matched in store: %d / %d" % (found.sum(), len(target)))
t = t[found].copy()
row = row[found]
for c in want:
    if c == "design_id":
        continue
    arr = np.load("%s/%s.npy" % (STORE, c), mmap_mode="r")
    v = np.asarray(arr[row])
    if c in vocab:
        v = np.array(vocab[c], dtype=object)[v]
    elif v.dtype.kind == "S":
        v = np.char.decode(v, "ascii")
    t[c] = v

t["in12k"] = t["stage"] >= 2
t = t[t["stack"] == "BG+LiMPNN"].copy()
print("BG+LiMPNN in the 12k with metrics: %d" % len(t))

# THE FIX: both backends must reproduce the design, not the better of them.
t["rmsd_worst"] = t[["rmsd_af3", "rmsd_esm"]].max(axis=1)
t["dockq_worst"] = t[["dockq_af3", "dockq_esm"]].min(axis=1)
t["iptm_worst"] = t[["iptm_af3", "iptm_esm"]].min(axis=1)
q = t[(t["rmsd_af3"] <= 3.0) & (t["rmsd_esm"] <= 3.0)
      & (t["dockq_worst"] >= 0.5) & (t["iptm_worst"] >= 0.7)].copy()
print("both backends reproduce it (rmsd_af3<=3 AND rmsd_esm<=3, "
      "min dockq>=0.5, min iptm>=0.7): %d" % len(q))


def z(s):
    return (s - s.mean()) / s.std() if len(s) > 1 and s.std() else s * 0


q["s_backbone_hbond"] = z(q["hb_af3"]) + z(q["f_backbone"])
q["s_head"] = z(q["f_loop"]) + 0.3 * z(q["c_loop"])
FLOOR = {"backbone_hbond": lambda d: d["f_backbone"] >= 0.85,
         "head": lambda d: d["f_loop"] >= 0.60}
WINDOW = {"backbone_hbond": ("bg_a / slot L2", lambda d: (d["x"] < 4.40) & (d["y"] > 2.60) & (d["y"] < 7.08)),
          "head": ("bg_b / slot R1", lambda d: (d["x"] > 4.40) & (d["y"] > 5.45))}
SHOW = ["design_id", "cluster_r10", "x", "y", "rmsd_af3", "rmsd_esm",
        "dockq_worst", "iptm_worst", "hb_af3",
        "f_backbone", "f_base", "f_loop", "c_total", "in12k", "in_window",
        "new_cluster", "have_af3"]

KEY = {"backbone_hbond": "f_backbone", "head": "f_loop"}
FLOORVAL = {"backbone_hbond": 0.85, "head": 0.60}
for mode in ("backbone_hbond", "head"):
    slot, win = WINDOW[mode]
    k = KEY[mode]
    print("\n--- %s: how close the both-backends pool gets on %s (floor %.2f)"
          % (mode, k, FLOORVAL[mode]))
    print("    best 8: " + ", ".join("%.3f" % v for v in
          q[k].sort_values(ascending=False).head(8)))
    inw = q[win(q)]
    print("    best 5 INSIDE the layout window (n=%d): " % len(inw)
          + (", ".join("%.3f" % v for v in inw[k].sort_values(ascending=False).head(5))
             if len(inw) else "none"))
    elig = q[FLOOR[mode](q)].copy()
    elig["in_window"] = win(elig)
    elig["new_cluster"] = ~elig["cluster_r10"].isin(SHOWN_CLUSTERS)
    elig["have_af3"] = elig["design_id"].isin(LOCAL)
    elig = elig.sort_values(["in_window", "new_cluster", "s_%s" % mode],
                            ascending=[False, False, False])
    print("\n=== %s  -> %s   %d clear the floor, %d of them in the layout window"
          % (mode, slot, len(elig), int(elig["in_window"].sum())))
    if len(elig):
        print(elig[SHOW].head(12).round(3).to_string(index=False))
    else:
        print("  NOTHING clears the floor under the both-backends filter.")

out = "/cluster/work/wicky/lvalsecchi/fig2b_repick_bg.csv"
q.sort_values("s_head", ascending=False)[
    ["design_id", "node_id", "coord_hash", "cluster_r10", "x", "y", "stage",
     "rmsd_af3", "rmsd_esm", "dockq_af3", "dockq_esm", "iptm_af3", "iptm_esm",
     "hb_af3", "pis_af3", "foot_af3", "f_backbone", "f_base", "f_loop",
     "f_stem_base", "c_total", "s_backbone_hbond", "s_head"]].to_csv(out, index=False)
print("\nwrote %s  (%d rows)" % (out, len(q)))
