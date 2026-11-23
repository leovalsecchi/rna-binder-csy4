#!/usr/bin/env python3
"""Shortlist binders that show specific RNA interaction modes, and locate them
on the campaign fold map (Figure 2b).

Runs on Euler. Candidate pool is every backbone that passed the design gate
(stage >= 1), because an example worth rendering should be a real candidate.

RNA frame, verified not assumed: the canonical 16-mer is ACUGCCGUAUAGGCAG, and
base_pairs() over all 50 browse50 AF3 refolds gives the SAME six pairs every
time: (2,16) (3,15) (4,14) (5,13) (6,12) (7,11). So

    nt 1            free 5' end
    nt 2-7, 11-16   stem (6 bp; 7-11 is the non-canonical closing pair)
    nt 8,9,10       loop / hairpin head  (U8 A9 U10)

Interaction classes scored here:
  backbone_hbond  many interface H-bonds AND contacts dominated by
                  phosphate+sugar rather than base
  pi_stack        high pi_stack_score AND a base-contact-rich interface
  head            contacts concentrated on the loop, nt 8-10
  groove          contacts concentrated on STEM BASES, which is what reading a
                  groove looks like in this binning. NOTE: the 48-dim
                  fingerprint bins phosphate/sugar/base only, so it CANNOT
                  separate major from minor groove. That needs atom-level
                  geometry (N7/O6/N4 vs O2/N3/O2') on the extracted structures.
                  Treat this column as "groove-like, VERIFY which groove".

Writes examples.csv with the identity fields needed for the structure registry
(design_id, uuid8, tpl_idx, campaign, stack, seq_run, length, foot) plus the
fold-map coordinates.
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


def col(i, m):
    return f"nt{i:02d}_{m}"


# ---------------------------------------------------------------- fold map
fp = pd.read_parquet(f"{D}/fp/canonical.parquet")
em = pd.read_parquet(f"{D}/embed/fold/embed.parquet")
nod = pd.read_csv(f"{D}/dedup/nodes.tsv", sep="\t")
can = pd.read_csv(f"{D}/dedup/canonical.tsv", sep="\t")
cl = pd.read_parquet(f"{D}/clusters/fold/clusters.parquet")[["node_id", "cluster_r10"]]

t = (nod.merge(can, on="struct_name").merge(em, on="node_id")
        .merge(cl, on="node_id").merge(fp, on="node_id"))

def hashes(tag):
    p = f"{FM}/results_bbmap_{tag}/design/dedup/canonical.tsv"
    return set(pd.read_csv(p, sep="\t")["coord_hash"])

hp, hl = hashes("pass"), hashes("lib")
t["stage"] = np.where(t["coord_hash"].isin(hl), 2,
             np.where(t["coord_hash"].isin(hp), 1, 0)).astype("int8")

# degree 0 in the THRESHOLDED clustering graph (not knn_layout)
e = pd.read_parquet(f"{D}/graph/fold/edges_cluster.parquet", columns=["i", "j"])
deg = pd.concat([e["i"], e["j"]]).value_counts()
t["deg0"] = ~t["node_id"].isin(deg.index)

t = t[t["stage"] >= 1].copy()
print(f"candidate pool: {len(t):,} backbones that passed the gate")

# ------------------------------------------------------- contact breakdown
for m in MOI:
    t[f"c_{m}"] = t[[col(i, m) for i in range(1, 17)]].sum(axis=1)
t["c_total"] = t[[f"c_{m}" for m in MOI]].sum(axis=1)
t = t[t["c_total"] >= 50].copy()          # a real interface, not a graze
print(f"with >=50 heavy-atom contacts: {len(t):,}")

per_nt = {i: t[[col(i, m) for m in MOI]].sum(axis=1) for i in range(1, 17)}
t["c_loop"] = sum(per_nt[i] for i in LOOP)
t["c_stem"] = sum(per_nt[i] for i in STEM)
t["c_nt1"] = per_nt[1]
t["c_stem_base"] = t[[col(i, "base") for i in STEM]].sum(axis=1)
t["c_loop_base"] = t[[col(i, "base") for i in LOOP]].sum(axis=1)

t["f_backbone"] = (t["c_phosphate"] + t["c_sugar"]) / t["c_total"]
t["f_base"] = t["c_base"] / t["c_total"]
t["f_loop"] = t["c_loop"] / t["c_total"]
t["f_stem"] = t["c_stem"] / t["c_total"]
t["f_stem_base"] = t["c_stem_base"] / t["c_total"]

# --------------------------------------------------------------- store join
want = ["design_id", "uuid8", "tpl_idx", "campaign", "stack", "seq_run",
        "length", "iptm_af3", "iptm_esm", "dockq_af3", "dockq_esm",
        "rmsd_af3", "rmsd_esm", "hb_af3", "hb_esm", "pis_af3", "pis_esm",
        "fbase_af3", "fbase_esm", "foot_af3", "foot_esm"]
meta = json.load(open(f"{STORE}/meta.json"))
vocab = meta["vocabs"]

ids = np.load(f"{STORE}/design_id.npy", mmap_mode="r")
# design_id.npy is |S64 bytes; the fingerprint table holds str. Compare as
# bytes or every lookup silently misses and the shortlist comes back empty.
assert ids.dtype.kind == "S", f"unexpected store id dtype {ids.dtype}"
target = pd.Index(t["design_id"].str.encode("ascii").values)
row = np.full(len(target), -1, dtype=np.int64)
CH = 1_000_000
for s in range(0, len(ids), CH):                 # chunked: 7M ids is 1.6 GB
    chunk = np.asarray(ids[s:s + CH])
    loc = target.get_indexer(chunk)
    hit = loc >= 0
    row[loc[hit]] = np.flatnonzero(hit) + s
found = row >= 0
print(f"matched in store: {found.sum():,} / {len(target):,}")
t = t[found].copy()
row = row[found]

for c in want:
    if c == "design_id":
        continue
    if c in ("uuid8",):
        arr = np.load(f"{STORE}/{c}.npy", mmap_mode="r")
        t[c] = np.char.decode(np.asarray(arr[row]), "ascii")
        continue
    arr = np.load(f"{STORE}/{c}.npy", mmap_mode="r")
    v = np.asarray(arr[row])
    if c in vocab:
        v = np.array(vocab[c], dtype=object)[v]
    elif v.dtype.kind == "S":
        v = np.char.decode(v, "ascii")
    t[c] = v

# Consensus = worst of the two refolds, matching the campaign summary
t["iptm"] = t[["iptm_af3", "iptm_esm"]].min(axis=1)
t["dockq"] = t[["dockq_af3", "dockq_esm"]].min(axis=1)
t["rmsd"] = t[["rmsd_af3", "rmsd_esm"]].max(axis=1)
t["hb"] = t[["hb_af3", "hb_esm"]].min(axis=1)
t["pis"] = t[["pis_af3", "pis_esm"]].min(axis=1)
t["foot"] = t[["foot_af3", "foot_esm"]].max(axis=1)

# Only show off binders that are credible to begin with.
q = t[(t["iptm"] >= 0.7) & (t["dockq"] >= 0.5) & (t["rmsd"] <= 3.0)].copy()
print(f"credible (iptm>=0.7, dockq>=0.5, rmsd<=3): {len(q):,}")


def z(s):
    return (s - s.mean()) / s.std()


q["s_backbone_hbond"] = z(q["hb"]) + z(q["f_backbone"])
q["s_pi_stack"] = z(q["pis"]) + z(q["f_base"])
q["s_head"] = z(q["f_loop"]) + 0.3 * z(q["c_loop"])
q["s_groove"] = z(q["f_stem_base"]) + 0.3 * z(q["c_stem"])

# A pick must actually EXEMPLIFY its mode, not merely rank top-n within a thin
# pool. Without these floors the groove and pi-stack lists fill up with head
# binders that happen to score least badly.
#
# The campaign cannot support all four modes equally. Over all 6,992,515
# designs: pi_stack_score is 0 for 98.8% (max 1.44), and interface H-bond count
# has median 1, p99 5, max 12. So "extensive H-bond network" is not available at
# any threshold; the honest backbone class is "contacts almost entirely to
# phosphate+sugar", which IS well populated.
FLOOR = {
    "backbone_hbond": lambda r: r["f_backbone"] >= 0.85,
    "pi_stack":       lambda r: r["pis"] >= 0.50,
    "head":           lambda r: r["f_loop"] >= 0.60,
    "groove":         lambda r: r["f_stem_base"] >= 0.40,
}
WANT = 4

OUT = ["design_id", "uuid8", "tpl_idx", "campaign", "stack", "seq_run",
       "length", "foot", "node_id", "coord_hash", "x", "y", "cluster_r10",
       "stage", "deg0", "iptm", "dockq", "rmsd", "hb", "pis",
       "c_total", "f_backbone", "f_base", "f_loop", "f_stem", "f_stem_base"]

# Scarcest mode first, so the thin classes get first claim on a shared design.
order = sorted(FLOOR, key=lambda m: q.apply(FLOOR[m], axis=1).sum())
used, picks = set(), []
for mode in order:
    elig = q[q.apply(FLOOR[mode], axis=1) & ~q["design_id"].isin(used)]
    elig = elig.sort_values(f"s_{mode}", ascending=False)
    # spread over fold clusters and stacks: four views of one cluster say
    # nothing about diversity
    sel = elig.drop_duplicates("cluster_r10").drop_duplicates("stack").head(WANT).copy()
    sel["mode"] = mode
    used.update(sel["design_id"])
    picks.append(sel)
    n_elig = int(q.apply(FLOOR[mode], axis=1).sum())
    flag = "" if len(sel) == WANT else f"   <-- SHORT: only {len(sel)} of {WANT}"
    print(f"\n=== {mode}  ({n_elig} eligible in pool){flag}")
    if len(sel):
        print(sel[["design_id", "stack", "stage", "cluster_r10", "iptm", "dockq",
                   "hb", "pis", "f_backbone", "f_base", "f_loop", "f_stem_base"]]
              .round(3).to_string(index=False))
    else:
        print("  NOTHING clears the floor. Do not fake this class.")

assert any(len(p) for p in picks), "no candidates survived; loosen the gate"
out = pd.concat([p for p in picks if len(p)])[["mode"] + OUT]
assert not out["design_id"].duplicated().any(), "a design was reused across modes"
out.to_csv("/cluster/work/wicky/lvalsecchi/fig02b_examples.csv", index=False)
print("\nwrote /cluster/work/wicky/lvalsecchi/fig02b_examples.csv", out.shape)
