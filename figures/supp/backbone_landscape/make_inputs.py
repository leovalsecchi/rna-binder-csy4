#!/usr/bin/env python3
"""Build backbone_landscape's inputs. RUNS ON EULER, not on the laptop.

    scp code/supp/backbone_landscape/make_inputs.py euler:/cluster/work/wicky/lvalsecchi/
    ssh euler '/cluster/work/wicky/lvalsecchi/binder_map_venv/bin/python3 \
        /cluster/work/wicky/lvalsecchi/make_inputs.py'          # union inputs
    ssh euler '/cluster/work/wicky/lvalsecchi/binder_map_venv/bin/python3 \
        /cluster/work/wicky/lvalsecchi/make_inputs.py --fig2b'  # the DRAWN table
    # then copy OUT/* into inputs/ (*.parquet into inputs/large/)

⚠️ `--fig2b` is the one the FIGURE reads. The rest builds the union-partition
tables, which the figure now only quotes (coverage, set quantities) and asserts.

A note on `cat file | ssh euler 'cat > f && ... &'`: that truncates `f` to zero
bytes, because the backgrounded job inherits the same stdin. Upload with scp, or
with `base64 | ssh euler 'base64 -d > f'`, and launch in a separate call.

Everything it reads lives only on the cluster: the UMAP embedding, the layer flags
and the k=30 cluster labels were never exported. Read-only apart from its own OUT.

THREE PARTITIONS, which are not the same object and do not share ids:
  union   results_union            526,467 nodes   the EMBEDDING + the stored layer flags
  k30     results_sw_k30_f50_mut    31,763 nodes   917 fold / 437 epitope clusters,
                                                   and both selection arms
  lib12k  results_lib12k             9,759 nodes   the delivered panel's own tree
A k30 or lib12k backbone is placed on the union embedding exactly the way
`scripts/stageU0_union.py` built `meta/layers.parquet`: that tree's canonical
struct_name -> coord_hash (the union's own hash shards) -> union node_id. The two
stored flags are then re-derived from scratch and asserted against the stored ones,
so a silent partition mismatch cannot pass.

k=30 r8.0 is the APPLIED map: 917 fold clusters. `results_v6` (= results_sw_k15_f50_mut)
gives 2,166, which is what the Methods text wrongly quotes.
"""
import glob
import json
import os
import sys

import numpy as np
import pandas as pd

FM = "/cluster/work/wicky/lvalsecchi/rna_binder/analyses/foldseek_map"
LS = "/cluster/work/wicky/lvalsecchi/rna_binder/analyses/library_selection"
U = f"{FM}/results_union/design"
K30 = f"{FM}/results_sw_k30_f50_mut/design"
BB = f"{FM}/results_bbmap/design"            # the MAIN-TEXT Fig 2b partition
OUT = "/cluster/work/wicky/lvalsecchi/_landscape_export"
ORDERED, COUNTERFACTUAL = "HIER fold->epitope", "top ipTM"
LIB12K = f"{LS}/csy4_binder_library_12k_v6_k30.csv"
APPLIED = f"{FM}/passing_v6_final.csv"       # the applied gate, 55,706 designs


def read_design_ids(csv):
    """A binder_filter export carries '#' comment lines before its header; the
    library-selection export does not. Count them rather than hardcoding 2."""
    skip = 0
    with open(csv) as fh:
        for line in fh:
            if not line.startswith("#"):
                break
            skip += 1
    return set(pd.read_csv(csv, skiprows=skip, usecols=["design_id"])["design_id"])


def fig2b():
    """Rebuild main-text Fig 2b's own scatter table, plus `in_counterfactual`.

    Columns and construction are `fig05_fold_clustering/export_fig.py`'s, which is
    what produced `inputs/large/fold_map_scatter.parquet`; the figure asserts that
    this table reproduces that one on all six of them. The added flag is the ipTM
    counterfactual - the `top ipTM` arm of hier_k30_max_selections.parquet - placed
    here by design_id, which is possible because `build_bbmap/meta/design_to_node`
    carries every one of the 6,992,515 designs.

    ⚠️ THE bbmap BACKBONE KEY IS NOT A COORD HASH. It is `<gen>_<uuid8>_<designN|
    recycleN>` with BG+BIF and BG+LiMPNN merged (build_bbmap/meta/build.json), so
    it is COARSER than the k30 partition's: the same applied gate gives 26,658
    backbones here against 31,763 there. Verified below, and the deficit is
    entirely BoltzGen's BIF merge plus HD's name collapse - RFD3 and BD1 agree to
    the unit on both layers.
    """
    os.makedirs(OUT, exist_ok=True)
    can = pd.read_csv(f"{BB}/dedup/canonical.tsv", sep="\t")
    nod = pd.read_csv(f"{BB}/dedup/nodes.tsv", sep="\t")
    em = pd.read_parquet(f"{BB}/embed/fold/embed.parquet")
    cl = pd.read_parquet(f"{BB}/clusters/fold/clusters.parquet")[["node_id",
                                                                  "cluster_r10"]]
    t = nod.merge(can, on="struct_name").merge(em, on="node_id").merge(cl, on="node_id")

    hp = set(pd.read_csv(f"{FM}/results_bbmap_pass/design/dedup/canonical.tsv",
                         sep="\t")["coord_hash"])
    hl = set(pd.read_csv(f"{FM}/results_bbmap_lib/design/dedup/canonical.tsv",
                         sep="\t")["coord_hash"])
    t["stage"] = np.where(t.coord_hash.isin(hl), 2,
                          np.where(t.coord_hash.isin(hp), 1, 0)).astype("int8")
    assert len(t) == 551_516 and (t.stage >= 1).sum() == 26_658 \
        and (t.stage == 2).sum() == 8_552, "bbmap cascade moved"

    # WHICH GATE? recomputed from the design ids, not trusted from the flags.
    meta = pd.read_parquet(f"{FM}/build_bbmap/meta/nodes.parquet")
    d2n = pd.read_parquet(f"{FM}/build_bbmap/meta/design_to_node.parquet")
    for ids, flag in ((read_design_ids(APPLIED), "in_passing"),
                      (read_design_ids(LIB12K), "in_12k")):
        assert set(d2n.loc[d2n.design_id.isin(ids), "node"]) == \
            set(meta.loc[meta[flag], "node"]), f"{flag} is not the gate it claims"
        print(f"{flag}: {len(ids):,d} designs reproduce the stored node set")

    sel = pd.read_parquet(f"{LS}/hier_k30_max_selections.parquet")
    cf_ids = set(sel.loc[sel.method == COUNTERFACTUAL, "design_id"])
    hit = d2n[d2n.design_id.isin(cf_ids)]
    assert len(hit) == len(cf_ids) == 12_000, (len(hit), len(cf_ids))
    cf_reps = set(meta.loc[meta.node.isin(set(hit.node)), "design_id"])
    t["in_counterfactual"] = t.struct_name.isin(cf_reps)
    print(f"ipTM arm: {len(hit):,d}/{len(cf_ids):,d} designs placed -> "
          f"{int(t.in_counterfactual.sum()):,d} drawn backbones")
    # the ordered arm placed the same way MUST land exactly on stage 2, which is
    # what proves the bridge rather than asserting it
    o_ids = set(sel.loc[sel.method == ORDERED, "design_id"])
    o_node = set(d2n.loc[d2n.design_id.isin(o_ids), "node"])
    o_rep = set(meta.loc[meta.node.isin(o_node), "design_id"])
    assert (t.struct_name.isin(o_rep) == (t.stage == 2)).all(), \
        "the ordered arm does not land on stage 2: the design_id bridge is wrong"
    print("bridge control: the ordered arm lands exactly on stage 2")

    STACK = {"RFD3_LiMPNN": "RFD3+LiMPNN", "BG_LiMPNN": "BG+LiMPNN",
             "BG_BIF": "BG+BIF", "BD1_LiMPNN": "BD1+LiMPNN",
             "HD_LiMPNN": "HD+LiMPNN"}
    pref = (t.struct_name.str.split("_len", n=1).str[0]
            .str.replace("csy4free_", "", regex=False))
    t["stack"] = pref.map(STACK).astype("category")
    e = pd.read_parquet(f"{BB}/graph/fold/edges_cluster.parquet", columns=["i", "j"])
    deg = pd.concat([e["i"], e["j"]]).value_counts()
    t["deg0"] = ~t.node_id.isin(deg.index)
    assert int(t.deg0.sum()) == 52_452, int(t.deg0.sum())

    out = t[["x", "y", "stage", "cluster_r10", "stack", "deg0",
             "in_counterfactual"]].copy()
    out["x"] = out.x.astype("float32")
    out["y"] = out.y.astype("float32")
    out["cluster_r10"] = out.cluster_r10.astype("int32")
    out.to_parquet(f"{OUT}/fig2b_layers.parquet", index=False, compression="zstd")
    print(pd.crosstab(out["stack"], out["stage"]).to_string())
    print(f"wrote {OUT}/fig2b_layers.parquet  {len(out):,d} rows")


def main():
    os.makedirs(OUT, exist_ok=True)
    nodes = pd.read_csv(f"{U}/dedup/nodes.tsv", sep="\t")
    emb = pd.read_parquet(f"{U}/embed/fold/embed.parquet")
    canon = pd.read_csv(f"{U}/dedup/canonical.tsv", sep="\t")
    layers = pd.read_parquet(f"{U}/meta/layers.parquet")

    sn2nid = pd.Series(nodes.node_id.values, index=nodes.struct_name)
    ch2nid = pd.Series(sn2nid.reindex(canon.struct_name).astype(int).values,
                       index=canon.coord_hash)
    layers["node_id"] = ch2nid.reindex(layers.coord_hash).astype(int).values

    # the union's own structure-name -> coord_hash map (stageU0's `hash_of`)
    hashes = pd.concat([pd.read_csv(f, sep="\t", usecols=["design_id", "coord_hash"])
                        for f in sorted(glob.glob(f"{U}/hashes/*.tsv"))],
                       ignore_index=True).drop_duplicates("design_id")
    name2ch = pd.Series(hashes.coord_hash.values, index=hashes.design_id)

    def to_union(struct_names, tag):
        ch = name2ch.reindex(pd.Index(struct_names))
        assert ch.notna().all(), f"{tag}: {int(ch.isna().sum())} names not hashed"
        nid = ch2nid.reindex(ch).astype(int)
        print(f"{tag}: {len(struct_names):,d} backbones -> {nid.nunique():,d} union nodes")
        return set(nid)

    k30n = pd.read_csv(f"{K30}/dedup/nodes.tsv", sep="\t")
    k2sn = pd.Series(k30n.struct_name.values, index=k30n.node_id)
    sel = pd.read_parquet(f"{LS}/hier_k30_max_selections.parquet")
    pool = pd.read_parquet(f"{LS}/v6k30_pool.parquet")
    arms = {m: set(sel.loc[sel.method == m, "node"].astype(int))
            for m in (ORDERED, COUNTERFACTUAL)}

    df = emb.merge(layers[["node_id", "in_sample500k", "in_passing60k", "in_12k"]],
                   on="node_id")
    df["in_pool_k30"] = df.node_id.isin(to_union(k30n.struct_name, "k30 gate pool"))
    df["in_ordered_k30"] = df.node_id.isin(
        to_union(k2sn.reindex(sorted(arms[ORDERED])), "ordered panel"))
    df["in_counterfactual"] = df.node_id.isin(
        to_union(k2sn.reindex(sorted(arms[COUNTERFACTUAL])), "ipTM counterfactual"))
    assert (df.in_passing60k == df.in_pool_k30).all(), "k30 pool != stored passing layer"
    assert (df.in_12k == df.in_ordered_k30).all(), "ordered arm != stored 12k layer"
    df.to_parquet(f"{OUT}/landscape_layers.parquet", index=False)

    # --- coverage and the set quantities, on the k30 map ---------------------
    fold = pd.read_parquet(f"{K30}/clusters/fold/clusters.parquet") \
        .set_index("node_id")["cluster_r80"]
    epi = pd.read_parquet(f"{K30}/clusters/epitope/clusters.parquet") \
        .set_index("node_id")["cluster_r80"]
    assert fold.nunique() == 917 and epi.nunique() == 437, "not the k=30 r8.0 map"

    sets, rows = {}, []
    for tag, nd, n_designs, iptm in (
            ("pool", set(pool.node.astype(int)), len(pool), pool.iptm),
            ("ordered", arms[ORDERED], 12_000, sel.loc[sel.method == ORDERED, "iptm"]),
            ("counterfactual", arms[COUNTERFACTUAL], 12_000,
             sel.loc[sel.method == COUNTERFACTUAL, "iptm"])):
        i = sorted(nd)
        f, e = fold.reindex(i).astype(int), epi.reindex(i).astype(int)
        sets[tag] = dict(nodes=nd, fold=set(f), epi=set(e), pair=set(zip(f, e)))
        rows.append(dict(arm=tag, designs=n_designs, backbones=len(nd),
                         folds=f.nunique(), epitopes=e.nunique(),
                         pairs=len(sets[tag]["pair"]),
                         mean_iptm=round(float(iptm.mean()), 4)))
    cov = pd.DataFrame(rows)
    cov.to_csv(f"{OUT}/landscape_coverage.csv", index=False)
    print(cov.to_string(index=False))

    q, o, c, pl = {}, sets["ordered"], sets["counterfactual"], sets["pool"]
    for key in ("nodes", "fold", "epi", "pair"):
        q.update({f"{key}_ordered": len(o[key]), f"{key}_cf": len(c[key]),
                  f"{key}_pool": len(pl[key]),
                  f"{key}_ordered_only": len(o[key] - c[key]),
                  f"{key}_cf_only": len(c[key] - o[key]),
                  f"{key}_both": len(o[key] & c[key]),
                  f"{key}_union": len(o[key] | c[key]),
                  f"{key}_net": len(o[key]) - len(c[key]),
                  f"{key}_pool_minus_ordered": len(pl[key] - o[key]),
                  f"{key}_pool_minus_cf": len(pl[key] - c[key])})
    q["fold_singletons_pool"] = int(
        (pool.groupby(pool.node.map(fold))["node"].nunique() == 1).sum())
    json.dump(q, open(f"{OUT}/landscape_quantities.json", "w"), indent=1, sort_keys=True)
    print(json.dumps(q, indent=1, sort_keys=True))
    print(f"wrote {OUT}")


if __name__ == "__main__":
    fig2b() if "--fig2b" in sys.argv else main()
