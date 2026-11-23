"""Helix/strand per cluster for the gated set, in the CAMPAIGN cluster frame.

The passing stage has its own map (results_bbmap_pass) but its cluster ids come from
clustering the induced subgraph, so they do not correspond to campaign ids at all
(see the non-nesting result). To answer "how did the clusters change" the frame has
to be held fixed: same cluster ids, same axes, recomputed over surviving backbones.
Backbones DO nest 100%, so this is well defined.
"""
import glob, csv, collections, json
import pyarrow.parquet as pq

D   = "/cluster/work/wicky/lvalsecchi/rna_binder/analyses/foldseek_map"
SSE = "/cluster/work/wicky/lvalsecchi/sse_bb"
OUT = SSE

def bbkey(nm):
    p = nm.split("_"); i = p.index("len80")
    return f"{p[1]}_{p[i+1]}_{p[i+2]}"

sse = {}
for f in glob.glob(f"{SSE}/shard_*.tsv"):
    with open(f) as fh:
        next(fh)
        for line in fh:
            p = line.rstrip("\n").split("\t")
            sse[p[0]] = (int(p[1]), int(p[2]), int(p[3]))

cl = pq.read_table(f"{D}/results_bbmap/design/clusters/fold/clusters.parquet",
                   columns=["node_id", "cluster_r05"]).to_pydict()
c_of = dict(zip(cl["node_id"], cl["cluster_r05"]))

full = {}                                  # bbkey -> (gen, cluster, struct_name)
for line in open(f"{D}/results_bbmap/design/dedup/nodes.tsv"):
    if line.startswith("node_id"):
        continue
    nid, nm, nd = line.rstrip("\n").split("\t")
    c = c_of.get(int(nid))
    if c is None:
        continue
    full[bbkey(nm)] = (nm.split("_")[1], c, nm)

passing = set()
for line in open(f"{D}/results_bbmap_pass/design/dedup/nodes.tsv"):
    if line.startswith("node_id"):
        continue
    _, nm, _ = line.rstrip("\n").split("\t")
    passing.add(bbkey(nm))
print(f"passing backbones: {len(passing):,d}   campaign backbones: {len(full):,d}")
print(f"passing not found in campaign map: {len(passing - set(full)):,d}")

def tally(keys):
    agg = collections.defaultdict(lambda: {"n": 0, "h": 0.0, "e": 0.0})
    for k in keys:
        v = full.get(k)
        if not v:
            continue
        gen, c, nm = v
        s = sse.get(nm)
        if not s:
            continue
        a = agg[(gen, c)]
        a["n"] += 1
        a["h"] += 100 * s[1] / s[0]
        a["e"] += 100 * s[2] / s[0]
    return agg

camp = tally(full.keys())
pas  = tally(passing)

gen_tot_c = collections.Counter(); gen_tot_p = collections.Counter()
for (g, c), a in camp.items(): gen_tot_c[g] += a["n"]
for (g, c), a in pas.items():  gen_tot_p[g] += a["n"]
print("\ngen   campaign bb   passing bb   survival")
for g in ("RFD3", "BG", "BD1", "HD"):
    print(f"  {g:5s} {gen_tot_c[g]:9,d} {gen_tot_p[g]:11,d}   {100*gen_tot_p[g]/gen_tot_c[g]:6.2f}%")

MIN = 15
rows = []
for (g, c), a in pas.items():
    if a["n"] < MIN:
        continue
    cmp_ = camp[(g, c)]
    rows.append({"gen": g, "cluster": c, "n_bb": a["n"],
                 "helix": round(a["h"] / a["n"], 2), "strand": round(a["e"] / a["n"], 2),
                 "n_bb_campaign": cmp_["n"],
                 "helix_campaign": round(cmp_["h"] / cmp_["n"], 2),
                 "strand_campaign": round(cmp_["e"] / cmp_["n"], 2),
                 "pass_rate": round(100 * a["n"] / cmp_["n"], 3),
                 "share_campaign": round(100 * cmp_["n"] / gen_tot_c[g], 2),
                 "share_passing": round(100 * a["n"] / gen_tot_p[g], 2)})
rows.sort(key=lambda r: -r["n_bb"])
print(f"\nclusters with >= {MIN} passing backbones: {len(rows)}")

with open(f"{OUT}/cluster_sse_pass.tsv", "w") as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t")
    w.writeheader(); w.writerows(rows)
json.dump({"min_bb": MIN, "gen_tot_campaign": dict(gen_tot_c),
           "gen_tot_passing": dict(gen_tot_p), "rows": rows},
          open(f"{OUT}/cluster_sse_pass.json", "w"))

print("\ntop 12 passing clusters (campaign frame):")
print(f"{'gen':5s} {'cl':>5s} {'pass bb':>8s} {'camp bb':>8s} {'rate':>7s} "
      f"{'share c->p':>12s}   H camp->pass   E camp->pass")
for r in rows[:12]:
    print(f"{r['gen']:5s} {r['cluster']:>5d} {r['n_bb']:>8,d} {r['n_bb_campaign']:>8,d} "
          f"{r['pass_rate']:>6.2f}% {r['share_campaign']:>5.1f}->{r['share_passing']:<5.1f} "
          f"  {r['helix_campaign']:5.1f}->{r['helix']:<5.1f}  {r['strand_campaign']:5.1f}->{r['strand']:<5.1f}")

# did any cluster that was small campaign-wide become prominent after the gate?
print("\nbiggest risers by share (passing share minus campaign share):")
for r in sorted(rows, key=lambda r: -(r["share_passing"] - r["share_campaign"]))[:6]:
    print(f"  {r['gen']:5s} cl{r['cluster']:<5d} {r['share_campaign']:5.1f}% -> {r['share_passing']:5.1f}%  "
          f"(pass rate {r['pass_rate']:.2f}%,  H {r['helix']:.1f} E {r['strand']:.1f})")
print("\nbiggest fallers:")
for r in sorted(rows, key=lambda r: (r["share_passing"] - r["share_campaign"]))[:6]:
    print(f"  {r['gen']:5s} cl{r['cluster']:<5d} {r['share_campaign']:5.1f}% -> {r['share_passing']:5.1f}%  "
          f"(pass rate {r['pass_rate']:.2f}%,  H {r['helix']:.1f} E {r['strand']:.1f})")
