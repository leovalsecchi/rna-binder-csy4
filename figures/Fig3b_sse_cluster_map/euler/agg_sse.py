"""Per-cluster secondary-structure composition, per generator.

Every backbone representative (551,516) carries P-SEA H/E counts from the shard
pass, so cluster means are exact, not resampled. Emits one JSON for the plot page
and one TSV for reuse.
"""
import glob, json, collections, csv
import pyarrow.parquet as pq

D   = "/cluster/work/wicky/lvalsecchi/rna_binder/analyses/foldseek_map/results_bbmap/design"
SSE = "/cluster/work/wicky/lvalsecchi/sse_bb"
LIB = "/cluster/work/wicky/lvalsecchi/rna_binder/analyses/library_selection/csy4_binder_library_12k_v6_k30.csv"
OUT = "/cluster/work/wicky/lvalsecchi/sse_bb"
MIN_BB = 200          # plot clusters with at least this many backbones in a generator

# struct_name -> (n_res, n_helix, n_strand)
sse = {}
for f in sorted(glob.glob(f"{SSE}/shard_*.tsv")):
    with open(f) as fh:
        next(fh)
        for line in fh:
            p = line.rstrip("\n").split("\t")
            sse[p[0]] = (int(p[1]), int(p[2]), int(p[3]))
print(f"sse rows: {len(sse):,d}")

cl = pq.read_table(f"{D}/clusters/fold/clusters.parquet",
                   columns=["node_id", "cluster_r05", "cluster_r80"]).to_pydict()
c05 = dict(zip(cl["node_id"], cl["cluster_r05"]))

nodes = {}                                     # node_id -> (gen, struct_name, bbkey)
for line in open(f"{D}/dedup/nodes.tsv"):
    if line.startswith("node_id"):
        continue
    nid, nm, nd = line.rstrip("\n").split("\t")
    p = nm.split("_"); i = p.index("len80")
    nodes[int(nid)] = (p[1], nm, f"{p[1]}_{p[i+1]}_{p[i+2]}", int(nd))

# library backbones, so each point can show how much of it was ordered
libkeys = collections.Counter()
for r in csv.DictReader(open(LIB)):
    p = r["design_id"].split("_"); i = p.index("len80")
    libkeys[f"{p[1]}_{p[i+1]}_{p[i+2]}"] += 1

agg = collections.defaultdict(lambda: {"n": 0, "h": 0.0, "e": 0.0, "hh": 0.0, "ee": 0.0,
                                       "designs": 0, "lib_bb": 0, "lib_designs": 0})
nomatch = 0
for nid, (gen, nm, bbkey, ndes) in nodes.items():
    v = sse.get(nm)
    c = c05.get(nid)
    if v is None or c is None:
        nomatch += 1; continue
    nres, nh, ne = v
    h = 100 * nh / nres; e = 100 * ne / nres
    a = agg[(gen, c)]
    a["n"] += 1; a["h"] += h; a["e"] += e; a["hh"] += h * h; a["ee"] += e * e
    a["designs"] += ndes
    if bbkey in libkeys:
        a["lib_bb"] += 1
        a["lib_designs"] += libkeys[bbkey]
print(f"nodes with no sse or no cluster: {nomatch:,d}")

rows = []
for (gen, c), a in agg.items():
    if a["n"] < MIN_BB:
        continue
    n = a["n"]
    mh = a["h"] / n; me = a["e"] / n
    sdh = max(0.0, a["hh"] / n - mh * mh) ** 0.5
    sde = max(0.0, a["ee"] / n - me * me) ** 0.5
    rows.append({"gen": gen, "cluster": c, "n_bb": n, "n_designs": a["designs"],
                 "helix": round(mh, 2), "strand": round(me, 2),
                 "sd_helix": round(sdh, 2), "sd_strand": round(sde, 2),
                 "lib_bb": a["lib_bb"], "lib_designs": a["lib_designs"]})
rows.sort(key=lambda r: -r["n_bb"])
print(f"clusters plotted (>= {MIN_BB} bb in a generator): {len(rows)}")

tot = collections.Counter()
for nid, (gen, nm, bbkey, ndes) in nodes.items():
    tot[gen] += 1
cov = collections.Counter()
for r in rows:
    cov[r["gen"]] += r["n_bb"]
print("\ngen  plotted-bb / all-bb  (coverage of the generator's backbones)")
for g in ("RFD3", "BG", "BD1", "HD"):
    print(f"  {g:5s} {cov[g]:7,d} / {tot[g]:7,d} = {100*cov[g]/tot[g]:5.1f}%   "
          f"clusters {sum(1 for r in rows if r['gen']==g)}")

with open(f"{OUT}/cluster_sse.json", "w") as fh:
    json.dump({"min_bb": MIN_BB, "totals": dict(tot), "rows": rows}, fh)
with open(f"{OUT}/cluster_sse.tsv", "w") as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t")
    w.writeheader(); w.writerows(rows)
print(f"\nwrote {OUT}/cluster_sse.{{json,tsv}}")
print("\ntop 12 by backbones:")
for r in rows[:12]:
    print(f"  {r['gen']:5s} cl{r['cluster']:<5d} bb={r['n_bb']:6,d}  "
          f"H {r['helix']:5.1f}±{r['sd_helix']:4.1f}  E {r['strand']:5.1f}±{r['sd_strand']:4.1f}  "
          f"lib_bb={r['lib_bb']}")
