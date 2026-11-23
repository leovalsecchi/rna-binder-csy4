"""Representatives that are THEMSELVES passing designs, not just backbones of one.

A backbone enters the gated set if ANY design on it passed, and the map elects one
representative name per node independently of that -- so only 16.7% of gated
representatives are passing designs. Here the candidate set is restricted to nodes
whose own elected name appears in passing_v6_final.csv.
"""
import glob, collections, os, subprocess, tarfile
import pyarrow.parquet as pq

D   = "/cluster/work/wicky/lvalsecchi/rna_binder/analyses/foldseek_map"
SSE = "/cluster/work/wicky/lvalsecchi/sse_bb"
OUT = "/cluster/work/wicky/lvalsecchi/cluster_reps_strict"
os.makedirs(OUT, exist_ok=True)
TOP = {"RFD3": 3, "BG": 9, "BD1": 0, "HD": 1}

passing_ids = set()
for line in open(f"{D}/passing_v6_final.csv"):
    if line.startswith("#"): continue
    did = line.split(",", 1)[0]
    if "len80" in did: passing_ids.add(did)

sse = {}
for f in glob.glob(f"{SSE}/shard_*.tsv"):
    with open(f) as fh:
        next(fh)
        for line in fh:
            p = line.rstrip("\n").split("\t")
            sse[p[0]] = (int(p[1]), int(p[2]), int(p[3]), int(p[4]), int(p[5]))

cl = pq.read_table(f"{D}/results_bbmap/design/clusters/fold/clusters.parquet",
                   columns=["node_id", "cluster_r05"]).to_pydict()
c_of = dict(zip(cl["node_id"], cl["cluster_r05"]))

cand = collections.defaultdict(list)
for line in open(f"{D}/results_bbmap/design/dedup/nodes.tsv"):
    if line.startswith("node_id"): continue
    nid, nm, _ = line.rstrip("\n").split("\t")
    g = nm.split("_")[1]
    if g in TOP and c_of.get(int(nid)) == TOP[g] and nm in passing_ids:
        cand[g].append(nm)

def he(nm):
    v = sse[nm]; return 100*v[1]/v[0], 100*v[2]/v[0], v[3], v[4]

CAMP_MEAN = {"RFD3": (87.8, 0.6), "BG": (37.0, 26.8), "BD1": (4.8, 38.4), "HD": (47.2, 2.7)}
picks = {}
print("candidates whose own design passed:")
for g, c in TOP.items():
    lst = cand[g]
    print(f"  {g:5s} cl{c}: {len(lst):,d}")
    if not lst: continue
    hs = [he(n)[0] for n in lst]; es = [he(n)[1] for n in lst]
    mh, me = sum(hs)/len(hs), sum(es)/len(es)
    nm = min(lst, key=lambda n: (he(n)[0]-mh)**2 + (he(n)[1]-me)**2)
    picks[g] = nm
    h, e, nh, ne = he(nm)
    ch, ce = CAMP_MEAN[g]
    print(f"        set mean H {mh:.1f} E {me:.1f}")
    print(f"        {nm}")
    print(f"        H {h:.1f}% E {e:.1f}%  {nh}H/{ne}E   vs campaign cluster H {ch:.1f} E {ce:.1f}"
          f"  (drift H {h-ch:+.1f} E {e-ce:+.1f})")
    assert nm in passing_ids, f"{nm} is not a passing design"

want = set(picks.values())
for f in glob.glob(f"{D}/results_bbmap/design/hashes/shard_*.tsv"):
    sh = os.path.basename(f)[:-4]
    names = [l.split("\t", 1)[0] for l in open(f) if l.split("\t", 1)[0] in want]
    if not names: continue
    open(f"{OUT}/.lst", "w").write("\n".join(f"{n}.pdb" for n in names) + "\n")
    subprocess.run(["tar", "-xzf", f"{D}/results_bbmap/design/shards/{sh}.tar.gz",
                    "-C", OUT, "-T", f"{OUT}/.lst"], check=True)
if os.path.exists(f"{OUT}/.lst"): os.remove(f"{OUT}/.lst")
with tarfile.open(f"{OUT}/reps_strict.tar", "w") as tf:
    for g in ("RFD3", "BG", "BD1", "HD"):
        if g in picks:
            tf.add(f"{OUT}/{picks[g]}.pdb", arcname=f"{g}_cl{TOP[g]}__{picks[g]}.pdb")
print("\nall four verified present in passing_v6_final.csv")
