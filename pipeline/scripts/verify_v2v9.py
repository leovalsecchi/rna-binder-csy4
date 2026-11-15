import sys, os, json
sys.path.insert(0, "/Users/leo/Documents/binder_filter")
import numpy as np
CC = "/cluster/work/wicky/lvalsecchi/rna_binder/analyses/cross_campaign"
os.environ["BF_STORE"] = CC + "/filter_store_v2v9"
from store import Store
s = Store(CC + "/filter_store_v2v9")
print(f"n = {s.n:,}  templates={len(s.templates)}  cols={len(s.numeric_cols())} numeric")
c, names = s.col("campaign"), s.vocabs["campaign"]
cnt = np.bincount(c, minlength=len(names))
print("\ncampaign breakdown:")
for i, nm in enumerate(names):
    if cnt[i]: print(f"  {nm:5s} {cnt[i]:>10,}")
print(f"  TOTAL {s.n:,}")

print("\nper-campaign medians (should be comparable across all 8):")
hdr = f"{'':5s}" + "".join(f"{k:>14s}" for k in ("dockq_af3","rmsd_af3","iptm_af3","pose_disagree"))
print(hdr)
for i, nm in enumerate(names):
    if not cnt[i]: continue
    m = c == i
    row = f"  {nm:3s}"
    for k in ("dockq_af3","rmsd_af3","iptm_af3","pose_disagree"):
        v = s.col(k)[m]; v = v[np.isfinite(v)]
        row += f"{np.median(v):>14.3f}" if v.size else f"{'—':>14s}"
    print(row)

# path reconstruction must work on BOTH sides of the join — template indices were remapped
print("\npath spot-check (2 rows per campaign):")
rng = np.random.default_rng(7)
bad = 0
for i, nm in enumerate(names):
    if not cnt[i]: continue
    idx = rng.choice(np.flatnonzero(c == i), 2, replace=False)
    for j in idx:
        p = s.paths(int(j))
        for key, ref in p.items():
            arc = ref.split("!")[0]
            full = "/cluster/work/wicky/lvalsecchi/rna_binder/" + arc
            if not os.path.exists(full):
                print(f"  MISSING {nm} row {j}: {arc}"); bad += 1
    print(f"  {nm}: 2 rows x 3 refs OK" if not bad else f"  {nm}: problems")
print(f"\npath check: {'ALL RESOLVED' if bad==0 else str(bad)+' MISSING'}")
