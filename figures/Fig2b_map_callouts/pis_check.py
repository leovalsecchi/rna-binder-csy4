"""How much pi-stacking exists at all? And is the quality gate starving the
groove and pi-stack classes? Diagnostic only, writes nothing."""
import json, numpy as np, pandas as pd

STORE = ("/cluster/work/wicky/lvalsecchi/rna_binder/analyses/cross_campaign/"
         "filter_store_v2v9_hdtop")
meta = json.load(open(f"{STORE}/meta.json"))

def arr(name):
    return np.load(f"{STORE}/{name}.npy", mmap_mode="r")

pis = np.minimum(np.asarray(arr("pis_af3")), np.asarray(arr("pis_esm")))
hb  = np.minimum(np.asarray(arr("hb_af3")),  np.asarray(arr("hb_esm")))
ipt = np.minimum(np.asarray(arr("iptm_af3")), np.asarray(arr("iptm_esm")))
dq  = np.minimum(np.asarray(arr("dockq_af3")), np.asarray(arr("dockq_esm")))
rm  = np.maximum(np.asarray(arr("rmsd_af3")), np.asarray(arr("rmsd_esm")))

print(f"all {len(pis):,} designs")
print("  pi_stack_score: zero in %.1f%%, >0.5 in %.2f%%, >1 in %.2f%%, max %.2f"
      % (100*(pis <= 0).mean(), 100*(pis > 0.5).mean(), 100*(pis > 1).mean(), pis.max()))
print("  hb_interface:   median %.0f, p99 %.0f, max %.0f"
      % (np.median(hb), np.percentile(hb, 99), hb.max()))

for name, m in [("iptm>=.7 dockq>=.5 rmsd<=3", (ipt>=.7)&(dq>=.5)&(rm<=3)),
                ("iptm>=.6 dockq>=.4 rmsd<=4", (ipt>=.6)&(dq>=.4)&(rm<=4)),
                ("iptm>=.5 dockq>=.3",         (ipt>=.5)&(dq>=.3))]:
    s = pis[m]; h = hb[m]
    print(f"\n{name}: {m.sum():,} designs")
    print("  pis>0.5: %d   pis>1: %d   max pis %.2f" % ((s>0.5).sum(), (s>1).sum(), s.max() if len(s) else -1))
    print("  hb>=10: %d   hb>=15: %d   max hb %.0f" % ((h>=10).sum(), (h>=15).sum(), h.max() if len(h) else -1))
