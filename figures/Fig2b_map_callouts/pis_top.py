"""Top pi-stacking designs across ALL 6,992,515 designs, not just the one
elected representative per backbone. The earlier search looked only at elected
representatives, which is why it found none."""
import json, numpy as np, pandas as pd

STORE = ("/cluster/work/wicky/lvalsecchi/rna_binder/analyses/cross_campaign/"
         "filter_store_v2v9_hdtop")
meta = json.load(open(f"{STORE}/meta.json")); V = meta["vocabs"]
def A(n): return np.asarray(np.load(f"{STORE}/{n}.npy", mmap_mode="r"))

pis = np.minimum(A("pis_af3"), A("pis_esm"))
ipt = np.minimum(A("iptm_af3"), A("iptm_esm"))
dq  = np.minimum(A("dockq_af3"), A("dockq_esm"))
rm  = np.maximum(A("rmsd_af3"), A("rmsd_esm"))
hb  = np.minimum(A("hb_af3"), A("hb_esm"))

print(f"designs with pi_stack_score > 0 : {int((pis>0).sum()):,} "
      f"({100*(pis>0).mean():.2f}%)")
for thr in (0.25, 0.5, 0.75, 1.0, 1.25):
    print(f"  pis > {thr:4}: {int((pis>thr).sum()):,}")

for label, m in [("credible (iptm>=.7, dockq>=.5, rmsd<=3)", (ipt>=.7)&(dq>=.5)&(rm<=3)),
                 ("relaxed  (iptm>=.6, dockq>=.4, rmsd<=4)", (ipt>=.6)&(dq>=.4)&(rm<=4))]:
    sel = m & (pis > 0.5)
    idx = np.flatnonzero(sel)
    idx = idx[np.argsort(-pis[idx])][:12]
    print(f"\n== {label}: {int(sel.sum())} designs with pis>0.5")
    if not len(idx):
        continue
    did = np.char.decode(A("design_id")[idx], "ascii")
    out = pd.DataFrame({
        "design_id": did,
        "stack": np.array(V["stack"], dtype=object)[A("stack")[idx]],
        "campaign": np.array(V["campaign"], dtype=object)[A("campaign")[idx]],
        "uuid8": np.char.decode(A("uuid8")[idx], "ascii"),
        "tpl_idx": A("tpl_idx")[idx],
        "pis": pis[idx].round(3), "iptm": ipt[idx].round(3),
        "dockq": dq[idx].round(3), "rmsd": rm[idx].round(2), "hb": hb[idx],
    })
    print(out.to_string(index=False))
    out.to_csv(f"/cluster/work/wicky/lvalsecchi/pis_top_{label.split()[0]}.csv", index=False)
