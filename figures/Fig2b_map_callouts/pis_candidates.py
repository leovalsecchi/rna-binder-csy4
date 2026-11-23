"""Candidate pi-stacking designs to test for stacking on the TOP flipped-out
nucleotide U10. Store pi_stack_score is a single number with no position, so
the only way to know WHERE a design stacks is to measure its refold. This just
produces the extraction list."""
import json, numpy as np, pandas as pd
STORE=("/cluster/work/wicky/lvalsecchi/rna_binder/analyses/cross_campaign/"
       "filter_store_v2v9_hdtop")
V=json.load(open(f"{STORE}/meta.json"))["vocabs"]
def A(n): return np.asarray(np.load(f"{STORE}/{n}.npy", mmap_mode="r"))
pis=np.minimum(A("pis_af3"),A("pis_esm"))
ipt=np.minimum(A("iptm_af3"),A("iptm_esm"))
dq =np.minimum(A("dockq_af3"),A("dockq_esm"))
rm =np.maximum(A("rmsd_af3"),A("rmsd_esm"))
m=(pis>0.35)&(ipt>=0.65)&(dq>=0.50)&(rm<=3.5)
idx=np.flatnonzero(m); idx=idx[np.argsort(-pis[idx])][:48]
print(f"{m.sum()} designs pass; taking top {len(idx)} by pis")
TPL=json.load(open(f"{STORE}/meta.json"))["templates"]
rows=[]
for i in idx:
    did=A("design_id")[i].decode(); u=A("uuid8")[i].decode(); t=int(A("tpl_idx")[i])
    p=TPL[t][0]                                  # af3 refold template
    sd=p.split("/")[2]; mr=p.split("/")[3]
    rows.append(dict(sd=sd,mr=mr,u=u,did=did,pis=round(float(pis[i]),3),
                     iptm=round(float(ipt[i]),3),dockq=round(float(dq[i]),3),
                     stack=np.array(V["stack"],dtype=object)[A("stack")[i]]))
d=pd.DataFrame(rows)
print(d.groupby(["sd","mr"]).size().to_string())
print(d.head(10).to_string(index=False))
with open("/cluster/work/wicky/lvalsecchi/joblist4.txt","w") as f:
    for r in rows: f.write(f"{r['sd']}|{r['mr']}|{r['u']}|{r['did']}\n")
d.to_csv("/cluster/work/wicky/lvalsecchi/pis_candidates.csv",index=False)
print("wrote joblist4.txt", len(rows))
