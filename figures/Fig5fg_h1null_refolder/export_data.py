#!/usr/bin/env python
"""Data for the two-panel figure: H1 is null on BOTH backends; interface measurements
are not reproducible between them.

Panel a uses a FIXED backend per mark - never the selected refold. Testing the better
of two after selecting on it is circular and is what produced the spurious ligand_RMSD
result in the first version."""
import json, sys, os
import numpy as np, pandas as pd
SP="/private/tmp/claude-501/-Users-leo/cbdfd1f4-189b-455f-b66a-20d7d6582507/scratchpad"
OUT=f"{SP}/fig2"
m=pd.read_pickle(f"{SP}/best/m32sel.pkl")
S=m[m.flow_FC.notna()].copy(); lab=(S.flow_FC>=1.2).to_numpy()
rng=np.random.default_rng(0)

def auc(v,lb):
    a,b=v[lb],v[~lb]; a,b=a[np.isfinite(a)],b[np.isfinite(b)]
    if len(a)<2 or len(b)<2: return np.nan
    return float(np.mean([(x>y)+0.5*(x==y) for x in a for y in b]))
def auc_ci(v,lb,B=5000):
    ia,ib=np.where(lb)[0],np.where(~lb)[0]
    out=[]
    for _ in range(B):
        idx=np.r_[rng.choice(ia,len(ia),True),rng.choice(ib,len(ib),True)]
        out.append(auc(v[idx],np.r_[np.ones(len(ia),bool),np.zeros(len(ib),bool)]))
    out=np.array([x for x in out if np.isfinite(x)])
    return float(np.percentile(out,2.5)),float(np.percentile(out,97.5))
def perm(v,lb,P=20000,seed=2):
    r=np.random.default_rng(seed)
    rk=pd.Series(v).rank().to_numpy(float); rk[~np.isfinite(v)]=np.nan
    obs=np.nanmean(rk[lb])-np.nanmean(rk[~lb])
    nul=np.array([(lambda p:np.nanmean(rk[p])-np.nanmean(rk[~p]))(r.permutation(lb)) for _ in range(P)])
    return ((np.abs(nul)>=abs(obs)).sum()+1)/(P+1)

# ---- panel a: ten flagship metrics, each on BOTH backends, no selection
H1=[("ipTM",            "af3_msaoff_iptm",            "esm_iptm"),
    ("DockQ vs design", "dockq_self_af3",             "dockq_self_esm"),
    ("ligand-RMSD",         "dockq_LRMSD_af3",            "dockq_LRMSD_esm"),
    ("Interface ΔSASA","dsasa_total_af3",        "dsasa_total_esm"),
    ("Heavy contacts",  "n_heavy_contacts_4A_af3",    "n_heavy_contacts_4A_esm"),
    ("Interface H-bonds","n_h_bonds_interface_af3",   "n_h_bonds_interface_esm"),
    ("Shape complementarity","sc_value_af3",          "sc_value_esm"),
    ("Packstat",        "packstat_af3",               "packstat_esm"),
    ("Buried unsat.",   "buns_vsasa_af3",             "buns_vsasa_esm"),
    ("π-stacking", "pi_stack_score_af3",         "pi_stack_score_esm")]
rows=[]
for nm,ca,ce in H1:
    for be,c in (("AF3",ca),("ESM",ce)):
        if c not in S: print("MISSING",c); continue
        v=S[c].to_numpy(float); lo,hi=auc_ci(v,lab)
        rows.append(dict(label=nm,backend=be,auc=auc(v,lab),lo=lo,hi=hi,p=perm(v,lab)))
A=pd.DataFrame(rows); A.to_csv(f"{OUT}/h1_both_backends.csv",index=False)
nsig=int((A.p<0.05).sum())
print(f"panel a: {len(A)} metric x backend marks, {nsig} at P<0.05 "
      f"(chance {0.05*len(A):.1f}); intervals covering 0.5: "
      f"{int(((A.lo<=0.5)&(A.hi>=0.5)).sum())}/{len(A)}")
print(A.round(3).to_string(index=False))

# ---- panel b: AF3 vs ESM reproducibility of each interface measurement, over the 12k
sp=lambda x,y:float(np.corrcoef(pd.Series(x).rank(),pd.Series(y).rank())[0,1])
FEAT=[("Arg at interface","if_R"),("Lys at interface","if_K"),
      ("Polar at interface","if_polar"),("Aromatic at interface","if_arom"),
      ("Tyr at interface","if_Y"),("Ala at interface","if_A"),
      ("Aliphatic at interface","if_aliph"),("Acidic at interface","if_neg")]
rows=[]
for nm,f in FEAT:
    ok=m[[f"af3_{f}",f"esm_{f}"]].dropna()
    rows.append(dict(label=nm,feat=f,rho=sp(ok[f"af3_{f}"],ok[f"esm_{f}"]),n=len(ok)))
B=pd.DataFrame(rows).sort_values("rho",ascending=False)
B.to_csv(f"{OUT}/backend_agreement.csv",index=False)
print("\npanel b: AF3 vs ESM agreement over the 12k")
print(B.round(3).to_string(index=False))

# the acid-phos flag concordance, for the annotation
a=(m.red_acid_phos_per_ifaceres_af3>0); e=(m.red_acid_phos_per_ifaceres_esm>0)
meta=dict(af3_rate=float(a.mean()),esm_rate=float(e.mean()),
          both=int((a&e).sum()),af3_only=int((a&~e).sum()),esm_only=int((e&~a).sum()),
          either=int((a|e).sum()),concord=float((a&e).sum()/(a|e).sum()),n=int(len(m)))
json.dump(meta,open(f"{OUT}/acidphos_concord.json","w"))
print(f"\nacid-phos flag: AF3 {100*meta['af3_rate']:.2f}%, ESM {100*meta['esm_rate']:.2f}%; "
      f"of {meta['either']} flagged by either, both {meta['both']} "
      f"({100*meta['concord']:.0f}%), AF3-only {meta['af3_only']}, ESM-only {meta['esm_only']}")
