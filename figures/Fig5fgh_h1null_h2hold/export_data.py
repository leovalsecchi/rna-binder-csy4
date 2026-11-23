#!/usr/bin/env python
"""Export the table behind the H1-null / H2-partial figure."""
import numpy as np, pandas as pd, json
SP="/private/tmp/claude-501/-Users-leo/02af95f9-88eb-48c3-927b-2ee10ca465de/scratchpad"
OUT=f"{SP}/h1h2"
m=pd.read_pickle(f"{SP}/best/m32best.pkl")
N=pd.read_csv(f"{SP}/nat/nat59_scored.csv")
win=(m.best_be=="af3").to_numpy()
m["iptm_best"]=np.where(win,m.af3_msaoff_iptm,m.esm_iptm)
S=m[m.flow_FC.notna()].copy(); lab=(S.flow_FC>=1.2).to_numpy()
rng=np.random.default_rng(0)

def perm(v,lb,P=20000,seed=2):
    r=np.random.default_rng(seed)
    rk=pd.Series(v).rank().to_numpy(float); rk[~np.isfinite(v)]=np.nan
    obs=np.nanmean(rk[lb])-np.nanmean(rk[~lb])
    nul=np.array([(lambda p:np.nanmean(rk[p])-np.nanmean(rk[~p]))(r.permutation(lb)) for _ in range(P)])
    return ((np.abs(nul)>=abs(obs)).sum()+1)/(P+1)
def auc(v,lb):
    a,b=v[lb],v[~lb]; a,b=a[np.isfinite(a)],b[np.isfinite(b)]
    if len(a)<2 or len(b)<2: return np.nan
    return float(np.mean([(x>y)+0.5*(x==y) for x in a for y in b]))
def auc_ci(v,lb,B=5000):
    ia,ib=np.where(lb)[0],np.where(~lb)[0]
    s=[auc(v[np.r_[rng.choice(ia,len(ia),True),rng.choice(ib,len(ib),True)]],
           np.r_[np.ones(len(ia),bool),np.zeros(len(ib),bool)]) for _ in range(B)]
    s=np.array([x for x in s if np.isfinite(x)])
    return float(np.percentile(s,2.5)),float(np.percentile(s,97.5))

# ---- panel a: the whole scan, both scorings
EXCL={"hit","n","has_af3_msaoff","has_b2","has_esm","has_consensus","removed_fraction",
 "avg_ON/OFF_ratio","intact_count","removed_count","backbone_id","fold_cluster",
 "epitope_cluster","flow_FC","binder","hts_ratio","tpl_idx","fid_af3","fid_esm","fid_best",
 "fid_worst","fidL_af3","fidL_esm","iptm_best"}
def scan(pick):
    cols=[c for c in m.columns if pd.api.types.is_numeric_dtype(m[c]) and c not in EXCL
          and "b2" not in c and pick(c) and m[c].notna().sum()>11000 and m[c].nunique()>2]
    ps=[]
    for c in cols:
        v=S[c].to_numpy(float)
        if np.isfinite(v).sum()<28 or np.unique(v[np.isfinite(v)]).size<3: continue
        ps.append(perm(v,lab))
    return np.sort(np.array(ps))
p_best=scan(lambda c: c.startswith("best_") or c.endswith("_best") or not ("af3" in c or "esm" in c))
p_af3 =scan(lambda c: "af3" in c or not ("af3" in c or "esm" in c))
np.save(f"{OUT}/p_best.npy",p_best); np.save(f"{OUT}/p_af3.npy",p_af3)
print(f"scan: best {len(p_best)} features, {int((p_best<.05).sum())} at p<0.05 "
      f"(chance {0.05*len(p_best):.1f}) | af3 {len(p_af3)}, {int((p_af3<.05).sum())} "
      f"(chance {0.05*len(p_af3):.1f})")

# ---- panel b: the flagship H1 metrics
H1=[("ipTM","iptm_best",+1),("DockQ vs design","dockq_self_best",+1),
    ("ligand_RMSD (A)","dockq_LRMSD_best",-1),("Interface dSASA","dsasa_total_best",+1),
    ("Heavy contacts","n_heavy_contacts_4A_best",+1),
    ("Interface H-bonds","n_h_bonds_interface_best",+1),
    ("Shape compl.","sc_value_best",+1),("Packstat","packstat_best",+1),
    ("Buried unsat.","buns_vsasa_best",-1),("pi-stacking","pi_stack_score_best",+1)]
rows=[]
for nm,c,d in H1:
    if c not in S: print("MISSING",c); continue
    v=S[c].to_numpy(float); a=auc(v,lab); lo,hi=auc_ci(v,lab)
    rows.append(dict(label=nm,col=c,auc=a,lo=lo,hi=hi,p=perm(v,lab),
                     binder=np.nanmean(v[lab]),non=np.nanmean(v[~lab])))
B=pd.DataFrame(rows); B.to_csv(f"{OUT}/h1_metrics.csv",index=False)
print("\nH1 flagship metrics, P(binder > non-binder):")
print(B[["label","auc","lo","hi","p"]].round(3).to_string(index=False))

# ---- panel c: H2 legs under the three stresses (p values already computed, recorded here)
LEGS=[("Aromatic fraction (seq)","seq",0.025,0.025,0.036),
      ("Alanine fraction (seq)","seq",0.788,0.788,0.537),
      ("Aromatic at interface","struct",0.131,0.230,0.195),
      ("Acidic at interface","struct",0.024,0.157,0.348),
      ("Lys at interface","struct",0.346,0.729,0.729),
      ("Csy4-direction score","mixed",0.010,0.028,0.247),
      ("Acid-phosphate flag","struct",0.0077,0.603,0.590)]
C=pd.DataFrame(LEGS,columns=["label","kind","p_af3","p_best","p_nobif"])
C.to_csv(f"{OUT}/h2_legs.csv",index=False)
print("\nH2 legs:"); print(C.to_string(index=False))

# ---- panel d: the surviving feature, with the natural population
d=dict(lib=m.arom.to_numpy().tolist(),
       binders=S.arom[lab].to_numpy().tolist(), non=S.arom[~lab].to_numpy().tolist(),
       nat=N.bind_arom.dropna().to_numpy().tolist(), csy4=0.085,
       lib_mean=float(m.arom.mean()), nat_med=float(N.bind_arom.median()))
nb=S[~S["stack"].eq("BG+BIF")]; lnb=(nb.flow_FC>=1.2).to_numpy()
d["binders_nobif"]=nb.arom[lnb].to_numpy().tolist(); d["non_nobif"]=nb.arom[~lnb].to_numpy().tolist()
json.dump(d,open(f"{OUT}/arom.json","w"))
print(f"\narom: 12k {m.arom.mean():.4f} | non {np.mean(d['non']):.4f} | "
      f"binders {np.mean(d['binders']):.4f} | naturals {N.bind_arom.median():.4f} | Csy4 0.085")
print(f"  no-BIF: non {np.mean(d['non_nobif']):.4f} binders {np.mean(d['binders_nobif']):.4f} "
      f"(n {len(d['non_nobif'])}/{len(d['binders_nobif'])})")
