#!/usr/bin/env python
"""Per-amino-acid test, sequence composition AND interface composition, binder vs non."""
import sys, numpy as np, pandas as pd, warnings
warnings.filterwarnings("ignore")
AN="/Users/leo/Downloads/design29_742db4d3/an32"; sys.path.insert(0,AN)
import geom, biotite.structure as struc
AAs="ARNDQEGHILKMFPSTWYV"   # C is excluded from these designs by construction
CUT=4.5
def comp(path):
    try: a=geom.load(path)
    except Exception: return None
    prot,rna=geom.split(a)
    if prot is None or rna is None: return None
    pres=[(r,prot[prot.res_id==r]) for r in np.unique(prot.res_id)]
    s=np.array([geom.AA3.get(x.res_name[0],"X") for _,x in pres])
    cell=struc.CellList(rna,cell_size=CUT)
    hit=np.array([sum(len(cell.get_atoms(c,radius=CUT)) for c in x.coord)>0 for _,x in pres])
    if not hit.any(): return None
    d={f"seq_{a_}":float((s==a_).mean()) for a_ in AAs}
    d.update({f"if_{a_}":float((s[hit]==a_).mean()) for a_ in AAs})
    d["n_if"]=int(hit.sum()); return d
m=pd.read_pickle("/private/tmp/claude-501/-Users-leo/cbdfd1f4-189b-455f-b66a-20d7d6582507/scratchpad/best/m32sel.pkl")
S=m[m.flow_FC.notna()].copy().reset_index(drop=True)
R=pd.DataFrame([comp(f"{AN}/cif/{d}.{be}.cif") or {} for d,be in zip(S.design_id,S.sel_be)])
T=pd.concat([S[["design_id","stack","flow_FC","scr_sel"]].reset_index(drop=True),R],axis=1)
T["binder"]=(T.flow_FC>=1.2).astype(int)
lab=(T.binder==1).to_numpy()
rng=np.random.default_rng(0)
def auc(v,lb):
    a,b=v[lb],v[~lb]; a,b=a[np.isfinite(a)],b[np.isfinite(b)]
    return float(np.mean([(x>y)+0.5*(x==y) for x in a for y in b]))
def ci(v,lb,B=4000):
    ia,ib=np.where(lb)[0],np.where(~lb)[0]
    s=[auc(v[np.r_[rng.choice(ia,len(ia),True),rng.choice(ib,len(ib),True)]],
        np.r_[np.ones(len(ia),bool),np.zeros(len(ib),bool)]) for _ in range(B)]
    s=np.array([x for x in s if np.isfinite(x)]); return float(np.percentile(s,2.5)),float(np.percentile(s,97.5))
def strat(D,col,lb,P=20000,seed=3):
    r=np.random.default_rng(seed)
    rk=D.groupby("stack")[col].rank(pct=True).to_numpy(float)
    idx={s:np.where(D["stack"].values==s)[0] for s in D["stack"].unique()}
    k={s:int(lb[i].sum()) for s,i in idx.items()}
    obs=np.nanmean(rk[lb])-np.nanmean(rk[~lb]);nul=np.empty(P)
    for i in range(P):
        pm=np.zeros(len(D),bool)
        for s,ix in idx.items():
            if k[s]: pm[r.choice(ix,k[s],replace=False)]=True
        nul[i]=np.nanmean(rk[pm])-np.nanmean(rk[~pm])
    return ((np.abs(nul)>=abs(obs)).sum()+1)/(P+1)
def bh(p):
    p=np.asarray(p,float);o=np.argsort(p);q=np.empty_like(p)
    q[o]=np.minimum.accumulate((p[o]*len(p)/np.arange(1,len(p)+1))[::-1])[::-1]
    return np.clip(q,0,1)
rows=[]
for a_ in AAs:
    vs=T[f"seq_{a_}"].to_numpy(float); vi=T[f"if_{a_}"].to_numpy(float)
    lo,hi=ci(vs,lab); loi,hii=ci(vi,lab)
    rows.append(dict(aa=a_, auc_seq=auc(vs,lab), lo_seq=lo, hi_seq=hi, p_seq=strat(T,f"seq_{a_}",lab),
                     auc_if=auc(vi,lab), lo_if=loi, hi_if=hii, p_if=strat(T,f"if_{a_}",lab),
                     b_seq=np.nanmean(vs[lab]), n_seq=np.nanmean(vs[~lab]),
                     b_if=np.nanmean(vi[lab]), n_if=np.nanmean(vi[~lab])))
A=pd.DataFrame(rows)
A["q_seq"]=bh(A.p_seq); A["q_if"]=bh(A.p_if)
A=A.sort_values("auc_seq")
A.to_csv("asnfig/aa20.csv",index=False)
print(A[["aa","auc_seq","p_seq","q_seq","auc_if","p_if","q_if"]].round(3).to_string(index=False))
print(f"\nsequence: {(A.p_seq<0.05).sum()} of {len(A)} at P<0.05 (chance {0.05*len(A):.1f}), min q {A.q_seq.min():.3f}")
print(f"interface: {(A.p_if<0.05).sum()} of {len(A)} at P<0.05, min q {A.q_if.min():.3f}")
T.to_csv("asnfig/aa20_perdesign.csv",index=False)
