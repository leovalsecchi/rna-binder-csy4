"""50 random designs from the k=30 hierarchical 12k: generator design + AF3 + ESM refold,
both refolds RNA-superposed onto the design so the target frame is shared."""
import os, sys, shutil
import numpy as np, pandas as pd
H="/cluster/work/wicky/lvalsecchi/rna_binder/code/scripts"
sys.path.insert(0,H); sys.path.insert(0,os.path.join(H,"..","src"))
import biotite.structure as struc
from biotite.structure.io import save_structure
from archive_read import resolve
from rna_binder.analysis.structural import load_structure, split_chains
from lrmsd_rna_aligned import RNA_BB, _bb, _match, _apply, MIN_REGISTER
FM="/cluster/work/wicky/lvalsecchi/rna_binder/analyses/foldseek_map"
LS="/cluster/work/wicky/lvalsecchi/rna_binder/analyses/library_selection"
OUT=FM+"/browse50"; N=50; SEED=7

raw=pd.read_csv(FM+"/passing_v6_final.csv",skiprows=2,low_memory=False)
raw["iptm_max"]=np.maximum(raw.iptm_af3,raw.iptm_esm)
raw["sc_rmsd"]=np.minimum(raw.rmsd_af3.fillna(np.inf),raw.rmsd_esm.fillna(np.inf))
lib=set(pd.read_parquet(LS+"/hier_k30_max_selections.parquet")
        .query('method == "HIER fold->epitope"').design_id)
L=raw[raw.design_id.isin(lib)]
rng=np.random.default_rng(SEED)
pick=L.iloc[rng.choice(len(L),N,replace=False)].sort_values("iptm_max",ascending=False)

def full(p):
    a=load_structure(resolve(p)); return a[a.element!="H"]

if os.path.isdir(OUT): shutil.rmtree(OUT)
os.makedirs(OUT)
man=[]
for i,r in enumerate(pick.itertuples()):
    tag=f"{i:02d}_{r.stack.replace('+','_')}_iptm{r.iptm_max:.3f}"
    d=os.path.join(OUT,tag); os.makedirs(d)
    try:
        D=full(r.design_cif)
        ref=_bb(split_chains(D)[1],RNA_BB)
        save_structure(os.path.join(d,f"1_design_{r.design_id}.cif"),D)
        oks=[]
        for lab,col in (("2_af3","af3_cif"),("3_esm","esm_cif")):
            M=full(getattr(r,col))
            rna=_bb(split_chains(M)[1],RNA_BB)
            ia,ib,reg=_match(ref,rna)
            if len(ia)<3 or reg<MIN_REGISTER:
                oks.append(f"{lab}:reg{reg:.2f}"); continue
            _,tr=struc.superimpose(ref[ia],rna[ib])
            save_structure(os.path.join(d,f"{lab}_{r.design_id}.cif"),_apply(tr,M))
            oks.append(f"{lab}:ok")
        man.append(dict(idx=i,folder=tag,design_id=r.design_id,stack=r.stack,campaign=r.campaign,
                        iptm_af3=round(r.iptm_af3,3),iptm_esm=round(r.iptm_esm,3),
                        iptm_max=round(r.iptm_max,3),sc_rmsd_best=round(r.sc_rmsd,2),
                        dockq_af3=round(r.dockq_af3,3),dockq_esm=round(r.dockq_esm,3),
                        status=";".join(oks)))
    except Exception as e:
        man.append(dict(idx=i,folder=tag,design_id=r.design_id,status=f"ERROR {e}"))
pd.DataFrame(man).to_csv(os.path.join(OUT,"manifest.csv"),index=False)
ok=sum(1 for m in man if m.get("status","").count("ok")==2)
print(f"{ok}/{N} complete triples -> {OUT}")
print(pd.DataFrame(man)["stack"].value_counts().to_string())
