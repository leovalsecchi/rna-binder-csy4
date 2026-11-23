#!/usr/bin/env python
"""Export the table behind the interface-Asn figure."""
import sys, glob, os, json, numpy as np, pandas as pd, warnings
warnings.filterwarnings("ignore")
AN="/Users/leo/Downloads/design29_742db4d3/an32"; sys.path.insert(0,AN)
import geom, biotite.structure as struc
CUT=4.5
def nat_split(a):
    prot=a[a.chain_id=="A"]; rna=a[a.chain_id=="B"]
    prot=prot[np.isin(prot.res_name,list(geom.AA3.keys()))]
    ok=[r for r in np.unique(prot.res_id) if (prot.atom_name[prot.res_id==r]=="CA").any()]
    prot=prot[np.isin(prot.res_id,ok)]
    rid=[r for r in np.unique(rna.res_id) if (rna.atom_name[rna.res_id==r]=="C1'").any()]
    return prot, rna[np.isin(rna.res_id,rid)]
def counts(path, natural=False):
    try: a=geom.load(path)
    except Exception: return None
    prot,rna=(nat_split(a) if natural else geom.split(a))
    if prot is None or rna is None or len(prot)==0: return None
    pres=[(r,prot[prot.res_id==r]) for r in np.unique(prot.res_id)]
    s=np.array([geom.AA3.get(x.res_name[0],"X") for _,x in pres])
    cell=struc.CellList(rna,cell_size=CUT)
    hit=np.array([sum(len(cell.get_atoms(c,radius=CUT)) for c in x.coord)>0 for _,x in pres])
    if not hit.any(): return None
    return dict(L=len(s), n_if=int(hit.sum()), if_N=float((s[hit]=="N").mean()),
                out_N=float((s[~hit]=="N").mean()) if (~hit).any() else np.nan,
                n_N_if=int((s[hit]=="N").sum()), all_N=float((s=="N").mean()))
OUT="/private/tmp/claude-501/-Users-leo/babadc09-87c1-4741-bc39-3ab6ba4f4349/scratchpad/asnfig"
os.makedirs(OUT,exist_ok=True)
M="/private/tmp/claude-501/-Users-leo/cbdfd1f4-189b-455f-b66a-20d7d6582507/scratchpad/best/m32sel.pkl"
m=pd.read_pickle(M)

# the 32, on the ligand_RMSD-selected refold AND on each backend
S=m[m.flow_FC.notna()].copy().reset_index(drop=True)
rows=[]
for d,be in zip(S.design_id,S.sel_be):
    sel=counts(f"{AN}/cif/{d}.{be}.cif") or {}
    a3=counts(f"{AN}/cif/{d}.af3.cif") or {}
    es=counts(f"{AN}/cif/{d}.esm.cif") or {}
    rows.append(dict(design_id=d, sel_be=be, **{f"sel_{k}":v for k,v in sel.items()},
                     **{f"af3_{k}":v for k,v in a3.items()}, **{f"esm_{k}":v for k,v in es.items()}))
R=pd.DataFrame(rows)
T=S[["design_id","stack","flow_FC","scr_sel"]].merge(R,on="design_id")
T["binder"]=(T.flow_FC>=1.2).astype(int)
T["uuid"]=T.design_id.str.extract(r"_([0-9a-f]{8})_")[0]
T.to_csv(f"{OUT}/asn32.csv",index=False)
print(f"the 32: {len(T)} rows, sel_if_N binders {T.sel_if_N[T.binder==1].mean():.4f} "
      f"non {T.sel_if_N[T.binder==0].mean():.4f}")

# a library reference sample on the selected backend
lib=m[m.flow_FC.isna()].sample(600,random_state=1)
rows=[counts(f"{AN}/cif/{d}.{be}.cif") for d,be in zip(lib.design_id,lib.sel_be)]
L=pd.DataFrame([r for r in rows if r]); L["stack"]=lib["stack"].to_numpy()[:len(L)]
L.to_csv(f"{OUT}/asn_lib600.csv",index=False)
print(f"library sample: {len(L)}  if_N {L.if_N.mean():.4f}  out_N {L.out_N.mean():.4f}  "
      f"enriched in {100*(L.if_N>L.out_N).mean():.0f}%")

# the 59 naturals + Csy4
rows=[]
for p in sorted(glob.glob("/Users/leo/rbp_rna_db/subset59/cif_clean/*.cif")):
    c=counts(p,natural=True)
    if c: rows.append({**c,"pdb":os.path.basename(p)[:-4]})
N=pd.DataFrame(rows); N.to_csv(f"{OUT}/asn_nat59.csv",index=False)
cs=counts("/Users/leo/hits32_browse_flowFC/csy4_ref.cif")
print(f"naturals: {len(N)}  if_N med {N.if_N.median():.4f}  out_N med {N.out_N.median():.4f}  "
      f"enriched in {100*(N.if_N>N.out_N).mean():.0f}%")
print(f"Csy4: {cs}")
json.dump({"csy4":cs,"lib_if_N":float(L.if_N.mean()),"lib_out_N":float(L.out_N.mean()),
           "lib_n":int(len(L))}, open(f"{OUT}/extra.json","w"))
