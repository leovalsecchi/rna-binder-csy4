import numpy as np, pandas as pd
FM="/cluster/work/wicky/lvalsecchi/rna_binder/analyses/foldseek_map"
LS="/cluster/work/wicky/lvalsecchi/rna_binder/analyses/library_selection"
RES=FM+"/results_sw_k30_f50_mut/design/"

raw=pd.read_csv(FM+"/passing_v6_final.csv",skiprows=2,low_memory=False)
seq=pd.read_csv(LS+"/v6_pool_seqs.csv").set_index("design_id")["sequence"]
lib=list(pd.read_parquet(LS+"/hier_k30_max_selections.parquet")
         .query('method == "HIER fold->epitope"').design_id)

nodes=pd.read_csv(RES+"dedup/nodes.tsv",sep="\t"); mem=pd.read_csv(RES+"dedup/members.tsv",sep="\t")
rep,m=[],[]
for sn,ids in zip(mem.struct_name,mem.design_ids):
    for i in str(ids).split(","): rep.append(sn); m.append(i)
nid=pd.Series(nodes.node_id.to_numpy(),index=nodes.struct_name)
raw["backbone_id"]=pd.Series(nid.reindex(pd.Index(rep)).to_numpy(),index=m).reindex(raw.design_id).to_numpy()
for ax,c in (("fold","fold_cluster"),("epitope","epitope_cluster")):
    cl=pd.read_parquet(RES+f"clusters/{ax}/clusters.parquet").set_index("node_id")
    raw[c]=cl["cluster_r80"].reindex(raw["backbone_id"]).to_numpy()

d=raw.set_index("design_id").reindex(lib)
out=pd.DataFrame({
    "design_id":lib,
    "sequence":seq.reindex(lib).values,
    "length":d.length.values,
    "stack":d["stack"].values,
    "generator":d.gen.values,
    "inverse_folder":d.invfold.values,
    "campaign":d.campaign.values,
    "backbone_id":d.backbone_id.values,
    "fold_cluster":d.fold_cluster.values,
    "epitope_cluster":d.epitope_cluster.values,
    "iptm_af3":d.iptm_af3.round(4).values,
    "iptm_esm":d.iptm_esm.round(4).values,
    "iptm_max":np.maximum(d.iptm_af3,d.iptm_esm).round(4).values,
    "sc_rmsd_af3":d.rmsd_af3.round(3).values,
    "sc_rmsd_esm":d.rmsd_esm.round(3).values,
    "dockq_af3":d.dockq_af3.round(4).values,
    "dockq_esm":d.dockq_esm.round(4).values,
    "dsasa_af3":d.dsasa_af3.round(1).values,
    "dsasa_esm":d.dsasa_esm.round(1).values,
    "n_hbonds_af3":d.hb_af3.values,
    "n_hbonds_esm":d.hb_esm.values,
}).sort_values("iptm_max",ascending=False).reset_index(drop=True)
out.insert(0,"n",range(1,len(out)+1))

print("=== VERIFICATION ===")
print("rows                    : %d"%len(out))
print("design_id unique        : %s (%d unique)"%(out.design_id.is_unique,out.design_id.nunique()))
print("design_id null          : %d"%out.design_id.isna().sum())
print("sequence null/empty     : %d"%(out.sequence.isna()|(out.sequence.astype(str).str.len()==0)).sum())
print("sequence unique         : %s (%d unique)"%(out.sequence.is_unique,out.sequence.nunique()))
print("sequence length         : min %d  max %d  (all==80: %s)"%(
      out.sequence.str.len().min(),out.sequence.str.len().max(),bool((out.sequence.str.len()==80).all())))
print("contains CYS            : %d"%out.sequence.str.contains("C").sum())
bad=set("".join(out.sequence.unique()))-set("ACDEFGHIKLMNPQRSTVWY")
print("non-standard AA chars   : %s"%(sorted(bad) if bad else "none"))
print("any other null cells    : %d"%out.isna().sum().sum())
p=LS+"/csy4_binder_library_12k_v6_k30.csv"
out.to_csv(p,index=False)
print("\nwrote %s  (%d cols)"%(p,out.shape[1]))
print(out.head(3).to_string(index=False,max_colwidth=22))
