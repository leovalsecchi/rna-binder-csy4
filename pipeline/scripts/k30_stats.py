import io, sys
import numpy as np, pandas as pd
FM="/cluster/work/wicky/lvalsecchi/rna_binder/analyses/foldseek_map"
LS="/cluster/work/wicky/lvalsecchi/rna_binder/analyses/library_selection"
RES=FM+"/results_sw_k30_f50_mut/design/"
buf=io.StringIO()
def P(s=""): print(s); print(s,file=buf)

raw=pd.read_csv(FM+"/passing_v6_final.csv",skiprows=2,low_memory=False)
raw["iptm_max"]=np.maximum(raw.iptm_af3,raw.iptm_esm)
raw["sc_rmsd"]=np.minimum(raw.rmsd_af3.fillna(np.inf),raw.rmsd_esm.fillna(np.inf))
raw["dockq"]=np.maximum(raw.dockq_af3,raw.dockq_esm)
raw["dsasa"]=np.maximum(raw.dsasa_af3,raw.dsasa_esm)
raw["hb"]=np.maximum(raw.hb_af3,raw.hb_esm)
raw["pis"]=np.maximum(raw.pis_af3,raw.pis_esm)
raw["foot"]=np.minimum(raw.foot_af3,raw.foot_esm)
raw["rog"]=np.minimum(raw.rog_af3,raw.rog_esm)
raw["fbase"]=np.maximum(raw.fbase_af3,raw.fbase_esm)

nodes=pd.read_csv(RES+"dedup/nodes.tsv",sep="\t"); mem=pd.read_csv(RES+"dedup/members.tsv",sep="\t")
rep,m=[],[]
for sn,ids in zip(mem.struct_name,mem.design_ids):
    for i in str(ids).split(","): rep.append(sn); m.append(i)
nid=pd.Series(nodes.node_id.to_numpy(),index=nodes.struct_name)
raw["node"]=pd.Series(nid.reindex(pd.Index(rep)).to_numpy(),index=m).reindex(raw.design_id).to_numpy()
for ax,c in (("fold","fold"),("epitope","epi")):
    cl=pd.read_parquet(RES+f"clusters/{ax}/clusters.parquet").set_index("node_id")
    raw[c]=cl["cluster_r80"].reindex(raw["node"]).to_numpy()

lib=set(pd.read_parquet(LS+"/hier_k30_max_selections.parquet")
        .query('method == "HIER fold->epitope"').design_id)
L=raw[raw.design_id.isin(lib)].copy(); PO=raw

P("="*74); P("LIBRARY STAT RUNDOWN  —  hierarchical fold→epitope, k=30, max-ipTM ranking"); P("="*74)
P(f"designs                 {len(L):>10,}   (from a pool of {len(PO):,})")
P(f"unique backbones        {L.node.nunique():>10,}   ({100*L.node.nunique()/len(L):.1f}% of designs)")
P(f"unique sequences        {L.sequence.nunique():>10,}" if "sequence" in L else "")
P(f"unique generation runs  {L.uuid8.nunique():>10,}   (pool {PO.uuid8.nunique():,})")
P(f"designs per backbone    {len(L)/L.node.nunique():>10.2f}   median {int(L.groupby('node').size().median())}, max {L.groupby('node').size().max()}")

P(); P("-"*74); P("COMPOSITION"); P("-"*74)
P(f"{'stack':16}{'library':>10}{'share':>9}{'pool share':>12}{'enrichment':>12}{'pool n':>9}")
pw=PO['stack'].value_counts(normalize=True)*100
for s,n in L['stack'].value_counts().items():
    P(f"{s:16}{n:>10,}{100*n/len(L):>8.1f}%{pw[s]:>11.1f}%{(100*n/len(L))/pw[s]:>11.2f}x{(PO['stack']==s).sum():>9,}")
P(); P(f"{'campaign':16}{'library':>10}{'share':>9}{'pool share':>12}")
cw=PO['campaign'].value_counts(normalize=True)*100
for c,n in sorted(L['campaign'].value_counts().items()):
    P(f"{c:16}{n:>10,}{100*n/len(L):>8.1f}%{cw[c]:>11.1f}%")

P(); P("-"*74); P("STRUCTURAL COVERAGE  (k=30 clustering, r80)"); P("-"*74)
NF,NE=PO.fold.nunique(),PO.epi.nunique(); NP=PO.groupby(["fold","epi"]).ngroups
bb=PO.groupby("fold")["node"].nunique(); sing=set(bb[bb==1].index)
P(f"fold clusters           {L.fold.nunique():>6,} / {NF:,}   ({100*L.fold.nunique()/NF:.1f}%)")
P(f"epitope clusters        {L.epi.nunique():>6,} / {NE:,}   ({100*L.epi.nunique()/NE:.1f}%)")
P(f"(fold,epitope) cells    {L.groupby(['fold','epi']).ngroups:>6,} / {NP:,}   ({100*L.groupby(['fold','epi']).ngroups/NP:.1f}%)")
P(f"fold singletons covered {len(set(L.fold)&sing):>6,} / {len(sing):,}")
occ=L.groupby(["fold","epi"]).size()
P(f"designs per occupied cell: median {occ.median():.0f}  p90 {occ.quantile(.9):.0f}  max {occ.max()}")

P(); P("-"*74); P("QUALITY  (library vs pool)"); P("-"*74)
P(f"{'metric':22}{'p05':>9}{'p25':>9}{'median':>9}{'p75':>9}{'p95':>9}{'  | pool median':>16}")
for lab,col,lo in (("ipTM max","iptm_max",0),("sc_RMSD best (Å)","sc_rmsd",1),
                   ("DockQ best","dockq",0),("ΔSASA max (Å²)","dsasa",0),
                   ("interface H-bonds","hb",0),("π-stack score","pis",0),
                   ("foot clearance (Å)","foot",0),("Rg (Å)","rog",0),("frac base contacts","fbase",0)):
    v=pd.to_numeric(L[col],errors="coerce").replace([np.inf,-np.inf],np.nan).dropna()
    pv=pd.to_numeric(PO[col],errors="coerce").replace([np.inf,-np.inf],np.nan).dropna()
    q=np.percentile(v,[5,25,50,75,95])
    fmt="%9.0f" if col in ("dsasa",) else "%9.2f" if col in ("sc_rmsd","hb","foot","rog") else "%9.3f"
    P(f"{lab:22}"+"".join(fmt%x for x in q)+f"{np.median(pv):>16.3f}")
P()
P(f"sc_RMSD <= 2 Å : {100*(L.sc_rmsd<=2).mean():.1f}%   <= 3 Å : {100*(L.sc_rmsd<=3).mean():.1f}%")
P(f"DockQ >= 0.7   : {100*(L.dockq>=0.7).mean():.1f}%   >= 0.8 : {100*(L.dockq>=0.8).mean():.1f}%")
P(f"ipTM  >= 0.8   : {100*(L.iptm_max>=0.8).mean():.1f}%   >= 0.9 : {100*(L.iptm_max>=0.9).mean():.1f}%")
P(f"ipTM  range    : {L.iptm_max.min():.3f} – {L.iptm_max.max():.3f}")
P(f"agreement: both refolders sc_RMSD<=3 : {100*((L.rmsd_af3<=3)&(L.rmsd_esm<=3)).mean():.1f}%")
open("/cluster/work/wicky/lvalsecchi/rna_binder/analyses/library_selection/k30_library_stats.txt","w").write(buf.getvalue())
