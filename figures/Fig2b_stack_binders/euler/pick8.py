import pandas as pd, numpy as np, json, os
B="/cluster/work/wicky/lvalsecchi/rna_binder/analyses/foldseek_map"
d=f"{B}/results_bbmap/design"
can=pd.read_csv(f"{d}/dedup/canonical.tsv", sep="\t")
nod=pd.read_csv(f"{d}/dedup/nodes.tsv", sep="\t")
em=pd.read_parquet(f"{d}/embed/fold/embed.parquet")
cl=pd.read_parquet(f"{d}/clusters/fold/clusters.parquet")[["node_id","cluster_r10"]]
t=nod.merge(can,on="struct_name").merge(em,on="node_id").merge(cl,on="node_id")
hl=set(pd.read_csv(f"{B}/results_bbmap_lib/design/dedup/canonical.tsv",sep="\t")["coord_hash"])
hp=set(pd.read_csv(f"{B}/results_bbmap_pass/design/dedup/canonical.tsv",sep="\t")["coord_hash"])
t["stage"]=np.where(t["coord_hash"].isin(hl),2,np.where(t["coord_hash"].isin(hp),1,0))
e=pd.read_parquet(f"{d}/graph/fold/edges_cluster.parquet",columns=["i","j"])
deg=pd.concat([e["i"],e["j"]]).value_counts()
t["deg0"]=~t["node_id"].isin(deg.index)
t=t[~t["deg0"]]
GEN={"RFD3_LiMPNN":"RFD3","BG_LiMPNN":"BG","BG_BIF":"BG","BD1_LiMPNN":"BD1","HD_LiMPNN":"HD"}
pref=t["struct_name"].str.split("_len",n=1).str[0].str.replace("csy4free_","",regex=False)
t["gen"]=pref.map(GEN)

S="/cluster/work/wicky/lvalsecchi/rna_binder/analyses/cross_campaign/filter_store_v2v9_hdtop"
ids=np.asarray(np.load(S+"/design_id.npy",mmap_mode="r"))
order=np.argsort(ids); sids=ids[order]
key=t["struct_name"].str.encode("ascii").values
pos=np.searchsorted(sids,key)
pos=np.clip(pos,0,len(sids)-1)
hit=sids[pos]==key
idx=order[pos]
print("store join:", hit.sum(), "/", len(t))
t=t[hit].copy(); idx=idx[hit]
def col(n): return np.load(S+f"/{n}.npy",mmap_mode="r")[idx]
for c in ["iptm_af3","iptm_esm","rmsd_af3","rmsd_esm","dockq_af3","dockq_esm",
          "seq_run","footnb_010_af3","footnb_010_esm","dsasa_af3","hb_af3",
          "packstat_af3","sc_value_af3","pis_af3"]:
    t[c]=np.asarray(col(c),dtype="float64")
t["uuid8"]=np.asarray(np.load(S+"/uuid8.npy",mmap_mode="r")[idx]).astype(str)
t["tpl_idx"]=np.asarray(np.load(S+"/tpl_idx.npy",mmap_mode="r")[idx])

# FIXED 2026-11-19. This gated on the BETTER of the two refolds - fmax on iptm
# and dockq, fmin on rmsd - which is the INVERSE of the campaign gate and of
# code/Fig2b_map_callouts/pick_examples.py, both of which take the worst of the
# two. A design that one backend placed well and the other placed 29 A away
# therefore passed, and two did: bg_a and bg_b are 29.34 A and 8.38 A from their
# design on AF3, the model the figure renders and measures on, against 1.60 A
# and 1.11 A on ESMFold2. Picking on best-of-two selects for backend DISAGREEMENT
# exactly where the figure needs agreement, because the call-out is a claim about
# the design, not about one predictor's opinion of it.
# A call-out must be reproduced by BOTH backends:
ip=np.fmin(t.iptm_af3,t.iptm_esm); dq=np.fmin(t.dockq_af3,t.dockq_esm)
rm=np.fmax(t.rmsd_af3,t.rmsd_esm); ft=np.fmax(t.footnb_010_af3,t.footnb_010_esm)
gate=(ip>=0.6)&(rm<=3)&(dq>=0.7)&(t.seq_run<=5)&(ft<=5)
t["gate"]=gate.fillna(False); t["ip"]=ip; t["dq"]=dq; t["rm"]=rm
lib=set(x.strip() for x in open("/cluster/work/wicky/lvalsecchi/rna_binder/analyses/library_selection/csy4_binder_library_12k_v6_k30.csv").read().split("\n")[1:] if x)
libids=set(l.split(",")[1] for l in lib if "," in l)
t["in12k"]=t["struct_name"].isin(libids)
g=t[t["gate"]].copy()
print("gate-passing representatives:", len(g), g.groupby("gen").size().to_dict())
# score: pose quality first, then interface
g["score"]=(g.dq-0.7)/0.3*2 + (g.ip-0.6)/0.4 + (3-g.rm)/3 + g.in12k*0.5 + g.dsasa_af3/3000
out=[]
for gen,sub in g.groupby("gen"):
    sub=sub.sort_values("score",ascending=False)
    picks=[sub.iloc[0]]
    for _,r in sub.iterrows():
        dmin=min(np.hypot(r.x-p.x, r.y-p.y) for p in picks)
        if dmin>4.0:
            picks.append(r); break
    else:
        for _,r in sub.iterrows():
            dmin=min(np.hypot(r.x-p.x, r.y-p.y) for p in picks)
            if dmin>2.0: picks.append(r); break
    out += picks[:2]
m=json.load(open(S+"/meta.json"))
res=[]
for r in out:
    p=m["templates"][int(r.tpl_idx)][2].replace("<U>",r.uuid8).replace("<ID>",r.struct_name)
    res.append(dict(gen=r.gen,design=r.struct_name,x=round(float(r.x),2),y=round(float(r.y),2),
                    cl=int(r.cluster_r10),in12k=bool(r.in12k),ip=round(float(r.ip),3),
                    dq=round(float(r.dq),3),rm=round(float(r.rm),3),
                    dsasa=round(float(r.dsasa_af3)),hb=float(r.hb_af3),
                    pis=round(float(r.pis_af3),3),pack=round(float(r.packstat_af3),3),
                    sc=round(float(r.sc_value_af3),3),path=p))
open("pick8.json","w").write(json.dumps(res,indent=1))
print("wrote", len(res), "picks")
