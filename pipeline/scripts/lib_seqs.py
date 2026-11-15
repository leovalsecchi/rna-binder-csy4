import csv, sys
import pandas as pd
csv.field_size_limit(10**7)
LS="/cluster/work/wicky/lvalsecchi/rna_binder/analyses/library_selection"
CC="/cluster/work/wicky/lvalsecchi/rna_binder/analyses/cross_campaign"
FM="/cluster/work/wicky/lvalsecchi/rna_binder/analyses/foldseek_map"
lib=set(pd.read_parquet(LS+"/hier_k30_max_selections.parquet")
        .query('method == "HIER fold->epitope"').design_id)
pool=set(pd.read_csv(FM+"/passing_v6_final.csv",skiprows=2,low_memory=False,
                     usecols=["design_id"]).design_id)
out={}
with open(CC+"/csy4_master_foot.csv",newline="") as f:
    r=csv.reader(f); h=[c.strip() for c in next(r)]
    i_id, i_seq = h.index("design_id"), h.index("sequence")
    for row in r:
        d=row[i_id]
        if d in pool: out[d]=row[i_seq]
print("sequences recovered: %d of %d pool designs"%(len(out),len(pool)))
pd.DataFrame({"design_id":list(out),"sequence":list(out.values())}).to_csv(LS+"/v6_pool_seqs.csv",index=False)
