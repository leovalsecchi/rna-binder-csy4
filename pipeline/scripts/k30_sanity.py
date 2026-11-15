import numpy as np, pandas as pd, collections, sys
FM="/cluster/work/wicky/lvalsecchi/rna_binder/analyses/foldseek_map"
LS="/cluster/work/wicky/lvalsecchi/rna_binder/analyses/library_selection"
raw=pd.read_csv(FM+"/passing_v6_final.csv",skiprows=2,low_memory=False)
lib=set(pd.read_parquet(LS+"/hier_k30_max_selections.parquet")
        .query('method == "HIER fold->epitope"').design_id)
L=raw[raw.design_id.isin(lib)]; PO=raw
print("cols:",len(raw.columns))
skip={"design_id","uuid8","stack","campaign","gen","invfold","design_cif","af3_cif","esm_cif"}
num=[c for c in raw.columns if c not in skip and pd.api.types.is_numeric_dtype(raw[c])]
print(f"\n{'metric':26}{'lib med':>10}{'pool med':>10}{'ratio':>8}{'lib p90':>10}  flag")
rows=[]
for c in num:
    a=pd.to_numeric(L[c],errors="coerce").replace([np.inf,-np.inf],np.nan).dropna()
    b=pd.to_numeric(PO[c],errors="coerce").replace([np.inf,-np.inf],np.nan).dropna()
    if len(a)<100 or len(b)<100: continue
    ma,mb=np.median(a),np.median(b)
    r=(ma/mb) if mb not in (0,) else np.nan
    rows.append((c,ma,mb,r,np.percentile(a,90)))
for c,ma,mb,r,p90 in sorted(rows,key=lambda x:-abs((x[3] or 1)-1) if x[3]==x[3] else 0):
    flag=""
    if r==r and (r<0.85 or r>1.18): flag="<-- shifted"
    if abs(ma)<1e-9 and abs(p90)<1e-9: flag="<-- all ~zero"
    print(f"{c:26}{ma:>10.3f}{mb:>10.3f}{(r if r==r else float('nan')):>8.2f}{p90:>10.3f}  {flag}")

print("\n"+"="*60); print("DUPLICATES & CYSTEINE"); print("="*60)
seqcol=[c for c in raw.columns if c.lower() in ("sequence","seq")]
if not seqcol:
    print("no sequence column in the export -> pulling from the v6 build metadata")
    b=pd.read_parquet(FM+"/build_v6/meta/designs.parquet")
    print("  build cols:",[c for c in b.columns][:12])
    sc=[c for c in b.columns if c.lower() in ("sequence","seq")]
    if sc:
        s=b.set_index("design_id")[sc[0]]
        L2=L.assign(sequence=s.reindex(L.design_id).values)
        PO2=PO.assign(sequence=s.reindex(PO.design_id).values)
    else: L2=PO2=None
else:
    L2,PO2=L.assign(sequence=L[seqcol[0]]),PO.assign(sequence=PO[seqcol[0]])
if L2 is not None and L2.sequence.notna().any():
    n=L2.sequence.notna().sum()
    print(f"  library sequences present : {n:,}")
    print(f"  unique sequences          : {L2.sequence.nunique():,}   duplicates: {n-L2.sequence.nunique():,}")
    d=L2.sequence.value_counts(); d=d[d>1]
    print(f"  sequences appearing >1x   : {len(d):,}   max copies: {d.max() if len(d) else 0}")
    cys=L2.sequence.dropna().str.contains("C")
    print(f"  designs containing CYS    : {cys.sum():,} ({100*cys.mean():.2f}%)")
    if cys.sum():
        by=L2.dropna(subset=["sequence"]).assign(c=cys.values).groupby("stack")["c"].agg(["sum","mean","size"])
        for s,r_ in by.iterrows():
            print(f"     {s:14} {int(r_['sum']):>5} / {int(r_['size']):>5}  ({100*r_['mean']:.1f}%)")
    pc=PO2.sequence.dropna().str.contains("C")
    print(f"  pool CYS rate             : {100*pc.mean():.2f}%")
