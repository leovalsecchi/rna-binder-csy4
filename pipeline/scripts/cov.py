import glob, csv, io, pandas as pd
PD="/cluster/work/wicky/lvalsecchi/rna_binder/analyses/library_selection/packed_designs"
LS="/cluster/work/wicky/lvalsecchi/rna_binder/analyses/library_selection"
fs=sorted(glob.glob(PD+"/validate_*.csv"))
print("files:",len(fs)," header:",open(fs[0]).readline().strip())
ok=set(); rows=0
for f in fs:
    txt=open(f,newline="").read().replace("\r\n","\n")
    for r in csv.DictReader(io.StringIO(txt)):
        rows+=1
        if (r.get("status") or "").strip()=="ok":
            k=r.get("design_id") or r.get("design") or r.get("name") or r.get("id")
            if k: ok.add(k.strip())
print("rows %d   ok-unique ids %d"%(rows,len(ok)))
L=pd.read_csv(LS+"/csy4_binder_library_12k_v6_k30.csv")
NEED={"RFD3+LiMPNN","BG+LiMPNN","BD1+LiMPNN"}
need=L[L["stack"].isin(NEED)]
have=need.design_id.isin(ok)
print("\n=== COVERAGE: 12k designs that need a packed file ===")
print("  need packed : %d"%len(need))
print("  covered     : %d (%.2f%%)"%(have.sum(),100*have.mean()))
print("  MISSING     : %d"%(~have).sum())
if (~have).sum():
    print(need[~have]["stack"].value_counts().to_string())
    need[~have][["design_id","stack"]].to_csv(LS+"/gap_v6_k30_missing.csv",index=False)
    print("  -> wrote gap_v6_k30_missing.csv")
nat=L[~L["stack"].isin(NEED)]
print("\n  native-sidechain designs (no packing needed): %d"%len(nat))
