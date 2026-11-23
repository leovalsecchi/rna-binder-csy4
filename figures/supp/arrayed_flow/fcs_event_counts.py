import glob,os,statistics as st
rows=[]
for f in sorted(glob.glob("/Volumes/wicky/expdata/flow/lvalsecchi/LV2_34_1/Matteo/**/*.fcs",recursive=True)):
    with open(f,"rb") as fh:
        h=fh.read(58)
        tb,te=int(h[10:18]),int(h[18:26])
        fh.seek(tb); txt=fh.read(te-tb+1).decode("latin-1")
    d=txt[0]; parts=txt[1:].split(d)
    kv={parts[i].upper():parts[i+1] for i in range(0,len(parts)-1,2)}
    rows.append((os.path.basename(f),int(kv.get("$TOT",-1)),kv.get("$DATE",""),kv.get("$CYT","")))
print("FCS files:",len(rows),"| instrument:",{r[3] for r in rows},"| dates:",sorted({r[2] for r in rows}))
tot=[r[1] for r in rows]
print(f"events per sample: min {min(tot):,}  median {int(st.median(tot)):,}  max {max(tot):,}")
for thr in (100000,):
    print(f"  at or above {thr:,}: {sum(1 for t in tot if t>=thr)}/{len(tot)}   below: {sum(1 for t in tot if t<thr)}")
print("  deciles:", [f"{int(x):,}" for x in st.quantiles(tot,n=10)])
low=sorted(rows,key=lambda r:r[1])[:8]
print("\nlowest eight:")
for n,t,_,_ in low: print(f"   {t:>8,}  {n[-18:]}")
print("\nexactly 100,000:",sum(1 for t in tot if t==100000))
