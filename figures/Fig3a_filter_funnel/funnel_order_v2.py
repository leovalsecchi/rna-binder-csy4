#!/usr/bin/env python3
"""The gate ACTUALLY applied to the full campaign, on all 6,992,515 designs.

Terms read verbatim from the filter JSON embedded in the header of
analyses/foldseek_map/passing_v6_final.csv, which records
"binder_filter export n=55706 store_n=6992515":

  iptm_max >= 0.6
  (rmsd_af3 <= 3)      OR  (rmsd_esm <= 3)          # ligand_RMSD
  (dockq_af3 >= 0.7)   OR  (dockq_esm >= 0.7)
  seq_run <= 5
  (footnb_010_af3 <= 5) AND (footnb_010_esm <= 5)   # hairpin-foot clash

All combined with AND, so the final count is order-independent; the funnel
order below is presentational. fig07c order: RMSD, DockQ, AA run, ipTM, foot clash.
"""
import numpy as np, os, json, pandas as pd
S="/cluster/work/wicky/lvalsecchi/rna_binder/analyses/cross_campaign/filter_store_v2v9_hdtop"
def col(n): return np.load(os.path.join(S,n+".npy"), allow_pickle=True)
VOC=json.load(open(os.path.join(S,"meta.json")))["vocabs"]["stack"]

stack=np.array([VOC[i] for i in col("stack")]); N=len(stack)
B=lambda m: np.nan_to_num(m, nan=False).astype(bool)

ipmax=np.fmax(col("iptm_af3"), col("iptm_esm"))   # filter used value="max"
TERMS=[("ligand_RMSD ≤ 3 Å",       B((col("rmsd_af3")<=3)|(col("rmsd_esm")<=3))),
       ("DockQ ≥ 0.7",         B((col("dockq_af3")>=0.7)|(col("dockq_esm")>=0.7))),
       ("Longest AA run ≤ 5",  B(col("seq_run")<=5)),
       ("ipTM (max) ≥ 0.6",    B(ipmax>=0.6)),
       ("Foot clash ≤ 5",      B((col("footnb_010_af3")<=5)&(col("footnb_010_esm")<=5)))]

print(f"store rows {N:,}\n")
print("MARGINAL (each term alone)")
for nm,m in TERMS: print(f"  {nm:22s} {m.sum():>9,}  {100*m.sum()/N:5.1f}%")

rows=[("All designs", np.ones(N,bool))]
cum=np.ones(N,bool)
for nm,m in TERMS:
    cum=cum&m; rows.append((nm,cum.copy()))

print("\nCUMULATIVE FUNNEL")
prev=N
for nm,m in rows:
    k=int(m.sum())
    print(f"  {nm:22s} {k:>9,}  {100*k/N:6.2f}%   kept {100*k/prev:5.1f}% of previous")
    prev=max(k,1)

final=int(rows[-1][1].sum())
print(f"\nFINAL {final:,}   expected 55,706   {'MATCH' if final==55706 else 'MISMATCH ('+str(final-55706)+')'}")

stacks=sorted(set(stack))
out=pd.DataFrame({"step":[r[0] for r in rows]})
for s in stacks: out[s]=[int((m&(stack==s)).sum()) for _,m in rows]
out["total"]=[int(m.sum()) for _,m in rows]
print("\nPER STACK (cumulative)"); print(out.to_string(index=False))
p="/cluster/work/wicky/lvalsecchi/funnel_real_gate_v2.csv"
out.to_csv(p,index=False); print("\nwrote",p)
