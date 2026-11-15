#!/bin/bash
#SBATCH --job-name=vrfv8top
#SBATCH --partition=normal.4h
#SBATCH --time=03:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=16G
#SBATCH --output=/cluster/work/wicky/lvalsecchi/rna_binder/ops/logs/vrfv8top_%j.log
set -euo pipefail
cd /cluster/work/wicky/lvalsecchi/rna_binder
python3 - <<'PY'
import csv, collections
csv.field_size_limit(10**7)
CC="/cluster/work/wicky/lvalsecchi/rna_binder/analyses/cross_campaign/"
OLD=CC+"csy4_master_foot.csv"; NEW=CC+"csy4_master_foot_v8top.csv"

# header must be byte-identical to the original
with open(OLD,"rb") as f: h_old=f.readline()
with open(NEW,"rb") as f: h_new=f.readline()
print("header byte-identical:", h_old==h_new)

PROBE=["dockq_self_af3","dockq_std_self_af3","dsasa_total_af3","foot_af3",
       "radius_of_gyration_ca_af3","n_h_bonds_interface_af3","dsasa_rna_total_af3",
       "iptm_min3","buns_vsasa_af3","dockq_self_b2"]
rows=0; ids=set(); cs=collections.Counter(); src=collections.Counter()
pop=collections.Counter(); pop_top=collections.Counter()
with open(NEW,newline="") as f:
    r=csv.reader(f); hdr=[c.strip() for c in next(r)]
    ix={c:i for i,c in enumerate(hdr)}
    i_st=[i for i,c in enumerate(hdr) if c=="source_table"]
    i_camp=ix["campaign"]; i_stack=ix["stack"]; i_did=ix["design_id"]
    ip={c:ix[c] for c in PROBE if c in ix}
    for row in r:
        rows+=1; ids.add(row[i_did])
        cs[(row[i_camp],row[i_stack])]+=1
        s=row[i_st[0]]; src[s]+=1
        istop = (s=="v8hdtopup")
        for c,i in ip.items():
            if row[i].strip() not in ("","nan"):
                pop[c]+=1
                if istop: pop_top[c]+=1
print(f"rows {rows:,}   unique design_id {len(ids):,}   duplicates {rows-len(ids):,}")
print("source_table:", dict(src))
ntop=src.get("v8hdtopup",0)
print()
print("appended block — populated fraction:")
for c in PROBE:
    if c in ip: print(f"  {c:28} {pop_top[c]:8,}/{ntop:,}  ({100*pop_top[c]/max(ntop,1):5.1f}%)")
print()
hd=sum(v for (cm,st),v in cs.items() if st=="HD+LiMPNN")
print(f"HD+LiMPNN total: {hd:,}")
print("per campaign x HD:")
for cm in sorted({k[0] for k in cs}):
    print(f"  {cm:6} total {sum(v for k,v in cs.items() if k[0]==cm):9,}   HD {cs.get((cm,'HD+LiMPNN'),0):9,}")
PY
