#!/bin/bash
#SBATCH --job-name=hdv8
#SBATCH --partition=normal.24h
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=8G
#SBATCH --output=/cluster/work/wicky/lvalsecchi/rna_binder/ops/logs/hdv8_%j.log
set -uo pipefail
cd /cluster/work/wicky/lvalsecchi/rna_binder/analyses/midrun_v8
python3 - <<'PY'
import csv, collections
csv.field_size_limit(10**7)
f=open('collected_raw.csv',newline=''); rd=csv.DictReader(f)
n=0; c=collections.Counter(); cys=0; seqs=collections.Counter()
for r in rd:
    if 'csy4free_HD_LiMPNN' not in r.get('design_id',''): continue
    n+=1
    a=r.get('has_af3_msaoff','') or r.get('af3_msaoff_cif_path','')
    b=r.get('has_b2','') or r.get('b2_cif_path','')
    e=r.get('has_esm','') or r.get('esm_cif_path','')
    ta=str(a).strip().lower() not in ('','0','false','nan')
    tb=str(b).strip().lower() not in ('','0','false','nan')
    te=str(e).strip().lower() not in ('','0','false','nan')
    c['af3' if ta else 'NO_af3']+=1
    c['b2'  if tb else 'NO_b2']+=1
    c['esm' if te else 'NO_esm']+=1
    if ta and tb and te: c['all3']+=1
    s=r.get('sequence','').strip()
    if 'C' in s: cys+=1
    if s: seqs[s]+=1
f.close()
print('midrun_v8 HD raw rows: %s\n'%format(n,','))
for k in ('af3','NO_af3','b2','NO_b2','esm','NO_esm','all3'):
    print('  %-8s %10s  %6.1f%%'%(k,format(c[k],','),100.0*c[k]/n if n else 0))
print('\n  cysteine-containing: %s (%.1f%%)'%(format(cys,','),100.0*cys/n if n else 0))
dup=sum(v-1 for v in seqs.values() if v>1)
print('  distinct sequences : %s   duplicate rows: %s'%(format(len(seqs),','),format(dup,',')))
print('\n  => all3 AND cys-free AND deduped is the clean gate; clean was 12,150')
PY
