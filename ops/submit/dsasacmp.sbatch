#!/bin/bash
#SBATCH --job-name=dsasacmp
#SBATCH --partition=normal.24h
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=4G
#SBATCH --output=/cluster/work/wicky/lvalsecchi/rna_binder/ops/logs/dsasacmp_%j.log
set -uo pipefail
cd /cluster/work/wicky/lvalsecchi/rna_binder/analyses/cross_campaign
python3 - csy4_master_foot.csv <<'PY'
import csv, sys, math, collections
csv.field_size_limit(10**7)
BE=['af3','b2','esm']
acc={b:{'n':0,'sx':0.,'sy':0.,'sxx':0.,'syy':0.,'sxy':0.,'r':[]} for b in BE}
bystack=collections.defaultdict(lambda: {'n':0,'r':[]})
f=open(sys.argv[1],newline=''); rd=csv.reader(f); hdr=next(rd)
ix={c:i for i,c in enumerate(hdr)}
n=0
for row in rd:
    n+=1
    if n%2000000==0: print('  %s rows'%format(n,','),flush=True)
    for b in BE:
        a=row[ix['dsasa_total_%s'%b]].strip(); c=row[ix['dsasa_rna_total_%s'%b]].strip()
        if not a or not c: continue
        try: a=float(a); c=float(c)
        except ValueError: continue
        if a!=a or c!=c or c<=0 or a<=0: continue
        d=acc[b]; d['n']+=1; d['sx']+=a; d['sy']+=c; d['sxx']+=a*a; d['syy']+=c*c; d['sxy']+=a*c
        if len(d['r'])<3000000: d['r'].append(a/c)
        if b=='af3':
            s=bystack[row[ix['stack']]]; s['n']+=1
            if len(s['r'])<2000000: s['r'].append(a/c)
f.close()
def pear(d):
    n=d['n']
    num=d['sxy']-d['sx']*d['sy']/n
    den=math.sqrt((d['sxx']-d['sx']**2/n)*(d['syy']-d['sy']**2/n))
    return num/den if den else float('nan')
def q(v,p):
    v=sorted(v); return v[min(len(v)-1,int(p*len(v)))]
print('\nFULL MASTER  rows=%s\n'%format(n,','))
print('%-5s %12s %10s %10s %8s %8s %8s %8s'%('be','n','med total','med rna','ratio p10','p50','p90','pearson r'))
for b in BE:
    d=acc[b]
    if not d['n']: continue
    print('%-5s %12s %10.1f %10.1f %8.3f %8.3f %8.3f %8.3f'%(
        b,format(d['n'],','),d['sx']/d['n'],d['sy']/d['n'],
        q(d['r'],.10),q(d['r'],.50),q(d['r'],.90),pear(d)))
print('\nratio total/rna by stack (af3):')
for s in sorted(bystack):
    v=bystack[s]
    if not v['r']: continue
    print('   %-14s n=%10s  p10 %6.3f  p50 %6.3f  p90 %6.3f'%(
        s,format(v['n'],','),q(v['r'],.10),q(v['r'],.50),q(v['r'],.90)))
PY
