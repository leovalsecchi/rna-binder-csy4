#!/bin/bash
#SBATCH --job-name=vfymerge
#SBATCH --partition=normal.4h
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=8G
#SBATCH --output=/cluster/work/wicky/lvalsecchi/rna_binder/ops/logs/vfymerge_%j.log
RB=/cluster/work/wicky/lvalsecchi/rna_binder
apptainer exec --bind /cluster/work/wicky --bind /tmp --env MPLCONFIGDIR=/tmp \
  /cluster/work/wicky/containers/boltzgen.sif python3 - <<'PY'
import pandas as pd, numpy as np
P='/cluster/work/wicky/lvalsecchi/rna_binder/analyses/cross_campaign/combined_clean_v2v8_fullmetrics.csv'
be=['af3_msaoff_iptm','b2_iptm','esm_iptm']
d=pd.read_csv(P,usecols=be+['iptm_min3','iptm_best2','mean_iptm3','consensus_iptm',
                            'campaign','source_table'],low_memory=False)
print('rows',f'{len(d):,}',flush=True)
print('by source  :',d.source_table.value_counts().to_dict())
print('by campaign:',{str(k):int(v) for k,v in d.campaign.value_counts().sort_index().items()})
mn=d[be].min(axis=1); md=d[be].median(axis=1); mu=d[be].mean(axis=1)
for c,ref,name in [('iptm_min3',mn,'min-of-3'),('iptm_best2',md,'best-2-of-3'),
                   ('mean_iptm3',mu,'mean-of-3'),('consensus_iptm',mn,'min-of-3')]:
    k=d[c].notna()&ref.notna()
    print(f'  {c:15s} == {name:12s} {100*np.isclose(d[c][k],ref[k],atol=1e-5).mean():6.2f}%'
          f'   median {d[c].median():.4f}   n={int(d[c].notna().sum()):,}')
print('violation consensus>best2:',int((d.consensus_iptm>d.iptm_best2+1e-9).sum()))
print('per-campaign consensus median:')
for k,g in d.groupby('campaign'):
    print(f'   {k:4s} n={len(g):>9,}  consensus {g.consensus_iptm.median():.4f}')
PY
