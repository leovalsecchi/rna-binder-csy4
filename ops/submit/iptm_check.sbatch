#!/bin/bash
#SBATCH --job-name=iptmchk
#SBATCH --partition=normal.4h
#SBATCH --time=00:40:00
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=8G
#SBATCH --output=/cluster/work/wicky/lvalsecchi/rna_binder/ops/logs/iptmchk_%j.log
RB=/cluster/work/wicky/lvalsecchi/rna_binder
apptainer exec --bind /cluster/work/wicky --bind /tmp --env MPLCONFIGDIR=/tmp \
  /cluster/work/wicky/containers/boltzgen.sif python3 - <<'PY'
import pandas as pd, numpy as np
P='/cluster/work/wicky/lvalsecchi/rna_binder/analyses/cross_campaign/combined_clean_v2v7_fullmetrics.csv'
be=['af3_msaoff_iptm','b2_iptm','esm_iptm']
d=pd.read_csv(P,usecols=be+["consensus_iptm","mean_iptm","campaign"],low_memory=False)
print('rows',len(d),flush=True)
mn3=d[be].min(axis=1); med3=d[be].median(axis=1); mean3=d[be].mean(axis=1)
ab=d[['af3_msaoff_iptm','b2_iptm']]
tests=[('min(all 3)',mn3),('median3=best2of3',med3),('mean(all 3)',mean3),
       ('min(AF3,B2)',ab.min(axis=1)),('mean(AF3,B2)',ab.mean(axis=1))]
for col in ('consensus_iptm','mean_iptm'):
    print('--',col,'non-null',int(d[col].notna().sum()))
    for name,ref in tests:
        k=d[col].notna()&ref.notna()
        if k.sum(): print(f'    == {name:18s} {100*np.isclose(d[col][k],ref[k],atol=1e-6).mean():6.2f}%')
    print(f'    median {d[col].median():.4f}  vs min-of-3 median {mn3.median():.4f}')
PY
