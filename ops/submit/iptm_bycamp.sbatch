#!/bin/bash
#SBATCH --job-name=iptmcamp
#SBATCH --partition=normal.4h
#SBATCH --time=00:40:00
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=8G
#SBATCH --output=/cluster/work/wicky/lvalsecchi/rna_binder/ops/logs/iptmcamp_%j.log
apptainer exec --bind /cluster/work/wicky --bind /tmp --env MPLCONFIGDIR=/tmp \
  /cluster/work/wicky/containers/boltzgen.sif python3 - <<'PY'
import pandas as pd, numpy as np
P='/cluster/work/wicky/lvalsecchi/rna_binder/analyses/cross_campaign/combined_clean_v2v7_fullmetrics.csv'
be=['af3_msaoff_iptm','b2_iptm','esm_iptm']
d=pd.read_csv(P,usecols=be+['consensus_iptm','campaign'],low_memory=False)
mn3=d[be].min(axis=1); med3=d[be].median(axis=1)
d['is_min3']=np.isclose(d.consensus_iptm,mn3,atol=1e-6)
d['is_med3']=np.isclose(d.consensus_iptm,med3,atol=1e-6)
d['n_be']=d[be].notna().sum(axis=1)
print(f"{'campaign':10s}{'n':>10s}{'min-of-3':>10s}{'best2of3':>10s}{'all3 present':>14s}")
for c,g in d.groupby('campaign'):
    print(f'{str(c):10s}{len(g):>10,}{100*g.is_min3.mean():>9.1f}%{100*g.is_med3.mean():>9.1f}%{100*(g.n_be==3).mean():>13.1f}%')
PY
