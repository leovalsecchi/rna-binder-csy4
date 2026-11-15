#!/bin/bash
#SBATCH --job-name=vfy2
#SBATCH --partition=normal.4h
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=10G
#SBATCH --output=/cluster/work/wicky/lvalsecchi/rna_binder/ops/logs/vfy2_%j.log
apptainer exec --bind /cluster/work/wicky --bind /tmp --env MPLCONFIGDIR=/tmp \
  /cluster/work/wicky/containers/boltzgen.sif python3 - <<'PY'
import pandas as pd
P='/cluster/work/wicky/lvalsecchi/rna_binder/analyses/cross_campaign/combined_clean_v2v8_fullmetrics.csv'
cols=['source_table','campaign','iptm_min3','iptm_best2','consensus_iptm',
      'n_h_bonds_interface_opth','delta_unsat_hbonds_opth','buns_vsasa','buns_vsasa_opth',
      'n_h_bonds_interface_opth_af3','n_polar_h_moved_af3','n_h_bonds_interface']
d=pd.read_csv(P,usecols=cols,low_memory=False)
print('rows',f'{len(d):,}')
print('\ncoverage of the optH block, by source (was 0% for v8 before):')
for c in ['n_h_bonds_interface_opth','delta_unsat_hbonds_opth','buns_vsasa','buns_vsasa_opth',
          'n_h_bonds_interface_opth_af3','n_polar_h_moved_af3']:
    g=d.groupby('source_table')[c].apply(lambda s:100*s.notna().mean())
    print(f'  {c:32s} ' + '  '.join(f'{k}={v:6.2f}%' for k,v in g.items()))
print('\nipTM aggregates present:',
      {k:f'{100*v:.2f}%' for k,v in d.groupby("source_table").iptm_min3.apply(lambda s:s.notna().mean()).items()})
print('\noptH effect on H-bonds (median, all rows):')
print(f'  n_h_bonds_interface      {d.n_h_bonds_interface.median():.3f}')
print(f'  n_h_bonds_interface_opth {d.n_h_bonds_interface_opth.median():.3f}')
print('\nper-campaign optH H-bond median:')
for k,g in d.groupby('campaign'):
    print(f'   {k:4s} n={len(g):>9,}  opth {g.n_h_bonds_interface_opth.median():.3f}  raw {g.n_h_bonds_interface.median():.3f}')
PY
