#!/bin/bash
#SBATCH --job-name=vrfchk4
#SBATCH --partition=normal.24h
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=4G
#SBATCH --output=/cluster/work/wicky/lvalsecchi/rna_binder/ops/logs/vrfchk4_%j.log
set -uo pipefail
cd /cluster/work/wicky/lvalsecchi/rna_binder/analyses/cross_campaign
python3 - <<'PY'
import csv
csv.field_size_limit(10**7)
IDS = {'csy4free_RFD3_LiMPNN_len80_893ee1bf_design10_seq4',
       'csy4free_RFD3_LiMPNN_len80_85e21c8b_design56_seq4',
       'csy4free_RFD3_LiMPNN_len80_e4994821_design50_seq6'}
COLS = ['dockq_self','dockq_fnat','dockq_std_self','dockq_consensus_n','fnat_factor',
        'dockq_self_af3','dockq_self_b2','dockq_self_esm']
for tag, path in (('PRE-MERGE ','combined_clean_v2v9_fullmetrics_bf_final.csv'),
                  ('POST-MERGE','csy4_master_foot.csv')):
    print('==== %s  %s' % (tag, path))
    with open(path, newline='') as f:
        rd = csv.reader(f); hdr = next(rd)
        ix = {c:i for i,c in enumerate(hdr)}
        did = ix['design_id']; hit = 0
        for row in rd:
            if row[did] in IDS:
                hit += 1
                print('  %s' % row[did])
                for c in COLS:
                    print('      %-20s %r' % (c, row[ix[c]] if c in ix else '<col absent>'))
                if hit == len(IDS): break
    print()
PY
