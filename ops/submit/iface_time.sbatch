#!/bin/bash
#SBATCH --job-name=ifacetime
#SBATCH --partition=normal.4h
#SBATCH --time=00:40:00
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=3G
#SBATCH --tmp=680000
#SBATCH --output=/cluster/work/wicky/lvalsecchi/rna_binder/ops/logs/ifacetime_%j.log
set -uo pipefail
RB=/cluster/work/wicky/lvalsecchi/rna_binder
cd $RB; export ARCH_READ_CACHE_MAX=3
# 400-design slice: enough to time, small enough to finish fast
head -401 $RB/analyses/cross_campaign/iface_v9/iface_manifest_af3.csv > /tmp/if400.csv
echo "rows: $(( $(wc -l < /tmp/if400.csv) - 1 ))"
t0=$(date +%s)
apptainer exec --bind /cluster/work/wicky --bind /tmp --env MPLCONFIGDIR=/tmp \
  /cluster/work/wicky/containers/boltzgen.sif \
  env PYTHONPATH=code/src:code/scripts python3 \
  code/scripts/contacts_sidecar.py --csv /tmp/if400.csv --path-col af3_msaoff_cif_path --jobs 8
t1=$(date +%s); echo "contacts: $((t1-t0))s"
apptainer exec --bind /cluster/work/wicky --bind /tmp --env MPLCONFIGDIR=/tmp \
  /cluster/work/wicky/containers/pyrosetta.sif \
  env PYTHONPATH=code/src:code/scripts python3 \
  code/scripts/interface_rosetta.py --csv /tmp/if400.csv --path-col af3_msaoff_cif_path --jobs 8 --with-hbonds
t2=$(date +%s); echo "rosetta: $((t2-t1))s  TOTAL: $((t2-t0))s for 400 designs"
echo "=> $(python3 -c "print(f'{400/max($t2-$t0,1):.1f}')") designs/s per 8-cpu task"
head -1 /tmp/if400.csv | tr ',' '\n' | tail -20
