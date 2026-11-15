#!/bin/bash
# Measure how Foldseek search actually scales on this node set before
# committing the full run. Builds two sub-DBs from the already-built stage-2
# database and times one query chunk against each.
set -uo pipefail
FM=/cluster/work/wicky/lvalsecchi/rna_binder/analyses/foldseek_map
RESULTS=${RESULTS:-results_full}
cd "$FM" || exit 1

for N in 50000 150000; do
  sbatch --job-name="probe_$N" --partition=normal.4h --time=04:00:00 \
    --cpus-per-task=16 --mem-per-cpu=4G \
    --output="$FM/logs/probe_${N}_%j.log" \
    --wrap "bash -lc '
      set -e
      module load stack/2024-06 gcc/12.2.0 openmpi/4.1.6 foldseek/8-ef4e960
      OUT=$FM/$RESULTS/design
      WORK=\${TMPDIR:-/tmp}/probe_$N
      mkdir -p \$WORK; trap \"rm -rf \$WORK\" EXIT
      head -n $N \$OUT/db/uniq.lookup | cut -f1 > \$WORK/keys.txt
      head -n 2000 \$WORK/keys.txt > \$WORK/qkeys.txt
      for suffix in \"\" _ss _ca; do
        foldseek createsubdb \$WORK/keys.txt  \$OUT/db/uniq\$suffix \$WORK/sub\$suffix
        foldseek createsubdb \$WORK/qkeys.txt \$OUT/db/uniq\$suffix \$WORK/q\$suffix
      done
      echo \"=== N=$N : 2000 queries against a $N-target DB ===\"
      /usr/bin/time -v foldseek search \$WORK/q \$WORK/sub \$WORK/aln \$WORK/tmp \
        -c 0.8 --cov-mode 0 -e 10 --max-seqs 300 --alignment-type 1 -a --threads 16
    '"
done
echo "Two probes submitted. When both finish:"
echo "  grep -H 'Elapsed (wall clock)' $FM/logs/probe_*.log"
echo "p = log(t150/t50)/log(3). Full map-1 search ~ t50 * (500000/50000)^p * (500000/2000)."
