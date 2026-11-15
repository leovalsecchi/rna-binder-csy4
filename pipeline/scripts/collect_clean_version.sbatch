#!/bin/bash
# collect_clean_version.sbatch — collect + clean ONE midrun version's 5 stacks into
# analyses/<version>/collected_clean.csv. NO metrics, NO combine (clean+aggregate only).
# Mirrors v3_reaggregate.sbatch's collect->clean, adds the /local archive-extract guard
# (collect reads design CIFs from .tar.zst via archive_read -> same leak risk as metrics).
#   sbatch collect_clean_version.sbatch <midrun_vN>
#SBATCH --job-name=cc_version
#SBATCH --partition=normal.24h
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu=12G
#SBATCH --tmp=680000
#SBATCH --output=/cluster/work/wicky/lvalsecchi/rna_binder/analyses/_cc_%x_%j.log
set -euo pipefail
RB=/cluster/work/wicky/lvalsecchi/rna_binder
cd "$RB"
export PATH="$RB/code/bin:$PATH"           # bundled zstd for packed (.tar.zst) runs
export ARCH_READ_CACHE_MAX=3               # bound node-local archive extract cache
S=code/scripts
C=campaigns/csy4_crrna_free
V="${1:?usage: sbatch collect_clean_version.sbatch <midrun_vN>}"
J="${SLURM_CPUS_PER_TASK:-16}"
OUT="analyses/$V"
mkdir -p "$OUT"

DIRS=$(ls -d "$C"/*/"$V" 2>/dev/null || true)
[ -n "$DIRS" ] || { echo "FATAL: no run dirs for $V under $C/*/"; exit 1; }
echo "=== collect $V (tmpdir=${TMPDIR:-UNSET}) $(date) ==="
echo "$DIRS" | sed 's/^/  /'

python3 "$S/collect.py" $DIRS -o "$OUT/collected_raw.csv" --jobs "$J" --progress
echo "=== collected: $(( $(wc -l < "$OUT/collected_raw.csv") - 1 )) rows $(date) ==="

python3 "$S/clean.py" --input "$OUT/collected_raw.csv" --outdir "$OUT"
echo "=== clean DONE: $(( $(wc -l < "$OUT/collected_clean.csv") - 1 )) clean rows $(date) ==="
