#!/bin/bash
# collect_stack.sbatch — collect ONE stack of a version into collected_raw_<stack>.csv
# (no clean). For resharding a large version's collect across parallel per-stack jobs so
# each single-process merge tail is ~5x smaller and the 5 run concurrently.
#   sbatch collect_stack.sbatch <midrun_vN> <stackdir>
#SBATCH --job-name=ccstk
#SBATCH --partition=normal.4h
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=48
#SBATCH --mem-per-cpu=2G
#SBATCH --tmp=40000
#SBATCH --output=/cluster/work/wicky/lvalsecchi/rna_binder/analyses/_cc_%x_%j.log
set -euo pipefail
RB=/cluster/work/wicky/lvalsecchi/rna_binder
cd "$RB"
export PATH="$RB/code/bin:$PATH"
export ARCH_READ_CACHE_MAX=3
S=code/scripts
C=campaigns/csy4_crrna_free
V="${1:?usage: collect_stack.sbatch <midrun_vN> <stackdir>}"
ST="${2:?need stackdir}"
J="${SLURM_CPUS_PER_TASK:-48}"
OUT="analyses/$V"
mkdir -p "$OUT"
D="$C/$ST/$V"
[ -d "$D" ] || { echo "FATAL: no run dir $D"; exit 1; }
echo "=== collect $V / $ST (tmpdir=${TMPDIR:-UNSET}, jobs=$J) $(date) ==="
python3 "$S/collect.py" "$D" -o "$OUT/collected_raw_${ST}.csv" --jobs "$J" --progress
echo "=== DONE $ST: $(( $(wc -l < "$OUT/collected_raw_${ST}.csv") - 1 )) rows $(date) ==="
