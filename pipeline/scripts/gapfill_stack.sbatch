#!/bin/bash
# gapfill_stack.sbatch — collect+clean ONE target stack fresh into its own gapfill dir.
# Re-aggregates the whole stack at current refold coverage so late-completing designs
# (Type-A gaps) are captured. Directives mirror collect_stack.sbatch (known-good).
#   sbatch gapfill_stack.sbatch <midrun_vN> <stackdir>
#SBATCH --job-name=gfstk
#SBATCH --partition=normal.4h
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=48
#SBATCH --mem-per-cpu=2G
#SBATCH --tmp=40000
#SBATCH --output=/cluster/work/wicky/lvalsecchi/rna_binder/analyses/_gf_%x_%j.log
set -euo pipefail
RB=/cluster/work/wicky/lvalsecchi/rna_binder
cd "$RB"
export PATH="$RB/code/bin:$PATH"
export ARCH_READ_CACHE_MAX=3
S=code/scripts
C=campaigns/csy4_crrna_free
V="${1:?usage: gapfill_stack.sbatch <midrun_vN> <stackdir>}"
ST="${2:?need stackdir}"
J="${SLURM_CPUS_PER_TASK:-48}"
D="$C/$ST/$V"
OUT="analyses/gapfill/${V}__${ST}"
mkdir -p "$OUT"
[ -d "$D" ] || { echo "FATAL: no run dir $D"; exit 1; }
echo "=== gapfill collect $V / $ST (tmp=${TMPDIR:-UNSET}, jobs=$J) $(date) ==="
python3 "$S/collect.py" "$D" -o "$OUT/collected_raw.csv" --jobs "$J" --progress
echo "=== raw: $(( $(wc -l < "$OUT/collected_raw.csv") - 1 )) rows; cleaning $(date) ==="
python3 "$S/clean.py" --input "$OUT/collected_raw.csv" --outdir "$OUT"
echo "=== clean: $(( $(wc -l < "$OUT/collected_clean.csv") - 1 )) rows $(date) ==="
echo "=== DONE $V/$ST $(date) ==="
