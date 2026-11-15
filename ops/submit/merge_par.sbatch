#!/bin/bash
# merge_par.sbatch — one shard of the chunk-parallel merge.
# MP_MODE=index | shard | recombine
#SBATCH --job-name=mergepar
#SBATCH --partition=normal.4h
#SBATCH --time=03:59:00
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=8G
#SBATCH --requeue
#SBATCH --open-mode=append
#SBATCH --output=/cluster/work/wicky/lvalsecchi/rna_binder/ops/logs/mergepar_%A_%a.log
set -uo pipefail
RB=/cluster/work/wicky/lvalsecchi/rna_binder
CC=$RB/analyses/cross_campaign
cd "$CC"
MASTER=combined_clean_v2v9_fullmetrics_bf_final.csv
OUT=csy4_master_complete_par.csv
IX=$CC/merge_offsets.json
PARTS=$CC/merge_parts
N=2729
NSH=${MP_NSH:-24}

case "${MP_MODE:?}" in
index)
  python3 merge_index.py --master "$MASTER" --chunk-dir lrfix_chunks --n-chunks $N --out "$IX"
  ;;
shard)
  mkdir -p "$PARTS"
  T=${SLURM_ARRAY_TASK_ID:?}
  A=$(( T * N / NSH )); B=$(( (T+1) * N / NSH ))
  echo "[$(date '+%F %T')] shard $T: chunks [$A,$B) host=$(hostname)"
  python3 merge_shard.py --master "$MASTER" --out "$PARTS/part_$(printf %03d $T).csv" \
    --chunk-dir lrfix_chunks --sasa-dir 'sasa_sidecars_{be}' --fnat-dir fnatf_sidecars \
    --shard-index "$IX" --chunk-start $A --chunk-end $B \
    --stats-out "$PARTS/stats_$(printf %03d $T).json"
  ;;
recombine)
  n=$(ls "$PARTS"/part_*.csv 2>/dev/null | wc -l)
  [ "$n" -eq "$NSH" ] || { echo "ABORT: $n parts, expected $NSH"; exit 3; }
  ls -1 "$PARTS"/part_*.csv | sort | xargs cat > "$OUT.tmp"
  RM=$(wc -l < "$MASTER"); RO=$(wc -l < "$OUT.tmp")
  CM=$(head -1 "$MASTER" | tr ',' '\n' | wc -l); CO=$(head -1 "$OUT.tmp" | tr ',' '\n' | wc -l)
  echo "rows $RM -> $RO ; cols $CM -> $CO"
  [ "$RM" -eq "$RO" ] || { echo "ABORT: row count changed"; exit 4; }
  [ "$CO" -eq $(( CM + 38 )) ] || { echo "ABORT: expected +38 cols"; exit 5; }
  mv "$OUT.tmp" "$OUT"
  echo "DONE -> $OUT"
  ;;
esac
