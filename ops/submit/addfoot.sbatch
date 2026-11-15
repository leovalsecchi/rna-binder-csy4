#!/bin/bash
# addfoot.sbatch - parallel append of foot_af3/foot_esm to the chunk-ordered master.
# AF_MODE=sidecar | index | shard | recombine
#SBATCH --job-name=addfoot
#SBATCH --partition=normal.4h
#SBATCH --time=03:59:00
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=16G
#SBATCH --requeue
#SBATCH --open-mode=append
#SBATCH --output=/cluster/work/wicky/lvalsecchi/rna_binder/ops/logs/addfoot_%A_%a.log
set -uo pipefail
RB=/cluster/work/wicky/lvalsecchi/rna_binder
CC=$RB/analyses/cross_campaign
cd "$CC"
MASTER=csy4_master_complete.csv
OUT=csy4_master_foot.csv
IX=$CC/foot_offsets.json
SC=$CC/foot_sidecars
PARTS=$CC/foot_parts
COLS=foot_af3,foot_esm
N=2729
NSH=${AF_NSH:-24}

case "${AF_MODE:?}" in
sidecar)
  python3 foot_sidecar_build.py --store "$CC/filter_store_v2v9" \
    --chunk-dir lrfix_chunks --out-dir "$SC" --n-chunks $N --cols "$COLS"
  ;;
index)
  python3 merge_index.py --master "$MASTER" --chunk-dir lrfix_chunks --n-chunks $N --out "$IX"
  ;;
shard)
  mkdir -p "$PARTS"
  T=${SLURM_ARRAY_TASK_ID:?}
  A=$(( T * N / NSH )); B=$(( (T+1) * N / NSH ))
  echo "[$(date '+%F %T')] shard $T: chunks [$A,$B) host=$(hostname)"
  python3 addcols_shard.py --master "$MASTER" --out "$PARTS/part_$(printf %03d $T).csv" \
    --chunk-dir lrfix_chunks --sidecar-dir "$SC" --cols "$COLS" \
    --index "$IX" --chunk-start $A --chunk-end $B
  ;;
recombine)
  n=$(ls "$PARTS"/part_*.csv 2>/dev/null | wc -l)
  [ "$n" -eq "$NSH" ] || { echo "ABORT: $n parts, expected $NSH"; exit 3; }
  ls -1 "$PARTS"/part_*.csv | sort | xargs cat > "$OUT.tmp"
  RM=$(wc -l < "$MASTER"); RO=$(wc -l < "$OUT.tmp")
  CM=$(head -1 "$MASTER" | tr ',' '\n' | wc -l); CO=$(head -1 "$OUT.tmp" | tr ',' '\n' | wc -l)
  echo "rows $RM -> $RO ; cols $CM -> $CO"
  [ "$RM" -eq "$RO" ] || { echo "ABORT: row count changed"; exit 4; }
  [ "$CO" -eq $(( CM + 2 )) ] || { echo "ABORT: expected +2 cols"; exit 5; }
  mv "$OUT.tmp" "$OUT"
  echo "DONE -> $OUT"
  ;;
esac
