#!/usr/bin/env bash
# basil_discussion — collect v3/v4/v5 SEPARATELY (3 parallel jobs), merge later.
set -uo pipefail
RB=/cluster/work/wicky/lvalsecchi/rna_binder
C=$RB/campaigns/csy4_crrna_free
OUT=$RB/analyses/basil_discussion
mkdir -p "$OUT"
SB="sbatch --partition=normal.24h --time=12:00:00 --cpus-per-task=32 --mem-per-cpu=3G"
COL="python3 $RB/code/scripts/collect.py"
ENV="cd $RB && export PATH=$RB/code/bin:\$PATH"

# v3 = HD + RFD3 only
$SB -J col_v3 --output=$OUT/_collect_v3_%j.log --wrap="$ENV && $COL \
  $C/rfd3__ligandmpnn/midrun_v3 \
  $C/halludesign__ligandmpnn/midrun_v3 \
  -o $OUT/collected_raw_v3.csv --jobs 32 --progress"

# v4 = all 5 stacks
$SB -J col_v4 --output=$OUT/_collect_v4_%j.log --wrap="$ENV && $COL \
  $C/rfd3__ligandmpnn/midrun_v4 \
  $C/boltzgen__boltzif/midrun_v4 \
  $C/boltzgen__ligandmpnn/midrun_v4 \
  $C/boltzdesign__ligandmpnn/midrun_v4 \
  $C/halludesign__ligandmpnn/midrun_v4 \
  -o $OUT/collected_raw_v4.csv --jobs 32 --progress"

# v5 = all 5 stacks
$SB -J col_v5 --output=$OUT/_collect_v5_%j.log --wrap="$ENV && $COL \
  $C/rfd3__ligandmpnn/midrun_v5 \
  $C/boltzgen__boltzif/midrun_v5 \
  $C/boltzgen__ligandmpnn/midrun_v5 \
  $C/boltzdesign__ligandmpnn/midrun_v5 \
  $C/halludesign__ligandmpnn/midrun_v5 \
  -o $OUT/collected_raw_v5.csv --jobs 32 --progress"
echo "submitted 3 per-campaign collect jobs"
