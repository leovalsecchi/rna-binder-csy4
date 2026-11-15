#!/bin/bash
#SBATCH --job-name=inodecensus
#SBATCH --partition=normal.24h
#SBATCH --time=08:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=4G
#SBATCH --output=/cluster/work/wicky/lvalsecchi/rna_binder/ops/logs/inode_census_%j.log
# Full inode census. Runs on a compute node because a whole-tree walk over ~3.6M
# inodes takes far longer than a login-node session tolerates.
set -uo pipefail
BASE=/cluster/work/wicky/lvalsecchi
OUT=$BASE/_inode_census
mkdir -p "$OUT"

echo "[$(date '+%F %T')] quota before"
lfs quota -p 510551 /cluster/work/wicky | tail -2

echo "[$(date '+%F %T')] === depth-2 census of $BASE ==="
du --inodes -d 2 "$BASE" 2>/dev/null | sort -rn > "$OUT/depth2.txt"
head -40 "$OUT/depth2.txt"

echo "[$(date '+%F %T')] === depth-4 census of campaigns ==="
du --inodes -d 4 "$BASE/rna_binder/campaigns" 2>/dev/null | sort -rn > "$OUT/campaigns_d4.txt"
head -40 "$OUT/campaigns_d4.txt"

echo "[$(date '+%F %T')] === depth-3 census of analyses ==="
du --inodes -d 3 "$BASE/rna_binder/analyses" 2>/dev/null | sort -rn > "$OUT/analyses_d3.txt"
head -30 "$OUT/analyses_d3.txt"

echo "[$(date '+%F %T')] DONE — full tables in $OUT/"
