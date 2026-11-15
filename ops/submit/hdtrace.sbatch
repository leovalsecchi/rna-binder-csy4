#!/bin/bash
#SBATCH --job-name=hdtrace
#SBATCH --partition=normal.24h
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=4G
#SBATCH --output=/cluster/work/wicky/lvalsecchi/rna_binder/ops/logs/hdtrace_%j.log
set -uo pipefail
RB=/cluster/work/wicky/lvalsecchi/rna_binder
printf "%-14s %14s %14s %14s %10s\n" campaign raw clean fullmetrics "clean/raw"
for v in midrun_v2_all midrun_v3 midrun_v4 midrun_v5 midrun_v6 midrun_v7 midrun_v8 midrun_v9; do
  d=$RB/analyses/$v
  [ -d "$d" ] || continue
  r=0; c=0; m=0
  [ -f "$d/collected_raw.csv" ]                 && r=$(grep -c "csy4free_HD_LiMPNN" "$d/collected_raw.csv" 2>/dev/null || echo 0)
  [ -f "$d/collected_clean.csv" ]               && c=$(grep -c "csy4free_HD_LiMPNN" "$d/collected_clean.csv" 2>/dev/null || echo 0)
  [ -f "$d/collected_clean_fullmetrics.csv" ]   && m=$(grep -c "csy4free_HD_LiMPNN" "$d/collected_clean_fullmetrics.csv" 2>/dev/null || echo 0)
  pct="-"
  [ "$r" -gt 0 ] && pct=$(awk -v a=$c -v b=$r 'BEGIN{printf "%.1f%%", 100*a/b}')
  printf "%-14s %14s %14s %14s %10s\n" "$v" "$r" "$c" "$m" "$pct"
done
echo
echo "=== unique HD design_ids in the MASTER, and in each clean table ==="
echo "master:"
grep -o "csy4free_HD_LiMPNN_len80_[0-9a-f]\{8\}" $RB/analyses/cross_campaign/csy4_master_foot.csv | sort -u | wc -l
echo "  ^ distinct HD uuids present in master"
