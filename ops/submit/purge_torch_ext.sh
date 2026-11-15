#!/bin/bash
#SBATCH --job-name=purge_torchext
#SBATCH --output=/cluster/work/wicky/lvalsecchi/rna_binder/ops/logs/hdclean/%x_%j.log
#SBATCH --error=/cluster/work/wicky/lvalsecchi/rna_binder/ops/logs/hdclean/%x_%j.log
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=2G
# Drop stale DeepSpeed evoformer_attn JIT build caches.
#
# containers/torch_ext_pertask/<slurm_jobid>/ is created per task by the HD JIT-lock fix
# so concurrent tasks don't fight over one ninja lock. Nothing ever reads another job's
# dir, so once the job is gone the dir is dead weight: ~23 inodes of .o/.so/.ninja files.
# Regenerable — the only cost of removing one is a recompile in a job that no longer exists.
#
# KEEP list = every JobId slurm currently knows for this user (running AND pending,
# per-array-task ids via scontrol, plus array parent ids), so a live task can never lose
# its build dir. Touches nothing outside torch_ext_pertask/<numeric>/.
set -uo pipefail
T=/cluster/work/wicky/lvalsecchi/rna_binder/containers/torch_ext_pertask
[ -d "$T" ] || { echo "no $T" >&2; exit 1; }
LIVE=$(mktemp)
scontrol show job -u lvalsecchi -o 2>/dev/null | grep -o 'JobId=[0-9]*' | cut -d= -f2 >> "$LIVE"
squeue -u lvalsecchi -h -o '%A %F' | tr ' ' '\n' | sed 's/_.*//' >> "$LIVE"
sort -u "$LIVE" -o "$LIVE"
echo "live job ids kept: $(wc -l < "$LIVE")"
n=0; k=0
for d in "$T"/*/; do
  b=$(basename "$d")
  case "$b" in ''|*[!0-9]*) echo "skip non-numeric $b"; continue ;; esac
  if grep -qx "$b" "$LIVE"; then k=$((k+1)); continue; fi
  rm -rf -- "$T/$b" && n=$((n+1))
  [ $((n % 2000)) -eq 0 ] && echo "  removed $n ..."
done
rm -f "$LIVE"
echo "removed $n stale dirs, kept $k live"
echo "remaining inodes under torch_ext_pertask: $(find "$T" | wc -l)"
