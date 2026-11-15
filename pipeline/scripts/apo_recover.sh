#!/bin/bash
# apo_recover.sh <tag>   (tag = v8 | v9)
# Resubmit the ORPHANED apo chunks listed in apo_chunks_<tag>/.recover_idx.txt.
# Orphans (neither done nor queued) do NOT self-heal: array_packer_loop only repacks
# STILL-PENDING tasks. Routing copied verbatim from apo_submit_v9.sh (known-good):
#   gpu.*   -> --nodelist=eu-g3-[062-077,083-084], NO gpumem
#   gpuhe.* -> --gres=gpumem:24g, NO nodelist
#   walltime ROUTES the tier: 4h=04:00:00 24h=06:00:00 120h=1-06:00:00; bulk submits
#   remap to .24h so they go in --hold then scontrol-move while pending.
set -uo pipefail
RB=/cluster/work/wicky/lvalsecchi/rna_binder; cd "$RB"
TAG="${1:?usage: apo_recover.sh v8|v9}"
S=apo_fold_array_${TAG}.sh
IDX=apo_chunks_${TAG}/.recover_idx.txt
PIN="eu-g3-[062-077,083-084]"
[ -s "$IDX" ] || { echo "FATAL: no $IDX"; exit 1; }
[ -f "$S" ]   || { echo "FATAL: no $S"; exit 1; }

mapfile -t ALL < "$IDX"
N=${#ALL[@]}
echo "recovering $N orphan chunks for $TAG"

# split the index list across 6 tiers, weighted to each tier is memory-cap concurrency
# (48 G/task): 4h=36, 24h=19, 120h=11 per pool.
declare -a NAMES=(gpu.4h gpuhe.4h gpu.24h gpuhe.24h gpu.120h gpuhe.120h)
declare -a QOSS=(es_wicky/gpu/4 es_wicky/gpuhe/4 es_wicky/gpu/24 es_wicky/gpuhe/24 es_wicky/gpu/120 es_wicky/gpuhe/120)
declare -a TLS=(04:00:00 04:00:00 06:00:00 06:00:00 1-06:00:00 1-06:00:00)
declare -a WTS=(36 36 19 19 11 11)
TOT=0; for w in "${WTS[@]}"; do TOT=$((TOT+w)); done

start=0
for i in "${!NAMES[@]}"; do
  if [ "$i" -eq $(( ${#NAMES[@]} - 1 )) ]; then take=$(( N - start ))
  else take=$(( N * ${WTS[$i]} / TOT )); fi
  [ "$take" -gt 0 ] || continue
  slice=$(IFS=,; echo "${ALL[*]:$start:$take}")
  [ -n "$slice" ] || continue
  if [[ "${NAMES[$i]}" == gpuhe* ]]; then EXTRA=(--gres=gpumem:24g); else EXTRA=(--nodelist="$PIN"); fi
  j=$(sbatch --parsable --array="$slice" --partition="${NAMES[$i]}" --qos="${QOSS[$i]}" \
        --time="${TLS[$i]}" "${EXTRA[@]}" "$S")
  echo "  ${NAMES[$i]}  job $j  ${take} chunks"
  start=$(( start + take ))
done
echo "submitted $start / $N"
