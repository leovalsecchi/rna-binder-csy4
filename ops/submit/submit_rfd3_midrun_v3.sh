#!/usr/bin/env bash
# Submit RFD3 midrun_v3 campaign: 282 uuids distributed across GPU pools x walltime tiers.
set -uo pipefail
RB=/cluster/work/wicky/lvalsecchi/rna_binder
cd "$RB/code"

# pool: name | nodelist | gres | partition_family | counts(4h 24h 120h)
POOLS=(
  "4090|eu-g6-[001-080]|gpumem:24g|gpuhe|40 40 40"
  "titanquadro|eu-g3-[062-077,083-084],eu-lo-g3-[049-061]|gpumem:24g|gpuhe|26 26 26"
  "3090|eu-g4-[001-031]|gpumem:24g|gpuhe|20 20 20"
  "a100|eu-a65-[02-07],eu-ts-[01-02]|gpumem:40g|gpupr|8 8 8"
)
TIERS=("4h:04:00:00" "24h:23:30:00" "120h:119:00:00")

total_sub=0
for prow in "${POOLS[@]}"; do
  IFS='|' read -r name nodelist gres fam counts <<< "$prow"
  read -r c4 c24 c120 <<< "$counts"
  declare -A CNT=( [4h]=$c4 [24h]=$c24 [120h]=$c120 )
  for trow in "${TIERS[@]}"; do
    tier="${trow%%:*}"; wall="${trow#*:}"
    n="${CNT[$tier]}"
    [ "$n" -gt 0 ] || continue
    part="${fam}.${tier}"
    echo ">>> POOL=$name  PART=$part  TIME=$wall  NODES=$nodelist  GRES=$gres  N=$n"
    python3 scripts/launch.py configs/rfd3__multif/default.yaml \
      --label midrun_v3 --num-runs "$n" \
      --override slurm.gpus=1 \
      --override slurm.gres="$gres" \
      --override slurm.partition="$part" \
      --override slurm.time="$wall" \
      --override slurm.nodelist="$nodelist" 2>&1 \
      | grep -iE "Submitted|jobid|error" | sed 's/^/    /'
    total_sub=$((total_sub + n))
  done
done
echo "=== requested $total_sub RFD3 midrun_v3 uuids ==="
