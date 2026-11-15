#!/usr/bin/env bash
# submit_bif_fill.sh — fresh BoltzGen (Route-B multif) to fill the BG+BIF gap (39,383
# cysteine designs excluded). Config cys-free (avoid-C + omit-C). NEW backbones.
# Spread: public gpu.24h (off-fairshare) + gpuhe.24h + gpuhe.4h(reduced num_designs).
set -uo pipefail
RB=/cluster/work/wicky/lvalsecchi/rna_binder
cd "$RB/code"
LABEL=bif_fill_2026-07-03
CFG=configs/boltzgen__multif/default.yaml

# name | nodelist | gres | partition | wall | num_designs | n
# 120h lanes use time>24h to force the .120h tier (router: wall<=24h -> .24h);
# gpuhe.4h uses num_designs=40 so BoltzGen (~2.3h) fits inside the 4h cap.
POOLS=(
  "public24|eu-g3-[062-070]|gpu:1|gpu.24h|23:30:00|60|18"
  "public120|eu-g3-[071-077,083-084]|gpu:1|gpu.120h|48:00:00|60|18"
  "owner24|eu-g6-[001-080]|gpumem:24g|gpuhe.24h|23:30:00|60|12"
  "owner120|eu-g6-[001-080]|gpumem:24g|gpuhe.120h|48:00:00|60|12"
  "owner4|eu-g6-[001-080]|gpumem:24g|gpuhe.4h|03:50:00|40|15"
)
ONLY="${ONLY:-}"   # optional: run only this lane name (e.g. ONLY=public24 for a probe)
CAP="${CAP:-}"     # optional: override n for the run (probe with CAP=1)

total=0
for prow in "${POOLS[@]}"; do
  IFS='|' read -r name nodelist gres part wall nd n <<< "$prow"
  [ -n "$ONLY" ] && [ "$name" != "$ONLY" ] && continue
  [ -n "$CAP" ] && n="$CAP"
  [ "$n" -gt 0 ] || continue
  echo ">>> LANE=$name  PART=$part  TIME=$wall  num_designs=$nd  NODES=$nodelist  GRES=$gres  N=$n"
  python3 scripts/launch.py "$CFG" \
    --label "$LABEL" --num-runs "$n" \
    --override slurm.gpus=1 \
    --override slurm.gres="$gres" \
    --override slurm.partition="$part" \
    --override slurm.time="$wall" \
    --override slurm.nodelist="$nodelist" \
    --override generator_params.bg_cli.num_designs="$nd" 2>&1 \
    | grep -iE "Submitted|jobid|error|uuid" | sed 's/^/    /'
  total=$((total + n))
done
echo "=== requested $total BoltzGen $LABEL gen-uuids ==="
