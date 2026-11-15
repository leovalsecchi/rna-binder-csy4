#!/usr/bin/env bash
# Usage: analyze_campaign.sh <campaign_label> [--structural] [--rerun-label LABEL]
set -euo pipefail

if [[ $# -lt 1 ]]; then
    echo "Usage: $0 <campaign_label> [--structural] [--rerun-label LABEL]" >&2
    exit 1
fi
CAMPAIGN="$1"; shift

STRUCTURAL=""
EXTRA_LABELS=()
while [[ $# -gt 0 ]]; do
    case "$1" in
        --structural) STRUCTURAL="--structural"; shift ;;
        --rerun-label) EXTRA_LABELS+=("--label" "$2"); shift 2 ;;
        *) echo "Unknown arg: $1" >&2; exit 1 ;;
    esac
done

BASE=/cluster/work/wicky/lvalsecchi/rna_binder/campaigns/csy4_crrna_free
OUT=/cluster/work/wicky/lvalsecchi/rna_binder/analyses/$CAMPAIGN
mkdir -p "$OUT"

module purge 2>/dev/null || true
module load stack/2024-06 python/3.12.8 2>/dev/null || true

echo "=== Discovering run dirs for '$CAMPAIGN' ==="
RUN_DIRS=()
for stack_dir in $BASE/*/; do
    for run_dir in "$stack_dir"${CAMPAIGN}_*/; do
        [[ -d "$run_dir" ]] && RUN_DIRS+=("$run_dir")
    done
done
if [[ ${#RUN_DIRS[@]} -eq 0 ]]; then
    echo "ERROR: No run dirs found for '$CAMPAIGN'" >&2; exit 1
fi
echo "Found ${#RUN_DIRS[@]}:"
for d in "${RUN_DIRS[@]}"; do echo "  $d"; done
echo

CSV_NAME="${CAMPAIGN}_all_stacks${STRUCTURAL:+_structural}.csv"
CSV_PATH="$OUT/$CSV_NAME"

echo "=== Aggregating ${STRUCTURAL:+(structural)} ==="
python3 /cluster/work/wicky/lvalsecchi/rna_binder/code/scripts/aggregate_metrics.py \
    "${RUN_DIRS[@]}" \
    --label "$CAMPAIGN" \
    "${EXTRA_LABELS[@]}" \
    $STRUCTURAL \
    --output "$CSV_PATH" \
    --progress \
    2>&1 | tee "$OUT/aggregate.log"

if [[ -n "$STRUCTURAL" ]]; then
    echo
    echo "=== Running full analysis ==="
    CSV_PATH="$CSV_PATH" OUT_DIR="$OUT" \
        python3 /cluster/work/wicky/lvalsecchi/rna_binder/analyses/Large_Gen_App_v1/full_analysis.py \
        2>&1 | tee "$OUT/full_analysis_output.txt"
fi

echo
echo "=== Done. Outputs in $OUT ==="
ls -lh "$OUT" | grep -v "^total"
