#!/bin/bash
# Chain one layer end to end.
#   bash slurm/submit_layer.sh <design|af3|esm> [n_shards] [n_chunks]
#
# Arrays are throttled at %16 on purpose: Euler evaluates only the top 250
# pending jobs, and flooding that window starves everything else in the
# account, including other people's work.
set -uo pipefail

LAYER="${1:?usage: submit_layer.sh <design|af3|esm> [n_shards] [n_chunks]}"
N_SHARDS="${2:-64}"
N_CHUNKS="${3:-32}"
TM_FLOOR="${TM_FLOOR:-0.5}"
EPI_FLOOR="${EPI_FLOOR:-0.5}"
ALN_TYPE="${ALN_TYPE:-2}"
# Overridable so the smoke test can run against build_smoke/results_smoke
# without colliding with a full run.
BUILD="${BUILD:-build}"
RESULTS="${RESULTS:-results}"

FM=/cluster/work/wicky/lvalsecchi/rna_binder/analyses/foldseek_map
cd "$FM" || exit 1
mkdir -p logs

echo "=== submitting layer=$LAYER shards=$N_SHARDS chunks=$N_CHUNKS"
echo "    TM_FLOOR=$TM_FLOOR EPI_FLOOR=$EPI_FLOOR ALN_TYPE=$ALN_TYPE"
echo "    BUILD=$BUILD RESULTS=$RESULTS"

# Stage 1 — extract + fingerprint
J1=$(sbatch --parsable --job-name="bm_s1_$LAYER" \
     --array="0-$((N_SHARDS - 1))%16" \
     --export=ALL,N_SHARDS="$N_SHARDS",BUILD="$BUILD",RESULTS="$RESULTS" \
     slurm/stage1_extract.sbatch "$LAYER") || exit 1
echo "  stage1 extract   $J1"

# Stage 2 — dedup + createdb + split queries
J2=$(sbatch --parsable --job-name="bm_s2_$LAYER" --dependency=afterok:"$J1" \
     --wrap "bash -lc 'cd $FM && RESULTS=$RESULTS bash scripts/stage2_builddb.sh $LAYER && \
       PYTHONPATH=src /cluster/work/wicky/lvalsecchi/binder_map_venv/bin/python \
       scripts/stage3_split_queries.py --layer $LAYER --n-chunks $N_CHUNKS \
       --results $FM/$RESULTS'" \
     --partition=normal.4h --time=04:00:00 --cpus-per-task=8 --mem-per-cpu=4G \
     --output="$FM/logs/bm_s2_${LAYER}_%j.log") || exit 1
echo "  stage2 builddb   $J2"

# Stage 3 — all-vs-all search
J3=$(sbatch --parsable --job-name="bm_s3_$LAYER" --dependency=afterok:"$J2" \
     --array="0-$((N_CHUNKS - 1))%16" \
     --export=ALL,ALN_TYPE="$ALN_TYPE",RESULTS="$RESULTS" \
     slurm/stage3_search.sbatch "$LAYER") || exit 1
echo "  stage3 search    $J3"

# Stage 4 — both graphs
J4=$(sbatch --parsable --job-name="bm_s4_$LAYER" --dependency=afterok:"$J3" \
     slurm/stage_single.sbatch scripts/stage4_graph.py \
       --layer "$LAYER" --tm-floor "$TM_FLOOR" --epi-floor "$EPI_FLOOR" \
       --results "$FM/$RESULTS") || exit 1
echo "  stage4 graph     $J4"

# Stages 5 and 6 — cluster and embed, per axis
PREV="$J4"
for AXIS in fold epitope; do
  J5=$(sbatch --parsable --job-name="bm_s5_${LAYER}_${AXIS}" \
       --dependency=afterok:"$PREV" \
       slurm/stage_single.sbatch scripts/stage5_cluster.py \
         --layer "$LAYER" --axis "$AXIS" --results "$FM/$RESULTS") || exit 1
  J6=$(sbatch --parsable --job-name="bm_s6_${LAYER}_${AXIS}" \
       --dependency=afterok:"$J5" \
       slurm/stage_single.sbatch scripts/stage6_embed.py \
         --layer "$LAYER" --axis "$AXIS" --results "$FM/$RESULTS") || exit 1
  echo "  stage5 $AXIS      $J5"
  echo "  stage6 $AXIS      $J6"
  PREV="$J6"
done

# Stage 7 — medoid traces
J7=$(sbatch --parsable --job-name="bm_s7_$LAYER" --dependency=afterok:"$PREV" \
     slurm/stage_single.sbatch scripts/stage7_medoids.py --layer "$LAYER" \
       --results "$FM/$RESULTS") || exit 1
echo "  stage7 medoids   $J7"

echo "$J7" > "$FM/logs/last_job_${LAYER}.txt"
echo "=== layer $LAYER submitted; final job $J7"
echo "Run stage 6b (comparisons) once ALL layers have finished:"
echo "  sbatch --dependency=afterok:<j7_design>:<j7_af3>:<j7_esm> \\"
echo "    slurm/stage_single.sbatch scripts/stage6b_compare.py"
