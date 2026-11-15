#!/bin/bash
# Three-map build: one shared 683 GB extraction, split per map, then stages 2-5
# per map, then the comparison. Map 2 (results_sw_k30_f50_mut) already exists
# and is never rebuilt.
# Parameters are frozen to match results_sw_k30_f50_mut exactly.
set -uo pipefail
FM=/cluster/work/wicky/lvalsecchi/rna_binder/analyses/foldseek_map
PY=/cluster/work/wicky/lvalsecchi/binder_map_venv/bin/python
cd "$FM" || exit 1
LAYER=design

# --- stage 1: one task per outer archive, over the UNION manifest -------------
J1=$(sbatch --parsable --job-name=cwfm_s1 --array=0-34%35 \
     --export=ALL,N_SHARDS=35,BUILD=build_shared,RESULTS=results_shared \
     slurm/stage1_extract_nested.sbatch) || exit 1
echo "stage1 extract (shared)  $J1"

# --- stage 1b: split into the two per-map trees ------------------------------
# lib12k is a census (9,759 nodes, not a sample of a larger population) --
# it MUST carry the :census tag or stage1b_split.py's shortfall check never runs.
J1B=$(sbatch --parsable --job-name=cwfm_s1b --dependency=afterok:$J1 \
      --partition=normal.4h --time=02:00:00 --cpus-per-task=4 --mem-per-cpu=4G \
      --output="$FM/logs/cwfm_s1b_%j.log" \
      --wrap "bash -lc 'cd $FM && PYTHONPATH=src $PY scripts/stage1b_split.py \
        --shared $FM/results_shared \
        --map full:build_full:results_full \
        --map lib12k:build_lib12k:results_lib12k:census'") || exit 1
echo "stage1b split            $J1B"

# --- stages 2-5, once per map ------------------------------------------------
submit_map () {
  local TAG=$1 RESULTS=$2 CHUNKS=$3 THROTTLE=$4 DEP=$5 GATE=${6:-0}
  local J2 J3 J4 J5 JGATE S3DEP
  # GATE is a switch, not a truthy value: `[ "$GATE" = "1" ]` alone would treat
  # a typo (`true`, `yes`, `01`) as "no gate" and launch the array unguarded,
  # silently. Anything but 0 or 1 is a bug at the call site, so refuse it.
  case "$GATE" in
    0|1) ;;
    *) echo "FATAL: submit_map $TAG: GATE must be 0 or 1, got '$GATE'" >&2; exit 1 ;;
  esac
  J2=$(sbatch --parsable --job-name=${TAG}_s2 --dependency=afterok:$DEP \
       --wrap "bash -lc 'cd $FM && RESULTS=$RESULTS bash scripts/stage2_builddb.sh $LAYER && \
         PYTHONPATH=src $PY scripts/stage3_split_queries.py --layer $LAYER \
         --n-chunks $CHUNKS --results $FM/$RESULTS'" \
       --partition=normal.4h --time=04:00:00 --cpus-per-task=8 --mem-per-cpu=8G \
       --output="$FM/logs/${TAG}_s2_%j.log") || exit 1
  # Optional hard gate between stage 2 and stage 3. Stage 3 for `full` is a
  # 512-task, ~2,100 CPU-h array, so its three correctness checks and the
  # scaling probe run as a Slurm job rather than by hand in the window between
  # the two -- afterok fails closed, so a failed gate leaves s3 in
  # DependencyNeverSatisfied instead of racing an operator with a scancel.
  # Only `full` is gated: lib12k is 4 chunks and cheap, and gating it would
  # make the gate depend on a tree it also validates.
  JGATE=""
  S3DEP=$J2
  if [ "$GATE" = "1" ]; then
    JGATE=$(sbatch --parsable --job-name=${TAG}_gate --dependency=afterok:$J2 \
            slurm/gate_full.sbatch) || exit 1
    S3DEP=$JGATE
  fi
  J3=$(sbatch --parsable --job-name=${TAG}_s3 --dependency=afterok:$S3DEP \
       --array=0-$((CHUNKS-1))%${THROTTLE} \
       --export=ALL,ALN_TYPE=1,RESULTS=$RESULTS \
       slurm/stage3_search.sbatch $LAYER) || exit 1
  J4=$(sbatch --parsable --job-name=${TAG}_s4 --dependency=afterok:$J3 \
       --partition=normal.4h --time=04:00:00 --cpus-per-task=16 --mem-per-cpu=8G \
       --output="$FM/logs/${TAG}_s4_%j.log" \
       --wrap "bash -lc 'cd $FM && PYTHONPATH=src $PY scripts/stage4_graph.py \
         --layer $LAYER --tm-floor 0.5 --epi-floor 0.6 --k 30 --mutual \
         --results $FM/$RESULTS'") || exit 1
  J5=$(sbatch --parsable --job-name=${TAG}_s5 --dependency=afterok:$J4 \
       slurm/stage_single.sbatch scripts/stage5_cluster.py \
         --layer $LAYER --axis fold --results $FM/$RESULTS) || exit 1
  echo "$TAG: s2=$J2${JGATE:+ gate=$JGATE} s3=$J3 s4=$J4 s5=$J5" >&2
  echo "$J5"
}

# The trailing 0/1 is GATE. The `|| exit 1` on each call site is load bearing:
# submit_map runs inside $( ), so its internal `|| exit 1` guards kill only the
# subshell -- without the outer guard a failed submission does not stop this
# script. Added 2026-09-24 after a sandbox repro. Do not remove or "simplify".
# `full` runs 512 chunks at 256 concurrent (was 128/64, widened 2026-09-24):
# chunking is a pure work partition and changes no result -- q$I naming widens
# past 99 identically in stage3_split_queries.py and stage3_search.sbatch, and
# stage4_graph.py globs hits_*.tsv.
J5_LIB=$(submit_map  lib12k results_lib12k 4   4   $J1B 0) || exit 1
J5_FULL=$(submit_map full   results_full   512 256 $J1B 1) || exit 1

# --- stage 9: compare all three ---------------------------------------------
J9=$(sbatch --parsable --job-name=cwfm_s9 \
     --dependency=afterok:$J5_FULL:$J5_LIB \
     --partition=normal.4h --time=01:00:00 --cpus-per-task=8 --mem-per-cpu=4G \
     --output="$FM/logs/cwfm_s9_%j.log" \
     --wrap "bash -lc 'cd $FM && PYTHONPATH=src $PY scripts/stage9_compare.py'") || exit 1
echo "stage9 compare           $J9"
echo "$J9" > $FM/logs/last_job_cwfm.txt
