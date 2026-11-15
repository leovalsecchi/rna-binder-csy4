#!/bin/bash
# Union fold map: stages 2-5 over results_union/design, FOLD AXIS ONLY.
#
# Task 1 (scripts/stageU0_union.py) already assembled results_union/design from
# results_shared + results_v6 -- 99 shards, 533,596 hash rows over 526,467
# distinct coord_hash, with 23,940 rows routed to missing/. So there is no
# stage 1 (extract) and no stage 1b (split) here: this driver starts at stage 2,
# which dedups those 533,596 rows down to the 526,467 nodes that get searched.
#
# Parameters are frozen to match results_sw_k30_f50_mut exactly, so the union
# clustering is comparable to the per-map ones: --alignment-type 1 and
# `-c 0.8 --cov-mode 0 -e 10 --max-seqs 300 -a` (both in
# slurm/stage3_search.sbatch), TM = min(qTM,tTM), then --tm-floor 0.5
# --epi-floor 0.6 --k 30 --mutual.
#
# Fold axis only (decided 2026-09-28). Stage 4 still writes graph/epitope/ --
# it builds both axes in one pass and needs fp/ to do it -- but stage 5 is run
# once, on --axis fold. results_union/design/fp/ is retained so the epitope axis
# can be clustered later without redoing the 512-task search.
#
# Error handling: `set -e` does NOT propagate out of a $( ) capture unless
# inherit_errexit is set, and this repo has been bitten by that twice. Every
# sbatch below therefore carries an explicit `|| exit 1` as well. Unlike
# submit_maps.sh this driver submits inline rather than through a helper
# function, so no $( ) here wraps more than a single command -- but the guards
# stay, because the next edit might reintroduce one.
set -euo pipefail
shopt -s inherit_errexit

# FM is overridable only so the tests can point the preflight at a fixture tree;
# nothing in production sets it, and the sbatch scripts hardcode their own copy.
FM="${FM:-/cluster/work/wicky/lvalsecchi/rna_binder/analyses/foldseek_map}"
PY=/cluster/work/wicky/lvalsecchi/binder_map_venv/bin/python
LAYER=design
RESULTS=results_union
# 512 chunks at 128 concurrent. Chunking is a pure work partition and changes no
# result: q$I naming widens past 99 identically in stage3_split_queries.py and
# stage3_search.sbatch, and stage4_graph.py globs hits_*.tsv.
CHUNKS=512
THROTTLE=128

cd "$FM" || exit 1
mkdir -p "$FM/logs" || exit 1

# --- preflight: task 1's tree must be on disk --------------------------------
# Fails closed. Without this, stage 2 would run for an hour over an empty shard
# glob and hand stage 3 a query DB with nothing in it.
UNION="$FM/$RESULTS/$LAYER"
shopt -s nullglob
UNION_SHARDS=("$UNION"/shards/shard_*.tar.gz)
UNION_HASHES=("$UNION"/hashes/shard_*.tsv)
shopt -u nullglob
if [ ${#UNION_SHARDS[@]} -eq 0 ] || [ ${#UNION_HASHES[@]} -eq 0 ]; then
  echo "FATAL: $UNION holds no shard_*.tar.gz / shard_*.tsv pair." >&2
  echo "Run scripts/stageU0_union.py (task 1) before this driver." >&2
  exit 1
fi
echo "union tree: ${#UNION_SHARDS[@]} shards, ${#UNION_HASHES[@]} hash tables in $UNION"

# --- stage 2: global dedup, Foldseek DB, then split the queries --------------
J2=$(sbatch --parsable --job-name=union_s2 \
     --partition=normal.4h --time=04:00:00 --cpus-per-task=8 --mem-per-cpu=8G \
     --output="$FM/logs/union_s2_%j.log" \
     --wrap "bash -lc 'cd $FM && RESULTS=$RESULTS bash scripts/stage2_builddb.sh $LAYER && \
       PYTHONPATH=src $PY scripts/stage3_split_queries.py --layer $LAYER \
       --n-chunks $CHUNKS --results $FM/$RESULTS'") || exit 1
echo "union s2 dedup+createdb+split   $J2"

# --- stage 3: all-vs-all search, one array task per query chunk --------------
# 00:20:00 and 1G/cpu are MEASURED, from the one-chunk test on 2026-09-25:
# 56 s elapsed, 927 MB peak RSS. The earlier 04:00:00 / 4G-per-cpu request
# throttled placement badly for no benefit -- do not inflate them back.
# slurm/stage3_search.sbatch is passed as the script, NOT re-wrapped: it stages
# the target DB to node-local scratch, which is what stopped 265 concurrent
# readers from collapsing the Lustre client. Keep it that way.
J3=$(sbatch --parsable --job-name=union_s3 --dependency=afterok:"$J2" \
     --array=0-$((CHUNKS - 1))%${THROTTLE} \
     --time=00:20:00 --mem-per-cpu=1G \
     --export=ALL,ALN_TYPE=1,RESULTS=$RESULTS \
     slurm/stage3_search.sbatch $LAYER) || exit 1
echo "union s3 all-vs-all search      $J3"

# --- stage 4: kNN graphs -----------------------------------------------------
J4=$(sbatch --parsable --job-name=union_s4 --dependency=afterok:"$J3" \
     --partition=normal.4h --time=04:00:00 --cpus-per-task=16 --mem-per-cpu=8G \
     --output="$FM/logs/union_s4_%j.log" \
     --wrap "bash -lc 'cd $FM && PYTHONPATH=src $PY scripts/stage4_graph.py \
       --layer $LAYER --tm-floor 0.5 --epi-floor 0.6 --k 30 --mutual \
       --results $FM/$RESULTS'") || exit 1
echo "union s4 knn graph              $J4"

# --- stage 5: Leiden clusters, fold axis only --------------------------------
J5=$(sbatch --parsable --job-name=union_s5_fold --dependency=afterok:"$J4" \
     slurm/stage_single.sbatch scripts/stage5_cluster.py \
       --layer $LAYER --axis fold --results $FM/$RESULTS) || exit 1
echo "union s5 fold clusters          $J5"

echo "chain: s2=$J2 -> s3=$J3 -> s4=$J4 -> s5 fold=$J5"
