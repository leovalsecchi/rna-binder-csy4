#!/bin/bash
# gpubench_v1 PHASE 2 -- fires automatically when every gen job has ended.
#
# Chained with --dependency=afterany on the gen jobs rather than watched from a
# login shell: long-running local watchers get killed in this environment, so
# detached finalisation belongs in a SLURM dependency.
#
# afterany, not afterok, on purpose: if a gen job fails we want this to RUN and
# report the shortfall, not hang forever. bench_refold.py refuses to submit
# unequal batch sizes, so a partial gen stops the refold phase loudly.
#
# The BG+BIF true-backbone rewrite is no longer a separate step against Lustre:
# a finished gen run has no designs/ at all (the finalize tail packs it into
# designs.tar.zst and deletes the tree), so bench_refold.py now applies the
# rewrite to the STAGED copy while it builds each pool -- production gen output
# is never mutated. What stays here is the self-check, which still GATES the
# pool build: if the transplant cannot reproduce BG+LiMPNN designs from their
# specs, no refold is submitted.
#SBATCH --job-name=gpubench_phase2
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=8G
#SBATCH --time=01:00:00
#SBATCH --tmp=20G
#SBATCH --output=/cluster/work/wicky/lvalsecchi/rna_binder/bench_gpueff/phase2_%j.log
#SBATCH --error=/cluster/work/wicky/lvalsecchi/rna_binder/bench_gpueff/phase2_%j.log
set -uo pipefail

BENCH=/cluster/work/wicky/lvalsecchi/rna_binder/bench_gpueff
cd "$BENCH" || exit 1
LABEL=gpubench_v1

echo "=== phase2 start $(date -Is) on $(hostname) ==="

echo
echo "### gen outcome (sacct -D, requeued attempts visible) ###"
JIDS=$(awk -F'\t' 'NR>1 && $1=="gen" {printf "%s%s", sep, $6; sep=","}' gen_manifest.tsv)
[ -n "$JIDS" ] && sacct -D -X -j "$JIDS" \
  --format=JobID,JobName%22,State,ElapsedRaw,NodeList%14 -P | column -s'|' -t

echo
echo "### step 1: BG+BIF true-backbone SELF-CHECK (gates the pool build) ###"
# Rebuilds BG+LiMPNN designs from their specs and requires exact agreement, so a
# broken transplant can never reach a pool. Reads the packed tarballs.
python3 bif_true_designs.py --label "$LABEL" --validate
rc1=$?
echo "bif_true_designs.py rc=$rc1"
if [ "$rc1" -ne 0 ]; then
  echo "!! true-backbone self-check failed -- NOT submitting refolds."
  echo "   BG+BIF would be scored against the re-prediction again, which is the"
  echo "   whole bug this bench was built to avoid. Fix, then run:"
  echo "     python3 bench_refold.py --label $LABEL --submit"
  exit 1
fi

echo
echo "### step 2: equal-N pools (extract -> BIF rewrite -> repack) + submit ###"
python3 bench_refold.py --label "$LABEL" --submit
echo "bench_refold.py rc=$?"

echo
echo "### current GPU-h table (gen only; refolds still to run) ###"
python3 collect_gpuh.py --label "$LABEL" 2>&1 | tail -40

# Chain the collector too, so the final table lands with no further action.
# Same reasoning as phase2 itself: a SLURM dependency, not a watcher.
echo
echo "### chaining phase3 (collector) on the refold jobs ###"
RJ=$(awk -F'\t' 'NR>1 && $1=="refold" {printf "%s%s", sep, $6; sep=":"}' \
     refold_manifest.tsv 2>/dev/null)
if [ -n "$RJ" ]; then
  cat > "$BENCH/_phase3.sbatch" <<EOF
#!/bin/bash
#SBATCH --job-name=gpubench_phase3
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=4G
#SBATCH --time=00:30:00
#SBATCH --output=$BENCH/phase3_%j.log
#SBATCH --error=$BENCH/phase3_%j.log
cd $BENCH
echo "=== FINAL gpubench_v1 RESULT \$(date -Is) ==="
python3 collect_gpuh.py --label $LABEL
EOF
  p3=$(sbatch --parsable --dependency=afterany:"$RJ" "$BENCH/_phase3.sbatch")
  echo "phase3 jobid=$p3 (afterany:$RJ)"
else
  echo "no refold jobs in refold_manifest.tsv -- run collect_gpuh.py by hand later"
fi

echo
echo "=== phase2 done $(date -Is) ==="
echo "Final table will be in $BENCH/phase3_<jobid>.log and gpuh_per_design.tsv"
