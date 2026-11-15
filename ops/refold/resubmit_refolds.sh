#!/bin/bash
# resubmit_refolds.sh — re-fold incomplete refold leaves on a NON-PREEMPT slot.
#
# For each target leaf it (1) skips if the uuid has a live job, (2) clears the
# finalize marker + my packed outputs/inputs so the refold actually re-runs from
# the intact source design slice, (3) resubmits the leaf's own run.sbatch with a
# non-preempt partition/walltime (af3 -> gpuhe.24h/23:30, b2 -> gpuhe.4h/4h).
#
# DRY-RUN by default. Usage:
#   resubmit_refolds.sh <af3|b2|all> [maxN] [--go]
#   TARGETS=<tsv> resubmit_refolds.sh ...      # override target list
# Examples:
#   resubmit_refolds.sh b2 0           # dry-run all B2
#   resubmit_refolds.sh b2 0 --go      # submit all B2 (cheap wave first)
#   resubmit_refolds.sh af3 20 --go    # submit a 20-leaf AF3 pilot
set -uo pipefail
RB=/cluster/work/wicky/lvalsecchi/rna_binder
CAMP=$RB/campaigns/csy4_crrna_free
TARGETS=${TARGETS:-$RB/_refold_resubmit_targets_2026-06-14.tsv}
TYPE=${1:-all}; MAXN=${2:-0}; GO=0; [ "${3:-}" = "--go" ] && GO=1
AF3_PART=gpuhe.24h; AF3_TIME=23:30:00; AF3_MEM=   # unset until an AF3 run is measured
B2_PART=gpuhe.4h;   B2_TIME=04:00:00; B2_MEM=3G   # per-cpu (Euler requires --mem-per-cpu); x4 cpus=12G vs 32G default; measured peak ~7.5G -> eases QOS-MaxMem at scale
LOG=$RB/_refold_resubmit_2026-06-14.log
[ -f "$TARGETS" ] || { echo "no target list: $TARGETS" >&2; exit 1; }

# refold jobs are named just refold_af3/refold_b2 (no uuid), so map active jobs to
# their LEAF via StdOut path (%o = <leaf>/slurm.log). Catches RUNNING and PENDING.
ACTIVE=$(mktemp)
# %o (StdOut) is (null) for stdin-submitted jobs -> resolve the leaf via scontrol StdOut. Catches RUNNING+PENDING.
for j in $(squeue -u lvalsecchi -h -o '%i %j' 2>/dev/null | grep -E 'refold_(af3|b2)' | awk '{print $1}'); do
  scontrol show job "$j" 2>/dev/null | grep -oE 'StdOut=\S+' | sed 's#StdOut=##; s#/[^/]*$##'
done | sort -u > "$ACTIVE"
is_active_leaf(){ grep -qxF "$1" "$ACTIVE"; }

n=0; skipped=0
while IFS=$'\t' read -r stack label rtype base uuid ndone nsrc pct; do
  [ "$TYPE" != all ] && [ "$rtype" != "refold_$TYPE" ] && continue
  [ "$MAXN" -gt 0 ] && [ "$n" -ge "$MAXN" ] && break
  leaf="$CAMP/$stack/$label/$rtype/$base"; sb="$leaf/run.sbatch"
  [ -f "$sb" ] || { echo "SKIP(no run.sbatch) $stack/$label/$rtype/$base"; skipped=$((skipped+1)); continue; }
  [ -f "$leaf/designs.tar.zst" ] && { echo "SKIP(finalized) $base"; skipped=$((skipped+1)); continue; }
  if is_active_leaf "$leaf"; then echo "SKIP(active) $base"; skipped=$((skipped+1)); continue; fi
  if find "$leaf/outputs" "$leaf"/slurm*.log -mmin -15 2>/dev/null | grep -q .; then echo "SKIP(recent) $base"; skipped=$((skipped+1)); continue; fi
  if [ "$rtype" = refold_af3 ]; then PART=$AF3_PART; TM=$AF3_TIME; MEM=$AF3_MEM; else PART=$B2_PART; TM=$B2_TIME; MEM=$B2_MEM; fi
  MEMOPT=; [ -n "$MEM" ] && MEMOPT="--mem-per-cpu=$MEM"
  if [ "$GO" = 1 ]; then
    rm -f  "$leaf"/designs.tar.zst "$leaf"/designs.tar.zst.* "$leaf"/designs.lock
    rm -f  "$leaf"/outputs.tar.zst "$leaf"/outputs.tar.zst.* "$leaf"/outputs.lock
    rm -f  "$leaf"/inputs.tar.zst  "$leaf"/inputs.tar.zst.*  "$leaf"/inputs.lock
    rm -rf "$leaf"/designs "$leaf"/outputs "$leaf"/inputs
    # stale cache dirs from a prior run break boltz/lightning ("Missing folder: lightning_logs/version_0")
    rm -rf "$leaf"/.cache "$leaf"/.mpl_cache "$leaf"/.numba_cache "$leaf"/.triton_cache "$leaf"/.cueq_triton_cache
    # strip stale #SBATCH --dependency=afterok:<gen-job> (gen jobs long gone -> "dependency problem"); submit via stdin
    jid=$(grep -vE '^#SBATCH[[:space:]]+--dependency' "$sb" | sbatch --parsable --partition="$PART" --time="$TM" $MEMOPT 2>&1)
    echo "$(date '+%F %T') SUBMIT $rtype $base (${pct}% done) $PART -> $jid" | tee -a "$LOG"
  else
    echo "DRY $rtype $base (${pct}% done $ndone/$nsrc) -> clean+submit $PART/$TM"
  fi
  n=$((n+1))
done < "$TARGETS"
echo "=== $( [ "$GO" = 1 ] && echo SUBMITTED || echo DRY-RUN ): $n leaves (type=$TYPE maxN=$MAXN), $skipped skipped ==="
rm -f "$ACTIVE"
