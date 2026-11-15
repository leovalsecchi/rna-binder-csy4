#!/usr/bin/env bash
# Move runnable B2 backlog (pending, source designs ready, not afterok-blocked)
# from gpuhe onto the public gpu.* queues (2080 Ti etc.). Round-robin 4h/24h/120h,
# gpumem:8g floor (11G 2080 Ti qualifies), walltime per queue. Cancels the gpuhe
# copy then resubmits (scontrol can't retype gres). DRY by default; --go to act.
set -u
GO=0; LIMIT=0
while [ $# -gt 0 ]; do case "$1" in --go) GO=1;; --limit) LIMIT=$2; shift;; esac; shift; done
TIERS=(gpu.4h gpu.24h gpu.120h); TW=(04:00:00 24:00:00 120:00:00)
# runnable pending b2 (NOT afterok Dependency) -> jobid + leaf
scontrol show job -o 2>/dev/null | grep "UserId=lvalsecchi" | grep "JobName=refold_b2" \
  | grep "JobState=PENDING" | grep -v "Reason=Dependency" | grep -E "Partition=[^[:space:]]*gpuhe" > /tmp/_pb2.txt
total=$(wc -l < /tmp/_pb2.txt); n=0; sub=0
while read -r line; do
  jid=$(echo "$line" | grep -oE 'JobId=[0-9]+' | cut -d= -f2)
  leaf=$(echo "$line" | grep -oE 'StdOut=\S+' | sed 's#StdOut=##; s#/[^/]*$##')
  [ -z "$jid" ] || [ -z "$leaf" ] && continue
  [ -f "$leaf/designs.tar.zst" ] && continue
  [ -f "$leaf/run.sbatch" ] || continue
  ti=$((n % 3)); part=${TIERS[$ti]}; wt=${TW[$ti]}
  n=$((n+1))
  if [ "$GO" = 1 ]; then
    scancel "$jid" 2>/dev/null
    out=$(grep -vE '^#SBATCH[[:space:]]+--dependency' "$leaf/run.sbatch" \
      | sbatch --parsable --partition="$part" --gres=gpumem:8g --time="$wt" 2>&1 | tail -1)
    case "$out" in [0-9]*) sub=$((sub+1));; *) echo "  SUBMIT FAIL ($part): $out"; esac
  fi
  [ "$LIMIT" -gt 0 ] && [ "$n" -ge "$LIMIT" ] && break
done < /tmp/_pb2.txt
echo "runnable b2 backlog: $total | $([ $GO = 1 ] && echo "moved=$sub to gpu.4h/24h/120h" || echo "would-move=$n")"
