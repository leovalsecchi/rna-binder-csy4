#!/bin/bash
# Resume all incomplete midrun_v2 refold leaves (af3 msaoff+msaon, b2) on bulk.
# No pre-clean: patched refold scripts self-unpack partial + clear caches at job start
# (queued leaves stay compact -> no inode spike; only running leaves go loose).
# bulk / 23:30 / --requeue(in sbatch) / mem-per-cpu 4G(af3) 3G(b2). Skips finalized+active.
# DRY by default; --go to submit.
set -uo pipefail
RB=/cluster/work/wicky/lvalsecchi/rna_binder
CAMP=$RB/campaigns/csy4_crrna_free
GO=0; [ "${1:-}" = "--go" ] && GO=1
LOG=$RB/_refold_resume_2026-06-16.log
STACKS="rfd3__ligandmpnn boltzdesign__ligandmpnn boltzgen__ligandmpnn boltzgen__boltzif"
ACT=$(mktemp)
# active = my running/pending refold leaves. ONE scontrol call (comma/space lists are "Invalid"),
# filter my refold jobs, map to leaf via StdOut. Reliable + scale-safe.
scontrol show job -o 2>/dev/null | grep 'UserId=lvalsecchi' | grep -E 'JobName=refold_(af3|b2|esm)' \
  | grep -oE 'StdOut=\S+' | sed 's#StdOut=##; s#/[^/]*$##' | sort -u > "$ACT"
echo "active leaves to skip: $(wc -l < "$ACT")"
submit() {  # $1=run.sbatch $2=mempercpu -> echoes jobid or ERR
  local out t=0
  while [ $t -lt 4 ]; do
    out=$(grep -vE '^#SBATCH[[:space:]]+--dependency' "$1" | sbatch --parsable --partition=bulk --time=23:30:00 --mem-per-cpu="$2" 2>&1 | tail -1)
    case "$out" in [0-9]*) echo "$out"; return 0;; *[Ss]ocket*|*timed*out*) t=$((t+1)); sleep 3;; *) echo "ERR:$out"; return 1;; esac
  done; echo ERR:timeoutx4; return 1
}
tot=0
for S in $STACKS; do
  for spec in "refold_af3 4G" "refold_b2 3G" "refold_esm 8G"; do
    set -- $spec; rt=$1; mem=$2; n=0
    for d in "$CAMP/$S/midrun_v2/$rt"/*/; do
      d=${d%/}; [ -d "$d" ] || continue
      [ -f "$d/designs.tar.zst" ] && continue
      grep -qxF "$d" "$ACT" 2>/dev/null && continue
      if [ "$GO" = 1 ]; then jid=$(submit "$d/run.sbatch" "$mem"); echo "$rt $(basename $d) -> $jid" >> "$LOG"; sleep 0.4; fi
      n=$((n+1))
    done
    printf "  %-26s %-12s %s %d\n" "$S" "$rt" "$([ "$GO" = 1 ] && echo submitted || echo would-submit)" "$n"
    tot=$((tot+n))
  done
done
echo "=== TOTAL $([ "$GO" = 1 ] && echo SUBMITTED || echo would-submit): $tot  (errors: $([ "$GO" = 1 ] && grep -c ERR: "$LOG" 2>/dev/null || echo 0)) ==="
rm -f "$ACT"
