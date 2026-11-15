#!/usr/bin/env bash
# Euler status digest: inode/disk quota + SLURM jobs + per-stack gen/refold funnel.
# Sends via Euler's /usr/sbin/sendmail. Must run ON Euler (needs squeue/sacct/lfs/sendmail).
# Intended trigger: a bs-gobi cron, e.g.
#   0 8,12,16,20 * * * ssh euler 'bash /cluster/work/wicky/lvalsecchi/rna_binder/code/scripts/status_digest.sh' >/dev/null 2>&1
# Flags:
#   --dry-run         print the message instead of emailing
#   --refresh-stacks  force-regenerate the (slow, ~6 min) stack funnel this run
# Env:
#   DIGEST_TO                recipient (default leo@javal.ch)
#   DIGEST_FUNNEL_MAX_AGE_MIN  reuse cached funnel if younger than this (default 300)
set +e

TO="${DIGEST_TO:-leo@javal.ch}"
FROM="lvalsecchi@euler.ethz.ch"
ACCT="lvalsecchi"
PROJ=510551
WORK=/cluster/work/wicky
RBROOT=/cluster/work/wicky/lvalsecchi/rna_binder
CACHE_DIR="$RBROOT/analyses/_digest"
CACHE="$CACHE_DIR/stacks.txt"
MAXAGE_MIN="${DIGEST_FUNNEL_MAX_AGE_MIN:-300}"

DRYRUN=0; FORCE=0
for a in "$@"; do
  [ "$a" = "--dry-run" ] && DRYRUN=1
  [ "$a" = "--refresh-stacks" ] && FORCE=1
done

HOST=$(hostname)
NOW=$(date '+%Y-%m-%d %H:%M %Z')

# ---- quota: inodes + disk (project 510551) ----
Q=$(lfs quota -p "$PROJ" "$WORK" 2>/dev/null | awk '/work\/wicky/{print $2,$4,$6,$7,$8; exit}')
KB=$(echo "$Q"|awk '{print $1}'); BLIMIT=$(echo "$Q"|awk '{print $2}')
FILES=$(echo "$Q"|awk '{print $3}'); ISOFT=$(echo "$Q"|awk '{print $4}'); IHARD=$(echo "$Q"|awk '{print $5}')
QUOTA=$(awk -v f="$FILES" -v s="$ISOFT" -v h="$IHARD" -v k="$KB" -v bl="$BLIMIT" 'BEGIN{
  printf "inodes: %d / soft %d / hard %d  (%.1f%% of hard, %.1f%% of soft; %d to hard)\n", f,s,h,100*f/h,100*f/s,h-f;
  printf "disk:   %.2f TB / %.2f TB  (%.1f%%)", k/1073741824, bl/1073741824, 100*k/bl;
  if (f>s) printf "\n*** OVER inode SOFT quota -- packing/cleanup needed ***";
}')

# ---- jobs ----
RUN=$(squeue -u "$ACCT" -h -t RUNNING 2>/dev/null | wc -l)
PEND=$(squeue -u "$ACCT" -h -t PENDING 2>/dev/null | wc -l)
TOT=$(squeue -u "$ACCT" -h 2>/dev/null | wc -l)
BYPART=$(squeue -u "$ACCT" -h -o '%P' 2>/dev/null | sort | uniq -c | sort -rn | awk '{printf "    %7s  %s\n",$1,$2}')
[ -z "$BYPART" ] && BYPART="    (no queued jobs)"
SINCE=$(date -d '24 hours ago' '+%Y-%m-%dT%H:%M' 2>/dev/null)
DONE=$(sacct -u "$ACCT" -X -n -S "$SINCE" -o State 2>/dev/null | awk '{print $1}' \
  | grep -vE '^(PENDING|RUNNING|REQUEUED|RESIZING|SUSPENDED)$' | sort | uniq -c | sort -rn \
  | awk '{printf "    %7s  %s\n",$1,$2}')
[ -z "$DONE" ] && DONE="    (none)"

# ---- stacks: per-stack gen/refold funnel (cached; ~6 min to regenerate) ----
mkdir -p "$CACHE_DIR" 2>/dev/null
NEED=1
if [ "$FORCE" = 0 ] && [ -f "$CACHE" ]; then
  AGE=$(( ( $(date +%s) - $(stat -c %Y "$CACHE") ) / 60 ))
  [ "$AGE" -lt "$MAXAGE_MIN" ] && NEED=0
fi
if [ "$NEED" = 1 ]; then
  module purge 2>/dev/null || true; module load stack/2024-06 python/3.12.8 2>/dev/null || true
  OUT=$(cd "$RBROOT" 2>/dev/null && timeout 600 python3 code/scripts/campaign_status.py --no-color 2>&1); RC=$?
  if [ "$RC" -eq 0 ] && [ -n "$OUT" ]; then printf '%s\n' "$OUT" > "$CACHE.tmp" && mv "$CACHE.tmp" "$CACHE"; fi
fi
if [ -f "$CACHE" ]; then
  CTS=$(date -d "@$(stat -c %Y "$CACHE")" '+%Y-%m-%d %H:%M')
  AGE=$(( ( $(date +%s) - $(stat -c %Y "$CACHE") ) / 60 ))
  STACKS="(funnel as of $CTS -- ${AGE} min ago)
$(cat "$CACHE")"
else
  STACKS="(stack funnel unavailable -- campaign_status.py failed and no cache yet)"
fi

# ---- compose + send ----
BODY=$(printf 'Euler status digest -- %s (%s)\n\n=== QUOTA (project %s, %s) ===\n%s\n\n=== JOBS (%s) ===\nrunning %s | pending %s | total in queue %s\nby partition:\n%s\nfinished last 24h:\n%s\n\n=== STACKS -- generation / refolding ===\n%s\n\n--\nstatus_digest.sh on %s\n' \
  "$NOW" "$HOST" "$PROJ" "$WORK" "$QUOTA" "$ACCT" "$RUN" "$PEND" "$TOT" "$BYPART" "$DONE" "$STACKS" "$HOST")
SUBJECT="Euler $(date '+%m-%d %H:%M'): inodes ${FILES}/${IHARD} | ${RUN} run/${PEND} pend"

if [ "$DRYRUN" = 1 ]; then
  echo "SUBJECT: $SUBJECT"; echo; echo "$BODY"; exit 0
fi
printf 'To: %s\nFrom: %s\nReply-To: %s\nSubject: %s\nContent-Type: text/plain; charset=UTF-8\n\n%s\n' \
  "$TO" "$FROM" "$TO" "$SUBJECT" "$BODY" | /usr/sbin/sendmail -t
echo "sent rc=$?"
