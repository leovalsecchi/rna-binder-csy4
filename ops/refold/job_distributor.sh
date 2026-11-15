#!/usr/bin/env bash
# job_distributor.sh — rebalance PENDING jobs across gpuhe queues by emptiness.
#
# Rules:
#   1. Generators NEVER on bulk.
#   2. Generators prioritised in 4h/24h/120h (placed in the emptiest tier they fit).
#   3. While ANY generator is pending -> refolds go to bulk only.
#      Once NO generator is pending -> refolds spread across ALL queues.
#   4. Walltime set per queue: 4h=04:00, 24h=24:00, 120h/bulk=120:00.
#
# Mechanism: post-submit `scontrol update partition=/TimeLimit=` (submit-time
# gpuhe.bulk gets remapped; scontrol moves STICK). TimeLimit can only DECREASE,
# so this assumes jobs are submitted with a generous (>=120h) ceiling.
#
# Emptiness (MVP): fractional free per-tier QOS GPU quota = (cap-running)/cap,
#   caps (es_wicky): 4h=96 24h=50 120h=29 bulk=240. Tie -> prefer faster tier.
#
# Job runtime needs (h, => shortest fitting tier): esm/b2=1, af3=15, gen=8.
#
#   job_distributor.sh            # dry-run: print the plan + would-be moves
#   job_distributor.sh --go       # actually scontrol-move
#   job_distributor.sh --user X   # default lvalsecchi
set -u
U=lvalsecchi; GO=0; MAXMOVES=300
while [ $# -gt 0 ]; do case "$1" in
  --go) GO=1;; --user) U=$2; shift;; --max) MAXMOVES=$2; shift;;
  *) echo "unknown arg: $1" >&2; exit 2;; esac; shift; done

SNAP="$(scontrol show job -o 2>/dev/null | grep "UserId=${U}")"
[ -z "$SNAP" ] && { echo "no jobs for $U"; exit 0; }

# awk emits: "PLAN|<text>" rows for the report, and "MOVE|jid|from|totier|settime"
PLAN="$(printf '%s\n' "$SNAP" | awk '
function tier_of(p){ if(p~/bulk/)return "bulk"; if(p~/120h/)return "120h"; if(p~/24h/)return "24h"; if(p~/4h/)return "4h"; return "other" }
BEGIN{
  split("4h 24h 120h bulk", ORD, " ");
  cap["4h"]=96; cap["24h"]=50; cap["120h"]=29; cap["bulk"]=240;
  sth["4h"]=4; sth["24h"]=24; sth["120h"]=120; sth["bulk"]=120;
  st["4h"]="04:00:00"; st["24h"]="24:00:00"; st["120h"]="120:00:00"; st["bulk"]="120:00:00";
  need["esm"]=1; need["b2"]=1; need["af3"]=15; need["gen"]=8;
}
{
  jid=jn=js=part=rs="";
  for(i=1;i<=NF;i++){ t=$i;
    if(t~/^JobId=/)jid=substr(t,7);
    else if(t~/^JobName=/)jn=substr(t,9);
    else if(t~/^JobState=/)js=substr(t,10);
    else if(t~/^Partition=/)part=substr(t,11);
    else if(t~/^Reason=/)rs=substr(t,8);
  }
  if(jn~/refold_af3/){kind="refold"; type="af3"}
  else if(jn~/refold_b2/){kind="refold"; type="b2"}
  else if(jn~/refold_esm/){kind="refold"; type="esm"}
  else {kind="gen"; type="gen"}
  tr=tier_of(part);
  if(js=="RUNNING"){ run[tr]++; next }
  if(js!="PENDING" || rs=="JobHeldUser") next;        # leave held/non-pending alone
  pn++; PJ[pn]=jid; PK[pn]=kind; PT[pn]=type; PP[pn]=part; PTR[pn]=tr;
  if(kind=="gen") gp++;
}
END{
  gens_pending = (gp>0);
  # process generators first, then refolds (stable: just two passes)
  for(pass=0; pass<2; pass++){
    for(k=1;k<=pn;k++){
      isgen=(PK[k]=="gen");
      if(pass==0 && !isgen) continue;
      if(pass==1 && isgen) continue;
      type=PT[k];
      # candidate tiers per rules
      if(isgen) nc=split("4h 24h 120h", CAND, " ");         # rule 1: never bulk
      else if(gens_pending) nc=split("bulk", CAND, " ");     # rule 3a
      else nc=split("4h 24h 120h bulk", CAND, " ");          # rule 3b
      # pick emptiest FITTING candidate (fractional free), tie -> faster (earlier) tier
      best=""; bestscore=-1e9;
      for(ci=1; ci<=nc; ci++){ T=CAND[ci];
        if(need[type] > sth[T]) continue;                   # walltime does not fit
        free=(cap[T]-run[T])/cap[T];
        if(free>bestscore+0.0001){ bestscore=free; best=T }
      }
      if(best==""){ printf "PLAN|  %-7s job %-9s  NO FITTING TIER (needs %sh walltime)\n", PT[k], PJ[k], need[type]; continue }
      run[best]++;                                          # reserve (so we spread)
      from=PTR[k];
      if(from==best){ stay[best,type]++; continue }
      mv[best,type]++; movecount++;
      # only the 4h tier needs a walltime change (decrease, to be accepted); other
      # tiers keep the job\47s submitted walltime (no scontrol increase needed).
      settime = (best=="4h") ? "04:00:00" : "-";
      printf "MOVE|%s|%s|%s|%s\n", PJ[k], from, best, settime;
    }
  }
  # summary
  printf "PLAN|gens_pending=%s  (refolds -> %s)\n", (gens_pending?"YES":"no"), (gens_pending?"bulk only":"all tiers");
  printf "PLAN|running by tier: 4h=%d 24h=%d 120h=%d bulk=%d\n", run["4h"]+0,run["24h"]+0,run["120h"]+0,run["bulk"]+0;
  printf "PLAN|moves planned: %d\n", movecount+0;
}
')"

echo "JOB DISTRIBUTOR — $U — $(date '+%Y-%m-%d %H:%M')   $([ $GO = 1 ] && echo '[--GO]' || echo '[dry-run]')"
printf '%s\n' "$PLAN" | grep '^PLAN|' | sed 's/^PLAN|//'
echo
echo "Planned moves (jobid: from -> tier @ walltime), by target+type:"
printf '%s\n' "$PLAN" | grep '^MOVE|' | awk -F'|' '{print $4"  "$3}' | sort | uniq -c | sort -rn | sed 's/^/  /'

if [ "$GO" = 1 ]; then
  echo; echo "Applying (cap=$MAXMOVES)..."
  n=0
  while IFS='|' read -r tag jid from to settime; do
    [ "$tag" = MOVE ] || continue
    n=$((n+1)); [ "$n" -gt "$MAXMOVES" ] && { echo "  hit move cap $MAXMOVES"; break; }
    # 4h placement needs the walltime lowered first so the move is accepted;
    # other tiers keep the job's existing (fitting) walltime.
    [ "$settime" != "-" ] && scontrol update jobid="$jid" TimeLimit="$settime" 2>/dev/null
    scontrol update jobid="$jid" partition="gpuhe.$to" 2>/dev/null
  done < <(printf '%s\n' "$PLAN" | grep '^MOVE|')
  echo "  applied $((n>MAXMOVES?MAXMOVES:n)) moves"
fi
