#!/usr/bin/env bash
# recent_jobs.sh — the last N *finished* SLURM jobs (any terminal state), newest first.
# Companion to tally.sh / job_summary.sh.  Alias: recent_jobs
#
#   recent_jobs [N] [-d DAYS] [-u USER]
#     N        how many jobs to show        (default 30)
#     -d DAYS  how far back to scan sacct   (default 7)
#     -u USER  which user                   (default $USER)
#
# Shows COMPLETED / CANCELLED / FAILED / TIMEOUT / REQUEUED / NODE_FAIL / PREEMPTED …
# i.e. everything that has an End time and is no longer RUNNING or PENDING.
set -u
N=30 ; DAYS=7 ; USER_Q="${USER:-$(whoami)}"
while [[ $# -gt 0 ]]; do
  case "$1" in
    -d) DAYS="${2:?-d needs a number}"; shift 2 ;;
    -u) USER_Q="${2:?-u needs a user}"; shift 2 ;;
    -h|--help) sed -n '2,9p' "$0" | sed 's/^# \{0,1\}//'; exit 0 ;;
    ''|*[!0-9]*) echo "recent_jobs: bad arg '$1' (want a number, or -d DAYS / -u USER)" >&2; exit 1 ;;
    *) N="$1"; shift ;;
  esac
done

USE_COLOR=0 ; [[ -t 1 ]] && USE_COLOR=1

# -X = one row per job (array tasks show as <id>_<idx>); -P = pipe-delimited.
# Keep only rows with a real End time that aren't still RUNNING/PENDING, sort by
# End descending (ISO timestamps sort chronologically), take the newest N.
sacct -X -u "$USER_Q" -S "now-${DAYS}days" -E now \
      -o JobID,State,Elapsed,End,ExitCode,Partition,NodeList,JobName%60 -P --noheader 2>/dev/null \
| awk -F'|' '$4!="Unknown" && $2!~/^(RUNNING|PENDING)/' \
| sort -t'|' -k4,4r \
| head -n "$N" \
| awk -F'|' -v color="$USE_COLOR" -v days="$DAYS" '
BEGIN{
  printf "%-12s  %-19s  %-20s  %-9s  %-6s  %-11s  %-14s  %s\n",
         "END","JOBID","STATE","ELAPSED","EXIT","PART","NODE","NAME";
}
{
  jobid=$1; state=$2; el=$3; end=$4; exitc=$5; part=$6; node=$7; name=$8;
  split(end,e,"T"); endf=substr(e[1],6) " " substr(e[2],1,5);   # 07-13 12:10
  split(state,sw," "); key=sw[1]; cnt[key]++;
  if(length(state)>20) state=substr(state,1,19) ".";
  if(length(name)>40)  name=substr(name,1,39) ".";
  if(length(node)>14)  node=substr(node,1,13) ".";
  c=""; r="";
  if(color=="1"){
    r="\033[0m";
    if      (key=="COMPLETED")                                             c="\033[32m";
    else if (key=="FAILED"||key=="NODE_FAIL"||key=="BOOT_FAIL"||key=="OUT_OF_MEMORY") c="\033[31m";
    else if (key=="TIMEOUT")                                               c="\033[35m";
    else if (key=="CANCELLED")                                             c="\033[33m";
    else if (key=="REQUEUED"||key=="PREEMPTED"||key=="SUSPENDED")          c="\033[36m";
  }
  stpad=sprintf("%-20s", state);
  printf "%-12s  %-19s  %s%s%s  %-9s  %-6s  %-11s  %-14s  %s\n",
         endf, jobid, c, stpad, r, el, exitc, part, node, name;
}
END{
  if(NR==0){ print "  (no finished jobs in the last " days " days — try a wider window: recent_jobs -d 30)"; exit }
  s="";
  for(k in cnt) s=s sprintf("%s %d \xc2\xb7 ", k, cnt[k]);
  sub(/ \xc2\xb7 $/,"",s);
  printf "\n  %d jobs shown \xc2\xb7 %s\n", NR, s;
}'
