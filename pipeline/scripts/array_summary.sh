#!/usr/bin/env bash
# array_summary.sh — per-array task tally for my Slurm job arrays.
#   One row per array (ARRAYJOB / NAME / PARTITION) with RUN / PEND / TOTAL = array TASKS
#   (squeue -r expands array elements; a pending range [0-209] counts as 210, not 1).
#   Companion to job_summary.sh (which now has an ARRAYS column + task-expanded counts).
# usage: bash array_summary.sh [user]      (default user: lvalsecchi)
set -u
U="${1:-lvalsecchi}"
echo "ARRAY SUMMARY — ${U} — $(date '+%Y-%m-%d %H:%M')"
echo
ARR="$(squeue -h -u "$U" -r -o '%F|%K|%T|%j|%P' 2>/dev/null | awk -F'|' '
  $2!="N/A" && $2!="" {
    key=$1"|"$4"|"$5; seen[key]=1; tot[key]++;
    if($3=="RUNNING")run[key]++; else if($3=="PENDING")pend[key]++;
    arr[$1]=1; TR+=($3=="RUNNING"); TP+=($3=="PENDING"); TT++;
  }
  END{
    for(k in seen){split(k,c,"|"); printf "%s|%s|%s|%d|%d|%d\n",c[1],c[2],c[3],run[k]+0,pend[k]+0,tot[k]+0}
    n=0;for(a in arr)n++;
    if(TT>0) printf "ZZTOTAL|%d array(s)||%d|%d|%d\n",n,TR+0,TP+0,TT+0
  }')"
if [ -n "$ARR" ]; then
  { echo "ARRAYJOB|NAME|PARTITION|RUN|PEND|TOTAL";
    printf '%s\n' "$ARR" | grep -v '^ZZTOTAL' | sort -t'|' -k1,1n;
    printf '%s\n' "$ARR" | grep '^ZZTOTAL' | sed 's/^ZZTOTAL/TOTAL/'; } | column -t -s'|'
else
  echo "  (no array jobs)"
fi
