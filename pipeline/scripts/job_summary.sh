#!/usr/bin/env bash
# job_summary.sh — summary of my current Slurm jobs by stack + role.
#   role = generation (rfd3/bd1/bg/hd) OR refolder (af3 / esm / b2)
#   (af3 msaon deactivated 2026-06-18; ESM is the 3rd refolder)
#   columns: ARRAYS (# distinct job arrays in the row) / RUN / PEND / HELD / OTH
#   counts are by TASK — array elements are expanded (a pending [0-209] = 210 PEND, not 1).
#   detailed per-array breakdown: see the `array_summary` command.
# usage: bash job_summary.sh [user]      (default user: lvalsecchi)
set -u
U="${1:-lvalsecchi}"
JOBS="$(scontrol show job -o 2>/dev/null | grep "UserId=${U}")"
if [ -z "$JOBS" ]; then echo "No jobs for ${U}."; exit 0; fi

echo "JOB SUMMARY — ${U} — $(date '+%Y-%m-%d %H:%M')"
echo

ROWS="$(printf '%s\n' "$JOBS" | awk '
function ss(s){
  if(s=="rfd3__ligandmpnn")return "rfd3";
  if(s=="boltzdesign__ligandmpnn")return "boltzdesign";
  if(s=="boltzgen__boltzif")return "bg_boltzif";
  if(s=="boltzgen__ligandmpnn")return "bg_limpnn";
  if(s=="halludesign__ligandmpnn")return "halludesign";
  return s;
}
function ntasks(ati,  n,parts,i,p,a,b,s,c,st,ab){   # tasks in an ArrayTaskId string (1 if not an array)
  if(ati=="") return 1;
  n=split(ati,parts,","); c=0;
  for(i=1;i<=n;i++){ p=parts[i]; sub(/%.*/,"",p);              # strip %concurrency
    if(p ~ /-/){ s=1; if(p ~ /:/){ split(p,st,":"); s=st[2]+0; p=st[1]; } split(p,ab,"-"); a=ab[1]+0; b=ab[2]+0; if(s<1)s=1; c+=int((b-a)/s)+1; }
    else c+=1; }
  return c;
}
{
  jn=js=rs=so=aji=ati="";
  for(i=1;i<=NF;i++){ t=$i;
    if(t ~ /^JobName=/) jn=substr(t,9);
    else if(t ~ /^JobState=/) js=substr(t,10);
    else if(t ~ /^Reason=/) rs=substr(t,8);
    else if(t ~ /^StdOut=/) so=substr(t,8);
    else if(t ~ /^ArrayJobId=/) aji=substr(t,12);
    else if(t ~ /^ArrayTaskId=/) ati=substr(t,13);
  }
  nt=ntasks(ati);
  stack="?";
  if(match(so,"csy4_crrna_free/[^/]+")) stack=ss(substr(so,RSTART+16,RLENGTH-16));
  if(jn ~ /refold_af3/){ kind="refold"; role="af3"; }
  else if(jn ~ /refold_b2/){ kind="refold"; role="b2"; }
  else if(jn ~ /refold_esm/){ kind="refold"; role="esm"; }
  else { kind="gen";
    if(jn ~ /rfd3/){role="rfd3_multif"; stack="rfd3";}
    else if(jn ~ /bd1/){role="bd1_multif"; stack="boltzdesign";}
    else if(jn ~ /bg_multif/){role="bg_multif"; stack="boltzgen";}
    else if(jn ~ /^hd/){role="hd_gen"; stack="halludesign";}
    else {role=jn;}
  }
  if(js=="RUNNING") b="RUN"; else if(js=="PENDING"&&rs=="JobHeldUser") b="HELD"; else if(js=="PENDING") b="PEND"; else b="OTH";
  g=kind"|"stack"|"role;
  R[g]+=(b=="RUN")*nt; P[g]+=(b=="PEND")*nt; H[g]+=(b=="HELD")*nt; O[g]+=(b=="OTH")*nt; seen[g]=1;
  if(ati!="" && aji!=""){ ak=g SUBSEP aji; if(!(ak in aseen)){ aseen[ak]=1; A[g]++; } }
}
END{ for(g in seen) print g"|"A[g]+0"|"R[g]+0"|"P[g]+0"|"H[g]+0"|"O[g]+0; }
' | sort -t'|' -k1,1 -k2,2 -k3,3)"

{ echo "KIND|STACK|ROLE|ARRAYS|RUN|PEND|HELD|OTH"; printf '%s\n' "$ROWS"; } | column -t -s'|'

echo
printf '%s\n' "$ROWS" | awk -F'|' '{a+=$4;r+=$5;p+=$6;h+=$7;o+=$8} END{printf "TOTAL: ARRAYS=%d  RUN=%d  PEND=%d  HELD=%d  OTH=%d  (tasks=%d)\n", a,r,p,h,o, r+p+h+o}'

echo
echo "By queue (partition) × type   [RUN / PEND / HELD]:"
printf '%s\n' "$JOBS" | awk '
{
  part="?"; state=""; jn=""; rs="";
  for(i=1;i<=NF;i++){ t=$i;
    if(t~/^Partition=/)part=substr(t,11);
    else if(t~/^JobState=/)state=substr(t,10);
    else if(t~/^JobName=/)jn=substr(t,9);
    else if(t~/^Reason=/)rs=substr(t,8);
  }
  if(jn~/refold_af3/)ty="refold:af3";
  else if(jn~/refold_b2/)ty="refold:b2";
  else if(jn~/refold_esm/)ty="refold:esm";
  else if(jn~/rfd3/)ty="gen:rfd3";
  else if(jn~/bd1/)ty="gen:bd1";
  else if(jn~/bg_multif/)ty="gen:bg";
  else if(jn~/^hd/)ty="gen:hd";
  else ty="other";
  rn=(state=="RUNNING"); hd=(state=="PENDING"&&rs=="JobHeldUser"); pd=(state=="PENDING"&&rs!="JobHeldUser");
  k=part"|"ty;
  R[k]+=rn; P[k]+=pd; H[k]+=hd; seen[k]=1;
  pr[part]+=rn; pp[part]+=pd; ph[part]+=hd; pseen[part]=1;
}
END{
  # col 2 = sort order (1=type rows, 2=per-partition TOTAL last); stripped before display
  for(k in seen){ split(k,a,"|"); print a[1]"|1|  "a[2]"|"R[k]+0"|"P[k]+0"|"H[k]+0; }
  for(p in pseen) print p"|2|  TOTAL|"pr[p]+0"|"pp[p]+0"|"ph[p]+0;
}' | sort -t'|' -k1,1 -k2,2n -k3,3 | cut -d'|' -f1,3,4,5,6 \
   | { echo "PARTITION|TYPE|RUN|PEND|HELD"; cat; } | column -t -s'|'

echo
echo "Per individual queue (multi-queue jobs counted in EVERY partition they can run in):"
printf '%s\n' "$JOBS" | awk '
{
  part=""; state=""; rs="";
  for(i=1;i<=NF;i++){ t=$i;
    if(t~/^Partition=/)part=substr(t,11);
    else if(t~/^JobState=/)state=substr(t,10);
    else if(t~/^Reason=/)rs=substr(t,8);
  }
  rn=(state=="RUNNING"); hd=(state=="PENDING"&&rs=="JobHeldUser"); pd=(state=="PENDING"&&rs!="JobHeldUser");
  n=split(part,qs,",");
  for(i=1;i<=n;i++){ q=qs[i]; R[q]+=rn; P[q]+=pd; H[q]+=hd; seen[q]=1; }
}
END{ for(q in seen) print q"|"R[q]+0"|"P[q]+0"|"H[q]+0; }' \
  | sort -t'|' -k1,1 | { echo "QUEUE|RUN|PEND|HELD"; cat; } | column -t -s'|'
echo "  (a job listing N partitions is counted in all N, so columns sum to >total jobs)"
