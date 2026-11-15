#!/bin/bash
# v2_reap_finalize.sbatch — self-chaining inode reaper for the v2 AF3/ESM gap
# refold campaign, which fires v2_finalize.sbatch once the refold arrays drain.
#
# Each pass (~minutes on normal.4h):
#   1. REAP: for every v2 refold_{af3,esm}/msaoff_*_gap leaf that is NOT live (its
#      run.sbatch is not inside a live array's inline SCRIPTS) and NOT recently
#      written (30 min), pack_unit its loose outputs/ + inputs/ -> *.tar.zst (rm dir).
#      Workers self-pack on success, so this only catches dead/leaked leaves.
#      pack_unit is flock'd + resume-safe (a later requeue restores outputs.tar.zst).
#   2. If the v2complete arrays are still in the queue  -> re-chain (sbatch --begin +30m self).
#      Else (drained) -> sbatch v2_finalize.sbatch (collect+combine+plot) and STOP.
#   Halt at any time: touch <CAMP>/_v2_reap.STOP
#SBATCH --job-name=v2_reap_finalize
#SBATCH --partition=normal.4h
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=4G
#SBATCH --requeue
#SBATCH --output=/cluster/work/wicky/lvalsecchi/rna_binder/analyses/_v2_reap_%j.log
set +e
RB=/cluster/work/wicky/lvalsecchi/rna_binder
CODE=$RB/code
export PATH="$CODE/bin:$PATH"
source "$CODE/scripts/stages/archive_pipeline.sh" || { echo "FATAL: no archive_pipeline.sh"; exit 1; }
declare -F pack_unit >/dev/null || { echo "FATAL: pack_unit missing"; exit 1; }
CAMP=$RB/campaigns/csy4_crrna_free
SELF=$CODE/scripts/v2_reap_finalize.sbatch
STOP=$CAMP/_v2_reap.STOP
NOW=$(date +%s)

# --- live leaf run.sbatch paths: direct jobs + tasks inside live arrays' SCRIPTS ---
ACT=$(mktemp)
python3 - > "$ACT" <<'PY'
import subprocess, re
from pathlib import Path
active = set()
for o in subprocess.run(["squeue","-u","lvalsecchi","-h","-o","%o"],
                        capture_output=True, text=True).stdout.splitlines():
    o = o.strip()
    if not o:
        continue
    active.add(o)
    if o.endswith(".sbatch"):
        try:
            m = re.search(r"SCRIPTS=\((.*?)\n\)", Path(o).read_text(), re.DOTALL)
        except OSError:
            m = None
        if m:
            for ln in m.group(1).splitlines():
                ln = ln.strip()
                if ln:
                    active.add(ln)
print("\n".join(sorted(active)))
PY

reaped=0; skip_live=0; skip_recent=0
for stack in boltzgen__boltzif boltzgen__ligandmpnn boltzdesign__ligandmpnn; do
  for be in refold_af3 refold_esm; do
    for leaf in "$CAMP/$stack/midrun_v2/$be"/msaoff_*_gap; do
      [ -d "$leaf" ] || continue
      if grep -qxF "$leaf/run.sbatch" "$ACT"; then skip_live=$((skip_live+1)); continue; fi
      recent=0
      for f in "$leaf"/outputs "$leaf"/inputs "$leaf"/slurm*.log "$leaf"/*.log; do
        [ -e "$f" ] || continue
        mt=$(stat -c %Y "$f" 2>/dev/null) || continue
        [ $((NOW - mt)) -lt 1800 ] && recent=1
      done
      [ "$recent" = 1 ] && { skip_recent=$((skip_recent+1)); continue; }
      for sub in outputs inputs; do
        [ -d "$leaf/$sub" ] && { pack_unit "$leaf/$sub" >/dev/null 2>&1 && reaped=$((reaped+1)); }
      done
    done
  done
done
rm -f "$ACT"
INODES=$(lfs quota -p 510551 /cluster/work/wicky 2>/dev/null | awk 'NR==3{print $6}')
echo "$(date '+%F %T') reap: packed=$reaped skip_live=$skip_live skip_recent=$skip_recent inodes=$INODES"

if [ -f "$STOP" ]; then echo "STOP sentinel -> exit (no re-chain, no finalize)"; exit 0; fi
NACT=$(squeue -u lvalsecchi -h -n rfaf3_v2complete,rfesm_v2complete -o '%A' 2>/dev/null | wc -l)
if [ "$NACT" -gt 0 ]; then
  echo "v2complete still active ($NACT queue lines) -> re-chain in 30min"
  sbatch --begin=now+30minutes "$SELF" >/dev/null 2>&1
else
  echo "v2complete DRAINED -> submit finalize (collect+combine+plot); loop stops"
  sbatch "$CODE/scripts/v2_finalize.sbatch" >/dev/null 2>&1
fi
