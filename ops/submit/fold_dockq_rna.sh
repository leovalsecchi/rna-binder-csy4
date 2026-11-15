#!/bin/bash
# fold_dockq_rna.sh — make the RNA-aligned DockQ the default in one master table.
#
# Safety order: write the replacement to <master>.new FIRST, verify its row count matches
# the original, only THEN move the original into the archive and swap the new one in. The
# original is moved (never deleted), so a failure at any step leaves a recoverable table.
#   fold_dockq_rna.sh <master.csv> <sidecar.csv> <archive_dir>
set -uo pipefail
RB=/cluster/work/wicky/lvalsecchi/rna_binder
M="$1"; S="$2"; A="$3"
SIF=/cluster/work/wicky/containers/boltzgen.sif
[ -s "$M" ] || { echo "ABORT: master missing $M"; exit 2; }
[ -s "$S" ] || { echo "ABORT: sidecar missing $S"; exit 2; }
mkdir -p "$A"

echo "[$(date '+%F %T')] folding $(basename "$M")"
BEFORE=$(wc -l < "$M")
apptainer exec --bind /cluster/work/wicky --bind /tmp --env MPLCONFIGDIR=/tmp \
  "$SIF" python3 "$RB/code/scripts/apply_dockq_rna.py" \
  --master "$M" --sidecar "$S" --out "$M.new"
rc=$?
[ $rc -eq 0 ] || { echo "ABORT: apply failed rc=$rc"; rm -f "$M.new"; exit $rc; }

AFTER=$(wc -l < "$M.new")
if [ "$BEFORE" -ne "$AFTER" ]; then
  echo "ABORT: row count changed ($BEFORE -> $AFTER); leaving master untouched"
  rm -f "$M.new"; exit 4
fi
NC_OLD=$(head -1 "$M" | tr ',' '\n' | wc -l)
NC_NEW=$(head -1 "$M.new" | tr ',' '\n' | wc -l)
echo "rows $BEFORE (unchanged); cols $NC_OLD -> $NC_NEW"

mv "$M" "$A/$(basename "$M")"          # archive the original (move, not copy)
mv "$M.new" "$M"
echo "[$(date '+%F %T')] DONE: original archived at $A/$(basename "$M")"
