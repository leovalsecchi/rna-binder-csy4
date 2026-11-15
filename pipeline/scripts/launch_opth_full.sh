#!/bin/bash
# launch_opth_full.sh — sizes and submits the full-campaign optH array once the
# manifest exists. ~2600 structures/task (~80 min at 8 workers, 4h walltime),
# UNTHROTTLED.
set -euo pipefail
RB=/cluster/work/wicky/lvalsecchi/rna_binder
cd "$RB"
MAN=analyses/cross_campaign/opth_manifest_full_af3.csv
[ -s "$MAN" ] || { echo "FATAL: manifest missing"; exit 1; }
ROWS=$(( $(wc -l < "$MAN") - 1 ))
N=$(( (ROWS + 2599) / 2600 ))
[ "$N" -gt 15000 ] && N=15000
[ "$N" -lt 1 ] && N=1
echo "manifest rows: $ROWS -> N=$N tasks (~$(( ROWS / N )) each), UNTHROTTLED"
THROTTLE=none bash code/scripts/opth_parallel.sh "$MAN" af3_msaoff_cif_path full_af3 "$N"
