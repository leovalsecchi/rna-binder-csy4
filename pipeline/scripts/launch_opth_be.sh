#!/bin/bash
# launch_opth_be.sh <be> — size + submit the full-campaign optH array for one backend.
set -euo pipefail
RB=/cluster/work/wicky/lvalsecchi/rna_binder
cd "$RB"
BE="${1:?need backend}"
case $BE in af3) COL=af3_msaoff_cif_path;; b2) COL=b2_cif_path;; esm) COL=esm_cif_path;; esac
MAN=analyses/cross_campaign/opth_manifest_full_${BE}.csv
[ -s "$MAN" ] || { echo "FATAL: manifest missing"; exit 1; }
ROWS=$(( $(wc -l < "$MAN") - 1 ))
N=$(( (ROWS + 2599) / 2600 )); [ "$N" -gt 15000 ] && N=15000; [ "$N" -lt 1 ] && N=1
echo "$BE: $ROWS rows -> N=$N tasks, UNTHROTTLED"
THROTTLE=none bash code/scripts/opth_parallel.sh "$MAN" "$COL" full_${BE} "$N"
