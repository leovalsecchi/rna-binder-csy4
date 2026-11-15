#!/bin/bash
# mkworklist.sh CHUNKDIR OUTDIR OUTPREFIX WORKLIST
# Emit the zero-padded indices of chunks whose sidecar is missing or short.
set -uo pipefail
CD="$1"; OD="$2"; PFX="$3"; WL="$4"
: > "$WL"
n_tot=0; n_need=0
for CIN in "$CD"/chunk_*.csv; do
  [ -s "$CIN" ] || continue
  IDX="${CIN##*/chunk_}"; IDX="${IDX%.csv}"
  n_tot=$((n_tot+1))
  OUT="$OD/${PFX}_${IDX}.csv"
  if [ -s "$OUT" ] && [ "$(wc -l < "$OUT")" -eq "$(wc -l < "$CIN")" ]; then continue; fi
  echo "$IDX" >> "$WL"; n_need=$((n_need+1))
done
echo "$WL: $n_need of $n_tot chunks need work"
