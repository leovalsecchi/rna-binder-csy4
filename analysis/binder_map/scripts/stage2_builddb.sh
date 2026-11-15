#!/bin/bash
# Stage 2: global dedup, then build the Foldseek DB from the tar.gz shards.
#   bash scripts/stage2_builddb.sh <design|af3|esm>
set -uo pipefail

LAYER="${1:?usage: stage2_builddb.sh <design|af3|esm>}"
FM=/cluster/work/wicky/lvalsecchi/rna_binder/analyses/foldseek_map
PY=/cluster/work/wicky/lvalsecchi/binder_map_venv/bin/python
OUT="$FM/${RESULTS:-results}/$LAYER"

module load stack/2024-06 gcc/12.2.0 openmpi/4.1.6 foldseek/8-ef4e960 || exit 1
command -v foldseek >/dev/null || { echo "FATAL: foldseek not on PATH"; exit 1; }

cd "$FM" || exit 1
PYTHONPATH=src "$PY" scripts/stage2_builddb.py --layer "$LAYER" \
  --results "$FM/${RESULTS:-results}" || exit 1

mkdir -p "$OUT/db"
rm -rf "$OUT/db/full"* "$OUT/db/uniq"*

echo "=== createdb from $(ls "$OUT"/shards/*.tar.gz | wc -l) shards ==="
foldseek createdb "$OUT"/shards/*.tar.gz "$OUT/db/full" || exit 1

[ -s "$OUT/dedup/keep_names.txt" ] || {
  echo "FATAL: $OUT/dedup/keep_names.txt missing or empty."
  echo "stage2_builddb.py wrote its dedup tables somewhere else -- check --results."
  exit 1
}
FULL=$(wc -l < "$OUT/db/full.lookup")
KEEP=$(wc -l < "$OUT/dedup/keep_names.txt")
echo "full db: $FULL entries; want $KEEP canonical"

# Map canonical names to DB keys via the lookup, then subset.
awk 'NR==FNR{want[$1]=1; next} ($2 in want){print $1}' \
  "$OUT/dedup/keep_names.txt" "$OUT/db/full.lookup" > "$OUT/db/keep_keys.txt"

FOUND=$(wc -l < "$OUT/db/keep_keys.txt")
echo "resolved $FOUND / $KEEP canonical names to DB keys"
if [ "$FOUND" -eq 0 ]; then
  echo "FATAL: no canonical name matched a DB key. createdb stores names as"
  echo "<design_id>.pdb; compare db/full.lookup column 2 against dedup/keep_names.txt."
  exit 1
fi
if [ "$FOUND" -ne "$KEEP" ]; then
  echo "FATAL: $((KEEP - FOUND)) canonical structures are not in the Foldseek DB."
  echo "That means createdb silently skipped them. Investigate before searching."
  exit 1
fi

foldseek createsubdb "$OUT/db/keep_keys.txt" "$OUT/db/full" "$OUT/db/uniq" || exit 1
foldseek createsubdb "$OUT/db/keep_keys.txt" "$OUT/db/full_ss" "$OUT/db/uniq_ss" 2>/dev/null
foldseek createsubdb "$OUT/db/keep_keys.txt" "$OUT/db/full_ca" "$OUT/db/uniq_ca" 2>/dev/null

echo "uniq db: $(wc -l < "$OUT/db/uniq.lookup") entries"
echo "=== stage 2 done for $LAYER ==="
