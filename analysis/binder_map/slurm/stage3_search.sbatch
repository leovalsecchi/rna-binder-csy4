#!/bin/bash
# Stage 3: all-vs-all Foldseek, one array task per query chunk.
#   sbatch --job-name=bm_s3_esm --array=0-31%16 slurm/stage3_search.sbatch esm
#SBATCH --partition=normal.4h
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu=4G
#SBATCH --tmp=10000
#SBATCH --output=/cluster/work/wicky/lvalsecchi/rna_binder/analyses/foldseek_map/logs/%x_%A_%a.log
set -uo pipefail

LAYER="${1:?usage: stage3_search.sbatch <design|af3|esm>}"
FM=/cluster/work/wicky/lvalsecchi/rna_binder/analyses/foldseek_map
OUT="$FM/${RESULTS:-results}/$LAYER"
I=$(printf "%02d" "$SLURM_ARRAY_TASK_ID")
ALN_TYPE="${ALN_TYPE:-2}"

module load stack/2024-06 gcc/12.2.0 openmpi/4.1.6 foldseek/8-ef4e960 || exit 1
command -v foldseek >/dev/null || { echo "FATAL: foldseek not on PATH"; exit 1; }

# Node-local scratch: Foldseek tmp dirs are inode-heavy and must never land on
# the project filesystem.
WORK="${TMPDIR:-/tmp}/bm_s3_${LAYER}_${I}"
mkdir -p "$WORK/db" "$OUT/hits"
trap 'rm -rf "$WORK"' EXIT

echo "=== stage3 layer=$LAYER chunk=$I alnType=$ALN_TYPE on $(hostname) $(date) ==="

# Stage the target DB onto node-local scratch before searching, and keep the
# query sub-DB there too. 2026-09-25: 265 array tasks reading these same files
# concurrently off the shared filesystem collapsed the Lustre client -- every
# task sat in cl_sync_io_wait at ~1.5% CPU and not one of 512 finished in 105
# minutes. The probe ran the identical search from local scratch in 34 s.
# This moves bytes only: no search parameter changes, results are unaffected.
T0=$(date +%s)
cp -L "$OUT"/db/uniq* "$WORK/db/" || exit 1
TDB="$WORK/db/uniq"
[ -s "$TDB.index" ]     || { echo "FATAL: target db not staged"; exit 1; }
[ -s "${TDB}_ss.index" ] || { echo "FATAL: target _ss not staged"; exit 1; }
[ -s "${TDB}_ca.index" ] || { echo "FATAL: target _ca not staged"; exit 1; }
echo "staged target db to $WORK/db in $(( $(date +%s) - T0 ))s"

QDB="$WORK/q$I"
# Foldseek needs the 3Di (_ss) and C-alpha (_ca) siblings alongside the main
# db; search reads q$I_ss directly and dies with "does not exist" without it.
foldseek createsubdb "$OUT/db/chunks/q$I.txt" "$TDB"      "$QDB"      || exit 1
foldseek createsubdb "$OUT/db/chunks/q$I.txt" "${TDB}_ss" "${QDB}_ss" || exit 1
foldseek createsubdb "$OUT/db/chunks/q$I.txt" "${TDB}_ca" "${QDB}_ca" || exit 1
[ -s "${QDB}_ss.index" ] || { echo "FATAL: ${QDB}_ss not built"; exit 1; }

foldseek search "$QDB" "$TDB" "$WORK/aln" "$WORK/tmp" \
  -c 0.8 --cov-mode 0 -e 10 --max-seqs 300 \
  --alignment-type "$ALN_TYPE" -a --threads "${SLURM_CPUS_PER_TASK:-16}" || exit 1

foldseek convertalis "$QDB" "$TDB" "$WORK/aln" "$OUT/hits/hits_$I.tsv" \
  --format-output "query,target,qtmscore,ttmscore,alntmscore,lddt,fident,alnlen,evalue,bits" \
  --threads "${SLURM_CPUS_PER_TASK:-16}" || exit 1

echo "hits: $(wc -l < "$OUT/hits/hits_$I.tsv")"
echo "=== stage3 chunk $I done $(date) ==="
