#!/bin/bash
# full_metrics_parallel.sh — compute WHOLE-SET structural metrics on an existing
# <dir>/collected_clean.csv IN PARALLEL (chunk -> array -> merge), producing
# <dir>/collected_clean_metrics.csv. No re-collect/clean (uses the clean CSV as-is).
#
#   sbatch full_metrics_parallel.sh <analysis_dir> [N_CHUNKS=32]
# On merge completion optionally chains a combine job if COMBINE_AFTER is set to the
# combine_and_rank args, e.g. COMBINE_AFTER="--v3 analyses/midrun_v3 --v2 <dir> --out analyses/cross_campaign".
#SBATCH --job-name=fmp_prep
#SBATCH --partition=normal.4h
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=4G
#SBATCH --output=/cluster/work/wicky/lvalsecchi/rna_binder/analyses/_fmp_prep_%j.log
set -euo pipefail
RB=/cluster/work/wicky/lvalsecchi/rna_binder
cd "$RB"; export PATH="$RB/code/bin:$PATH"
DIR="${1:?usage: full_metrics_parallel.sh <analysis_dir> [N]}"
N="${2:-32}"
CLEAN="$DIR/collected_clean.csv"
CK="$DIR/_metrics_chunks"
[ -f "$CLEAN" ] || { echo "FATAL: no $CLEAN"; exit 1; }

echo "=== split $CLEAN into $N chunks $(date) ==="
rm -rf "$CK"; mkdir -p "$CK/logs"
python3 - "$CLEAN" "$CK" "$N" <<'PY'
import sys, math, pandas as pd
inp, ck, N = sys.argv[1], sys.argv[2], int(sys.argv[3])
df = pd.read_csv(inp, low_memory=False)
per = math.ceil(len(df) / N)
for i in range(N):
    df.iloc[i*per:(i+1)*per].to_csv(f"{ck}/chunk_{i:03d}.csv", index=False)
print(f"  {len(df):,} rows -> {N} chunks (~{per} each)")
PY

# --- generate the compute array (DIR baked in) ---
cat > "$CK/_compute.sbatch" <<EOF
#!/bin/bash
#SBATCH --job-name=fmp_compute
#SBATCH --partition=normal.4h
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=3G
#SBATCH --requeue
#SBATCH --output=$CK/logs/%A_%a.log
set -euo pipefail
cd $RB; export PATH="$RB/code/bin:\$PATH"
IDX=\$(printf "%03d" "\$SLURM_ARRAY_TASK_ID")
IN="$CK/chunk_\${IDX}.csv"; OUT="$CK/chunk_\${IDX}_metrics.csv"
[ -s "\$IN" ] || { echo "empty \$IN"; exit 0; }
echo "chunk \$IDX: \$(( \$(wc -l < "\$IN") - 1 )) designs on \$(hostname)"
python3 $RB/code/scripts/compute_metrics.py --input "\$IN" --output "\$OUT" --jobs "\${SLURM_CPUS_PER_TASK:-8}" --progress
EOF

# --- generate the merge (concat -> collected_clean_metrics.csv, + optional combine) ---
cat > "$CK/_merge.sbatch" <<EOF
#!/bin/bash
#SBATCH --job-name=fmp_merge
#SBATCH --partition=normal.4h
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=4G
#SBATCH --output=$DIR/_fmp_merge_%j.log
set -euo pipefail
cd $RB; export PATH="$RB/code/bin:\$PATH"
FULL="$DIR/collected_clean_metrics.csv"
# glob into an array, NOT `ls | head -1`: head exits after one line, ls takes
# SIGPIPE, and under set -euo pipefail that killed the merge with ExitCode 13
# (128+13) before printing anything. Raced on the pipe buffer, so it only bit
# once the chunk count got large (96/128). Fixed 2026-08-14.
shopt -s nullglob; _F=( $CK/chunk_*_metrics.csv ); FIRST="\${_F[0]:-}"
[ -n "\$FIRST" ] || { echo "FATAL: no chunk metrics"; exit 1; }
head -1 "\$FIRST" > "\$FULL"
n=0; for f in $CK/chunk_*_metrics.csv; do tail -n +2 "\$f" >> "\$FULL"; n=\$((n+1)); done
echo "merged \$n chunks -> \$(( \$(wc -l < "\$FULL") - 1 )) rows -> \$FULL"
rm -rf "$CK"
if [ -n "${COMBINE_AFTER:-}" ]; then
  echo "=== combine_and_rank $COMBINE_AFTER ==="
  python3 $RB/code/scripts/combine_and_rank.py $COMBINE_AFTER --top 10 --per-stack 10
fi
echo "=== fmp merge DONE \$(date) ==="
EOF

ARR=$(sbatch --parsable --array=0-$((N-1)) "$CK/_compute.sbatch")
echo "compute array: $ARR (0-$((N-1)))"
MRG=$(sbatch --parsable --dependency=afterok:"$ARR" "$CK/_merge.sbatch")
echo "merge: $MRG (afterok:$ARR)"
echo "=== fmp prep DONE $(date) ==="
