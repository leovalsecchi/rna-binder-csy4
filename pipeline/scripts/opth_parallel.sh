#!/bin/bash
# opth_parallel.sh — shard the paired H-optimisation experiment across a SLURM array.
#
#   bash opth_parallel.sh <manifest.csv> <path_col> <tag> [N=64]
#
# Output: analyses/cross_campaign/opth_<tag>.csv
#
# Directives copied from groove_sasa_parallel.sh (known-good): normal.4h, 8 cpus,
# --mem-per-cpu=3G, --tmp=680000 + ARCH_READ_CACHE_MAX=3 (node-local /local extract-cache
# guard), --requeue. CPU-only, pyrosetta.sif — no GPU.
set -euo pipefail
RB=/cluster/work/wicky/lvalsecchi/rna_binder
cd "$RB"
MAN="${1:?usage: opth_parallel.sh <manifest.csv> <path_col> [N]}"
COL="${2:?need path-col name}"
TAG="${3:?need output tag}"
N="${4:-64}"
# THROTTLE=none (or empty) submits the array UNTHROTTLED (no %N cap)
THROTTLE="${THROTTLE-64}"
CK="${MAN%.csv}_chunks"
OUTCSV="$RB/analyses/cross_campaign/opth_${TAG}.csv"
[ -f "$MAN" ] || { echo "FATAL: no $MAN"; exit 1; }

rm -rf "$CK"; mkdir -p "$CK/logs"
python3 - "$MAN" "$CK" "$N" <<'PY'
import sys, math, pandas as pd
inp, ck, N = sys.argv[1], sys.argv[2], int(sys.argv[3])
df = pd.read_csv(inp, low_memory=False)
per = math.ceil(len(df) / N)
for i in range(N):
    df.iloc[i*per:(i+1)*per].to_csv(f"{ck}/chunk_{i:03d}.csv", index=False)
print(f"{len(df):,} rows -> {N} chunks (~{per} each)")
PY

cat > "$CK/_compute.sbatch" <<EOF
#!/bin/bash
#SBATCH --job-name=opth_${TAG}
#SBATCH --partition=normal.4h
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=3G
#SBATCH --tmp=680000
#SBATCH --requeue
#SBATCH --output=$CK/logs/%A_%a.log
set -uo pipefail
cd $RB
export ARCH_READ_CACHE_MAX=3
IDX=\$(printf "%03d" "\$SLURM_ARRAY_TASK_ID")
IN="$CK/chunk_\${IDX}.csv"
OUT="$CK/out_\${IDX}.csv"
[ -s "\$IN" ] || { echo "empty \$IN"; exit 0; }
J="\${SLURM_CPUS_PER_TASK:-8}"
echo "task \$IDX tmpdir=\${TMPDIR:-UNSET} on \$(hostname) \$(date)"
apptainer exec --bind /cluster/work/wicky --bind /tmp --env MPLCONFIGDIR=/tmp \\
  /cluster/work/wicky/containers/pyrosetta.sif \\
  python3 code/scripts/opth_metrics.py --csv "\$IN" --out "\$OUT" \\
  --path-col "$COL" --jobs "\$J" --progress
RC=\$?
[ \$RC -eq 0 ] || { echo "FATAL chunk \$IDX: opth_metrics exited \$RC"; exit 1; }
head -1 "\$OUT" | grep -q "buns_ia_opth" || {
  echo "FATAL chunk \$IDX: output columns missing"; exit 1; }
echo "chunk \$IDX DONE \$(date)"
EOF

cat > "$CK/_merge.sbatch" <<EOF
#!/bin/bash
#SBATCH --job-name=opth_${TAG}_merge
#SBATCH --partition=normal.4h
#SBATCH --time=00:20:00
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=4G
#SBATCH --output=$CK/_merge_%j.log
set -euo pipefail
cd $RB
OUT="$OUTCSV"
# bash globs expand PRE-SORTED; ls|sort|head races on SIGPIPE under pipefail
CHUNKS=( $CK/out_*.csv )
FIRST="\${CHUNKS[0]}"
REF=\$(head -1 "\$FIRST")
for f in "\${CHUNKS[@]}"; do
  [ "\$(head -1 "\$f")" = "\$REF" ] || { echo "FATAL: header mismatch in \$f"; exit 1; }
done
echo "\${#CHUNKS[@]} chunks, headers consistent"
echo "\$REF" > "\$OUT"
for f in "\${CHUNKS[@]}"; do tail -n +2 "\$f" >> "\$OUT"; done
echo "merged -> \$OUT (\$(( \$(wc -l < "\$OUT") - 1 )) rows)"
EOF

if [ -z "$THROTTLE" ] || [ "$THROTTLE" = "none" ]; then
  ARRSPEC="0-$((N-1))"
else
  ARRSPEC="0-$((N-1))%${THROTTLE}"
fi
ARR=$(sbatch --parsable --array="$ARRSPEC" "$CK/_compute.sbatch")
echo "compute array: $ARR ($ARRSPEC)"
MRG=$(sbatch --parsable --dependency=afterok:"$ARR" "$CK/_merge.sbatch")
echo "merge: $MRG (afterok:$ARR)"
