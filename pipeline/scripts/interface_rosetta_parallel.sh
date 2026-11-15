#!/bin/bash
# interface_rosetta_parallel.sh — fill Rosetta interface dSASA/sc/packstat in a metrics
# CSV IN PLACE, sharded across a SLURM array. Each task gets a FRESH node-local /tmp, so
# the cumulative /tmp exhaustion that failed the single 48-worker run cannot happen.
# --resume means only blank dsasa_total rows are (re)computed; filled rows are untouched.
#   sbatch interface_rosetta_parallel.sh <metrics.csv> [N_CHUNKS=32]
#SBATCH --job-name=ifrose_prep
#SBATCH --partition=normal.4h
#SBATCH --time=00:30:00
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=4G
#SBATCH --output=/cluster/work/wicky/lvalsecchi/rna_binder/analyses/_ifrose_prep_%j.log
set -euo pipefail
RB=/cluster/work/wicky/lvalsecchi/rna_binder
cd "$RB"
CSV="${1:?usage: interface_rosetta_parallel.sh <metrics.csv> [N]}"
N="${2:-32}"
CK="${CSV%.csv}_ifrose_chunks"
[ -f "$CSV" ] || { echo "FATAL: no $CSV"; exit 1; }
BK="$CSV.bak_pre_ifrose"; [ -f "$BK" ] || cp "$CSV" "$BK"
echo "backup: $BK  ($(( $(wc -l < "$CSV") - 1 )) rows)"

rm -rf "$CK"; mkdir -p "$CK/logs"
python3 - "$CSV" "$CK" "$N" <<'PY'
import sys, math, pandas as pd
inp, ck, N = sys.argv[1], sys.argv[2], int(sys.argv[3])
df = pd.read_csv(inp, low_memory=False)
per = math.ceil(len(df) / N)
for i in range(N):
    df.iloc[i*per:(i+1)*per].to_csv(f"{ck}/chunk_{i:03d}.csv", index=False)
blank = df["dsasa_total"].isna().sum() if "dsasa_total" in df.columns else len(df)
print(f"{len(df):,} rows -> {N} chunks (~{per} each); {blank:,} blank dsasa to compute")
PY

cat > "$CK/_compute.sbatch" <<EOF
#!/bin/bash
#SBATCH --job-name=ifrose
#SBATCH --partition=normal.4h
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=3G
#SBATCH --requeue
#SBATCH --output=$CK/logs/%A_%a.log
set -uo pipefail
cd $RB
IDX=\$(printf "%03d" "\$SLURM_ARRAY_TASK_ID")
IN="$CK/chunk_\${IDX}.csv"
[ -s "\$IN" ] || { echo "empty \$IN"; exit 0; }
export TMPDIR=/tmp            # node-local; sharding bounds per-task usage
apptainer exec --bind /cluster/work/wicky --bind /tmp --env MPLCONFIGDIR=/tmp \
  /cluster/work/wicky/containers/pyrosetta.sif \
  python3 $RB/code/scripts/interface_rosetta.py --csv "\$IN" --resume \
  --jobs "\${SLURM_CPUS_PER_TASK:-8}" --progress
EOF

cat > "$CK/_merge.sbatch" <<EOF
#!/bin/bash
#SBATCH --job-name=ifrose_merge
#SBATCH --partition=normal.4h
#SBATCH --time=00:30:00
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=4G
#SBATCH --output=${CSV%.csv}_ifrose_merge_%j.log
set -euo pipefail
cd $RB
FIRST=\$(ls $CK/chunk_*.csv 2>/dev/null | head -1)
[ -n "\$FIRST" ] || { echo "FATAL: no chunks"; exit 1; }
TMP="$CSV.ifrose_merged"
head -1 "\$FIRST" > "\$TMP"
for f in \$(ls $CK/chunk_*.csv | sort); do tail -n +2 "\$f" >> "\$TMP"; done
mv "\$TMP" "$CSV"
echo "merged \$(ls $CK/chunk_*.csv | wc -l) chunks -> \$(( \$(wc -l < "$CSV") - 1 )) rows"
python3 - "$CSV" <<'PYX'
import csv, sys
csv.field_size_limit(10**7)
rows = list(csv.DictReader(open(sys.argv[1])))
for c in ["dsasa_total","dsasa_hphobic","sc_value","packstat"]:
    n = sum(1 for r in rows if str(r.get(c,"")).strip() not in ("","nan"))
    print(f"  {c:16s} {n}/{len(rows)}")
PYX
rm -rf "$CK"
echo "=== ifrose merge DONE \$(date) ==="
EOF

ARR=$(sbatch --parsable --array=0-$((N-1)) "$CK/_compute.sbatch")
echo "compute array: $ARR (0-$((N-1)))"
MRG=$(sbatch --parsable --dependency=afterok:"$ARR" "$CK/_merge.sbatch")
echo "merge: $MRG (afterok:$ARR)"
echo "=== ifrose prep DONE $(date) ==="
