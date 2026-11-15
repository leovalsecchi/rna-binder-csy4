#!/bin/bash
# basil_full_dockq_driver.sh — after the 3 per-campaign collects finish: merge -> clean
# -> chunk -> fan out FULL-SET DockQ (af3/b2/esm) as a throttled Slurm array -> merge into
# a reusable sidecar analyses/basil_discussion/collected_clean_dockq.csv.
# Submit with a dependency on the collect jobs:
#   sbatch --dependency=afterany:<v3>:<v4>:<v5> basil_full_dockq_driver.sh
#SBATCH --job-name=basil_dq_prep
#SBATCH --partition=normal.4h
#SBATCH --time=03:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu=4G
#SBATCH --output=/cluster/work/wicky/lvalsecchi/rna_binder/analyses/basil_discussion/_dq_driver_%j.log
set -euo pipefail
RB=/cluster/work/wicky/lvalsecchi/rna_binder
cd "$RB"; export PATH="$RB/code/bin:$PATH"
DIR="$RB/analyses/basil_discussion"; S="$RB/code/scripts"

echo "=== (1) merge per-campaign raw $(date) ==="
python3 - <<PY
import pandas as pd, glob
fs = sorted(glob.glob("$DIR/collected_raw_v*.csv"))
print("merging:", fs)
df = pd.concat([pd.read_csv(f, low_memory=False) for f in fs], ignore_index=True)
df.to_csv("$DIR/collected_raw.csv", index=False)
print(f"  merged {len(df):,} designs -> collected_raw.csv")
print("  per campaign:", {f.split('_')[-1].replace('.csv',''): len(pd.read_csv(f, low_memory=False)) for f in fs})
PY

echo "=== (2) clean (full, no downsampling) $(date) ==="
python3 "$S/clean.py" --input "$DIR/collected_raw.csv" --outdir "$DIR"
CLEAN="$DIR/collected_clean.csv"
ROWS=$(( $(wc -l < "$CLEAN") - 1 ))
echo "  clean set: $ROWS designs"
echo "  per-stack clean counts:"
python3 - <<PY
import pandas as pd
d = pd.read_csv("$CLEAN", low_memory=False)
print(d.groupby("stack").size().sort_values().to_string())
PY

echo "=== (3) chunk + fan out DockQ $(date) ==="
N=$(( (ROWS + 2999) / 3000 )); [ "$N" -lt 16 ] && N=16; [ "$N" -gt 256 ] && N=256
CK="$DIR/_dockq_chunks"; rm -rf "$CK"; mkdir -p "$CK/logs"
python3 - "$CLEAN" "$CK" "$N" <<'PY'
import sys, math, pandas as pd
inp, ck, N = sys.argv[1], sys.argv[2], int(sys.argv[3])
df = pd.read_csv(inp, low_memory=False)
per = math.ceil(len(df) / N)
w = 0
for i in range(N):
    sub = df.iloc[i*per:(i+1)*per]
    if len(sub):
        sub.to_csv(f"{ck}/chunk_{i:03d}.csv", index=False); w += 1
print(f"  {len(df):,} rows -> {w} chunks (~{per} each)")
PY

cat > "$CK/_compute.sbatch" <<EOF
#!/bin/bash
#SBATCH --job-name=basil_dq
#SBATCH --partition=normal.4h
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=3G
#SBATCH --requeue
#SBATCH --open-mode=append
#SBATCH --output=$CK/logs/%A_%a.log
set -euo pipefail
cd $RB; export PATH="$RB/code/bin:\$PATH"
IDX=\$(printf "%03d" "\$SLURM_ARRAY_TASK_ID")
IN="$CK/chunk_\${IDX}.csv"; OUT="$CK/chunk_\${IDX}_dockq.csv"
[ -s "\$IN" ] || { echo "no \$IN"; exit 0; }
[ -s "\$OUT" ] && { echo "\$OUT exists -> skip"; exit 0; }
echo "chunk \$IDX: \$(( \$(wc -l < "\$IN") - 1 )) designs on \$(hostname) \$(date)"
python3 $RB/analyses/plots_tools/dockq_refolds.py --input "\$IN" --out "\$OUT" --jobs "\${SLURM_CPUS_PER_TASK:-8}"
EOF

cat > "$CK/_merge.sbatch" <<EOF
#!/bin/bash
#SBATCH --job-name=basil_dq_merge
#SBATCH --partition=normal.4h
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=4G
#SBATCH --output=$DIR/_dq_merge_%j.log
set -euo pipefail
FULL="$DIR/collected_clean_dockq.csv"
FIRST=\$(ls $CK/chunk_*_dockq.csv 2>/dev/null | head -1)
[ -n "\$FIRST" ] || { echo "FATAL: no chunk dockq outputs"; exit 1; }
head -1 "\$FIRST" > "\$FULL"
have=0; for f in $CK/chunk_*_dockq.csv; do tail -n +2 "\$f" >> "\$FULL"; have=\$((have+1)); done
done_chunks=\$(ls $CK/chunk_*.csv 2>/dev/null | grep -vc _dockq || true)
echo "merged \$have chunk-dockq files -> \$(( \$(wc -l < "\$FULL") - 1 )) rows -> \$FULL"
echo "(source chunks kept in $CK for any re-run of missing ones)"
echo "=== basil DockQ merge DONE \$(date) ==="
EOF

ARR=$(sbatch --parsable --array=0-$((N-1))%40 "$CK/_compute.sbatch")
echo "DockQ compute array: $ARR ($N chunks, 40 concurrent)"
MRG=$(sbatch --parsable --dependency=afterany:"$ARR" "$CK/_merge.sbatch")
echo "DockQ merge: $MRG (afterany:$ARR)"
echo "=== basil DockQ driver prep DONE $(date) — array=$ARR merge=$MRG ==="
