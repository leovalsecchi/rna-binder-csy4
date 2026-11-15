#!/bin/bash
# v3_reagg_prep.sbatch — stage 1 of the parallel full v3 aggregation.
# collect(rfd3+hd) -> clean -> split collected_clean.csv into N chunks, then dispatch
# the compute array (stage 2) + merge (stage 3). Submitted with --dependency=afterany
# on the RFD3 v3 ESM/AF3 gap-fill arrays so it fires once those land.
#SBATCH --job-name=v3_reagg_prep
#SBATCH --partition=normal.24h
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=48
#SBATCH --mem-per-cpu=2G
#SBATCH --output=/cluster/work/wicky/lvalsecchi/rna_binder/analyses/_v3_reagg_prep_%j.log
set -euo pipefail
RB=/cluster/work/wicky/lvalsecchi/rna_binder
cd "$RB"
export PATH="$RB/code/bin:$PATH"
S=code/scripts
C=campaigns/csy4_crrna_free
O=analyses/midrun_v3
J="${SLURM_CPUS_PER_TASK:-48}"
N=32                                   # array width for the structural compute

echo "=== (1) re-collect midrun_v3 (rfd3 + halludesign) $(date) ==="
python3 "$S/collect.py" \
    "$C/rfd3__ligandmpnn/midrun_v3" \
    "$C/halludesign__ligandmpnn/midrun_v3" \
    -o "$O/collected_raw.csv" --jobs "$J" --progress

echo "=== (2) clean + sample $(date) ==="
python3 "$S/clean.py" --input "$O/collected_raw.csv" --outdir "$O"

echo "=== (3) split collected_clean.csv into $N chunks (quoting-safe) $(date) ==="
CK="$O/_metrics_chunks"
rm -rf "$CK"; mkdir -p "$CK/logs"
python3 - "$O/collected_clean.csv" "$CK" "$N" <<'PY'
import sys, math
import pandas as pd
inp, ck, N = sys.argv[1], sys.argv[2], int(sys.argv[3])
df = pd.read_csv(inp, low_memory=False)
per = math.ceil(len(df) / N)
n = 0
for i in range(N):
    part = df.iloc[i*per:(i+1)*per]
    part.to_csv(f"{ck}/chunk_{i:03d}.csv", index=False)   # header on every chunk
    if len(part): n += 1
print(f"  {len(df):,} designs -> {N} chunks (~{per} each; {n} non-empty)")
PY

echo "=== (4) dispatch compute array + merge $(date) ==="
ARR=$(sbatch --parsable --array=0-$((N-1)) "$S/v3_reagg_compute.sbatch")
echo "  compute array: $ARR (0-$((N-1)))"
MRG=$(sbatch --parsable --dependency=afterok:"$ARR" "$S/v3_reagg_merge.sbatch")
echo "  merge job: $MRG (afterok:$ARR)"
echo "=== prep DONE $(date) ==="
