#!/bin/bash
# interface_metrics_parallel.sh — add BOTH interface passes (Biotite contacts+pairwise,
# then Rosetta ΔSASA/sc/packstat) to a metrics CSV IN PLACE, sharded across a SLURM array.
# Each task processes one chunk with a fresh node-local /tmp, so the cumulative /tmp
# exhaustion that broke the single-job 30k run cannot happen. Rosetta uses --resume.
#   sbatch interface_metrics_parallel.sh <metrics.csv> [N_CHUNKS=48]
#SBATCH --job-name=ifmet_prep
#SBATCH --partition=normal.4h
#SBATCH --time=00:40:00
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=6G
#SBATCH --output=/cluster/work/wicky/lvalsecchi/rna_binder/analyses/_ifmet_prep_%j.log
set -euo pipefail
RB=/cluster/work/wicky/lvalsecchi/rna_binder
cd "$RB"
CSV="${1:?usage: interface_metrics_parallel.sh <metrics.csv> [N]}"
N="${2:-48}"
CK="${CSV%.csv}_ifmet_chunks"
[ -f "$CSV" ] || { echo "FATAL: no $CSV"; exit 1; }
BK="$CSV.bak_pre_ifmet"; [ -f "$BK" ] || cp "$CSV" "$BK"
echo "backup: $BK  ($(( $(wc -l < "$CSV") - 1 )) rows)"

rm -rf "$CK"; mkdir -p "$CK/logs"
python3 - "$CSV" "$CK" "$N" <<'PY'
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
#SBATCH --job-name=ifmet
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
export TMPDIR=/tmp
J="\${SLURM_CPUS_PER_TASK:-8}"
echo "chunk \$IDX contacts \$(date)"
apptainer exec --bind /cluster/work/wicky --bind /tmp --env MPLCONFIGDIR=/tmp \
  /cluster/work/wicky/containers/boltzgen.sif \
  env PYTHONPATH=code/src:code/scripts python3 \
  code/scripts/contacts_sidecar.py --csv "\$IN" --jobs "\$J" --progress
echo "chunk \$IDX rosetta \$(date)"
apptainer exec --bind /cluster/work/wicky --bind /tmp --env MPLCONFIGDIR=/tmp \
  /cluster/work/wicky/containers/pyrosetta.sif \
  python3 code/scripts/interface_rosetta.py --csv "\$IN" --resume --jobs "\$J" --progress
echo "chunk \$IDX DONE \$(date)"
EOF

cat > "$CK/_merge.sbatch" <<EOF
#!/bin/bash
#SBATCH --job-name=ifmet_merge
#SBATCH --partition=normal.4h
#SBATCH --time=00:30:00
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=6G
#SBATCH --output=${CSV%.csv}_ifmet_merge_%j.log
set -euo pipefail
cd $RB
FIRST=\$(ls $CK/chunk_*.csv 2>/dev/null | sort | head -1)
[ -n "\$FIRST" ] || { echo "FATAL: no chunks"; exit 1; }
TMP="$CSV.ifmet_merged"
head -1 "\$FIRST" > "\$TMP"
for f in \$(ls $CK/chunk_*.csv | sort); do tail -n +2 "\$f" >> "\$TMP"; done
mv "\$TMP" "$CSV"
echo "merged \$(ls $CK/chunk_*.csv | wc -l) chunks -> \$(( \$(wc -l < "$CSV") - 1 )) rows"
python3 - "$CSV" <<'PYX'
import csv, sys
csv.field_size_limit(10**7)
rows = list(csv.DictReader(open(sys.argv[1])))
for c in ["dsasa_total","dsasa_hphobic","sc_value","packstat","n_heavy_contacts_4A","frac_base_contacts"]:
    n = sum(1 for r in rows if str(r.get(c,"")).strip() not in ("","nan"))
    print(f"  {c:20s} {n}/{len(rows)}")
PYX
rm -rf "$CK"
echo "=== ifmet merge DONE \$(date) ==="
EOF

ARR=$(sbatch --parsable --array=0-$((N-1)) "$CK/_compute.sbatch")
echo "compute array: $ARR (0-$((N-1)))"
MRG=$(sbatch --parsable --dependency=afterok:"$ARR" "$CK/_merge.sbatch")
echo "merge: $MRG (afterok:$ARR)"
echo "=== ifmet prep DONE $(date) ==="
