#!/bin/bash
# clash_farep_parallel.sh — backfill the steric-clash metrics onto a metrics CSV that
# already has the other interface columns, sharded (fresh /tmp per task). Per chunk:
#   1) contacts_sidecar (new interface_geom) -> adds n_interface_clashes/worst_clash/
#      clash_overlap_sum (recomputes the exact contact cols too; identical values)
#   2) interface_rosetta --farep-only -> adds interface_fa_rep (fast, no packstat/IA)
#   sbatch clash_farep_parallel.sh <metrics.csv> [N_CHUNKS=48]
#SBATCH --job-name=clfr_prep
#SBATCH --partition=normal.4h
#SBATCH --time=00:40:00
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=6G
#SBATCH --output=/cluster/work/wicky/lvalsecchi/rna_binder/analyses/_clfr_prep_%j.log
set -euo pipefail
RB=/cluster/work/wicky/lvalsecchi/rna_binder
cd "$RB"
CSV="${1:?usage: clash_farep_parallel.sh <metrics.csv> [N]}"
N="${2:-48}"
CK="${CSV%.csv}_clfr_chunks"
[ -f "$CSV" ] || { echo "FATAL: no $CSV"; exit 1; }
BK="$CSV.bak_pre_clfr"; [ -f "$BK" ] || cp "$CSV" "$BK"
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
#SBATCH --job-name=clfr
#SBATCH --partition=normal.4h
#SBATCH --time=02:00:00
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
apptainer exec --bind /cluster/work/wicky --bind /tmp --env MPLCONFIGDIR=/tmp \
  /cluster/work/wicky/containers/boltzgen.sif \
  env PYTHONPATH=code/src:code/scripts python3 \
  code/scripts/contacts_sidecar.py --csv "\$IN" --jobs "\$J" --progress
apptainer exec --bind /cluster/work/wicky --bind /tmp --env MPLCONFIGDIR=/tmp \
  /cluster/work/wicky/containers/pyrosetta.sif \
  python3 code/scripts/interface_rosetta.py --csv "\$IN" --resume --jobs "\$J" --progress
apptainer exec --bind /cluster/work/wicky --bind /tmp --env MPLCONFIGDIR=/tmp \
  /cluster/work/wicky/containers/pyrosetta.sif \
  python3 code/scripts/interface_rosetta.py --csv "\$IN" --farep-only --jobs "\$J" --progress
echo "chunk \$IDX DONE \$(date)"
EOF

cat > "$CK/_merge.sbatch" <<EOF
#!/bin/bash
#SBATCH --job-name=clfr_merge
#SBATCH --partition=normal.4h
#SBATCH --time=00:30:00
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=6G
#SBATCH --output=${CSV%.csv}_clfr_merge_%j.log
set -euo pipefail
cd $RB
FIRST=\$(ls $CK/chunk_*.csv 2>/dev/null | sort | head -1)
[ -n "\$FIRST" ] || { echo "FATAL: no chunks"; exit 1; }
TMP="$CSV.clfr_merged"
head -1 "\$FIRST" > "\$TMP"
for f in \$(ls $CK/chunk_*.csv | sort); do tail -n +2 "\$f" >> "\$TMP"; done
mv "\$TMP" "$CSV"
echo "merged \$(ls $CK/chunk_*.csv | wc -l) chunks -> \$(( \$(wc -l < "$CSV") - 1 )) rows"
python3 - "$CSV" <<'PYX'
import csv, sys
csv.field_size_limit(10**7)
rows = list(csv.DictReader(open(sys.argv[1])))
for c in ["n_interface_clashes","worst_clash","interface_fa_rep","dsasa_total"]:
    n = sum(1 for r in rows if str(r.get(c,"")).strip() not in ("","nan"))
    print(f"  {c:20s} {n}/{len(rows)}")
PYX
rm -rf "$CK"
echo "=== clfr merge DONE \$(date) ==="
EOF

ARR=$(sbatch --parsable --array=0-$((N-1)) "$CK/_compute.sbatch")
echo "compute array: $ARR (0-$((N-1)))"
MRG=$(sbatch --parsable --dependency=afterok:"$ARR" "$CK/_merge.sbatch")
echo "merge: $MRG (afterok:$ARR)"
echo "=== clfr prep DONE $(date) ==="
