#!/bin/bash
# interface_be_parallel.sh — recompute the interface suite (Biotite contacts +
# Rosetta dSASA/sc/packstat) on ONE backend's CIFs into a per-backend sidecar,
# sharded across a SLURM array with fresh node-local /tmp per task. Does NOT touch
# the master. Mirrors interface_metrics_parallel.sh but threads --path-col and writes
# a backend sidecar (<manifest>_<be>_chunks/iface_be_<be>.csv).
#   sbatch interface_be_parallel.sh <iface_manifest_<be>.csv> <path_col> <be_tag> [N=40]
#SBATCH --job-name=ifbe_prep
#SBATCH --partition=normal.4h
#SBATCH --time=00:40:00
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=6G
#SBATCH --output=/cluster/work/wicky/lvalsecchi/rna_binder/analyses/_ifbe_prep_%j.log
set -euo pipefail
RB=/cluster/work/wicky/lvalsecchi/rna_binder
cd "$RB"
MAN="${1:?usage: interface_be_parallel.sh <manifest.csv> <path_col> <be_tag> [N]}"
COL="${2:?need path-col name}"
BE="${3:?need backend tag}"
N="${4:-40}"
CK="${MAN%.csv}_${BE}_chunks"
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
#SBATCH --job-name=ifbe_${BE}
#SBATCH --partition=normal.4h
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=3G
#SBATCH --tmp=680000
#SBATCH --requeue
#SBATCH --output=$CK/logs/%A_%a.log
set -uo pipefail
cd $RB
export ARCH_READ_CACHE_MAX=3          # bound node-local .tar.zst extract cache (/local leak fix)
IDX=\$(printf "%03d" "\$SLURM_ARRAY_TASK_ID")
IN="$CK/chunk_\${IDX}.csv"
[ -s "\$IN" ] || { echo "empty \$IN"; exit 0; }
J="\${SLURM_CPUS_PER_TASK:-8}"
echo "task \$IDX tmpdir=\${TMPDIR:-UNSET} ($BE) \$(date)"
echo "chunk \$IDX contacts ($BE) \$(date)"
apptainer exec --bind /cluster/work/wicky --bind /tmp --env MPLCONFIGDIR=/tmp \\
  /cluster/work/wicky/containers/boltzgen.sif \\
  env PYTHONPATH=code/src:code/scripts python3 \\
  code/scripts/contacts_sidecar.py --csv "\$IN" --path-col "$COL" --jobs "\$J" --progress
echo "chunk \$IDX rosetta ($BE) \$(date)"
apptainer exec --bind /cluster/work/wicky --bind /tmp --env MPLCONFIGDIR=/tmp \\
  /cluster/work/wicky/containers/pyrosetta.sif \\
  python3 code/scripts/interface_rosetta.py --csv "\$IN" --path-col "$COL" --resume --jobs "\$J" --progress
echo "chunk \$IDX DONE \$(date)"
EOF

cat > "$CK/_merge.sbatch" <<EOF
#!/bin/bash
#SBATCH --job-name=ifbe_${BE}_merge
#SBATCH --partition=normal.4h
#SBATCH --time=00:20:00
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=6G
#SBATCH --output=$CK/_merge_%j.log
set -euo pipefail
cd $RB
OUT="$CK/iface_be_${BE}.csv"
FIRST=\$(ls $CK/chunk_*.csv | sort | head -1)
head -1 "\$FIRST" > "\$OUT"
for f in \$(ls $CK/chunk_*.csv | sort); do tail -n +2 "\$f" >> "\$OUT"; done
echo "merged -> \$OUT (\$(( \$(wc -l < "\$OUT") - 1 )) rows)"
python3 - "\$OUT" <<'PYX'
import csv, sys
csv.field_size_limit(10**7)
rows = list(csv.DictReader(open(sys.argv[1])))
for c in ["dsasa_total","sc_value","n_heavy_contacts_4A","frac_base_contacts","worst_clash"]:
    n = sum(1 for r in rows if str(r.get(c,"")).strip() not in ("","nan"))
    print(f"  {c:20s} {n}/{len(rows)}")
PYX
EOF

ARR=$(sbatch --parsable --array=0-$((N-1)) "$CK/_compute.sbatch")
echo "compute array ($BE): $ARR (0-$((N-1)))"
MRG=$(sbatch --parsable --dependency=afterok:"$ARR" "$CK/_merge.sbatch")
echo "merge ($BE): $MRG (afterok:$ARR)"
echo "=== ifbe prep DONE ($BE) $(date) ==="
