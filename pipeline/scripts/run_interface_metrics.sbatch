#!/bin/bash
# run_interface_metrics.sbatch — add the field-standard interface metrics to a
# metrics CSV IN PLACE: exact contacts + pairwise AA table (Biotite) then Rosetta
# InterfaceAnalyzer dSASA/sc/packstat. Backs the CSV up first.
#
#   sbatch run_interface_metrics.sbatch <metrics.csv> [--no-packstat]
#SBATCH --job-name=iface_metrics
#SBATCH --partition=normal.4h
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=48
#SBATCH --mem-per-cpu=3G
#SBATCH --output=/cluster/work/wicky/lvalsecchi/rna_binder/analyses/_iface_%j.log
set -euo pipefail
RB=/cluster/work/wicky/lvalsecchi/rna_binder
cd "$RB"
CSV="${1:?usage: sbatch run_interface_metrics.sbatch <metrics.csv> [--no-packstat]}"
shift || true
EXTRA="${*:-}"
J="${SLURM_CPUS_PER_TASK:-48}"
[ -f "$CSV" ] || { echo "FATAL: no $CSV"; exit 1; }

BK="$CSV.bak_pre_interface"
[ -f "$BK" ] || cp "$CSV" "$BK"
echo "backup: $BK  ($(( $(wc -l < "$CSV") - 1 )) rows)"

echo "=== contacts_sidecar (boltzgen.sif) $(date) ==="
apptainer exec --bind /cluster/work/wicky --bind /tmp --env MPLCONFIGDIR=/tmp \
  /cluster/work/wicky/containers/boltzgen.sif \
  env PYTHONPATH=code/src:code/scripts python3 \
  code/scripts/contacts_sidecar.py --csv "$CSV" --jobs "$J" --progress

echo "=== interface_rosetta (pyrosetta.sif) $EXTRA $(date) ==="
apptainer exec --bind /cluster/work/wicky --bind /tmp --env MPLCONFIGDIR=/tmp \
  /cluster/work/wicky/containers/pyrosetta.sif \
  python3 code/scripts/interface_rosetta.py --csv "$CSV" --jobs "$J" --progress $EXTRA

echo "=== DONE $(date) ==="
python3 - "$CSV" <<'PY'
import csv, sys
csv.field_size_limit(10**7)
rows=list(csv.DictReader(open(sys.argv[1])))
def frac(c):
    n=sum(1 for r in rows if r.get(c,"") not in ("",None))
    return f"{n}/{len(rows)}"
for c in ["dsasa_total","dsasa_hphobic","sc_value","packstat","n_heavy_contacts_4A","frac_base_contacts","aa_base_ARG"]:
    print(f"  filled {c:22s} {frac(c)}")
PY
