#!/bin/bash
# pistack_pernt_array.sbatch — per-NUCLEOTIDE pi-stacking for all three refolders.
#
# Clone of pistack_consensus_array.sbatch (same chunking, same node-local hard guard), but
# running pi_stack_pernt.py so the stack can be attributed to a specific RNA base. Needed to
# drop nt 1 — the dangling 5' base that is only exposed because the hairpin was cropped.
# Each backend is a SEPARATE process so the extraction scratch frees between them.
#SBATCH --job-name=pistack_pernt
#SBATCH --array=0-19
#SBATCH --partition=normal.4h
#SBATCH --time=03:59:00
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=3G
#SBATCH --output=/cluster/work/wicky/lvalsecchi/rna_binder/ops/logs/pistack_pernt_%A_%a.log
set -uo pipefail
RB=/cluster/work/wicky/lvalsecchi/rna_binder
cd "$RB"

SCRATCH="${TMPDIR:-}"
case "$(readlink -f "$SCRATCH" 2>/dev/null)" in /cluster/*|"") SCRATCH="" ;; esac
if [ -z "$SCRATCH" ] && [ -d /scratch ]; then
  case "$(readlink -f /scratch)" in /cluster/*) : ;; *) SCRATCH=$(mktemp -d /scratch/pisnt.XXXXXX);; esac
fi
if [ -z "$SCRATCH" ] || [ ! -w "$SCRATCH" ]; then echo "ABORT: no node-local scratch"; exit 3; fi
export TMPDIR="$SCRATCH"
export ARCH_READ_CACHE_MAX=4

T="$SLURM_ARRAY_TASK_ID"
CH=$(printf "%s/analyses/cross_campaign/pernt_chunks/chunk_%02d.csv" "$RB" "$T")
OUT="$RB/analyses/cross_campaign/pernt_out"
mkdir -p "$OUT"
[ -s "$CH" ] || { echo "no chunk $CH"; exit 0; }
echo "[$(date '+%F %T')] task $T -> $CH  TMPDIR=$TMPDIR  host=$(hostname)"

SIF=/cluster/work/wicky/containers/boltzgen.sif
run() {  # $1 = path column, $2 = backend tag
  local o
  o=$(printf "%s/chunk_%02d_%s.csv" "$OUT" "$T" "$2")
  if [ -s "$o" ]; then echo "  $o exists — skip"; return 0; fi
  apptainer exec --bind /cluster/work/wicky --bind "$TMPDIR" --env MPLCONFIGDIR="$TMPDIR" \
    --env TMPDIR="$TMPDIR" --env ARCH_READ_CACHE_MAX=4 --env PYTHONWARNINGS=ignore \
    "$SIF" env PYTHONPATH=code/src:code/scripts python3 code/scripts/pi_stack_pernt.py \
    --csv "$CH" --output "$o" --path-col "$1" --jobs 8 || echo "WARN rc=$? $1"
}
run af3_msaoff_cif_path af3
run b2_cif_path         b2
run esm_cif_path        esm
echo "[$(date '+%F %T')] DONE task $T"
