#!/bin/bash
# aant_chunk.sbatch - protein<->RNA interface at (amino-acid x RNA nt position) for ONE
# shard, via plots_tools/interface_aa_nt_map.py (PyRosetta: H-bonds + 4.5 A heavy-atom
# contacts in one pose pass). The per-structure long CSV stays on NODE-LOCAL scratch and
# only a (stack, aa, nt) AGGREGATE is kept: the full long form would be ~177M rows / ~27 GB.
#
#   sbatch --array=0-<N>%<throttle> aant_chunk.sbatch
#SBATCH --job-name=aant
#SBATCH --partition=normal.4h
#SBATCH --time=03:59:00
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=1900M
#SBATCH --requeue
#SBATCH --open-mode=append
#SBATCH --output=/cluster/work/wicky/lvalsecchi/rna_binder/ops/logs/aant_%A_%a.log
set -uo pipefail
RB=/cluster/work/wicky/lvalsecchi/rna_binder
cd "$RB"
SIF_PR=/cluster/work/wicky/containers/pyrosetta.sif
J="${SLURM_CPUS_PER_TASK:-8}"
SH="$RB/analyses/contact_makeup/aant_shards"
OUT="$RB/analyses/contact_makeup/aant_agg"
mkdir -p "$OUT"

# ---- HARD GUARD: archive extraction MUST be node-local, never /cluster Lustre --------
SCRATCH="${TMPDIR:-}"
case "$(readlink -f "$SCRATCH" 2>/dev/null)" in /cluster/*|"") SCRATCH="" ;; esac
if [ -z "$SCRATCH" ] && [ -d /scratch ]; then
  case "$(readlink -f /scratch)" in /cluster/*) : ;; *) SCRATCH=$(mktemp -d /scratch/aant.XXXXXX);; esac
fi
if [ -z "$SCRATCH" ] || [ ! -w "$SCRATCH" ] || { case "$(readlink -f "$SCRATCH")" in /cluster/*) true;; *) false;; esac; }; then
  echo "ABORT: no node-local scratch (SCRATCH='$SCRATCH')"; exit 3
fi
export TMPDIR="$SCRATCH"
export ARCH_READ_CACHE_MAX=3

IDX=$(printf "%05d" "${SLURM_ARRAY_TASK_ID:-0}")
IN="$SH/shard_${IDX}.csv"
AGG="$OUT/agg_${IDX}.csv"
[ -s "$IN" ] || { echo "empty/missing $IN"; exit 0; }

LONG="$TMPDIR/map_${IDX}.csv"
echo "[$(date '+%F %T')] shard=$IN host=$(hostname) TMPDIR=$TMPDIR"

apptainer exec --bind /cluster/work/wicky --bind "$TMPDIR" \
  --env MPLCONFIGDIR="$TMPDIR" --env TMPDIR="$TMPDIR" --env ARCH_READ_CACHE_MAX=3 \
  "$SIF_PR" python3 analyses/plots_tools/interface_aa_nt_map.py \
  --csv "$IN" --output "$LONG" --jobs "$J" --progress
rc=$?
[ -s "$LONG" ] || { echo "FAIL: no output (rc=$rc)"; exit 1; }

# Designs attempted per stack (the honest denominator) + (stack,aa,nt) sums.
awk -F, -v OFS=, '
  FNR==NR { if (FNR>1) nin[$2]++; next }
  FNR==1  { next }
  { c[$1 OFS $5 OFS $6] += $8; h[$1 OFS $5 OFS $6] += $7; if (!( $4 in seenp )) { seenp[$4]=1; nout[$1]++ } }
  END {
    print "stack,aa,nt_pos,n_contacts,n_hbonds,n_struct_in,n_struct_scored"
    for (k in c) { split(k, a, SUBSEP); print a[1], a[2], a[3], c[k], h[k], nin[a[1]]+0, nout[a[1]]+0 }
  }
' "$IN" "$LONG" > "$AGG.tmp" && mv -f "$AGG.tmp" "$AGG"

echo "[$(date '+%F %T')] DONE $AGG rows=$(( $(wc -l < "$AGG") - 1 )) mapper_rc=$rc"
