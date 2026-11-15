#!/bin/bash
# fullmetrics_chunk.sbatch — run the FULL interface-metric suite (minus ipsae) on ONE
# chunk of designs for ONE refold backend, writing bare metric columns into that chunk
# CSV in place. Sequential family passes share ONE node-local archive-extract cache, so
# each packed cif is unpacked once per chunk (not 5-6x). Consensus across af3/b2/esm is
# built later by the per-family *_consensus_integrate.py from the 3 per-backend chunk sets.
#
# Families run (all CPU): contacts (Biotite) + Rosetta dSASA/sc/packstat + H-bonds (FUSED,
# one pose) + structural/SSE + tip + red + DockQ.
#
# Driven by env (set at submit): FM_CHUNKDIR, FM_BACKEND (af3|b2|esm), FM_PATHCOL
# (<be>_cif_path). Array index -> $FM_CHUNKDIR/chunk_${FM_BACKEND}_<NN>.csv.
#   sbatch --array=0-<N>%<throttle> --export=ALL,FM_CHUNKDIR=...,FM_BACKEND=b2,FM_PATHCOL=b2_cif_path fullmetrics_chunk.sbatch
#SBATCH --job-name=fullmetrics
#SBATCH --partition=normal.4h
#SBATCH --time=03:59:00
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=1900M
#SBATCH --requeue
#SBATCH --output=/cluster/work/wicky/lvalsecchi/rna_binder/ops/logs/fullmetrics_%A_%a.log
set -uo pipefail
RB=/cluster/work/wicky/lvalsecchi/rna_binder
cd "$RB"
SIF_BG=/cluster/work/wicky/containers/boltzgen.sif
SIF_PR=/cluster/work/wicky/containers/pyrosetta.sif
J="${SLURM_CPUS_PER_TASK:-8}"
: "${FM_CHUNKDIR:?set FM_CHUNKDIR}"; : "${FM_BACKEND:?set FM_BACKEND}"; : "${FM_PATHCOL:?set FM_PATHCOL}"

# ---- HARD GUARD: archive extraction MUST be node-local, never /cluster Lustre --------
SCRATCH="${TMPDIR:-}"
case "$(readlink -f "$SCRATCH" 2>/dev/null)" in /cluster/*|"") SCRATCH="" ;; esac
if [ -z "$SCRATCH" ] && [ -d /scratch ]; then
  case "$(readlink -f /scratch)" in /cluster/*) : ;; *) SCRATCH=$(mktemp -d /scratch/fm.XXXXXX);; esac
fi
if [ -z "$SCRATCH" ] || [ ! -w "$SCRATCH" ] || { case "$(readlink -f "$SCRATCH")" in /cluster/*) true;; *) false;; esac; }; then
  echo "ABORT: no node-local scratch (SCRATCH='$SCRATCH')"; exit 3
fi
export TMPDIR="$SCRATCH"
export ARCH_READ_CACHE_MAX=3      # bound the node-local .tar.zst extract cache (/local leak fix)

IDX=$(printf "%05d" "${SLURM_ARRAY_TASK_ID:-0}")
CIN="$FM_CHUNKDIR/chunk_${IDX}.csv"                   # shared input (backend-agnostic path cols)
CH="$FM_CHUNKDIR/chunk_${FM_BACKEND}_${IDX}.csv"      # per-backend writable output
[ -s "$CIN" ] || { echo "empty/missing $CIN"; exit 0; }
cp -f "$CIN" "$CH"                                    # fresh per-backend copy (idempotent on requeue)
echo "[$(date '+%F %T')] backend=$FM_BACKEND col=$FM_PATHCOL chunk=$CH host=$(hostname) TMPDIR=$TMPDIR"

bg(){ apptainer exec --bind /cluster/work/wicky --bind "$TMPDIR" \
        --env MPLCONFIGDIR="$TMPDIR" --env TMPDIR="$TMPDIR" --env ARCH_READ_CACHE_MAX=3 \
        "$SIF_BG" env PYTHONPATH=code/src:code/scripts python3 "$@"; }
pr(){ apptainer exec --bind /cluster/work/wicky --bind "$TMPDIR" \
        --env MPLCONFIGDIR="$TMPDIR" --env TMPDIR="$TMPDIR" --env ARCH_READ_CACHE_MAX=3 \
        "$SIF_PR" python3 "$@"; }
stamp(){ echo "[$(date '+%F %T')] --- $1 ---"; }

stamp "contacts (biotite)"
bg code/scripts/contacts_sidecar.py --csv "$CH" --path-col "$FM_PATHCOL" --jobs "$J" --progress || echo "WARN contacts rc=$?"
stamp "structural + tip + red (compute_metrics --no-hbonds)"
bg code/scripts/compute_metrics.py --input "$CH" --output "$CH" --path-col "$FM_PATHCOL" --no-hbonds --jobs "$J" --progress || echo "WARN compute_metrics rc=$?"
stamp "pi_stack (PLIP)"
bg code/scripts/pi_stack_sidecar.py --csv "$CH" --path-col "$FM_PATHCOL" --jobs "$J" --progress || echo "WARN pistack rc=$?"
stamp "rosetta dSASA/sc/packstat + H-bonds (FUSED)"
pr code/scripts/interface_rosetta.py --csv "$CH" --path-col "$FM_PATHCOL" --with-hbonds --jobs "$J" --progress || echo "WARN rosetta rc=$?"
stamp "DockQ (host)"
python3 code/scripts/dockq_self.py --csv "$CH" --path-col "$FM_PATHCOL" --ref-col designed_cif_path --jobs "$J" --progress || echo "WARN dockq rc=$?"
echo "[$(date '+%F %T')] DONE $CH"
