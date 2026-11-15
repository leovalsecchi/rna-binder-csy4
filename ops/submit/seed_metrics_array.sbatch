#!/bin/bash
# seed_metrics_array.sbatch — SASA(base/backbone), H-bonds(base/backbone), BUNS and
# pi-stacking on the SEED replicates, to measure how reliably each metric reproduces.
#
# The seed self-consistency run gives 5 independent draws per (design, backend). Every
# interface metric we quote is computed on ONE production refold; running the same metrics
# across the replicates says how much of each number is signal and how much is the draw.
#
# One array task = one (backend, chunk) unit of vsdesign_manifest.tsv (~1000 structures).
# The chunk archive is extracted ONCE and all four tools run over the same loose CIFs —
# extraction, not compute, is what made the earlier per-design runs slow.
#
#   sbatch --array=0-28%29 seed_metrics_array.sbatch
#SBATCH --job-name=seedmetrics
#SBATCH --partition=normal.4h
#SBATCH --time=03:50:00
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=3G
#SBATCH --requeue
#SBATCH --output=/cluster/work/wicky/lvalsecchi/rna_binder/ops/logs/seedmetrics_%A_%a.log
set -uo pipefail
RB=/cluster/work/wicky/lvalsecchi/rna_binder
D=$RB/seedtest_2026-07-30
BG=/cluster/work/wicky/containers/boltzgen.sif
PR=/cluster/work/wicky/containers/pyrosetta.sif
cd "$RB"
OUT=$D/seed_metrics
mkdir -p "$OUT"

SCRATCH="${TMPDIR:-}"
case "$(readlink -f "$SCRATCH" 2>/dev/null)" in /cluster/*|"") SCRATCH="" ;; esac
if [ -z "$SCRATCH" ] && [ -d /scratch ]; then
  case "$(readlink -f /scratch)" in /cluster/*) : ;; *) SCRATCH=$(mktemp -d /scratch/sm.XXXXXX);; esac
fi
[ -n "$SCRATCH" ] && [ -w "$SCRATCH" ] || { echo "ABORT: no node-local scratch"; exit 3; }
export TMPDIR="$SCRATCH"

T="${SLURM_ARRAY_TASK_ID:-0}"
MAN="${SM_MANIFEST:-$D/vsdesign_manifest.tsv}"   # override to score the seed-6 top-up
LINE=$(sed -n "$((T+1))p" "$MAN")
[ -n "$LINE" ] || { echo "no manifest line $T"; exit 0; }
BE=$(echo "$LINE" | cut -f1)
CHUNK=$(echo "$LINE" | cut -f2)
TAG="${BE}_${CHUNK%.tar.zst}"

WORK="$TMPDIR/$TAG"
mkdir -p "$WORK/designs"
echo "[$(date '+%F %T')] $TAG host=$(hostname) TMPDIR=$TMPDIR"
tar --zstd -xf "$D/outputs/$BE/$CHUNK" -C "$WORK" || { echo "extract failed"; exit 1; }

IN="$WORK/in.csv"
python3 - "$WORK" "$BE" "$IN" <<'PY'
import csv, os, re, sys
work, be, out = sys.argv[1:4]
# af3 replicates are seeds (_AF3_seed<N>_sample0); b2/esm are diffusion samples.
pat = re.compile(r"^(?P<did>.+?)_(?:AF3|B2|ESM)(?:_seed(?P<seed>\d+))?_sample(?P<samp>\d+)\.cif$")
rows = []
root = os.path.join(work, "designs")
for fn in sorted(os.listdir(root)):
    m = pat.match(fn) if fn.endswith(".cif") else None
    if not m:
        continue
    rows.append({"design_id": f"{m.group('did')}|{int(m.group('seed') or m.group('samp'))}",
                 "base_design_id": m.group("did"), "backend": be,
                 "seed": int(m.group("seed") or m.group("samp")),
                 "model_path": os.path.join(root, fn)})
with open(out, "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=["design_id", "base_design_id", "backend", "seed",
                                      "model_path"])
    w.writeheader(); w.writerows(rows)
print(f"{len(rows)} structures", flush=True)
PY

run_bg(){ apptainer exec --bind /cluster/work/wicky --bind "$TMPDIR" --env MPLCONFIGDIR="$TMPDIR" \
  --env TMPDIR="$TMPDIR" --env PYTHONWARNINGS=ignore "$BG" \
  env PYTHONPATH=code/src:code/scripts python3 "$@"; }
run_pr(){ apptainer exec --bind /cluster/work/wicky --bind "$TMPDIR" --env MPLCONFIGDIR="$TMPDIR" \
  --env TMPDIR="$TMPDIR" --env PYTHONWARNINGS=ignore "$PR" \
  env PYTHONPATH=code/src:code/scripts python3 "$@"; }

O="$OUT/${TAG}"
[ -s "${O}_pistack.csv" ] || run_bg code/scripts/pi_stack_pernt.py \
    --csv "$IN" --output "${O}_pistack.csv" --path-col model_path --jobs 8 || echo "WARN pistack"
[ -s "${O}_sasa.csv" ]    || run_bg code/scripts/contacts_sidecar.py \
    --csv "$IN" --output "${O}_sasa.csv" --path-col model_path --jobs 8 --sasa || echo "WARN sasa"
[ -s "${O}_hbbb.csv" ]    || run_pr code/scripts/hbond_base_backbone.py \
    --csv "$IN" --output "${O}_hbbb.csv" --path-col model_path --jobs 8 || echo "WARN hbbb"
[ -s "${O}_buns.csv" ]    || run_pr code/scripts/opth_metrics.py \
    --csv "$IN" --out "${O}_buns.csv" --path-col model_path --jobs 8 || echo "WARN buns"

cp "$IN" "${O}_index.csv"
rm -rf "$WORK"
echo "[$(date '+%F %T')] DONE $TAG"
