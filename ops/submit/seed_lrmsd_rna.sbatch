#!/bin/bash
# seed_lrmsd_rna.sbatch — RNA-aligned LRMSD for the per-SEED refolds vs the designed backbone.
#
# Same array shape and staging as seed_vs_design.sbatch (one task = one backend x chunk of
# vsdesign_manifest.tsv), but computes BOTH LRMSD directions instead of DockQ, so the stored
# per-seed DockQ can be rescored target-fixed by ratio anchoring.
#
#   sbatch --array=0-28%29 seed_lrmsd_rna.sbatch
#SBATCH --job-name=seedlrmsdrna
#SBATCH --partition=normal.4h
#SBATCH --time=00:45:00
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=600M
#SBATCH --requeue
#SBATCH --output=/cluster/work/wicky/lvalsecchi/rna_binder/ops/logs/seedlrmsdrna_%A_%a.log
set -uo pipefail
RB=/cluster/work/wicky/lvalsecchi/rna_binder
D=$RB/seedtest_2026-07-30
SIF=/cluster/work/wicky/containers/boltzgen.sif
cd "$RB"
J="${SLURM_CPUS_PER_TASK:-8}"
OUTDIR="${VSD_OUTDIR:-$D/vsdesign_lrmsdrna}"
mkdir -p "$OUTDIR"

# node-local scratch (same hard guard as the metrics jobs; never /cluster Lustre)
SCRATCH="${TMPDIR:-}"
case "$(readlink -f "$SCRATCH" 2>/dev/null)" in /cluster/*|"") SCRATCH="" ;; esac
if [ -z "$SCRATCH" ] && [ -d /scratch ]; then
  case "$(readlink -f /scratch)" in /cluster/*) : ;; *) SCRATCH=$(mktemp -d /scratch/sl.XXXXXX);; esac
fi
[ -n "$SCRATCH" ] && [ -w "$SCRATCH" ] || { echo "ABORT: no node-local scratch"; exit 3; }
export TMPDIR="$SCRATCH"

T="${SLURM_ARRAY_TASK_ID:-0}"
LINE=$(sed -n "$((T+1))p" "${VSD_MANIFEST:-$D/vsdesign_manifest.tsv}")
[ -n "$LINE" ] || { echo "no manifest line $T"; exit 0; }
BE=$(echo "$LINE" | cut -f1)
CHUNK=$(echo "$LINE" | cut -f2)
TAG="${BE}_${CHUNK%.tar.zst}"
OUT="$OUTDIR/${TAG}.csv"
[ -s "$OUT" ] && { echo "[$(date '+%F %T')] $OUT exists — skip"; exit 0; }

WORK="$TMPDIR/$TAG"
mkdir -p "$WORK/designs"
echo "[$(date '+%F %T')] $TAG host=$(hostname) TMPDIR=$TMPDIR"
tar --zstd -xf "$D/outputs/$BE/$CHUNK" -C "$WORK" || { echo "extract failed"; exit 1; }

IN="$WORK/in.csv"
python3 - "$WORK" "$D/design_refs.csv" "$BE" "$IN" <<'PY'
import csv, os, re, sys
work, refs_csv, be, out = sys.argv[1:5]
refs = {r["design_id"]: r["designed_cif_path"] for r in csv.DictReader(open(refs_csv))}
# af3 replicates are seeds (_AF3_seed<N>_sample0); b2/esm are diffusion samples
# (_B2_sample<N> / _ESM_sample<N>). Both are the 5 replicates -> "seed" index.
pat = re.compile(r"^(?P<did>.+?)_(?:AF3|B2|ESM)(?:_seed(?P<seed>\d+))?_sample(?P<samp>\d+)\.cif$")
rows, miss = [], 0
root = os.path.join(work, "designs")
for fn in sorted(os.listdir(root)):
    if not fn.endswith(".cif"):
        continue
    m = pat.match(fn)
    if not m:
        continue
    did = m.group("did")
    ref = refs.get(did)
    if not ref:
        miss += 1
        continue
    rows.append({"design_id": did, "backend": be, "seed": int(m.group("seed") or m.group("samp")),
                 "model_path": os.path.join(root, fn), "designed_cif_path": ref})
with open(out, "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=["design_id", "backend", "seed",
                                      "model_path", "designed_cif_path"])
    w.writeheader()
    w.writerows(rows)
print(f"{len(rows)} structures, {miss} without a design ref", flush=True)
PY

apptainer exec --bind /cluster/work/wicky --bind "$TMPDIR" \
  --env MPLCONFIGDIR="$TMPDIR" --env TMPDIR="$TMPDIR" --env PYTHONWARNINGS=ignore \
  "$SIF" env PYTHONPATH=code/src:code/scripts python3 \
  code/scripts/seed_lrmsd_rna.py --csv "$IN" --output "$OUT" --jobs "$J"
rc=$?
rm -rf "$WORK"
echo "[$(date '+%F %T')] DONE rc=$rc $OUT"
exit $rc
