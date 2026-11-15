#!/bin/bash
# seed_vs_design.sbatch — per-SEED DockQ vs the DESIGNED backbone for the seed
# self-consistency experiment (2000 designs x 5 seeds x 3 backends = 30000 structures).
#
# seedtest_results.csv only has seed-vs-seed pairwise DockQ (self-consistency). The
# threshold-flip question ("how many of the 5 cross DockQ 0.4?") needs each seed scored
# against the design, which is what this adds.
#
# One array task = one (backend, chunk) unit from vsdesign_manifest.tsv.
#   sbatch --array=0-28%29 seed_vs_design.sbatch
#SBATCH --job-name=seedvsdesign
#SBATCH --partition=normal.4h
#SBATCH --time=00:45:00
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=1200M
#SBATCH --requeue
#SBATCH --output=/cluster/work/wicky/lvalsecchi/rna_binder/ops/logs/seedvsdesign_%A_%a.log
set -uo pipefail
RB=/cluster/work/wicky/lvalsecchi/rna_binder
D=$RB/seedtest_2026-07-30
cd "$RB"
J="${SLURM_CPUS_PER_TASK:-8}"
OUTDIR="${VSD_OUTDIR:-$D/vsdesign}"
mkdir -p "$OUTDIR"

# node-local scratch (same hard guard as the metrics jobs; never /cluster Lustre)
SCRATCH="${TMPDIR:-}"
case "$(readlink -f "$SCRATCH" 2>/dev/null)" in /cluster/*|"") SCRATCH="" ;; esac
if [ -z "$SCRATCH" ] && [ -d /scratch ]; then
  case "$(readlink -f /scratch)" in /cluster/*) : ;; *) SCRATCH=$(mktemp -d /scratch/sv.XXXXXX);; esac
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
# <design_id>_<BE>_seed<N>_sample<M>.cif  -> design_id, seed
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

PATH="$RB/code/bin:$PATH" python3 code/scripts/dockq_self.py --csv "$IN" --output "$OUT" \
  --path-col model_path --ref-col designed_cif_path --jobs "$J" --progress
rc=$?
rm -rf "$WORK"
echo "[$(date '+%F %T')] DONE rc=$rc $OUT"
exit $rc
