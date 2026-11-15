#!/bin/bash
# Stratified ~100k sample from the 1.75M combined_clean_metrics.csv (min(PER, n) per
# stack), shuffled, -> analyses/cross_campaign/sample100k.csv. Base metrics + AF3 refold
# paths carried through; interface metrics + DockQ added by later steps.
#   sbatch sample100k.sbatch [PER_STACK=20000]
#SBATCH --job-name=sample100k
#SBATCH --partition=normal.4h
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=24G
#SBATCH --output=/cluster/work/wicky/lvalsecchi/rna_binder/analyses/_sample100k_%j.log
set -euo pipefail
RB=/cluster/work/wicky/lvalsecchi/rna_binder
cd "$RB"
PER="${1:-20000}"
python3 - analyses/cross_campaign/combined_clean_metrics.csv \
          analyses/cross_campaign/sample100k.csv "$PER" <<'PY'
import sys, pandas as pd
C, OUT, PER = sys.argv[1], sys.argv[2], int(sys.argv[3])
print("reading", C, flush=True)
df = pd.read_csv(C, low_memory=False)
print("full:", df.shape, "| stacks:", df["stack"].value_counts().to_dict(), flush=True)
parts = [g.sample(n=min(PER, len(g)), random_state=0) for _, g in df.groupby("stack")]
samp = pd.concat(parts).sample(frac=1, random_state=0).reset_index(drop=True)
samp.to_csv(OUT, index=False)
print("wrote", OUT, samp.shape, "| per stack:", samp["stack"].value_counts().to_dict())
PY
echo "=== sample100k DONE $(date) ==="
