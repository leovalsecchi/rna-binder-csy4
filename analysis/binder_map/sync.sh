#!/bin/bash
# Push local source/tests/web assets to Euler, where the venv with pytest/pyarrow/biotite
# lives and where the pages are baked. web/ and web_union/ carry the templates and the
# vendored fflate + 3Dmol the bake reads, so a page built on Euler needs them there too.
# Laptop authors, Euler tests. Never touches result trees or build dirs.
set -euo pipefail
FM=/cluster/work/wicky/lvalsecchi/rna_binder/analyses/foldseek_map
cd "$(dirname "$0")"
for d in src scripts tests slurm web web_union; do
  [ -d "$d" ] || continue
  rsync -az --exclude='__pycache__' --exclude='*.pyc' "$d"/ euler:"$FM/$d"/
done
for f in pyproject.toml sync.sh probe_scaling.sh submit_maps.sh submit_union.sh; do
  [ -f "$f" ] && rsync -az "$f" euler:"$FM"/
done
echo "synced -> euler:$FM"
