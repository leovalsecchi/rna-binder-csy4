#!/usr/bin/env bash
# Rebuild the code archive from the as-run tree.
#
# SRC defaults to the NAS mirror of Euler (works without VPN). Point it at
# /cluster/work/wicky/lvalsecchi/rna_binder on Euler itself for the authoritative
# copy. Re-running overwrites tracked files in place; git shows what moved.
#
# Only source-like files are copied. ops/ and analyses/plots_tools/ hold GB of run
# logs and figure output next to their scripts, so this filters by extension rather
# than copying directories whole.
set -euo pipefail

SRC="${SRC:-/Volumes/wicky/users/lvalsecchi/euler_work/rna_binder}"
DEST="${DEST:-$(cd "$(dirname "$0")/.." && pwd)}"
LOCAL_DOCS="${LOCAL_DOCS:-$HOME/docs}"

[ -d "$SRC" ] || { echo "SRC not found: $SRC" >&2; exit 1; }

# Code extensions only. No data, no logs, no images, no checkpoints.
FILTERS=(
  --include='*/'
  --include='*.py'   --include='*.sh'    --include='*.sbatch'
  --include='*.yaml' --include='*.yml'   --include='*.toml'
  --include='*.patch' --include='*.diff' --include='*.orig'
  --include='*.md'   --include='*.cfg'
  # Keep .bak_pre_* variants of SOURCE files: several are the only record of pre-fix
  # behaviour. Qualified by base extension on purpose -- a bare '*.bak_pre_*' also
  # matches data backups, and one of those is a 167 MB csv.
  --include='*.py.bak_pre_*'    --include='*.sh.bak_pre_*'
  --include='*.yaml.bak_pre_*'  --include='*.yml.bak_pre_*'
  --include='*.sbatch.bak_pre_*' --include='*.md.bak_pre_*'
  --exclude='__pycache__/' --exclude='*.pyc' --exclude='.DS_Store'
  --exclude='*'
)

# Run output sits interleaved with the scripts that produced it. ops/logs alone holds
# enough files that a bare `ls` over the NAS mirror does not return, so these are
# pruned by name BEFORE the include filters above are considered.
PRUNE=(
  --exclude='logs/' --exclude='launch_logs/' --exclude='slurm_logs/'
  --exclude='*.log' --exclude='*.slurmout' --exclude='*.err' --exclude='*.out'
  --exclude='output/' --exclude='outputs/' --exclude='results/' --exclude='_work/'
  --exclude='.pytest_cache/' --exclude='*.csv' --exclude='*.csv.gz' --exclude='*.tsv'
  # Virtualenvs, vendored web assets, and the per-tool result/data/export stores that
  # the laptop repos keep next to their source.
  --exclude='.venv/' --exclude='venv/' --exclude='node_modules/' --exclude='.git/'
  --exclude='local_results*/' --exclude='store/' --exclude='data/' --exclude='exports/'
  --exclude='*.npy' --exclude='*.npz' --exclude='*.parquet'
)

copy() {  # copy <src-subdir> <dest-subdir>
  local from="$SRC/$1" to="$DEST/$2"
  [ -d "$from" ] || { echo "  skip (absent): $1"; return 0; }
  mkdir -p "$to"
  rsync -a --prune-empty-dirs "${PRUNE[@]}" "${FILTERS[@]}" "$from/" "$to/"
  echo "  $1 -> $2  ($(find "$to" -type f | wc -l | tr -d ' ') files)"
}

echo "SRC=$SRC"
echo "DEST=$DEST"

echo "pipeline/"
copy code/src      pipeline/src
copy code/scripts  pipeline/scripts
copy code/pylibs   pipeline/pylibs

echo "configs/, patches/, ops/, analysis/, docs/"
copy code/configs            configs

# Target inputs are DATA, so copy() would filter them out (its include list is source
# extensions only). ~316 KB: 4AL5-derived coordinates, the 16-nt MSA as run, the AF3 and
# Protenix templates. Small enough to archive, and the configs are meaningless without it.
if [ -d "$SRC/targets" ]; then
  rsync -a --exclude='.DS_Store' --exclude='._*' "$SRC/targets/" "$DEST/targets/"
  echo "  targets -> targets  ($(find "$DEST/targets" -type f | wc -l | tr -d ' ') files)"
else
  echo "  skip (absent): targets"
fi
copy code/patches            patches/rna_binder
copy boltzgen_patch          patches/boltzgen_routeb

# Named explicitly rather than copying ops/ whole: ops/logs holds so many files that
# rsync's source-side enumeration does not return over the NAS mirror, even with the
# directory excluded.
for sub in agg loops refold submit watchers; do
  copy "ops/$sub" "ops/$sub"
done

copy analyses/plots_tools    analysis/plots_tools
copy campaign_kit            analysis/campaign_kit
copy docs                    docs/as_run

# Third-party patch sets live beside the project, not inside it.
for p in halludesign_upgrade_patches; do
  if [ -d "$SRC/../$p" ]; then
    mkdir -p "$DEST/patches/$p"
    rsync -a --prune-empty-dirs "${PRUNE[@]}" "${FILTERS[@]}" "$SRC/../$p/" "$DEST/patches/$p/"
    echo "  ../$p -> patches/$p"
  fi
done

# The clustering, filtering and browsing tools were developed on the laptop as their own
# git repositories. They are copied in rather than added as submodules: none has a
# remote, so a submodule reference would not resolve for anyone else. Each one's HEAD is
# recorded in provenance/local_repos.tsv so the copy can be traced back.
LOCAL_REPOS="${LOCAL_REPOS:-$HOME/Documents}"
if [ -d "$LOCAL_REPOS" ]; then
  : > "$DEST/provenance/local_repos.tsv"
  printf 'repo\tsource\thead\tdate\n' >> "$DEST/provenance/local_repos.tsv"
  for r in binder_map binder_filter metric_explorer rna_binder_gui; do
    [ -d "$LOCAL_REPOS/$r" ] || continue
    mkdir -p "$DEST/analysis/$r"
    rsync -a --prune-empty-dirs "${PRUNE[@]}" "${FILTERS[@]}" \
      --exclude='build/' --exclude='dist/' --exclude='*.egg-info/' \
      "$LOCAL_REPOS/$r/" "$DEST/analysis/$r/"
    head=$(git -C "$LOCAL_REPOS/$r" rev-parse --short HEAD 2>/dev/null || echo "untracked")
    when=$(git -C "$LOCAL_REPOS/$r" log -1 --format=%cs 2>/dev/null || echo "-")
    printf '%s\t%s\t%s\t%s\n' "$r" "$LOCAL_REPOS/$r" "$head" "$when" \
      >> "$DEST/provenance/local_repos.tsv"
    echo "  $r -> analysis/$r  (HEAD $head)"
  done
fi

# Methods cites several design/findings notes that only ever lived on the laptop.
if [ -d "$LOCAL_DOCS" ]; then
  mkdir -p "$DEST/docs/design_notes"
  rsync -a --include='*.md' --exclude='*' "$LOCAL_DOCS/" "$DEST/docs/design_notes/"
  echo "  $LOCAL_DOCS -> docs/design_notes"
fi

echo
echo "total tracked files: $(find "$DEST" -type f -not -path '*/.git/*' | wc -l | tr -d ' ')"
