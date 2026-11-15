#!/bin/bash
# cleanup_campaigns.sh — archive historical runs from campaigns/
#
# Categorizes everything in campaigns/ into:
#   ARCHIVE: smoketests + historical production runs (move to archive/)
#   DELETE: yamltest runs (already superseded by launchertest runs)
#   KEEP: launchertest runs + README + currently-pending jobs
#
# Usage:
#   ./cleanup_campaigns.sh --dry-run    # show what would happen
#   ./cleanup_campaigns.sh              # actually do it

set -euo pipefail

DRY=false
if [ "${1:-}" = "--dry-run" ]; then DRY=true; fi

PROJ=/cluster/work/wicky/lvalsecchi/rna_binder
ARCHIVE=/cluster/work/wicky/lvalsecchi/archive/rna_binder_2026-05

run() {
    if $DRY; then
        echo "DRY: $*"
    else
        echo "RUN: $*"
        eval "$@"
    fi
}

# Bail if there are pending or running jobs that touch these dirs.
PENDING=$(squeue -u "$USER" -h -o "%i %j" | wc -l)
if [ "$PENDING" -gt 0 ] && ! $DRY; then
    echo "WARNING: $PENDING active job(s) for user $USER."
    squeue -u "$USER"
    read -p "Proceed anyway? [y/N] " ans
    case "$ans" in y|Y|yes) ;; *) echo "Aborted."; exit 1 ;; esac
fi

# ------------------------------------------------------------------------
# 1. Create archive subdirs
# ------------------------------------------------------------------------
run "mkdir -p $ARCHIVE/smoketests"
run "mkdir -p $ARCHIVE/campaigns/halludesign__ligandmpnn"
run "mkdir -p $ARCHIVE/campaigns/rfd3__ligandmpnn"

# ------------------------------------------------------------------------
# 2. Archive _smoketests/ (whole dir)
# ------------------------------------------------------------------------
# Note: keeping boltzdesign1/scripts and configs intact in archive
#   because we still reference them until BD1 has a launcher path.
if [ -d "$PROJ/campaigns/_smoketests" ]; then
    for d in $PROJ/campaigns/_smoketests/*; do
        name=$(basename "$d")
        run "mv $d $ARCHIVE/smoketests/$name"
    done
    run "rmdir $PROJ/campaigns/_smoketests"
fi

# ------------------------------------------------------------------------
# 3. Archive historical halludesign__ligandmpnn runs
# ------------------------------------------------------------------------
HD_ARCHIVE=(
    "poc_66952183"
    "small_66957224"
    "medium_66957226"
    "large_66957230"
    "large_66961351"
    "large_66961351_af3_refold"
    "_logs"
)

for name in "${HD_ARCHIVE[@]}"; do
    src="$PROJ/campaigns/csy4_crrna_free/halludesign__ligandmpnn/$name"
    if [ -e "$src" ]; then
        run "mv $src $ARCHIVE/campaigns/halludesign__ligandmpnn/$name"
    fi
done

# ------------------------------------------------------------------------
# 4. Delete yamltest dirs (no historical value, superseded by launchertest)
# ------------------------------------------------------------------------
YAMLTEST_DELETE=(
    "$PROJ/campaigns/csy4_crrna_free/halludesign__ligandmpnn/yamltest_164943"
    "$PROJ/campaigns/csy4_crrna_free/halludesign__ligandmpnn/yamltest_172732"
    "$PROJ/campaigns/csy4_crrna_free/boltzgen__boltzif/yamltest_164946"
    "$PROJ/campaigns/csy4_crrna_free/boltzgen__boltzif/yamltest_171820"
    "$PROJ/campaigns/csy4_crrna_free/boltzgen__boltzif/test_spec.yaml"
    "$PROJ/campaigns/csy4_crrna_free/rfd3__ligandmpnn/yamltest_173110"
)

for path in "${YAMLTEST_DELETE[@]}"; do
    if [ -e "$path" ]; then
        run "rm -rf $path"
    fi
done

# ------------------------------------------------------------------------
# 5. Archive rfd3 launchertest as canonical reference
# ------------------------------------------------------------------------
# This was the first proven launcher run; keeping a copy in archive.
RFD3_REF="$PROJ/campaigns/csy4_crrna_free/rfd3__ligandmpnn/launchertest_b4ed9efe"
if [ -d "$RFD3_REF" ]; then
    run "cp -r $RFD3_REF $ARCHIVE/campaigns/rfd3__ligandmpnn/launchertest_b4ed9efe_REFERENCE"
fi

# ------------------------------------------------------------------------
# 6. Summary
# ------------------------------------------------------------------------
echo ""
echo "=== Cleanup summary ==="
if $DRY; then
    echo "DRY-RUN: nothing was changed. Re-run without --dry-run to apply."
else
    echo "Cleanup complete."
    echo ""
    echo "Archive:  $ARCHIVE"
    echo "Campaigns now:"
    ls -la $PROJ/campaigns/csy4_crrna_free/ 2>/dev/null | head
    for stack in $PROJ/campaigns/csy4_crrna_free/*/; do
        echo "--- $(basename $stack) ---"
        ls "$stack" 2>/dev/null | head
    done
fi
