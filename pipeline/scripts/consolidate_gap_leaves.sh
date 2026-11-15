#!/bin/bash
# consolidate_gap_leaves.sh — fold each refold_*/msaoff_<uuid>_gap leaf into the
# canonical msaoff_<uuid> leaf so leaf names match their uuid (no _gap suffix).
#   * original has a pack -> MERGE (unpack orig + gap cifs, verify count, repack)
#   * original empty/absent -> RENAME gap -> msaoff_<uuid>
# SAFE: verifies merged member count == orig+gap BEFORE removing anything; a failed
# leaf is left untouched (both dirs intact). DRYRUN=1 (default) only prints the plan.
#SBATCH --job-name=consolidate_gap
#SBATCH --partition=normal.24h
#SBATCH --time=08:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=4G
#SBATCH --output=/cluster/work/wicky/lvalsecchi/rna_binder/analyses/_consolidate_gap_%j.log
set -uo pipefail
RB=/cluster/work/wicky/lvalsecchi/rna_binder
CODE=$RB/code; export PATH="$CODE/bin:$PATH"
source "$CODE/scripts/stages/archive_pipeline.sh" || { echo "no archive_pipeline"; exit 1; }
declare -F pack_unit >/dev/null || { echo "no pack_unit"; exit 1; }
DRYRUN="${DRYRUN:-1}"
ONLY="${ONLY:-}"                 # optional: restrict to leaves whose path matches this grep
CAMP=$RB/campaigns/csy4_crrna_free
SCR="${TMPDIR:-/tmp}/consolidate.$$"; mkdir -p "$SCR"
merged=0; renamed=0; fail=0; skipped=0

# count cifs from the .list sidecar (instant, no decompress); fall back to the archive
cifcount(){ if [ -f "$1.list" ]; then grep -c "\.cif$" "$1.list"; else zstd -dc "$1" 2>/dev/null | tar -t 2>/dev/null | grep -c "\.cif$"; fi; }

for D in "$CAMP"/*/midrun_v*/refold_af3 "$CAMP"/*/midrun_v*/refold_esm; do
  [ -d "$D" ] || continue
  for L in "$D"/msaoff_*_gap; do
    [ -d "$L" ] || continue
    [ -n "$ONLY" ] && ! echo "$L" | grep -q "$ONLY" && continue
    u=$(basename "$L" | sed "s/^msaoff_//; s/_gap$//")
    O="$D/msaoff_$u"
    [ -f "$L/designs.tar.zst" ] || { echo "SKIP $L (gap has no pack)"; skipped=$((skipped+1)); continue; }
    gg=$(cifcount "$L/designs.tar.zst")

    if [ -f "$O/designs.tar.zst" ]; then
      go=$(cifcount "$O/designs.tar.zst")
      exp=$((go+gg))
      if [ "$DRYRUN" = 1 ]; then echo "MERGE  ${O#$CAMP/}  orig=$go gap=$gg -> $exp"; merged=$((merged+1)); continue; fi
      T="$SCR/m_${u}_$RANDOM"; mkdir -p "$T"
      zstd -dc "$O/designs.tar.zst" 2>/dev/null | tar -x -C "$T" 2>/dev/null
      zstd -dc "$L/designs.tar.zst" 2>/dev/null | tar -x -C "$T" 2>/dev/null
      got=$(find "$T/designs" -name "*.cif" 2>/dev/null | wc -l)
      # union count is in [max(orig,gap), orig+gap]; overlap (dup cifs) makes it < sum.
      # only a real loss (got below the larger source) is a failure.
      mx=$go; [ "$gg" -gt "$mx" ] && mx=$gg
      if [ "$got" -lt "$mx" ] || [ "$got" -gt "$exp" ]; then echo "FAIL  $u  orig=$go gap=$gg union=$got (LEFT UNTOUCHED)"; rm -rf "$T"; fail=$((fail+1)); continue; fi
      rm -f "$O/designs.tar.zst" "$O/designs.tar.zst.list"
      rm -rf "$O/designs"; mv "$T/designs" "$O/designs"
      if pack_unit "$O/designs" >/dev/null 2>&1; then rm -rf "$L" "$T"; merged=$((merged+1)); else echo "FAIL pack $u"; fail=$((fail+1)); fi
    else
      if [ "$DRYRUN" = 1 ]; then echo "RENAME ${O#$CAMP/}  gap=$gg (orig empty/absent)"; renamed=$((renamed+1)); continue; fi
      rm -rf "$O"; mv "$L" "$O"; renamed=$((renamed+1))
    fi
  done
done
rm -rf "$SCR"
echo "=== merged=$merged renamed=$renamed skipped=$skipped fail=$fail (DRYRUN=$DRYRUN) ==="
