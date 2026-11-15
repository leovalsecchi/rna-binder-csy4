#!/bin/bash
# groove_all.sh — full 3-backend RNA moiety/groove dSASA pass for a master metrics CSV.
# Computes the Biotite dSASA + groove decomposition on the af3, b2 and esm refolds, then
# integrates the worst-case consensus into <master>_groove.csv (master left untouched).
#
#   bash groove_all.sh <master.csv> [N=96]
#   THROTTLE=64 bash groove_all.sh ...   # max concurrent tasks per backend array
#
# Sibling of interface_consensus_all.sh, but SASA-only (no Rosetta) and af3 is computed as
# a sidecar too, since the master has no dsasa_rna_* columns yet.
set -euo pipefail
RB=/cluster/work/wicky/lvalsecchi/rna_binder
cd "$RB"
MASTER="${1:?usage: groove_all.sh <master.csv> [N]}"
N="${2:-96}"
export THROTTLE="${THROTTLE:-64}"
[ -f "$MASTER" ] || { echo "FATAL: no $MASTER"; exit 1; }
BASE="${MASTER%.csv}_groove"
mkdir -p "$BASE"

echo "=== build af3+b2+esm manifests ==="
python3 code/scripts/build_groove_manifests.py "$MASTER" "$BASE"

run_be() {   # $1=be tag, $2=path-col -> prints the merge job id on stdout
  local be="$1" col="$2" out
  out=$(bash code/scripts/groove_sasa_parallel.sh "$BASE/groove_manifest_${be}.csv" "$col" "$be" "$N")
  echo "$out" >&2
  echo "$out" | sed -nE "s/^merge \(${be}\): ([0-9]+).*/\1/p"
}

echo "=== submit af3 + b2 + esm arrays (N=$N, throttle=$THROTTLE) ==="
MAF=$(run_be af3 af3_msaoff_cif_path)
MB2=$(run_be b2  b2_cif_path)
MES=$(run_be esm esm_cif_path)
[ -n "$MAF" ] && [ -n "$MB2" ] && [ -n "$MES" ] || {
  echo "FATAL: could not capture all merge job ids (af3=$MAF b2=$MB2 esm=$MES)"; exit 1; }
echo "merge jobs: af3=$MAF b2=$MB2 esm=$MES"

echo "=== queue integrate (afterok all three merges) ==="
cat > "$BASE/_integrate.sbatch" <<EOF
#!/bin/bash
#SBATCH --job-name=groove_integrate
#SBATCH --partition=normal.4h
#SBATCH --time=00:30:00
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=8G
#SBATCH --output=${BASE}/_integrate_%j.log
set -euo pipefail
cd $RB
mkdir -p $BASE/side
cp $BASE/groove_manifest_af3_af3_chunks/groove_af3.csv $BASE/side/
cp $BASE/groove_manifest_b2_b2_chunks/groove_b2.csv    $BASE/side/
cp $BASE/groove_manifest_esm_esm_chunks/groove_esm.csv $BASE/side/
python3 code/scripts/groove_consensus_integrate.py "$MASTER" "$BASE/side"
EOF
INT=$(sbatch --parsable --dependency=afterok:${MAF}:${MB2}:${MES} "$BASE/_integrate.sbatch")
echo "integrate job: $INT (afterok:$MAF:$MB2:$MES)"
echo "=== groove_all submitted: af3=$MAF b2=$MB2 esm=$MES integrate=$INT ==="
echo "    output -> ${MASTER%.csv}_groove.csv"
