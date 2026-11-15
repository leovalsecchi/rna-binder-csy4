#!/bin/bash
# interface_consensus_all.sh — full 3-backend interface consensus for a master metrics CSV.
# Recompute the interface suite (Biotite contacts + Rosetta dSASA/sc) on the B2 and ESM
# refolds, then integrate the worst-case consensus (bare col = min/max across af3/b2/esm,
# AF3 preserved as <col>_af3). One command; use on any new campaign master to keep the
# consensus columns up to date (go-forward infra). The integrate backs the master up first.
#   bash interface_consensus_all.sh <master.csv> [N=48]
#
# All compute tasks carry --tmp=680000 + ARCH_READ_CACHE_MAX=3 (node-local /local guard),
# inherited from interface_be_parallel.sh.
set -euo pipefail
RB=/cluster/work/wicky/lvalsecchi/rna_binder
cd "$RB"
MASTER="${1:?usage: interface_consensus_all.sh <master.csv> [N]}"
N="${2:-48}"
[ -f "$MASTER" ] || { echo "FATAL: no $MASTER"; exit 1; }
BASE="${MASTER%.csv}_ifbe"
mkdir -p "$BASE"

echo "=== build B2+ESM manifests ==="
python3 code/scripts/build_iface_be_manifests.py "$MASTER" "$BASE"

run_be() {   # $1=be tag, $2=path-col -> prints the merge job id on stdout
  local be="$1" col="$2" out
  out=$(bash code/scripts/interface_be_parallel.sh "$BASE/iface_manifest_${be}.csv" "$col" "$be" "$N")
  echo "$out" >&2
  echo "$out" | sed -nE "s/^merge \(${be}\): ([0-9]+).*/\1/p"
}
echo "=== submit B2 + ESM recompute arrays (N=$N) ==="
MB2=$(run_be b2  b2_cif_path)
MES=$(run_be esm esm_cif_path)
[ -n "$MB2" ] && [ -n "$MES" ] || { echo "FATAL: could not capture merge job ids"; exit 1; }
echo "merge jobs: b2=$MB2 esm=$MES"

echo "=== queue integrate (afterok both merges) ==="
cat > "$BASE/_integrate.sbatch" <<EOF
#!/bin/bash
#SBATCH --job-name=ifbe_integrate
#SBATCH --partition=normal.4h
#SBATCH --time=00:30:00
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=8G
#SBATCH --output=${BASE}/_integrate_%j.log
set -euo pipefail
cd $RB
mkdir -p $BASE/be
cp $BASE/iface_manifest_b2_b2_chunks/iface_be_b2.csv    $BASE/be/
cp $BASE/iface_manifest_esm_esm_chunks/iface_be_esm.csv $BASE/be/
python3 code/scripts/interface_consensus_integrate.py "$MASTER" "$BASE/be"
python3 code/scripts/validate_consensus.py "$MASTER"
EOF
INT=$(sbatch --parsable --dependency=afterok:${MB2}:${MES} "$BASE/_integrate.sbatch")
echo "integrate job: $INT (afterok:$MB2:$MES)"
echo "=== interface_consensus_all submitted: b2=$MB2 esm=$MES integrate=$INT ==="
