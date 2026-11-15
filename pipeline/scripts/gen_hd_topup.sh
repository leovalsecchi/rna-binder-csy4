#!/bin/bash
# gen_hd_topup.sh — stage N MORE midrun_v3 HalluDesign gen jobs (de-novo, 18nt, 150 designs each)
# to close the gap to the 180k-design goal. Clones the live run.sbatch template, assigns each
# a FRESH uuid8 + FRESH unique seed + a (partition,time,nodelist) tuple from the SPREAD table.
#
# Modes:
#   --plan         (default) compute uuids/seeds/spread, PRINT the breakdown, write ONE sample
#                  run.sbatch to $OUT/_sample_run.sbatch. Creates NO run dirs, submits NOTHING.
#   --stage        actually create the N run dirs + run.sbatch + run_meta.json, write the
#                  submit list to $OUT/<stamp>.list. Submits NOTHING.
# Never calls sbatch. Submission is a separate, explicit step.
set -uo pipefail

RB=/cluster/work/wicky/lvalsecchi/rna_binder
CDIR=$RB/campaigns/csy4_crrna_free/halludesign__ligandmpnn
TEMPLATE=$CDIR/midrun_v3_001a4983              # known: gpuhe.4h / 04:00:00 / eu-g6 / seed 10135
OUT=$RB/.cache/hd_topup
N="${HD_TOPUP_N:-173}"
MODE="${1:---plan}"

mkdir -p "$OUT"

# ---- SPREAD table: weight | partition | time | nodelist  (recommended: throughput-weighted to gpupr) ----
# Weights are RELATIVE (normalized to N below) so N is the only knob. gpupr 50 / gpuhe 50.
# gpupr has NO QOS mem cap and good fairshare (0.79) vs gpuhe's ~0.197 wall -> new throughput lives there.
# walltimes staggered so jobs don't all expire together. All tuples are proven-valid (already running).
SPREAD=(
  "30|gpupr.24h|23:30:00|eu-a65-[02-07],eu-ts-[01-02]"
  "20|gpupr.120h|119:00:00|eu-a65-[02-07],eu-ts-[01-02]"
  "23|gpuhe.24h|23:30:00|eu-g6-[001-080]"
  "14|gpuhe.120h|119:00:00|eu-g6-[001-080]"
  "13|gpuhe.24h|16:00:00|eu-g4-[001-011,013-031]"
)

# ---- existing uuids (collision set) and max seed ----
mapfile -t EXIST < <(ls -d "$CDIR"/midrun_v3_* 2>/dev/null | sed 's#.*/midrun_v3_##')
declare -A HAVE; for u in "${EXIST[@]}"; do HAVE[$u]=1; done
# seeds live in run.sbatch as the bare line after "--seed \"; take the global max, start above it
maxseed=$(grep -A1 -h '^      --seed' "$CDIR"/midrun_v3_*/run.sbatch 2>/dev/null \
          | grep -E '^      [0-9]+' | tr -d ' \\' | sort -n | tail -1)
[ -z "$maxseed" ] && maxseed=11340
seed=$((maxseed+1))

# ---- normalize weights -> integer counts summing EXACTLY to N (largest-remainder) ----
Wsum=0; for row in "${SPREAD[@]}"; do Wsum=$((Wsum + ${row%%|*})); done
COUNTS=(); base=0
for row in "${SPREAD[@]}"; do w=${row%%|*}; c=$(( w * N / Wsum )); COUNTS+=("$c"); base=$((base+c)); done
rem=$((N - base))                                   # leftover from integer division
for ((j=0; j<rem; j++)); do COUNTS[$j]=$(( COUNTS[$j] + 1 )); done   # hand to highest-weight buckets first

# ---- expand into an ordered tuple list of length N ----
TUPLES=()
for idx in "${!SPREAD[@]}"; do
  rest=${SPREAD[$idx]#*|}
  for ((i=0; i<COUNTS[$idx]; i++)); do TUPLES+=("$rest"); done
done

# ---- generate N fresh (uuid,seed) and pair with tuples ----
UUIDS=(); SEEDS=()
made=0
while [ "$made" -lt "$N" ]; do
  u=$(openssl rand -hex 4)
  [ -n "${HAVE[$u]:-}" ] && continue
  HAVE[$u]=1
  UUIDS+=("$u"); SEEDS+=("$seed"); seed=$((seed+1)); made=$((made+1))
done

render(){  # $1=uuid $2=seed $3=partition $4=time $5=nodelist  -> run.sbatch on stdout
  local u=$1 s=$2 p=$3 t=$4 n=$5
  # ORDER MATTERS: seed (10135) first, THEN uuid (001a4983) — avoids clobbering a hex uuid that
  # might contain "10135". Template's own partition/time/nodelist are the gpuhe.4h/04:00:00/eu-g6 lines.
  sed -e "s/10135/$s/g" \
      -e "s/001a4983/$u/g" \
      -e "s|#SBATCH --partition=gpuhe.4h|#SBATCH --partition=$p|" \
      -e "s|#SBATCH --time=04:00:00|#SBATCH --time=$t|" \
      -e "s|#SBATCH --nodelist=eu-g6-\[001-080\]|#SBATCH --nodelist=$n|" \
      "$TEMPLATE/run.sbatch"
}

echo "=== HD top-up plan ==="
echo "target new uuids N=$N  | each = 150 designs = $((N*150)) designs"
echo "fresh seed range: ${SEEDS[0]} .. ${SEEDS[$((N-1))]}  (existing max was $maxseed)"
echo "template: $TEMPLATE/run.sbatch"
echo
echo "spread (count | partition | time | nodelist):"
gpupr=0; gpuhe=0
for idx in "${!SPREAD[@]}"; do
  c=${COUNTS[$idx]}; rest=${SPREAD[$idx]#*|}
  printf '  %3d | %s\n' "$c" "$rest"
  case "${SPREAD[$idx]}" in *gpupr*) gpupr=$((gpupr+c));; *gpuhe*) gpuhe=$((gpuhe+c));; esac
done
echo
echo "totals: gpupr=$gpupr ($((gpupr*100/N))%)  gpuhe=$gpuhe ($((gpuhe*100/N))%)  sum=$((gpupr+gpuhe))"

if [ "$MODE" = "--plan" ]; then
  IFS='|' read -r p t n <<< "${TUPLES[0]}"
  render "${UUIDS[0]}" "${SEEDS[0]}" "$p" "$t" "$n" > "$OUT/_sample_run.sbatch"
  echo
  echo "PLAN ONLY: no run dirs created, nothing submitted."
  echo "sample run.sbatch written to: $OUT/_sample_run.sbatch  (uuid ${UUIDS[0]}, seed ${SEEDS[0]}, $p)"
  exit 0
fi

if [ "$MODE" = "--stage" ]; then
  stamp=$(date '+%Y%m%d_%H%M%S')
  LIST="$OUT/topup_$stamp.list"
  : > "$LIST"
  for ((k=0;k<N;k++)); do
    u=${UUIDS[$k]}; s=${SEEDS[$k]}
    IFS='|' read -r p t n <<< "${TUPLES[$k]}"
    rd="$CDIR/midrun_v3_$u"
    mkdir -p "$rd"
    render "$u" "$s" "$p" "$t" "$n" > "$rd/run.sbatch"; chmod +x "$rd/run.sbatch"
    [ -f "$TEMPLATE/resolved_config.yaml" ] && cp "$TEMPLATE/resolved_config.yaml" "$rd/resolved_config.yaml"
    cat > "$rd/run_meta.json" <<JSON
{
  "uuid8": "$u",
  "label": "midrun_v3",
  "submitted": "$stamp",
  "generator": "halludesign",
  "invfold": "ligandmpnn",
  "target": "csy4_crrna_free",
  "source": "gen_hd_topup.sh",
  "overrides": { "slurm": { "gpus": 1, "gres": "gpumem:24g", "partition": "$p", "time": "$t", "nodelist": "$n" },
                 "generator_params": { "params": { "seed": $s } } }
}
JSON
    echo "$rd/run.sbatch" >> "$LIST"
  done
  echo
  echo "STAGED $N run dirs (no submit). Submit list: $LIST"
  echo "To submit (after review):  while read f; do sbatch \"\$f\"; done < $LIST"
  exit 0
fi

echo "unknown mode '$MODE' (use --plan or --stage)" >&2; exit 1
