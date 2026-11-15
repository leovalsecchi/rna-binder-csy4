set -u
mapfile -t IDS < <(squeue -u lvalsecchi -h -t PENDING -n refold_esm -o "%i %P" | awk '$2 ~ /gpuhe\.24h/ {print $1}')
echo "esm jobs on gpuhe.24h to move: ${#IDS[@]}"
[ "${#IDS[@]}" -eq 0 ] && { echo "none to move"; exit 0; }
for ((i=0;i<${#IDS[@]};i+=100)); do
  chunk=$(IFS=,; echo "${IDS[*]:i:100}")
  scontrol update jobid="$chunk" partition=gpuhe.bulk 2>/dev/null
done
echo "=== esm partition split after move ==="
squeue -u lvalsecchi -h -n refold_esm -o "%P" | sort | uniq -c
