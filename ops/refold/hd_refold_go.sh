#!/usr/bin/env bash
# hd_refold_go.sh — submit HD AF3 (-> gpuhe.bulk) and HD B2 (-> gpu.bulk) refolds
# for all HalluDesign packed designs lacking a leaf. B2 lands on gpu.24h at submit
# (cluster remaps *.bulk -> *.24h), then we scontrol-move each to gpu.bulk
# (scontrol bypasses the submit-time remap).
set -u
RB=/cluster/work/wicky/lvalsecchi/rna_binder
echo "===== HD AF3 -> gpuhe.bulk ($(date)) ====="
PART=bulk bash "$RB/ops/submit/submit_af3_backfill.sh" --stack halludesign__ligandmpnn --go
echo "===== HD B2 submit (lands gpu.24h via remap) ====="
b2out=$(PART=gpu.bulk bash "$RB/ops/submit/submit_b2_backfill.sh" --stack halludesign__ligandmpnn --go)
echo "$b2out"
echo "===== scontrol-move new HD B2 jobs onto gpu.bulk ====="
jids=$(echo "$b2out" | grep -oE "> [0-9]+" | grep -oE "[0-9]+")
moved=0
for j in $jids; do
  scontrol update job "$j" partition=gpu.bulk 2>/dev/null && moved=$((moved+1))
done
echo "moved $moved B2 jobs to gpu.bulk"
echo "DONE_HD_REFOLD $(date)"
