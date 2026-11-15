#!/bin/bash
# Re-pull the enriched cross-campaign metrics CSV from Euler into data/metrics.csv.
# Run this after the interface-metrics job finishes to gain the ΔSASA/sc/packstat
# columns, then restart the app (Ctrl-C the server, `python3 app.py`).
set -e
cd "$(dirname "$0")"
SRC=/cluster/work/wicky/lvalsecchi/rna_binder/analyses/cross_campaign/app_metrics_100k.csv
ssh -o ConnectTimeout=20 euler "gzip -c $SRC" > data/metrics.csv.gz
gunzip -f data/metrics.csv.gz
echo "refreshed: $(( $(wc -l < data/metrics.csv) - 1 )) rows"
echo "live metric columns now present:"
head -1 data/metrics.csv | tr ',' '\n' | grep -E \
  "consensus_iptm|dsasa_total|sc_value|sc_value_af3|frac_base_contacts|frac_base_contacts_af3|worst_clash|consensus_n_structures" \
  | sed 's/^/  /'
echo "restart the app to pick them up:  python3 app.py"
