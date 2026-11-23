# Fig5pr: precision-recall of the Fig 5 metrics, plus AF3 vs ESM agreement

Emits `output/Fig5pr_precision_recall/pr_curves.{png,pdf,svg}` and `backend_agreement.{png,pdf,svg}`.
Run `export_data.py` (writes `inputs/fig5pr_v2_*.csv`, refuses to overwrite), then `make_figure.py`.
Python: `~/.venvs/binder_map/bin/python` (numpy only; scipy/sklearn fail to import in that venv).

## pr_curves
Left, middle: step precision-recall curves per metric, AF3 and ESM scored separately (never the better of
the two); dotted line = baseline precision (prevalence). Right: average precision (AP) with 95% bootstrap
interval (within class) per metric and backend, dotted line = baseline. Log axes. Line style + colour = metric,
marker shape = backend (circle AF3, square ESM, diamond sequence).

Classes (not rebalanced, not subsampled), v2:
- **Positives, n = 15**: ordered designs with flow FC (`ONOFF_avg`) >= 1.2 (`inputs/lv2_34_1_flow_on_off.csv`, 2 null-design_id control rows dropped).
- **Negatives, n = 8,849**: all HTS non-hits with > 75 reads. Pool = `inputs/poolcraft_gates/gate3_positive.csv` rows with on-target and negative-control reads both > 75 (8,881 designs = the 32 hits + 8,849 non-hits; HITS in that file are exactly the 32 ordered, so classes are disjoint, 0 removed). Same pool with the narrower splicing < 0.07 rule = 8,836 (21 ordered removed); the broad rule adds 13. Sensitivity pools (v1 screen_points rule 9,272; screen_points broad 9,293) in `fig5pr_v2_ap_sensitivity.csv`.
- Baseline precision = 15 / 8,864 = 0.169%.

Ten metrics, a-priori direction (opposite-direction AP stored): ipTM, DockQ, interface dSASA, H-bonds, heavy-atom contacts (4 A), shape complementarity (sc_value), packstat, pi-stacking score higher = binder; ligand-RMSD and BUNS (delta_unsat_hbonds) lower = binder; plus Asn fraction of sequence. AF3 and ESM scored separately. Secondary variants (heavy contacts 4.5 A, BUNS optH, BUNS VSASA) have `primary = 0` and are not drawn. P = 3000 label permutations; designs with a missing value rank last (3 negatives lack ligand-RMSD/DockQ).
BG+BIF ligand-RMSD/DockQ use the true-backbone columns (`inputs/large/bif_true_per_backend.csv.gz`); the stored library values are wrong for those 2,400 designs.
The five recovered metrics: `inputs/fig5pr_recovered_metrics_12k.csv.gz`, pulled from Euler `analyses/cross_campaign/csy4_master_foot.csv` (per-backend cols `*_af3/_b2/_esm`); dockq/dsasa/H-bonds in that file match the 12k library table.

## backend_agreement
AF3-vs-ESM Spearman rho per metric for the 15 proposed binders, the 32 ordered, and the 12,000 library,
95% bootstrap intervals (the library's is too narrow to see). Source: `inputs/csy4_binder_library_12k_v6_k30.csv`.
Percentile of the subset rho among 2,000 random same-size library subsets is in `fig5pr_backend_agreement.csv`.

Data: csy4_binder_library_12k_v6_k30.csv, fig5pr_recovered_metrics_12k.csv.gz, lv2_34_1_flow_on_off.csv, poolcraft_gates/gate3_positive.csv, screen_points.csv. Outputs: fig5pr_v2_{pr_curves,ap_summary,ap_sensitivity,coverage,backend_agreement,class_ids}.csv (v1 files kept). pi-stacking is zero for all 15 binders in ESM (AF3: 1 of 15), so its agreement rho is undefined there.
