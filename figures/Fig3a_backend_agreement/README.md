# Fig3a_backend_agreement (new Fig 3 panel A)

AF3 (x) vs ESMFold2 (y) for ipTM, DockQ, ligand_RMSD over the whole campaign.
Run: `~/.venvs/binder_map/bin/python code/Fig3a_backend_agreement/make_figure.py` (~6 min; reads a 259 MB CSV; set FIG3A_CACHE=<file.pkl> to cache while iterating).
Output: `output/Fig3a_backend_agreement/Fig3a_backend_agreement.{svg,pdf,png}`, 178 x 37 mm (504 x 104.4 pt), bare (no letter), live Arial text in the SVG.

## Data
`inputs/large/campaign_metrics_full.csv.gz`, per-backend `_af3`/`_esm` columns only (unsuffixed iptm/sc_dockq/sc_rmsd are consensus).
BG+BIF DockQ and RMSD: TRUE-backbone per-backend values from `inputs/large/bif_true_per_backend.csv.gz`, joined in `fix_bif()`
(1,436,724 of 1,480,140 rows). The other 43,416 BG+BIF rows could not be rescored on Euler (42,447 no_native_interface,
899 unpaired / no_limpnn_uuid, 70 ref_error) and are excluded from the DockQ and RMSD sub-panels. Nothing is provisional.
ipTM is backbone-independent: all rows (minus 74 NaN).

## Reading the panel (caption text)
Dashed lines are the campaign thresholds (ipTM 0.6, DockQ 0.7, ligand_RMSD 3 A) on both axes; dotted = y = x. The % in each
quadrant is the share of drawn points (the four sum to 100.0; largest-remainder rounding, ipTM AF3-only is 30.66% true,
shown 30.6). Top-left of each panel: Pearson r. Density = log counts (~7 M points cannot be drawn singly); ligand_RMSD on log
axes, values beyond 0.5-100 A clipped into the edge bins.
**The campaign gate is NOT the "both pass" corner.** A design passed if EITHER backend passed (ipTM: the higher of the two),
so the gate's pass region is the L-shaped union of three quadrants (both + AF3 only + ESMFold2 only); only the quadrant where
both backends fail is rejected. Gate-passing designs: ipTM 49.1%, DockQ 3.0%, ligand_RMSD 2.4% of points. This is deliberately not
drawn on the panel; the caption must say it.

## Numbers
| metric | points drawn | Pearson | Spearman | opposite sides (not drawn) | both | AF3 only | ESM only | neither |
|---|---|---|---|---|---|---|---|---|
| ipTM >= 0.6 | 6,992,441 | 0.251 | 0.242 | 2,680,763 (38.34%) | 753,998 (10.78%) | 2,143,693 (30.66%) | 537,070 (7.68%) | 3,557,680 (50.88%) |
| DockQ >= 0.7 | 6,898,579 | 0.495 | 0.535 | 186,274 (2.70%) | 22,849 (0.33%) | 65,395 (0.95%) | 120,879 (1.75%) | 6,689,456 (96.97%) |
| ligand_RMSD <= 3 A | 6,898,579 | 0.568 (log10 0.539) | 0.592 | 149,808 (2.17%) | 16,211 (0.23%) | 59,740 (0.87%) | 90,068 (1.31%) | 6,732,560 (97.59%) |

Gate-passing (any backend): ipTM 3,434,761; DockQ 209,123; ligand_RMSD 166,019. Of those, only 11% (DockQ) and 10% (ligand_RMSD) pass in both backends.
Sub-panels are wider than tall (43 x 23 mm) to fit the 37 mm ceiling, so y = x is not at 45 degrees.
