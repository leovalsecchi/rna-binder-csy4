# Csy4 RNA-binder — plotting suite

Single-command suite: 10 plots, each an independently-runnable script producing
**both** a static PNG (matplotlib) and an interactive HTML (plotly). All thresholds,
colors, and paths live in `config.py` and nowhere else.

## Layout
```
analyses/plots_tools/        # this code + sidecars/
analyses/plots/<DATASET>/    # outputs: png/ + html/  (DATASET set in config.py; here combined_v2_v3)
analyses/archive/            # retired one-off analyses + logs
```

## Run
```bash
cd analyses/plots_tools
python3 run_all.py                       # all 10 → ../plots/<DATASET>/{png,html}/ + index.html
python3 plot_04_scatter_iptm_dockq.py    # any single plot, standalone
python3 run_all.py --stacks rfd3+ligandmpnn boltzgen+boltzif   # subset
python3 run_all.py --master <other.csv> --outdir <dir>         # different data / output
```

## Files
| file | role |
|---|---|
| `config.py` | **single source of truth** — THRESHOLDS, STACK_COLORS (Okabe-Ito, CVD-safe, 2 blues for the BoltzGen pair), REFOLD_CONDITIONS (MSA-off only), paths |
| `data_loader.py` | wide CSV → long (design × refold_condition); MSA-off filter; column mapping; per-nt interface contacts; sidecar merges; `passing_designs()`/`funnel()` |
| `plot_helpers.py` | shared house style, `save_pair`, threshold lines, `pending_pair`, `funnel_figure` |
| `plot_01..10_*.py` | one per spec plot; each exposes `generate(df, cfg) → (mpl_fig, plotly_fig)` |
| `run_all.py` | entry point (`--stacks`, `--outdir`, `--master`) + writes `output/html/index.html` |
| `dockq_refolds.py` | computes DockQ(refold vs designed) per refold condition → `sidecars/dockq_10k.csv` |

## Data model
The aggregated CSV is wide (one row/design). `data_loader` melts it to **one row per
(design × refold_condition)** with conditions `af3_msa_off` and `b2_msa_off` (MSA-off
only — ESM and any MSA-on are dropped here, never per-plot). Per-refold metrics
(iptm, ptm, dockq) differ by row; design-level metrics (protein_rmsd, helix/sheet,
n_hbonds, interface_contacts) repeat on both rows.

## Metric status (2026-07-03, on the 3k sample)
- **live:** iptm, ptm (af3), **dockq** (3k sidecar), protein_rmsd, n_hbonds, helix_frac,
  sheet_frac, interface_contacts
- **pending → graceful placeholders / skipped filters:** pae, fraction_disordered,
  n_clashes, per-nt hbond_contacts (plot 6), b2 ptm. Each auto-populates when its
  column/sidecar appears — no plot code change.

## DockQ: 3k now → 900k later
`dockq_refolds.py` wraps `code/scripts/dockq_self.py` once per refold condition. It's
run on the 3k sample today. To scale to the full ~900k set: point `--input` at chunks
of `combined_clean_metrics.csv` and fan out with the same array pattern as
`full_metrics_parallel.sh`, then flip `config.MASTER_CSV` to the full table.
