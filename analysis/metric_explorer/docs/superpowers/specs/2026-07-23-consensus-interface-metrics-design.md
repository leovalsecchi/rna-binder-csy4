# Consensus interface metrics — 3-backend worst-case aggregation

**Date:** 2026-07-23
**Status:** design approved (v2 — aligned to house convention), ready for plan
**Primary home:** `rna_binder` metric-calc pipeline on Euler (source of truth).
**Consumer:** `~/Documents/metric_explorer` (inherits columns via `refresh_data.sh`).

## Goal

Extend the interface/structural design metrics that are currently AF3-only to a
**worst-case consensus across the three refold backends (AF3 / Boltz-2 / ESMFold2)**, using
the **same convention already in place** for `n_h_bonds_interface`, `pi_stack_*`, and
`red_*`. This must be part of the metric-calc infrastructure (every future campaign gets
it), plus a one-time backfill of the current master.

## House convention (verified — MUST match)

From `hbond_consensus_integrate.py` (the script that gave `n_h_bonds_interface` its
backends) and `red_consensus_integrate.py`:

```python
m[f"{COL}_af3"] = m[COL]                       # preserve original AF3 as _af3
m[COL] = m[[af3, b2, esm]].min(axis=1)         # BARE column becomes the consensus
```

So the convention is:
- **bare `M`** = the worst-case consensus (min for higher-is-better, max for lower-is-better)
- **`M_af3`** = the preserved original AF3 value
- **`M_b2`, `M_esm`** = per-backend values
- there is **no `M_consensus` suffix** — the bare column *is* the consensus.

This means existing consumers (app, plots, collections) that read bare `dsasa_total`,
`sc_value`, etc. will, after this change, be reading the **consensus** — exactly as they
already read bare `n_h_bonds_interface` as a consensus today. `_af3` remains available for
anyone who needs the single-structure AF3 value.

## Where metrics are written — integration point (verified)

```
run_analysis.sh
  → collect.py → clean.py                         (per-design rows; 3 refold CIF cols)
  → compute_metrics.py                            (structural / tip / red / hbonds)
  → interface_metrics_parallel.sh                 (contacts_sidecar.py + interface_rosetta.py, AF3)
  → aggregate → analyses/cross_campaign/app_metrics_100k.csv   ← MASTER (source of truth)
  → refresh_data.sh  →  metric_explorer/data/metrics.csv       (SNAPSHOT the app reads)
```

- **Master the explorer actually uses:** `analyses/cross_campaign/app_metrics_100k.csv`
  (169,342,007 B — exact match to the local `metrics_100k.csv`).
- **`refresh_data.sh` bug:** its `SRC` still names `combined_sample10k_metrics.csv`. Fix
  `SRC` to `app_metrics_100k.csv` as part of this work (otherwise a refresh silently
  downgrades the app to a stale 10k file without the new columns).
- Consensus columns are produced **upstream** (in the pipeline that writes the master); the
  explorer inherits them via the corrected `refresh_data.sh`. No local merge.

## Current state (verified)

- All 100k rows have all three refold CIFs (`has_af3_msaoff`=`has_b2`=`has_esm`=100%).
- **Already 3-backend + consensus (do NOT recompute):** `n_h_bonds_interface`,
  `pi_stack_*`, `red_*` (bare=consensus, `_af3/_b2/_esm` present); `iptm`/`dockq` per-backend
  present, consensus derived on load in `app.py`.
- **AF3-only today (this work targets these):** the `interface_metrics_parallel.sh` outputs —
  Rosetta (`interface_rosetta.py`) `dsasa_total/polar/hphobic`, `dsasa_hphobic_frac`,
  `sc_value`, `packstat`, `delta_unsat_hbonds`; Biotite (`contacts_sidecar.py`)
  `n_heavy_contacts_4A/4p5A`, `n_contact_res_prot/rna`, `n_contacts_base/backbone/
  phosphate/sugar`, `frac_base_contacts`, `n_interface_clashes`, `clash_overlap_sum`,
  `worst_clash`, `aa_base_*`, `aa_bb_*`; plus transferable structural
  `n_interface_contacts`, `interface_bsa_approx`.
- **Producers accept `--path-col`** (default `af3_msaoff_cif_path`); they fill fixed column
  names in place, so per-backend runs write to per-backend sidecar CSVs (do NOT overwrite the
  AF3 master in place).

## Architecture (mirror the existing per-backend pattern)

Exactly how `red`/`hbond`/`pi_stack` got their backends:

1. **Per-backend recompute driver** — a variant of `interface_metrics_parallel.sh` that
   takes a backend + path-col and writes a **backend sidecar** (`design_id` + the recomputed
   COLS) instead of overwriting in place. Run it twice: `b2` (`--path-col b2_cif_path`) and
   `esm` (`--path-col esm_cif_path`). Sharded SLURM array with fresh node-local `/tmp` per
   task (the `/local`-scratch-leak fix), `interface_rosetta.py --resume`, merge per backend →
   `interface_be_chunks/…_b2.csv`, `…_esm.csv`. AF3 already exists in the master; not rerun.

2. **`interface_consensus_integrate.py`** (new) — modeled on `hbond_consensus_integrate.py`.
   Joins the b2/esm sidecars onto the master by `design_id`; for each metric in the direction
   map: set `M_af3 = M` (bare), add `M_b2`/`M_esm`, overwrite bare `M` with
   `min`/`max`(af3,b2,esm) (NaN-tolerant). Adds `consensus_n_structures` (0–3). SKIP metrics
   are left untouched (bare stays AF3; no per-backend columns written for them). Backs up the
   master to `app_metrics_100k.bak_pre_iface_consensus` first; idempotent (drops any existing
   `_af3/_b2/_esm` for its cols before recomputing).

3. **Wire into `run_analysis.sh`** (go-forward): after `interface_metrics_parallel.sh`, run
   the b2/esm passes + `interface_consensus_integrate.py` so future campaigns' masters carry
   consensus automatically.

## Direction map

Scoped to exactly what the `interface_metrics_parallel.sh` driver produces
(`interface_rosetta.py` COLS + `contacts_sidecar.py` `CONTACT_COLUMNS`):

**MIN (higher = better → bare becomes worst-case floor):**
`dsasa_total`, `dsasa_polar`, `dsasa_hphobic`, `sc_value`, `packstat`,
`n_heavy_contacts_4A`, `n_heavy_contacts_4p5A`, `n_contact_res_prot`, `n_contact_res_rna`,
`n_contacts_base`, `n_contacts_backbone`, `n_contacts_phosphate`, `n_contacts_sugar`,
`frac_base_contacts`

**MAX (lower = better → bare becomes worst-case ceiling):**
`delta_unsat_hbonds`, `n_interface_clashes`, `clash_overlap_sum`, `worst_clash`

**SKIP (no direction → bare stays AF3, no per-backend columns):**
`dsasa_hphobic_frac` (ratio), `aa_base_*` (20), `aa_bb_*` (20)

**Out of this driver's scope (structural pass, left AF3-only):** `n_interface_contacts`,
`interface_bsa_approx`, `radius_of_gyration_ca`, `helix/sheet/loop_fraction`, `sse_*` —
produced by `compute_all_structural`, not by `interface_metrics_parallel.sh`; they are fold
descriptors (near backend-invariant) and are not worth a second structural pass here.

Already handled elsewhere (untouched here): `n_h_bonds_interface` (min), `red_*` (max),
`pi_stack_*` (min), `iptm`/`dockq` (min, on load).

## Two deliverables

1. **Go-forward infra** — per-backend driver + `interface_consensus_integrate.py` wired into
   `run_analysis.sh`.
2. **One-time backfill** — run the b2/esm passes over the current master's CIF refs +
   `interface_consensus_integrate.py` → `app_metrics_100k.csv` gains consensus now. Keep
   `app_metrics_100k.bak_pre_iface_consensus` for rollback. Explorer picks it up via the
   corrected `refresh_data.sh`.

## Metric explorer integration (`app.py`) — minimal

Because bare columns become the consensus (house convention), the app needs **no new
`_consensus` metric entries** — the metrics it already exposes (`dsasa_total`, `sc_value`,
`frac_base_contacts`, …) automatically become consensus, like `n_h_bonds_interface` already
is. Update only the METRIC **descriptions** to state "worst-of-3 consensus across AF3/B2/ESM
(AF3-only value in `_af3`)". Fix `refresh_data.sh` `SRC`. Viewer/3Dmol serving unchanged.

## Validation

- **Invariants (all rows):** MIN metrics → bare `== min(_af3,_b2,_esm)` and `<= _af3`; MAX
  metrics → bare `== max(...)` and `>= _af3`. Assert programmatically post-backfill.
- **Spot check:** manually recompute the suite on 5 designs' B2 + ESM CIFs; diff vs sidecars
  (exact for counts, float tol for dSASA/sc).
- **Integrate fixture test:** synthetic 3-row master + b2/esm sidecars with known values →
  assert bare consensus, `_af3` preservation, `consensus_n_structures`, NaN handling, and
  idempotency (run twice = same result).
- **Regression:** re-run the natural-vs-design comparison using the `_af3` columns (now the
  single-structure design values) so the comparison stays like-for-like against the naturals.

## Interpretation caveat (carry into the natural comparison)

After this change, bare interface columns are **consensus** (worst-of-3), not AF3. The fair
column for a like-for-like comparison against the naturals (single crystal structure, no
refolds) is **`M_af3`**, not bare `M`. Consensus makes designs look worse than naturals by
construction; the design↔natural comparison must read `_af3`.

## Out of scope (YAGNI)

- No recompute of `n_h_bonds_interface`/`red_*`/`pi_stack_*`/`iptm`/`dockq` (already done).
- No consensus for descriptors (composition, SSE, Rg, hphobic_frac) — no sensible direction;
  and no per-backend columns written for them (avoid ~80 junk `aa_*_b2/_esm` cols).
- No re-refolding; uses existing B2/ESM CIFs.
- No surgery inside `interface_rosetta.py` / `contacts_sidecar.py` — they accept `--path-col`.
- No change to structure serving, collections, or `sync_euler_to_nas.sh`.

## Open items — RESOLVED

- Master file = `analyses/cross_campaign/app_metrics_100k.csv`. ✔
- Interface driver = `interface_metrics_parallel.sh`; producers `interface_rosetta.py`
  (Rosetta) + `contacts_sidecar.py` (Biotite), both `--path-col`. ✔
- Consensus template = `hbond_consensus_integrate.py` / `red_consensus_integrate.py`. ✔
- `refresh_data.sh` `SRC` is stale (`combined_sample10k`) — fix to `app_metrics_100k.csv`. ✔
