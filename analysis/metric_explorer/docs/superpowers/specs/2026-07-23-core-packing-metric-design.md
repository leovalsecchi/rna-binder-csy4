# Core-packing metric — `rnacon` port + geometric packing quality

**Date:** 2026-07-23
**Status:** ✅ IMPLEMENTED 2026-07-23 — backfilled on 100k (0 errors, 98.5% scored), live in the
explorer (21 metrics). By-stack signal: BD1 worst cores (hyd/pack 0.50/0.50), BG best (0.80/0.70),
RFD3 clean chemistry but looser pack (0.82/0.51). `core_packing_score` CN bounds calibrated on real
cores (p5=8/p50=12/p95=15 → LO=8, HI=15). Refinable from the full-100k distribution.
**Scope chosen:** faithful `Bbindercore` port **+ geometric packing-density** (no Rosetta).
**Spans two repos:** the `rna_binder` analysis pipeline (metric computation + 100k backfill,
authoritative code at Euler `$RB/code/src/...`; working copy `~/scratch_scoring/src/...`) and
`~/Documents/metric_explorer` (drivable metric + 3D viz + data refresh).

## Goal

Turn the ChimeraX `rnacon`/`Bbindercore.cxc` core view into computed per-design metrics on the
100k dataset + metric explorer: (1) a **faithful** reproduction of "how hydrophobic/yellow is
the buried core", and (2) an **improved, continuous core-packing-quality** score based on
packing density — a new *fold-quality* axis (the explorer is otherwise almost all interface
metrics).

## What the preset does (ported exactly)

On the folded binder **in isolation**: per-residue SASA → **core = residues with SASA < 30 Å²**
→ split **hydrophobic {ALA VAL LEU ILE MET PHE TRP TYR PRO CYS} (yellow/good)** vs
**polar/charged (red/bad)**; **good** = C–C packing contacts ≤ 3.8 Å among hydrophobic core
sidechains; **bad** = buried polar N/O/S contacts ≤ 3.5 Å + steric clashes. Tunables in the
preset: core threshold 30 Å² (15–25 stricter), pack 3.8 Å, polar 3.5 Å.

Note vs `red_contacts.py`: that metric's hydrophobic set **excludes Tyr**; `Bbindercore`
**includes Tyr** in the core set — the port follows `Bbindercore` (Tyr in).

## Module: `analysis/core_packing.py` (biotite, no Rosetta)

Sibling of `red_contacts.py`, same contract:
- `from .structural import load_structure, split_chains` → protein/RNA auto-detected (no
  hardcoded chain); **drop the RNA, keep the protein binder chain** (monomer core, matching the
  preset's model-tiling isolation).
- SASA via `biotite.structure.sasa(protein, point_number=200)` summed residue-wise
  (`apply_residue_wise`) — the same tool `structural.py:interface_bsa_approx` already uses.
- Geometry via `scipy.spatial.cKDTree` (as `red_contacts`/`tip_contact` do).
- `compute_core_packing(path, sasa_core_cutoff=30.0, pack_dist=3.8, polar_dist=3.5,
  cn_radius=10.0, cn_ref=10.0, underpack_cn=6.0) -> dict`; **never raises** (failure →
  `{"core_packing_error": "<Type>: <msg>"}`); flat scalar dict; `CORE_COLUMNS` tuple; `main()`
  argparse→JSON CLI (mirror `red_contacts.main()`).

### Columns — faithful port (reproduces the ChimeraX picture)

| column | definition |
|---|---|
| `core_n_res` | # protein residues with monomer SASA < 30 Å² |
| `core_frac` | `core_n_res` / total protein residues |
| **`core_hydrophobic_frac`** | **THE faithful metric** — hydrophobic core ÷ all core (the "yellow fraction"); higher = better |
| `core_n_buried_polar` | polar/charged core residues (the "red" count) |
| `core_buried_polar_frac` | `core_n_buried_polar` / `core_n_res` (lower = better) |
| `core_good_pack_contacts` | C–C atom pairs ≤ 3.8 Å between hydrophobic-core sidechains (yellow) |
| `core_bad_polar_contacts` | buried polar N/O/S atom pairs ≤ 3.5 Å involving polar-core residues (red) |
| `core_clashes` | intra-core sidechain heavy-atom VdW overlaps ≥ 0.6 Å (Bondi radii, as `red_contacts`) |

### Columns — improved: geometric packing density (the "quality" metric)

Packing quality = how densely the buried core coordinates, independent of composition.
For each **core residue** (SASA < 30): its Cβ coordination `cn_i` = # of *other* protein Cβ
atoms within `cn_radius` (default 10 Å); its heavy-atom contact number = # protein heavy atoms
within 5 Å of the residue's heavy atoms.

| column | definition |
|---|---|
| `core_mean_cb_cn` | mean Cβ coordination number over core residues (higher = tighter) |
| `core_mean_heavy_cn` | mean heavy-atom contact number over core residues |
| `core_underpacked_frac` | fraction of core residues with `cn_i < underpack_cn` (void proxy; lower = better) |
| **`core_packing_score`** | **the improved quality scalar** — mean over core residues of `min(cn_i / cn_ref, 1)` ∈ [0,1]; higher = better packed |
| `core_packing_error` | error string when a design can't be scored |

`cn_ref`, `cn_radius`, `underpack_cn` are **tunables to calibrate** once on real data (and vs
the [[rbp-rna-reference-set]] naturals) so `core_packing_score` sits in a sensible 0–1 band;
defaults above are first-pass. `core_packing_score` correlates only weakly with
`core_hydrophobic_frac` (composition) — that's the point: it catches right-chemistry-but-loose
cores.

## Pipeline integration (`compute_metrics.py`)

Per the mapped contract (all reproduced from `red_contacts`):
1. `from rna_binder.analysis.core_packing import compute_core_packing, CORE_COLUMNS` (beside the
   others, ~L70). No `analysis/__init__.py` change.
2. Fork-inherited `_CORE_CFG = {"enabled": True}` near `_RED_CFG`; CLI `--no-core`.
3. In `_compute_one`, after the red block: `if _CORE_CFG["enabled"]: try: m.update(compute_core_packing(Path(cif)))
   except Exception as e: m["core_packing_error"] = ...` — on the **af3_msaoff refold cif**.
4. Columns merge automatically (every returned key becomes a column; all scalars, no filter
   needed).

## 100k backfill (Euler, CPU — clone the proven `red` path)

- Bare spack/host `python3` (zstd-capable `tar`; **not** a container), sharded `sbatch` array
  0–19 (~5k/chunk), `--account=es_wicky`, `--time=01:00:00` (→ `normal.4h`), like
  `red_100k_array.sh`. One extra SASA pass on the ~1k-atom protein monomer → **~1 h, no GPU**.
  Cheapest: run it as a `core`-only fast path in `_compute_one` (like `--red-only`), or simply
  extend the existing red backfill to also emit `core_*`.
- Merge with a `core_integrate.py` cloned from `red_integrate.py`: concat chunks →
  `drop_duplicates("design_id")` → LEFT-merge `core_*` onto
  `$RB/analyses/cross_campaign/app_metrics_100k.csv` (idempotent; drops any pre-existing
  `core_*`).

## Metric explorer wiring

**The deliverable is the numeric metric columns** (step 1); the 3D coloring (step 2) is
OPTIONAL eye-candy that mimics the ChimeraX look and can be dropped without affecting the
metrics.

1. `app.py METRICS_SPEC`: add rows — `core_hydrophobic_frac` ("Core hydrophobic frac", higher
   better) and **`core_packing_score`** ("Core packing quality", higher better), each with an
   ⓘ description; optionally `core_underpacked_frac` (lower better). `load_one` picks them up
   automatically once the column exists. **This is what makes them drivable metrics.**
2. **3D viz (OPTIONAL)** — new `VIZ_MODE` key `"core"` + `VIZ_DESC["core"]="buried core sidechains —
   hydrophobic (yellow) vs polar (red)"`; new `applyViz` `else if(mode==="core")` branch that
   shows core sidechains as sticks coloured yellow(hydrophobic)/red(polar), mirroring ChimeraX.
   **Caveat:** the frontend has no SASA, only coordinates → burial is approximated **client-side
   by neighbour-atom density** (a core residue = one whose sidechain is surrounded); the *metric
   columns* use real SASA on Euler. Map both new metric keys to `"core"`.
3. **Data:** `refresh_data.sh` currently writes `data/metrics.csv`, which the app does **not**
   load — must also refresh `data/metrics_100k.csv` (and optionally `metrics_30k.csv`) from the
   enriched `app_metrics_100k.csv`. Restart the app (kill via `lsof -tiTCP:8000`, not
   `pkill -f app.py`).

## Testing

- **Unit** (`test_core_packing.py`): a tiny synthetic structure with a known buried
  hydrophobic residue vs an exposed polar one → assert core selection, `core_hydrophobic_frac`,
  contact/clash counts, and monotonicity of `core_packing_score` (a tightly-packed cluster
  scores higher than a spread-out one); error path returns `core_packing_error`, never raises.
- **Live**: run the CLI on a handful of real AF3 refold CIFs; sanity-check `core_*` ranges and
  eyeball one against ChimeraX `Bbindercore` (same core residues coloured).
- **Calibration**: histogram `core_packing_score` across stacks + the naturals subset to set
  `cn_ref` so the scale is meaningful before wiring the explorer slider.

## Out of scope (per user's choice)

- Rosetta `packstat`/RosettaHoles2, Voronoi void-volume, buried-unsat-polar — proposed but not
  selected. `core_packing.py` is structured so a `packstat` column could be added later via a
  separate `pyrosetta.sif` array without touching the geometric path.

## Decisions (locked 2026-07-23)

1. Core definition: **absolute SASA < 30 Å²** (faithful to the preset, shared by both metrics).
2. Backfill: **standalone `--core-only` array**; the `core_integrate.py` LEFT-merge is
   idempotent so `core_*` can be added to `app_metrics_100k.csv` at any time without re-running
   the rest.
3. **100k only** (not 30k).
4. **Include** the optional yellow/red 3D core viz in the explorer.
