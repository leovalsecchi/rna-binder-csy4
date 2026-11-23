# Data availability

Three tiers: what is deposited publicly, what is retained on institutional storage, and
what is not kept at all.

## Deposited

| Deposit | Contents | Approximate size |
|---|---|---|
| Reduced dataset | The gated design set (55,706 designs passing the five-criterion cascade), the 12,000-member ordered library, and the characterised hits with their measured values | tens of MB |
| Container images | The Apptainer images listed in `software_versions.md`, cited by sha256 digest | ~60 GB |
| This code archive | `README.md` describes the contents | ~5 MB |

The reduced dataset is the useful unit: it is what every figure in the chapter is drawn
from, and it is small enough to deposit whole. The deposit carries the metric columns as
computed, including the corrections in `deviations.md` — in particular the `*_true`
columns for BG+BIF, since the stored `sc_RMSD` and `DockQ` values for that stack are
measured against a re-prediction rather than the designed backbone.

## Retained on institutional storage

The full campaign lives on the ETH Zurich Euler cluster under
`/cluster/work/wicky/lvalsecchi/rna_binder/`, mirrored nightly to group storage:

- Design coordinates for all 6,992,515 designs, packed as per-campaign `tar.zst`
  archives (35 archives, ~683 GB).
- Structure re-predictions from all three backends.
- The assembled metric masters, 12–15 GB each: `csy4_master_foot.csv` (218 columns) and
  `combined_clean_v2v9_fullmetrics_bf_final.csv` (178 columns).

These have been handed to the supervisor for retention under HST guidelines 2.5. They
are not deposited publicly: the coordinate set is too large to be useful as a download,
and the masters are wide intermediate tables rather than results.

## Not retained

- **Per-design random state** for RFdiffusion3, BoltzGen and BoltzDesign1 — never
  seeded, never logged. See the reproducibility note in `README.md`.
- **Full L×L PAE matrices** from re-prediction, which is why the ipSAE family exists
  only for designs re-predicted after 2026-07-02.
- **The build recipe for `esmc.sif`**, whose source `docker-archive` tar was deleted.
- **Counts of designs predating the cysteine exclusion and BoltzIF Route B**, derivable
  only from run timestamps rather than recorded directly.

## Archived in this repository

The target inputs are small enough to ship with the code and are in `targets/` (316 KB):
the 4AL5-derived coordinates the stacks were run against, the 16-nt MSA used for
re-prediction, the AlphaFold3 and Protenix template JSONs, and the empty hotspot file.
The filter that defines the gated set is `provenance/gate_v6_final.json`.

## Licence

The deposited data is released under **CC-BY-4.0** (`LICENSE-DATA` in the repository
root). The pipeline code is released separately under the **MIT License** (`LICENSE`).
