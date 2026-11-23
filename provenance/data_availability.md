# Data availability

Three tiers: what is deposited publicly, what is retained on institutional storage, and
what is not kept at all.

## Deposited

| Deposit | Contents | Approximate size |
|---|---|---|
| This code archive | The pipeline, the figure and supplementary-table code, and this provenance record; `README.md` describes the contents | ~9 MB |
| Reduced dataset | The gated design set (55,706 designs passing the five-criterion cascade), the 12,000-member ordered library, and the characterised hits with their measured values | tens of MB |

The code archive is at https://github.com/leovalsecchi/rna-binder-csy4, release
`v1.0-thesis-freeze`, archived on Zenodo under **10.5281/zenodo.24485695**
(<https://doi.org/10.5281/zenodo.24485695>); the concept DOI 10.5281/zenodo.24485694 resolves to the latest
version. The deposit is an upload of the tagged tree rather than an automatic
capture, so a later version is uploaded the same way. The deposited snapshot does
not contain its own DOI: a published Zenodo record is immutable, so the DOI is
written into the repository afterwards and appears in the files of any later
version.

The reduced dataset is the useful unit: it is what every figure in the chapter is drawn
from, and it is small enough to deposit whole. The deposit carries the metric columns as
computed, including the corrections in `deviations.md` — in particular the `*_true`
columns for BG+BIF, since the stored `sc_RMSD` and `DockQ` values for that stack are
measured against a re-prediction rather than the designed backbone.

## Not deposited: the container images

The nine production images are **not** deposited. Together they are about 60 GB, which
is past what a single repository record takes, and several carry model weights whose
licences do not permit redistribution — the AlphaFold3 parameters in particular are used
under Google DeepMind's terms and are neither included nor redistributable.

They are identified instead by **sha256 digest** in `software_versions.md`, which also
records the tool version, the commit where one exists, and the weights each image
loaded. `container_recipes/` holds a definition file for each, but be clear about what those
are worth: only `foundry_0.1.9.def` carries real build steps, and the rest are two-line
pointers at sources that in most cases no longer resolve — local Docker daemon tags, a
deleted `docker-archive` tar, a sandbox directory that is gone. **No image is rebuildable
from its recipe**, which is why the digest, not the deffile, is the identifier.

A digest identifies an image; it does not reconstitute one. Anyone needing a byte-exact
environment should ask for the images directly; they are held with the rest of the
campaign on institutional storage.

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
