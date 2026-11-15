# Foldseek Binder Map — Design

Date: 2026-08-21

## Purpose

Build an exploratory landscape map of the passing binder designs along **two
independent axes**, and compare them:

1. **Fold** — backbone similarity via Foldseek. *What shape is the binder?*
2. **Epitope** — which part of the RNA target it contacts. *Where does it bind?*

Both are shipped in one self-contained HTML explorer as side-by-side linked maps,
coloured by generation stack.

The goal is descriptive — understand what fold space the campaign covers, whether
the five stacks occupy different regions of it, and whether fold predicts epitope.
This is not a selection gate and produces no ordering decision.

Inspiration: `2000targets_Protein_Library_Explorer` (Matteo Lampis). Conventions
borrowed from it are the conservative `min(qTM, tTM)` similarity, a sparse
thresholded neighbour graph rather than a dense one, community detection on that
graph, and a single portable HTML with no CDN or server dependency.

## Inputs

`/Users/leo/Documents/binder_filter/exports/passing.csv` — 71,866 rows.

| Property | Value |
|---|---|
| Unique uuid8 | 6,780 |
| Unique design backbones | 32,848 |
| Binder length | 80 aa (all rows) |
| Stacks | RFD3+LiMPNN 37,951 · BG+LiMPNN 19,162 · BG+BIF 9,881 · BD1+LiMPNN 3,344 · HD+LiMPNN 1,528 |
| Campaigns | v2–v9 |

Each row carries three structure paths in `<archive>.tar.zst!<member>` form,
relative to `/cluster/work/wicky/lvalsecchi/rna_binder/`:
`design_cif`, `af3_cif`, `esm_cif`.

Each CIF holds an 80-residue protein binder and the RNA target.

### Chain identity — do not assume chain A

Verified by probe on 2026-08-21 across all five stacks and all three layers:

| Layer | Protein chain | Notes |
|---|---|---|
| `af3`, `esm` | A (RNA = B) | consistent across all five stacks |
| `design` — BD1+LiMPNN, BG+LiMPNN, HD+LiMPNN | A (RNA = B) | |
| `design` — **RFD3+LiMPNN** | **B** (RNA = **A**) | **inverted** |
| `design` — BG+BIF | A | different mmCIF dialect; needs a real parser |
| `design` — BD1+LiMPNN | A | carries the 18-nt RNA variant |

RFD3+LiMPNN is 37,951 of 71,866 rows. Hardcoding chain A would feed the RNA
target into Foldseek as the binder for **53% of the design layer**, producing a
map that looks plausible and is meaningless.

**Rule: chains are selected by composition, never by chain ID.** Reuse
`rna_binder.analysis.structural.split_chains`, which defines RNA as the chain
whose residues subset `{A, C, G, U}` and protein as the other — already proven
across every dialect here. Structures are read with
`rna_binder.analysis.structural.load_structure` (biotite), not a hand-rolled
parser; a throwaway parser written during design silently failed on the BG+BIF
dialect, which is precisely the failure mode to avoid.

**Assertion:** the selected protein chain must be 80 residues. Violations are
written to `missing/` with reason `chain_select_failed` and excluded.

Design CIFs for the three LiMPNN stacks carry no sidechains. Irrelevant for
Foldseek (backbone geometry) and handled for contacts by counting whatever heavy
atoms are present, with sidechain availability recorded per layer.

## RNA reference frame

`rna_binder.analysis.structural` already defines the canonical target as the
16-nt csy4 RNA `ACUGCCGUAUAGGCAG` with `RNA_LEN = 16`, and computes
`rna_contact_nt_1..16`. That canonical numbering is the epitope frame.

Length variants exist. The stage-0b census over 105 sampled structures
(2026-08-23, real data) found exactly two, both aligning:

| Variant | Length | Offset | Where |
|---|---|---|---|
| `ACUGCCGUAUAGGCAG` | 16 | 0 | every af3 and esm refold; design layer for BG+BIF, BG+LiMPNN, RFD3+LiMPNN |
| `CACUGCCGUAUAGGCAGC` | 18 | 1 | design layer only, and only BD1+LiMPNN and HD+LiMPNN |

So the 18-nt form is a property of two design pipelines' output; refolding
normalises everything to the canonical 16-mer. (Two earlier drafts got this
wrong — first attributing it to BG+BIF from an `entity_poly` header, then to
BD1 alone from a single fixture. The census governs.) Variants are mapped to
canonical indices by exact substring alignment of the RNA sequence; flanking
nucleotides outside the canonical frame are dropped from the fingerprint and
their contact mass is reported separately as `contacts_offframe`.

**Stage 0 emits an RNA-variant census** (sequence → count, per stack and
campaign) before anything else runs. If any variant fails to align to the
canonical 16-mer, the pipeline stops rather than silently misregistering
epitopes — a wrong frame would make every epitope comparison invalid while
looking entirely normal.

## Layers

Three independent layers, toggleable in the explorer:

| Layer | Source | Structures after dedup |
|---|---|---|
| `design` | `design_cif` | ~32,848 |
| `af3` | `af3_cif` | ~71,866 |
| `esm` | `esm_cif` | ~71,866 |

Each layer is a fully independent compute chain with its own DB, graph, clusters
and embedding. They run concurrently on Euler.

Within each layer, **fold** and **epitope** are two independent partitions and two
independent embeddings over the *same* node set, so every point carries both a
fold cluster and an epitope cluster.

## Deduplication

**One point per unique backbone.** Dedup is by coordinate hash, not by filename.

The hash is MD5 over the chain-A backbone atom coordinates (N, CA, C, O) rounded
to 3 decimal places, in file order. Verified on 2026-08-21: sequence variants of
the same design (`design45_seq4/seq7/seq9`) are coordinate-identical, so this
collapses the `_seqN` family automatically, and also catches duplicates across
uuids or stacks that naming would not reveal.

Dedup runs twice: locally within each extraction task (cheap, removes most of the
design-layer redundancy since seq variants share an archive), then globally in
stage 2 to catch cross-task collisions.

Each surviving point carries the full list of `design_id`s behind it. The point
count changes when layers are toggled; the UI displays it.

## Placement

| Side | Path |
|---|---|
| Local | `~/Documents/binder_map/` |
| Euler | `/cluster/work/wicky/lvalsecchi/rna_binder/analyses/foldseek_map/` |
| Results | `analyses/foldseek_map/results/<layer>/` |

Compute runs on Euler because the structures and Foldseek live there. Only small
tables cross the wire. Foldseek is the Euler module `foldseek/8-ef4e960`, loaded
under `stack/2024-06 gcc/12.2.0 openmpi/4.1.6`.

## Pipeline

Eight stages, each a separate script with a file contract, so any stage can be
rerun alone. Stages 1–7 are chained with `--dependency=afterok`.

### Stage 0 — `select.py` (local)

Reads `passing.csv`. Emits per layer `manifests/<layer>.tsv`:
`row_id, design_id, uuid8, stack, campaign, gen, invfold, archive_path, member_path`.

Also emits `meta/designs.parquet` — a curated subset of the CSV columns keyed by
`design_id`, for the explorer. Curated rather than all 65 columns to keep the HTML
payload bounded: `dockq_af3/esm`, `iptm_af3/esm/mean/min`, `rmsd_af3/esm`,
`dsasa_af3/esm`, `hb_af3/esm`, `rog_af3/esm`, `foot_af3/esm`, `pose_disagree`,
`fold_disagree`.

Finally, emits the **RNA-variant census** (`meta/rna_variants.tsv`) described
under *RNA reference frame*, on a sample of structures per stack × campaign. Any
variant that fails to align to the canonical 16-mer stops the pipeline here,
before any compute is spent.

### Stage 1 — `extract.py` (Euler, array 64 × 4 cores, `%16`)

This stage does **two** jobs in one archive pass: backbone extraction for Foldseek
and epitope fingerprinting. The 6,780 archives per layer are opened exactly once;
fingerprints are therefore near-free in wall-clock.

Manifest grouped by `.tar.zst` archive, batched ~106 archives per task. Each
archive is opened once. For each needed member:

1. `load_structure` (biotite), then `split_chains` → protein, RNA **by
   composition**. Assert the protein chain is 80 residues.
2. **Backbone** — keep N/CA/C/O/CB of the protein chain, write a minimal PDB into
   the shard.
3. **Epitope fingerprint** — protein↔RNA heavy-atom contacts at **4.5 Å**, the
   same `RNA_CONTACT_CUTOFF` used by `rna_contact_counts`, so the recomputed
   `rna_contact_nt_*` totals are directly comparable to the existing tables.
   Counts are binned per canonical nucleotide (1..16) × moiety
   (phosphate / sugar / base) → a 48-dim vector, plus scalar
   `contacts_offframe`, `contacts_total`.
4. **Approach geometry** — superpose the complex onto the reference RNA
   (target-fixed, the same RNA-aligned convention adopted for DockQ on
   2026-08-04), then record the binder centroid `(x, y, z)` and its first
   principal axis `(u, v, w)` in that frame.

All work goes through `$TMPDIR` with a bounded cache and streams straight into the
task's `.tar.gz` shard — nothing lands as a loose file. This pipeline has hit the
project inode cap before; the whole point of shards is to keep the footprint in
the tens of inodes rather than 215k.

`struct_name` is the PDB member name inside the shard and equals the source
`design_id`; after dedup, one `struct_name` stands for a set of `design_id`s.

Outputs per layer: `shards/shard_%02d.tar.gz`, `hashes/shard_%02d.tsv`
(`design_id, struct_name, coord_hash, kept`), `fp/shard_%02d.parquet`
(`design_id` + 48 fingerprint columns + `contacts_offframe`, `contacts_total`,
centroid and axis), `missing/shard_%02d.tsv`.

Missing structures (archived, failed, or absent refolds) are recorded with a
reason and excluded. The count is reported and carried into the manifest; the
layer simply has fewer points.

### Stage 2 — `builddb.sh` (Euler, 1 × 8 cores)

Concatenates the hash tables, elects one canonical `struct_name` per `coord_hash`
globally, and writes `dedup/<layer>.tsv` (hash → canonical name) and
`members/<layer>.tsv` (canonical name → design_id list).

`foldseek createdb shards/*.tar.gz db/full` (createdb reads tarballs directly),
then maps canonical names to DB keys via `db/full.lookup` and
`foldseek createsubdb` down to `db/uniq`.

### Stage 3 — `search.sh` (Euler, array 32 × 16 cores, `%16`)

Query keys split into 32 chunks via `createsubdb` → `db/q00..q31`. Task *i* runs:

```
foldseek search db/q$i db/uniq aln/aln_$i tmp_$i \
  -c 0.8 --cov-mode 0 -e 10 --max-seqs 300 -a --threads 16
foldseek convertalis db/q$i db/uniq aln/aln_$i out/hits_$i.tsv \
  --format-output query,target,qtmscore,ttmscore,alntmscore,lddt,fident,alnlen,evalue,bits
```

Each task searches its chunk against the **full** target DB. The 72k × 80aa target
DB is well under 1 GB, so every task holds it in RAM.

Query-chunk arrays are used rather than the module's MPI support: a failed chunk
is a one-line resubmit, whereas a failed MPI rank kills the run.

`--alignment-type` is fixed by the smoke test (see Validation) — default 2
(3Di+AA) versus 1 (true TMalign) on 80-mers.

### Stage 4 — `graph.py` (Euler, 1 × 16 cores, 64 GB)

Run **twice per layer**, once per axis. Stages 4–6 are source-agnostic: they take
a neighbour source and emit graph, clusters and coordinates, so fold and epitope
share one code path and one set of parameters.

**Fold axis** — neighbours from the Foldseek hits:

Similarity is `TM = min(qtmscore, ttmscore)` — conservative, so a pair only scores
high if it is high in both normalisations. Self-hits dropped; the two directions
of a pair reconciled by taking the max of the observed `min` values.

Edges kept where `TM >= TM_FLOOR` (default 0.50, fixed by the smoke test), then
top *k* = 15 neighbours per node, unioned into an undirected weighted graph.

`redundancy` — the count of neighbours at `TM >= 0.95` — is computed from the
**full** hit list before top-*k* truncation, otherwise it saturates at *k*. It is
still bounded above by `--max-seqs` (300); the explorer labels the scale
`redundancy (capped at 300)` rather than implying an exact count.

Two graphs come out of this stage, and they are deliberately different:

- **Clustering graph** — thresholded at `TM_FLOOR`, then top-*k*. A node with no
  hit above the floor becomes a genuine isolate and therefore its own singleton
  cluster. Absence of evidence stays absence of evidence; it is never turned into
  a weak edge.
- **Layout graph** — raw top-*k* with no floor, so every node has *k* neighbours.
  UMAP's `precomputed_knn` requires a full *k* per node, and padding the
  thresholded graph to satisfy it would fabricate exactly the edges the floor was
  meant to reject.

A point can therefore sit inside a visual cluster while being a singleton by
cluster id. That is informative, not a bug, and the explorer surfaces it.

**Epitope axis** — neighbours from the 48-dim fingerprints. Vectors are
L1-normalised (so a promiscuous binder and a sparse one that touch the same place
are neighbours), similarity is cosine, and the kNN is exact via a brute-force
blocked matrix product: 72k × 48 floats is ~14 MB, so the full similarity matrix
is computed in blocks on one node in a couple of minutes. No approximate index,
no extra dependency.

Designs with `contacts_total == 0` — no interface at all — are excluded from the
epitope axis and flagged `no_interface`. They keep their fold coordinates and
render greyed out on the epitope map. Zero contacts is not an epitope, and
imputing one would invent the very signal being measured.

`EPI_FLOOR` (default cosine 0.5, fixed by the smoke test) plays the role `TM_FLOOR`
plays for fold.

Outputs, per axis, under `graph/<axis>/`: `edges_cluster.parquet` (`i, j, sim`),
`knn_layout.parquet` (`i, j, sim`), `nodes.parquet` (`struct_name, degree,
redundancy, is_isolate`).

### Stage 5 — `cluster.py` (Euler, 1 × 16 cores, 64 GB)

Leiden on `graph/<axis>/edges_cluster.parquet` at three resolutions (0.5, 1.0,
2.0), fixed seed → `cluster_r05`, `cluster_r10`, `cluster_r20`, per axis.

Medoid per cluster = the node maximising summed TM to its cluster members over the
sparse graph.

Output: `clusters/<axis>/clusters.parquet`.

### Stage 6 — `embed.py` (Euler, 1 × 16 cores, 64 GB)

UMAP on `graph/<axis>/knn_layout.parquet` (the unthresholded top-*k*), distance
`1 - sim`, `n_neighbors=15`, `min_dist=0.1`, fixed seed. Run per axis.

Output: `embed/<axis>/embed.parquet` (`struct_name, x, y`).

Screen distance is illustrative only. Cluster membership and similarity are
defined by the graph and the exact TM scores, never inferred from coordinates —
the explorer states this.

### Stage 6b — `compare.py` (Euler, 1 × 8 cores)

The point of having two axes. Four comparisons, all cheap table operations:

1. **Fold × epitope contingency**, per layer and per resolution, plus **Adjusted
   Mutual Information** and **Adjusted Rand Index** between the two partitions.
   Answers: does backbone fold predict where a binder lands? A low AMI is a real
   result, not a failure.
2. **Cross-layer epitope agreement.** For each design present in more than one
   layer, cosine similarity between its fingerprints, and whether its epitope
   cluster label matches. Emitted per design as `epi_agree_af3_esm`,
   `epi_agree_design_af3`, `epi_agree_design_esm`.
3. **Intended vs achieved.** Design-layer epitope versus refold epitope — did the
   binder dock where it was drawn to? This is comparison 2 read with the design
   layer as ground truth.
4. **Stack profiles.** Per stack, the mean fingerprint over the 16 canonical
   nucleotides, so "RFD3+LiMPNN prefers nt 3–11" becomes a directly readable
   curve rather than an inference from cluster colours.

Outputs: `compare/contingency_<layer>_<res>.parquet`,
`compare/partition_scores.parquet`, `compare/epitope_agreement.parquet`,
`compare/stack_profiles.parquet`.

**Interpretation caveat, carried into the explorer UI.** Epitope assignment
inherits pose reliability. This campaign's own results show ipTM does not certify
pose and that 81% of fold-agreeing designs dock more than 10 Å apart, so a single
refold's epitope is a noisy estimate. Comparison 2 is the mitigation: designs
whose af3 and esm epitopes agree are the trustworthy ones, and `epi_agree_af3_esm`
is surfaced as a colour scale and a filter rather than buried in a table.
`passing.csv` already gates on DockQ ≥ 0.7 in at least one backend, so this set is
better-posed than the campaign at large — but not uniformly so.

### Stage 7 — `medoids.py` (Euler, 1 × 4 cores)

For the union of medoids across all three resolutions **and both axes**, pull the
CA trace from the shards and bundle every medoid into a single gzipped file
(`medoids/medoids.json.gz`) rather than one file per medoid — again an inode
decision.

### Stage 8 — `build_html.py` (local)

Pulls `embed/<axis>`, `clusters/<axis>`, `graph/<axis>/nodes`, `medoids` and
`compare/*` for **both axes** of every layer, plus `meta/designs.parquet`, and
emits `dist/binder_map.html`.

Coordinates are baked as base64 Float32 arrays, categoricals as integer code
arrays with a lookup, medoid PDBs as gzip+base64 inflated in-browser. Two axes
means two coordinate pairs and two cluster columns per point; the 48-dim
fingerprints are baked as int16 for the per-point sparkline. Estimated size
~25–40 MB across all three layers.

## The explorer

One self-contained HTML. No CDN, no server, no NAS dependency. 3Dmol.js and fflate
are the only inlined libraries.

**Two maps, side by side, linked** — fold on the left, epitope on the right, each
a plain 2D canvas scatter with manual pan/zoom. 72k points redraw in a few
milliseconds, and avoiding a charting library keeps the file small and
dependency-free.

Hover, selection and lasso are **shared**: highlighting a fold cluster lights up
its members in epitope space and vice versa. This is the whole reason both maps
are on screen at once — "does this fold family share an epitope" should be one
glance, not a remembered comparison.

Either map can be collapsed to give the other full width.

**Controls**

- *Layer*: design / af3 / esm. Swaps coordinates, clusters and point count.
- *Colour by* (applies to **both** maps, so either partition can be viewed in the
  other's space): stack (default), campaign, gen, invfold, **fold cluster**,
  **epitope cluster**, **approach geometry** (binder centroid and principal axis
  in the RNA frame, mapped to a colour wheel — this is how a same-epitope,
  opposite-face pair becomes visible), **epitope agreement**
  (`epi_agree_af3_esm`), redundancy, or any curated numeric column.
- *Cluster resolution*: the three Leiden resolutions, applied to both axes.
- *Stack visibility*: per-stack show/hide, so "where does HD+LiMPNN sit relative to
  RFD3+LiMPNN" is a direct question.
- Search by `design_id`; lasso select.

**Right panel**, three tabs

- *Point* — the backbone's design_ids, stack, campaign, metrics, its fold and
  epitope cluster, its per-nucleotide contact profile as a 16-bar sparkline, and
  its cross-layer epitope agreement.
- *Cluster* — id, size, stack-composition bar, and the medoid rendered in 3Dmol.
  For a fold cluster, the mean epitope profile of its members; for an epitope
  cluster, the fold-cluster composition. Each side answers "what does the *other*
  axis look like in here", which is what turns coloured blobs into an atlas.
- *Compare* — the fold × epitope contingency heatmap, AMI/ARI for the current
  layer and resolution, and the per-stack mean nucleotide profiles.

**Provenance footer** — Foldseek version, search parameters, TM floor, epitope
cutoff (4.5 Å) and `EPI_FLOOR`, Leiden resolutions and seed, per-layer and
per-axis counts, and the input CSV hash, so any figure taken from the explorer can
be traced back.

## Reuse of existing contact tables

`rna_contact_nt_1..16` already exists in some campaign tables. Checked
2026-08-21:

| Table | `rna_contact_nt_*` |
|---|---|
| `cross_campaign/combined_clean_v2v8_fullmetrics.csv` (9.7 GB master) | absent |
| `midrun_v9/collected_clean_fullmetrics.csv` | absent |
| `midrun_v9/collected_clean_metrics.csv` (Aug 14) | present |
| `midrun_v3`, `v4`, `v5` `collected_clean_metrics.csv` | present |
| `midrun_v2`, `v6`, `v7`, `v8` | no such file |

**Decision: recompute, do not reuse.** Coverage is v3/v4/v5/v9 only — roughly 42%
of `passing.csv`, missing v7 and v8, the two largest campaigns — and the columns
carry no backend suffix, so they are one unlabelled set rather than per-layer.
That single set cannot support cross-layer epitope agreement, which is the
mitigation for pose noise. Recomputing rides along with an archive pass already
being made, covers all three layers uniformly, and keeps the pipeline
self-contained.

The existing columns are used instead as a **validation cross-check**: on the ~42%
overlap, recomputed per-nucleotide totals must reproduce the stored values (same
4.5 Å heavy-atom cutoff, same `RNA_LEN = 16` frame). Any systematic disagreement
stops the pipeline.

## Validation

**Smoke test first.** The full pipeline runs end-to-end on a 5,000-design
stratified subsample (1,000 per stack) before any full-scale job. It takes ~15
minutes, exercises all eight stages, and settles two parameters that should not be
guessed at scale:

1. `TM_FLOOR` for graph edges, and `EPI_FLOOR` for the epitope graph — chosen
   from the observed similarity distributions so each graph is sparse but not
   fragmented. **If the smoke test is inconclusive, 0.50 stands for both.**
2. `--alignment-type` — whether true TMalign (1) is worth its cost over 3Di+AA
   (2) on 80-residue structures. **If inconclusive, 2 stands**, since it is the
   cheaper option and the one Foldseek is tuned for.

**Checks**

- Known-identical backbones score TM = 1.0.
- ~20 random pairs verified against an independent TMalign run.
- Every `design_id` in `passing.csv` resolves to exactly one point per layer, or
  appears in `missing/` with a reason. Counts must reconcile.
- **Chain selection**: the protein chain is 80 residues for every structure, and
  the RNA chain is 16 nt after canonical-frame mapping. Per-stack tallies of which
  chain ID turned out to be protein are logged, and RFD3+LiMPNN design structures
  must come out inverted (protein = B) — if that tally ever flips, something
  upstream changed and the run stops.
- **Epitope cross-check** against the stored `rna_contact_nt_*` on the ~42%
  overlap, as above.
- **RNA-variant census** from stage 0: every variant aligns to the canonical
  16-mer, or the pipeline stops.
- `results/manifest.json` records tool versions, parameters, input hashes, and all
  stage counts.

## Wall-clock

Roughly 2–3 hours for all three layers, versus about a day serial. Per layer:
~96 array tasks plus 6 singletons (stages 2, 4, 5, 6, 6b, 7); arrays throttled at
`%16` so the ~310 total tasks never flood Euler's top-250 pending evaluation
window. The epitope axis adds negligible wall-clock: fingerprints ride along with
stage 1, and its kNN is a two-minute blocked matrix product.

## Out of scope

- Structures for non-medoid points. Embedding 71,866 × 3 would produce a
  50–150 MB HTML per layer; medoids give the same interpretive value for a few MB.
- A sequence-similarity (MMseqs2) layer. Structure only.
- Cross-layer point linking beyond shared `design_id`.
- A third independent map for approach geometry. It ships as a colour overlay on
  both maps instead — enough to expose same-epitope-opposite-face pairs without a
  third partition to interpret.
- Any selection, ranking, or ordering decision. This is a map, not a filter.
