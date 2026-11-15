# Campaign-Wide Fold Map — Design

Date: 2026-09-14

## Amendments — 2026-09-24 (first real launch)

The design below is the 2026-09-14 record and is left unedited. These five points
supersede it where they conflict.

1. **Stage 3 for map 1 is 512 chunks × 256 concurrent**, not 128 × 64. Chunking is a
   pure partition of the query set (`stage3_split_queries.py` strides the keys, stage 4
   globs `hits_*.tsv`), so it changes no result — only wall-clock and backfill
   placement. Verified safe at this width: `q{i:02d}` and `printf "%02d"` widen
   identically past 99, and Euler's `MaxArraySize` is 15000. The "128" in the stage
   table and the cost table below is the historical figure.

2. **A gate job now sits between `full_s2` and `full_s3`** (`slurm/gate_full.sbatch`,
   `scripts/gate_checks.py`, `scripts/probe_project.py`). It runs the hash-identity
   check against `results_v6`, the manifest-id ↔ shard-stem check, the split-non-empty
   check, then both probe arms, and fails if the projected cost exceeds 6000 CPU-h.
   `full_s3` depends `afterok` on it. Validation steps 4, 4b and the scaling probe are
   therefore no longer a manual race to `scancel` — they are a scheduler dependency
   that fails closed.

3. **`map3_n` is 9,760, not 9,759.** 9,760 is the stage-0 name-rule row count; 9,759 is
   the post-stage-2 coord-hash node count. They are different stages and comparing them
   is an error. The task-8 brief's Step 3 expects 9,759 from `sample.json` and is wrong
   — the live `build_shared/meta/sample.json` correctly holds 9,760.

4. **`MIN_SHARED_HASHES` is 3000, not 5000.** Measured 2026-09-24: the shared manifest
   holds 508,207 design_ids, `results_v6` holds 55,706, and the intersection is **6,377**
   — a hard upper bound on the gate's `n_shared`. A 5000 floor left only 27% of margin,
   so a modest extraction loss would have blocked a good run after `full_s2` had burned.

5. **First launch: 2026-09-24**, jobs 15057505 (extract) → 15057604 (compare), with the
   gate at 15057555 and the 512-task search array at 15057569. Nothing in this project
   had ever been run before that date.

## Purpose

Produce **three fold clusterings under identical conditions, differing only in
sample size**, so the fold space of the campaign can be compared against the fold
space of what was filtered and of what was actually ordered.

| # | map | designs | nodes | status |
|---|---|---|---|---|
| 1 | uniform random sample of the whole campaign | ~1.18M | **500,000** | to build |
| 2 | the filtered pool | 55,706 | **31,763** | **exists** — `results_sw_k30_f50_mut` |
| 3 | the ordered library | 12,000 | **9,759** | to build |
| 4 | fold-representative census (deferred, see below) | — | **702,045** | decision after the probe |

Map 2 is not rebuilt. Maps 1 and 3 are built with the same pipeline, the same
parameters, and the same dedup key that produced map 2.

This is descriptive. It produces no selection, no gate, and no ordering decision,
and it does not revisit the library.

## The reference clustering being reproduced

The operative map is `results_sw_k30_f50_mut`, and the resolution the 12k
selection consumed is **`cluster_r80`** — verified at
`analyses/library_selection/v6hier_q.py:61`, which reads `cluster_r80` for both
the fold and the epitope axis. 917 fold clusters on that pool.

| stage | parameters |
|---|---|
| Foldseek | v8-ef4e960, `module load stack/2024-06 gcc/12.2.0 openmpi/4.1.6 foldseek/8-ef4e960` |
| search | `--alignment-type 1` (TMalign), `-c 0.8 --cov-mode 0 -e 10 --max-seqs 300 -a` |
| similarity | `TM = min(qTM, tTM)`, self-hits dropped, both directions reconciled by max-of-min |
| fold graph | floor **0.50**, k = **30**, **mutual** kNN |
| epitope | 48-dim contact fingerprint at 4.5 Å, L1-normalised, cosine, floor **0.60**, same k and mutual flag |
| clusters | Leiden at r05/r10/r20/r40/r80, fixed seed; **r80 is operative** |
| dedup key | MD5 over the protein chain's N/CA/C/O coordinates, rounded to 3 dp |

Every one of these is carried over unchanged, **including the dedup key and its
known defects** (below). Matching map 2 is the requirement; improving on it would
break the comparison the three maps exist to support.

## The dedup key is not a fold identity — measured, and it matters

The coord-hash key is sensitive to rigid-body placement and to sub-0.1 Å noise, so
"distinct hashes" is **not** "distinct backbones". Measured 2026-09-14 on real
pairs drawn from `results_v6`:

- **BG+BIF** `design8_seq2` vs `_seq4`: **superposed RMSD 0.94 Å**, but **mean raw
  per-atom offset 26.6 Å** (max 47 Å). The same fold, placed elsewhere in the
  coordinate frame. BoltzIF re-predicts the complex per sequence; it does not
  rebuild the backbone.
- **HD** `recycle7_sample1` vs `_sample8`: **max 0.022 Å over all 395 backbone
  atoms**, identical atom order and chain id. The same structure; the LigandMPNN
  packer rewrites coordinates with sub-0.1 Å jitter that 3-dp rounding does not
  absorb. `sample0..9` are the ten sequences on one hallucinated backbone.

So the campaign has two different counts, and they differ by 4.2×:

| stack | designs | coord-hash nodes | distinct folds | inflation |
|---|---|---|---|---|
| BD1+LiMPNN | 1,444,614 | 144,468 | 144,468 | 1.0× |
| BG+LiMPNN | 1,519,958 | 151,999 | 151,999 | 1.0× |
| RFD3+LiMPNN | 1,524,796 | 152,512 | 152,512 | 1.0× |
| **BG+BIF** | 1,480,140 | **1,480,140** | 150,526 | **9.8×** |
| **HD+LiMPNN** | 1,023,007 | **1,023,007** | 102,540 | **10.0×** |
| **total** | **6,992,515** | **2,952,126** | **702,045** | **4.2×** |

**Decision: sample over coord-hash nodes, not folds.** Map 2 was built on this key
and carries the same inflation, so matching it is what keeps the three maps
comparable. The cost is stated rather than fixed: BG+BIF and HD folds carry ~10×
the weight per fold in map 1, and any statement about "fold space coverage" reads
an inflated denominator. Foldseek sees these near-duplicates at TM ≈ 0.95–1.00, so
they raise the `redundancy` column and merge into their own cluster rather than
inventing new ones.

## Map 4 — the fold-level view (DEFERRED, decision point named)

Maps 1-3 sample and cluster **coord-hash nodes**, which over-represent BG+BIF and
HD folds ~10x each (§ above). Map 4 would repeat map 1 on **fold
representatives** — one structure per distinct fold, 702,045 of them, a census
rather than a sample — with the same frozen parameters.

**Why it cannot be replaced by re-weighting.** Node-level cluster sizes can be
converted to fold counts after the fact, and stage 9 does exactly that (below).
But the inflation changes **the clustering itself**, not only the tallies: Leiden
sees ~10 near-identical copies of each BG+BIF and HD backbone as a dense blob, so
those stacks form tighter and more numerous clusters than their fold content
warrants. No post-hoc weight repairs a community structure that was built on
duplicated nodes. Map 4 is therefore a genuinely different experiment, not a
different summary of the same one.

**Cost.** Search scales ~N², so 702,045 vs 500,000 nodes is ~2x the expensive
stage, plus a second extraction pass. Roughly double the project.

**Decision point.** Do not build map 4 up front. Build maps 1-3, read stage 9's
node-vs-fold columns, and decide with two things in hand that do not exist yet:
the scaling probe's measured cost exponent (replacing the N² guess), and evidence
of whether the inflation actually moved any conclusion. If the node-level and
fold-level readings agree, map 4 buys nothing.

## Node sets

Backbones are enumerated from `design_id` with a per-stack rule, verified against
the map-2 ground truth (`results_v6/design/dedup/members.tsv`, 55,706 → 31,763):

| stack | rule | evidence |
|---|---|---|
| RFD3+LiMPNN | strip `_seq<N>` | 14,455 prefixes, **0** split |
| BG+LiMPNN | strip `_seq<N>` | 5,694 prefixes, **0** split (4 backbones span 2 prefixes) |
| BD1+LiMPNN | strip `_seq<N>` | 1,774 prefixes, **0** split |
| BG+BIF | one per design | 6,082 prefixes, 1,629 (27%) split |
| HD+LiMPNN | one per design | 802 designs → 802 hashes across 447 `(uuid, recycle)` groups |

The rule only decides **which member to extract**; the coord-hash remains
authoritative. Applying it to the map-2 pool yields 31,767 keys against 31,763 true
nodes — an over-count of exactly the 4 BG+LiMPNN backbones spanning two prefixes.

**Map 1** — uniform random sample of **500,000** of the 2,952,126 nodes,
`SEED = 42`, without replacement. Uniform, not stack-stratified: cluster size as
generated is what the map is read against, and stratifying would rescale exactly
that.

**Map 3** — every backbone behind
`analyses/library_selection/csy4_binder_library_12k_v6_k30.csv`: 12,000 designs.
**Stage 0 yields 9,760 name-rule keys; stage 2's coord-hash dedup collapses them
to 9,759 graph nodes.** The two differ by exactly one because
`BG_LiMPNN_..._04bdeaa9_design7` and `BG_LiMPNN_..._fadec0de_design67` are
coordinate-identical backbones generated under different uuids (verified) — one
of the four cross-uuid duplicates measured on the full v6 pool. Not a sample;
the whole set.

The two are built as **independent runs**. Map 3 is not injected into map 1 — that
would give its nodes inclusion probability 1 against 0.169 for the rest and make
map 1's cluster sizes unreadable.

## Labels — joined from existing exports, never recomputed

The filtering is already done and audited. Nothing here re-derives a gate. Each
node carries `n_designs` plus two membership flags obtained by a straight join on
`design_id` against files that already exist:

| flag | source | designs | nodes |
|---|---|---|---|
| `in_passing` | `foldseek_map/passing_v6_final.csv` | 55,706 | 31,763 |
| `in_12k` | `library_selection/csy4_binder_library_12k_v6_k30.csv` | 12,000 | 9,760 name-rule / 9,759 coord-hash |

A node is flagged if **any** design on it is in that file.

## The repack constraint

On **2026-09-11** every midrun directory was collapsed into one solid
`midrun_vN.tar.zst`. Design structures are now double-nested:

```
midrun_v8.tar.zst  →  midrun_v8/designs/<uuid>.tar.zst  →  <design_id>.cif
```

- **683 GB** across **35** outer archives (8–36 GB each), of which only the first
  8-10% of each is design output — ~68 GB is actually read
- measured throughput ~163 MB/s per stream; parse+hash+fingerprint is
  **0.031 cpu-s/structure**, so parsing all 508,100 costs ~4.4 CPU-h
- `.list` sidecars index members, but a solid zstd stream has **no random access**
- designs are a slice of each archive — 1,653 of 20,442 entries in `midrun_v8`
- **the store's path templates are stale**: 40 of 40 sampled design archives absent
- `code/scripts/archive_read.py` resolves only single-level `<dir>.tar.zst!<member>`

The access pattern inverts: **one pass per outer archive, pulling everything
needed from it**, stopping once the wanted members are in hand — design members
sit in the first 8-10% of an archive (measured), so the pass never touches the
refold bulk behind them. Map 1 and map 3 share this pass — map 3's 9,759 members are extracted alongside map 1's
500,000 in the same sweep, then split into two builds. Streaming 683 GB twice for
the sake of 9,759 structures would be indefensible.

## Pipeline

Stages 2–5 are unchanged and run **once per map**. Runs on Euler CPU partitions.
Python is `/cluster/work/wicky/lvalsecchi/binder_map_venv/bin/python`.

### Stage 0 — `stage0_sample.py` (new, local)

Enumerates the node universe, draws map 1, collects map 3, and emits **two**
manifests plus the shared extraction plan:

- `build_full/manifests/design.tsv` — map 1, 500,000 rows
- `build_lib12k/manifests/design.tsv` — map 3, 9,759 rows
- `build_shared/manifests/design.tsv` — the union (≈508,100 rows), which is what
  stage 1 actually extracts
- `meta/nodes.parquet` per build — `design_id, stack, campaign, n_designs,
  in_passing, in_12k, map1, map3`
- `meta/sample.json` — seed, universe, draw size, per-stack counts

Manifest columns are the stage-1 contract plus the nested triple:
`row_id, design_id, uuid8, stack, campaign, gen, invfold, outer_path, inner_path,
member_path`.

### Stage 1 — `stage1_extract_nested.py` (rewritten, 35 tasks × 2 cpu)

**One array task per outer archive**, over the union manifest. Streams the archive
once and **stops as soon as every wanted inner archive is consumed** — measured
2026-09-14, design members occupy only the first 8-10% of an archive and the rest
is refold output, so early exit turns 683 GB of reads into ~68 GB. Streams the
archive, expands needed inner blobs **in memory**, and per member: `load_structure`
→ `split_chains` **by composition, never by chain ID** (RFD3+LiMPNN design CIFs
have chain A = RNA) → assert 80 residues → coord-hash → backbone PDB into the
task's `.tar.gz` shard → 48-dim epitope fingerprint at 4.5 Å + approach geometry.

Scratch discipline is a correctness requirement: nothing expands to a loose file,
`--tmp` is requested explicitly, and each task asserts parsed **content**, not file
existence. A blank shard with the right row count is the failure mode that defeated
the resume check in July 2026.

### Stage 1b — `split_shards.py` (new, 1 × 4 cpu)

Splits the shared shards into `results_full/` and `results_lib12k/` by manifest
membership, so stages 2–5 see two ordinary result trees and need no changes.

### Stages 2–5 — unchanged, run twice

`stage2_builddb.sh` (coord-hash dedup + `createdb`/`createsubdb`),
`stage3_search.sbatch` (`ALN_TYPE=1`, 128 chunks for map 1, 4 for map 3),
`stage4_graph.py`, `stage5_cluster.py`.

### Stage 4 — two scale fixes (map 1 only, but applied once)

- `stage4_graph.py:77` `pd.concat` of every hit TSV holds ~150M rows with string
  columns (~20 GB at map-1 scale) and exhausts a 64 GB node → stream per file,
  map names to node ids on read, keep `(int32, int32, float32)` ≈ 1.8 GB.
- `src/binder_map/knn.py:167` `block = 4096` allocates 4096 × 508k float32 =
  **8.3 GB per block** in the epitope kNN → drop to 1024.

### Stage 9 — `compare_maps.py` (new, 1 × 8 cpu)

Reads all three completed maps and emits `compare/`:

- `map_summary.parquet` — per map: nodes, clusters at each resolution, singletons,
  cluster-size quantiles, stack composition, median redundancy
- `cluster_table_<map>.parquet` — per cluster: **node count AND distinct-fold
  count**, stack composition, `n_in_passing`, `n_in_12k`. Both size columns are
  mandatory and must appear side by side: the node count is inflated ~10x for
  BG+BIF and HD, so a cluster size or a "fold space coverage" figure quoted from
  it alone reads an inflated denominator. The fold key is `_seq<N>` stripped for
  BG+BIF and the three LiMPNN stacks, `_sample<S>_pack<P>` stripped for HD.
- `crosswalk.parquet` — for designs present in more than one map, the cluster id
  in each, with ARI and AMI between the partitions on the shared designs

## Validation

1. **Hash identity.** Every design present in both a new run and `results_v6` must
   reproduce its `coord_hash` exactly. A hard assert — it proves the repack rewrite
   did not change what gets extracted. Gate the search on it.
2. **TM identity.** Spot-check `qtmscore` on query/target pairs shared with the
   cached v6 hits; max |ΔTM| < 1e-3.
3. **Chain selection.** Protein chain is 80 residues everywhere; the per-stack tally
   must show RFD3+LiMPNN inverted (protein = B). If that flips, stop.
4. **Counts reconcile.** Every manifest row resolves to exactly one node or appears
   in `missing/` with a reason.
5. **Map 3 is complete.** It is a census, not a sample. Stage 0 must emit
   **9,760** rows and stage 2 must dedup them to **9,759** graph nodes; any
   shortfall must be fully accounted for in `missing/`. Do not compare a
   name-rule count against a coord-hash count — they differ by construction.
6. **Sample is uniform.** Map 1's per-stack draw counts match the multinomial
   expectation from the universe table within tolerance.
7. **Scaling probe before the big run.** Stage 3 runs on 50k and 150k node subsets
   first; fit the cost exponent and re-derive the estimate. The figure below is an
   N² extrapolation from the v6 run (31,763 nodes → 8.5 CPU-h), not a measurement.

## Cost

| stage | shape | estimate |
|---|---|---|
| 0 sample | local | minutes |
| 1 extract (shared) | 35 × 2 cpu, ~10 min | **~5 CPU-h**, reads ~68 GB |
| 1b split | 1 × 4 cpu | minutes |
| 2–5 map 3 (9,759 nodes) | as the v6 run | **~20 CPU-h** |
| 3 search map 1 | 128 × 16 cpu | **~2,100 CPU-h** (probe replaces this) |
| 4 graph map 1 | 1 × 16 cpu, 128 GB | ~2 h |
| 5 cluster map 1 | 1 × 16 cpu, 64 GB | ~1 h |
| 9 compare | 1 × 8 cpu | minutes |

Hit table ~20 GB (map 1), shards ~2.5 GB.

## Known limits, stated up front

- **Map 1 is a 16.9% sample of nodes.** Clusters below ~30 true nodes are not
  reliably resolved, and singletons are sampled at 16.9% like everything else.
- **Cluster counts are not comparable across maps without care.** Leiden cluster
  count grows with node count; 500k nodes will yield far more clusters than 31,763.
  Compare size *distributions* and composition, not raw cluster counts.
- **The three partitions are independent.** A cluster id in map 1 has no relation
  to a cluster id in map 2. The crosswalk reports ARI/AMI on shared designs as a
  result, not as a correspondence.
- **BG+BIF and HD are ~10× over-represented at fold level** in every map, map 2
  included. Consistent, therefore comparable — but not a neutral view of fold space.
- Epitope labels on map 1 are noisier than on maps 2–3: map 1 is the whole campaign,
  not a DockQ-gated slice. The fold axis is the deliverable; epitope rides along
  because it is nearly free at stage 1.

## Out of scope

- Rebuilding map 2. It exists and is the reference.
- The af3 and esm layers. Design layer only.
- UMAP embedding and the HTML explorer.
- Any change to the 12k library.
- **Repairing the campaign-wide stale path templates.** Stage 0 rewrites refs for
  its own manifests; the store, `archive_read.py`, binder_filter's ChimeraX
  click-through and `opendesign` remain broken by the repack. Tracked separately.
