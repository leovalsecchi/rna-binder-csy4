# Chosen parameters

Measured on the esm smoke shard (614 unique backbones, 52,061 Foldseek hits), 2026-08-23.

## TM_FLOOR = 0.40  (was 0.50)

Fold-graph edge threshold, on `TM = min(qTM, tTM)`.

| floor | edges | nodes connected |
|---|---|---|
| 0.30 | 30,599 | 600/614 (97.7%) |
| **0.40** | **14,994** | **563/614 (91.7%)** |
| 0.50 | 8,706 | 491/614 (80.0%) |
| 0.60 | 5,104 | 391/614 (63.7%) |
| 0.70 | 2,672 | 291/614 (47.4%) |

0.50 was the spec default and is wrong for this data: it strands 20% of designs as
singletons, so a fifth of the map carries no neighbourhood information. 0.30 keeps
almost everything but at twice the edge count, risking one undifferentiated blob.
0.40 is the knee — 92% connected, still sparse.

Also defensible on the science: TM >= 0.5 is the conventional same-fold threshold,
but TM-score is length-sensitive and these are 80-residue binders, which score lower
than full-length domains at equivalent structural similarity.

Observed TM distribution (excluding self-hits, n=51,447):
q50 0.324 | q75 0.425 | q90 0.599 | q99 0.847 | max 0.997

## ALN_TYPE = 2 (3Di+AA, the Foldseek default)

Kept. All 614 self-hits score exactly 1.000, so the alignment is behaving correctly
on 80-mers; no evidence TMalign (type 1) would buy anything for its cost.

## EPI_FLOOR = 0.50

Unchanged pending the epitope distribution from stage 4.

## Confirmed reachable: zero-edge nodes

Minimum hits per query is 1 — at least one structure matched only itself. After
`symmetrise_min` drops self-hits it has zero edges, which is exactly the case that
produced NaN coordinates and crashed UMAP's spectral init. Not hypothetical.
