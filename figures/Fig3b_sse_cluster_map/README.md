# Fig3b_sse_cluster_map

Fold-cluster secondary-structure map: helix content against strand content, one
mark per fold cluster per backbone generator, for the full campaign and for the
gated 55k.

**Adopted into this repo on 2026-11-09** from `~/rna_binder_plots/sse_cluster_map_2026-11-06/`
and from loose scripts in the Euler work root. That directory's own README said
this repo was where the work belonged; it had no version control. The originals
are still in place and can be deleted once this copy is confirmed.

| File | What it is |
|---|---|
| `make_panel.py` | single panel, campaign only |
| `make_panel_compare.py` | two panels, campaign against gated |
| `make_plot.py` | self-contained interactive HTML with hover |
| `euler/sse_worker.py`, `euler/sse.sbatch` | 29-task array, P-SEA on every backbone |
| `euler/agg_sse.py` | per (generator, cluster) aggregation, campaign |
| `euler/agg_pass.py` | same for the gated set, campaign cluster ids |
| `euler/reps_strict.py` | representatives that are themselves passing designs |

Inputs are `inputs/cluster_sse.tsv` (78 clusters) and
`inputs/cluster_sse_pass.tsv` (60 clusters), plus `inputs/cluster_sse_meta.json`.
Outputs land in `output/Fig3b_sse_cluster_map/`.

Interpreter: `~/.venvs/binder_map/bin/python3`. The scripts import this repo's
style module and fall back to inline values when the NAS is unmounted.

## Rebuild
```
python code/Fig3b_sse_cluster_map/make_panel.py
python code/Fig3b_sse_cluster_map/make_panel_compare.py
```
The Euler stage runs first and only needs rerunning if the map changes:
```
sbatch euler/sse.sbatch          # -> sse_bb/shard_NN.tsv, 551,516 rows
python euler/agg_sse.py          # -> cluster_sse.{tsv,json}
python euler/agg_pass.py         # -> cluster_sse_pass.{tsv,json}
```

## Six things that bit during this work
1. **A cluster medoid can misdescribe its own cluster's dominant generator.**
   4 of 11 top clusters are medoid'd by a BG+LiMPNN structure while dominated by
   another generator. BD1 cl0: the medoid reads 2H/3E, helix 45%; its own BD1
   members are 0H/4E, helix 5%. Sample the generator's own members.
2. **A gated backbone is usually not a passing design.** A backbone enters the
   gated set if ANY of its ~10-20 designs passed, and the map elects one
   representative per node independently of that. Only **16.7%** of gated
   representatives are themselves passing designs. Filter on
   `passing_v6_final.csv` design ids, not on backbones.
3. **Cluster ids do not nest across stages.** Each stage clusters its own
   induced subgraph and counts are not even monotonic. Hold the campaign frame
   fixed and recompute over survivors; backbones nest 100%. Never put two
   independent clusterings side by side.
4. **The stored `sse_helix_fraction` / `sse_strand_fraction` columns are the
   wrong data.** They were computed on refold structures, not design backbones,
   and cover only ~42% of the campaign. Recompute on the bbmap shards.
5. **Shard tars are gzipped**, so `tar -xzOf <shard> <one>.pdb` re-decompresses
   the whole archive. Group by shard and use `tar -xzf -T <list>`.
6. **BG+BIF and BG+LiMPNN share backbones** (150,419 of 150,526). The key is
   `<gen>_<uuid8>_<designN|recycleN>`, ignoring the inverse folder. Never count
   both. The representative is the LiMPNN copy.

## Known flaw, not yet fixed
The tinted region per generator in `make_panel_compare.py` covers "the clusters
holding 80% of that generator's backbones". That threshold resolves to a
different number of clusters per panel (BoltzGen: 23 campaign, 9 gated), so the
two regions do not summarise the same thing and some of the apparent shrinkage
is the rule rather than the data. A fixed share threshold, for example every
cluster at or above 2% of its generator, applied to both panels, would fix it.

`VERIFY:` the hard-coded `GEN_TOT_CAMP` / `GEN_TOT_PASS` totals in
`make_panel_compare.py` are inline constants, not read from the data. They were
not re-checked during adoption.

## Related
`fig06b_selection_counterfactual` reuses the same per-backbone P-SEA pass to
show how helix and strand content shift between selection rules.
