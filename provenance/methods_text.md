# Drafted Methods text

Two passages for the Notion chapter. Past tense, passive, per HST guidelines 4.5.
Placeholders in `«»` need a value before submission.

---

## §2.3.2 Orchestration — replacement for the version-control sentence

The current draft ends with:

> Source code was version controlled; the exact commit used for the work reported here
> is given in section 2.13.5.

This is not accurate. The pipeline was developed directly on the cluster filesystem and
was not under version control while the campaigns ran; there is no commit series and no
commit to quote. Replace that sentence with:

> The pipeline was developed directly on the cluster filesystem and was not placed under
> version control during the campaigns. The production tree was archived verbatim at the
> computational data freeze and is cited in section 2.13.5; configuration files quoted
> throughout this chapter are the archived production files rather than reconstructions.

That last clause is worth keeping, because it is the fact that makes the rest of the
chapter citable: campaign submission overrode only scheduler fields — partition,
nodelist, GRES and walltime — together with the per-job HalluDesign seed, so every
`configs/*/default.yaml` is the parameter set that ran.

---

## §2.13.5 Code and data availability — full replacement

> Pipeline source code, configuration files and the modifications applied to third-party
> design tools were archived at the computational data freeze on «freeze date» and are
> available at «repository URL», release «tag», archived under «Zenodo DOI». The archive
> is an as-run snapshot of the production tree rather than a development history: the
> pipeline was not under version control while the campaigns ran, and the snapshot
> records the state of the code at the freeze. Container images, which range from 3.5 to
> 15.7 GB, are deposited separately under «Zenodo DOI» and are cited by digest in
> Supplementary Table S1. Mid-campaign changes to generation and metric code that affect
> the interpretation of stored values are documented in the archive's deviation record
> and summarised in «Supplementary Table S2».
>
> Design coordinates, metric tables and selection outputs are held on the ETH Zurich
> Euler cluster under `/cluster/work/wicky/lvalsecchi/rna_binder/`, mirrored nightly to
> the group storage system, and have been handed to the supervisor for retention in
> accordance with HST guidelines 2.5. The assembled metric tables are 12 to 15 GB and are
> not deposited in full; the deposited reduction comprises the gated design set, the
> 12,000-member ordered library and the characterised hits, at «Zenodo DOI».
>
> Individual designs cannot be regenerated from this code. Of the five generative stacks,
> RFdiffusion3 ran with its seed unset, BoltzGen and BoltzDesign1 expose no seed
> parameter, and although HalluDesign was seeded per job, its vendored inverse-folding
> wrapper reseeds the global random state from operating-system entropy before sequence
> sampling. Aggregate distributions are reproducible; individual designs are not. The
> metric and selection code is deterministic given the archived structures.

---

## Note on the three open numeric conflicts

**Foldseek / Leiden resolution — resolved, 8.0.** The chapter flags a conflict between
8.0 and 0.80 for the chosen fine partition. It is **8.0**. `RESOLUTIONS` in
`src/binder_map/clustering.py` lines 10–11 maps `cluster_r80` to `8.0`; the `r80` label
is shorthand for that resolution and is not a decimal. No change to any reported
clustering result follows — only the stated parameter.

**The delivered library was clustered at k = 30, not k = 15, and the coverage figures in
§2.9 are wrong.** §2.9 states k = 15 and reports "917 fold and 341 epitope clusters out
of 2,166 and 1,190 over the gated 31,763 backbones (42% / 29%)". That sentence mixes two
different maps.

The delivered library is `csy4_binder_library_12k_v6_k30.csv`, selected on
`results_sw_k30_f50_mut`; `submit_maps.sh` records the sweep parameters as "frozen to
match `results_sw_k30_f50_mut` exactly". The k = 15 map is a sweep arm, not the delivered
selection. Read directly from that map's cluster tables
(`design/clusters/{fold,epitope}/clusters.parquet`):

| | gated backbones | clusters at resolution 8.0 |
|---|---|---|
| fold axis | 31,763 | **917** |
| epitope axis | 31,763 | **437** |

So 917 is the k = 30 fold count, 2,166 is the k = 15 fold count, and 1,190 is not the
cluster count of either map. Against the correct denominators the panel covers **917 of
917 fold clusters (100%)** and **341 of 437 epitope clusters (78.0%)** — not 42% and 29%.

The library itself is unaffected; only the reported coverage changes. Parameters for the
delivered map, for the chapter: Foldseek v8-ef4e960 all-against-all,
`--alignment-type 1` (TMalign), `-c 0.8`, `--max-seqs 300`, `TM = min(qTM, tTM)`, floor
0.50, **mutual** kNN at **k = 30**; epitope axis a 48-dimensional fingerprint
(16 nucleotides × phosphate/sugar/base) over protein–RNA heavy-atom contacts within
4.5 Å, cosine kNN on L1-then-L2-normalised rows, floor 0.6, mutual, k = 30; Leiden under
the RBConfiguration quality function on both axes.

One further correction while this sentence is being rewritten: the chapter describes the
selection as "run under a fixed seed of 42". `SEED = 42` is consumed only by the
`random` comparison arm. The hierarchical method that produced the delivered library is
deterministic by construction (`rank(method="first")` and stable sorts), so calling it
seeded is misleading — it needs no seed.

**HalluDesign refinement time steps.** State the production value, 60, as the setting
used. No further commentary.
