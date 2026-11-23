# Supplementary tables T1-T6

Working labels only. **Supplementary numbering is not settled**, so no S-number
appears in any filename here. Renumbering = rename the files and the caption
entries; nothing depends on the number.

⚠️ **One convention conflict to settle.** The slugs `T1_`-`T4_` were specified for
this task, but `RUNLOG.md` (2026-11-19, "SM plan superseded: slugs, not
S-numbers") sets a repo-wide rule that **numbers must not appear in any path**,
and the sibling supplementary figures were renamed to bare slugs for exactly that
reason (`FigS2_metric_reliability` -> `seed_reliability` -> `metric_reliability`, etc.). The T-numbers
here are stable working labels for table *identity*, not citation order, so they
do not force a renumber when the main text changes - which is what that rule was
protecting against. If you want strict consistency anyway, rename to
`construct_inventory`, `primers`, `route_statistics`, `container_inventory`:
only this README and the two generators' output filenames refer to them, and
nothing reads them programmatically.

Built 2026-11-19, all six complete 2026-11-23. Each table is written as both a
machine-readable `.csv` and a `.md` rendering, by slug. T3, T4, T5 and T6 are
**generated** from primary files by the scripts beside them; T1 and T2 are
transcriptions measured from primary records, rendered by
`make_T1_T2_renderings.py`.

**Each `<slug>.md` is a rendering of its `.csv` and nothing else.** Provenance -
where a value came from, what contradicted it, what is still unverified - lives in
`T2_findings.md`, which is not part of this archive, never in the table files. T3-T6 need no
findings file: their generators read the primary records directly, so the script
*is* the provenance.

| slug | state | generator |
|---|---|---|
| `T1_construct_inventory` | complete - **87 plasmids**, eight of them redacted (see below); ordered synthesis fragments excluded by decision | `make_T1_T2_renderings.py` (`.md` only; the `.csv` is measured from the DNA records) |
| `T2_primers` | complete - **9 rows, 9 sequences**; the last two came from the group's primer register rather than from this project tree | `make_T1_T2_renderings.py` (`.md` only) |
| `T3_route_statistics` | complete and reproducing | `make_T3_route_statistics.py` |
| `T4_container_inventory` | complete - two versions are not recorded **inside the images** and say so; the digest is their identifier | `make_T4_container_inventory.py` |
| `T5_candidates` | **complete, no `NEED:` column** - it was listed as blocked on screen data and is not | `make_T5_candidates.py` |
| `T6_compute_requests` | complete - 7 production stages, parsed from the code archive | `make_T6_compute_requests.py` |

Run the generators from the repo root:

```bash
~/.venvs/binder_map/bin/python3 code/supp/tables/make_T3_route_statistics.py
~/.venvs/binder_map/bin/python3 code/supp/tables/make_T4_container_inventory.py
~/.venvs/binder_map/bin/python3 code/supp/tables/make_T5_candidates.py
```

Both carry `assert`-based self-checks that fail loudly if a source changes; the
numbers in this README are what those checks pin.

Terminology throughout: a generator + inverse-folder combination is a **design
route**, never a "stack". Literal CSV column names and the `stack` column inside
primary files are left as they are. Fixed route order everywhere:
RFdiffusion3+LigandMPNN, BoltzGen+LigandMPNN, BoltzGen+BoltzIF,
BoltzDesign1+LigandMPNN, HalluDesign+LigandMPNN.

---

## Sources, with mtimes

Copied into `sources/` so the generators do not need `ssh`:

| file in `sources/` | origin | origin mtime |
|---|---|---|
| `result.txt` | `euler:/cluster/work/wicky/lvalsecchi/gpu_hours/result.txt` | 2026-10-06 10:41 (4,271 B) |
| `rollup.py` | `euler:/cluster/work/wicky/lvalsecchi/gpu_hours/rollup.py` | 2026-10-06 10:42 (1,236 B) |
| `cascade_by_stack.csv` | `code/Fig3e_gpu_efficiency/methods/sources/gate/cascade_by_stack.csv` | 2026-11-12 12:51 (6,247 B) |
| `container_digests.txt` | `/Volumes/wicky/users/lvalsecchi/euler_work/container_digests.txt` | 2026-11-15 15:22 (1,034 B) |

Read in place, not copied:

| file | mtime |
|---|---|
| `inputs/Fig3e_gpu_efficiency.csv` | 2026-11-12 12:16 (989 B) |
| `inputs/csy4_binder_library_12k_v6_k30_CORRECTED.csv` | 2026-11-17 22:05 |
| `code/Fig3e_gpu_efficiency/methods/METHODS.md` | 2026-11-12 12:54 |
| `euler:.../rna_binder/analyses/cross_campaign/csy4_master_foot.csv` | 2026-09-01 13:35 (14,553,905,627 B = 14.55 GB) |
| `euler:.../rna_binder/analyses/foldseek_map/passing_v6_final.csv` | 2026-09-01 17:12 (50,007,971 B) |
| `/Users/leo/LIBRARY_SELECTION_REPORT.md` | 2026-09-04 11:24 |

`sources/container_digests.txt` is byte-identical to the Euler-side copy
`/cluster/work/wicky/lvalsecchi/container_digests.txt` (both 1,034 B, same
mtime).

---

## T3 - per-route campaign statistics

### The GPU-hour trap, and what the column headers mean

GPU-hours per route exist in **two quantities that must never be mixed**. Both
are in the table, each in its own labelled column, because the figure uses one
and the campaign record is the other:

**`GPU-h (campaign ledger, sacct)`** and
**`Designs/GPU-h (campaign ledger, sacct)`** - the real campaign spend.
GPU-h = `ElapsedRaw x gres/gpu / 3600`, summed over **every allocation,
including requeues and failures**. From `sources/result.txt`; re-roll with
`ssh euler 'cd /cluster/work/wicky/lvalsecchi/gpu_hours && python3 rollup.py < result.txt'`
(`rollup.py` reads the result block on **stdin** - running it with no input
silently prints all zeros).

**`Designs/GPU-h (pinned-RTX-4090 benchmark)`** - `gpubench_v1`, campaign-identical
jobs re-run on pinned RTX 4090s, computed as `1 / (gpuh_gen + gpuh_refold)` from
`inputs/Fig3e_gpu_efficiency.csv` (upstream
`euler:.../rna_binder/bench_gpueff/gpuh_per_design.tsv`). **This is what the
repo figure `code/Fig3e_gpu_efficiency` plots.**

Its own `methods/METHODS.md` section 1 argues the ledger **cannot rank routes**,
for four reasons, all confirmed present in that file: only 68% of campaign GPU-h
ended `COMPLETED`; routes ran on different GPU generations, confounding GPU with
route; BoltzGen generation is one job serving two routes with no recorded split;
and the ledger has no per-route design denominator. The ledger column is here as
the **campaign record**, not as a ranking.

### The shared-BoltzGen convention - two different conventions, do not merge them

The two BoltzGen routes share one generation allocation, and the two columns
split it **differently**:

- **Benchmark** (header of `inputs/Fig3e_gpu_efficiency.csv`): each BoltzGen
  route is charged the **full shared backbone stage (16.7% of the job)** plus
  **only its own inverse folding** (BoltzIF 72.6%, LigandMPNN 10.7%). The two
  shares therefore **deliberately exceed 100%** of the job - the question being
  answered is "what would this route cost run alone". Reproduced from the
  benchmark TSV: shared job `20640689` = 1.2939 GPU-h for 1,200 designs
  (600 + 600), so BoltzGen+LigandMPNN = `1.2939 x 0.274 / 600` = 0.000591
  (file: 0.000592) and BoltzGen+BoltzIF = `1.2939 x 0.893 / 600` = 0.001926
  (file: 0.001930).
- **Ledger** (this table's sacct columns): each BoltzGen route is charged the
  **entire shared job**, all 6,671.4 GPU-h of it - a cruder over-attribution
  than the benchmark's percentage split, and **not** the 16.7/72.6/10.7 rule.
  `rollup.py` itself does **not** do this: it reports BoltzGen+BoltzIF 17,728.9
  and BoltzGen+LigandMPNN 17,959.7 with the shared 6,671.4 listed separately.
  The 24,400.3 / 24,631.1 figures are that convention applied on top of
  `rollup.py`'s output.

**Consequence, stated in the table:** the ledger GPU-h column sums to
**190,724.6**, which exceeds the five-route total by exactly the shared BoltzGen
job (6,671.4). That is not an error in the column; it is the convention.

### Totals - both correct for their own scope, never in conflict

Re-derived by running `rollup.py < result.txt`:

- **184,053.2 GPU-h** - the five design routes, shared BoltzGen counted **once**.
- **197,796.5 GPU-h** - total project GPU time, which adds apo self-consistency
  13,603.2, seed self-consistency 26.1, dev/smoke tests 72.1, and 41.9
  unattributable refold GPU-h (`refold_?` 34.4 + `refold_AF3` 2.7 +
  `refold_B2` 3.4 + `refold_ESM` 1.4).
- `184,053.2 + 13,603.2 + 26.1 + 72.1 + 41.9 = 197,796.5` exactly. The apparent
  "one off" against a quoted 197,797 is a rounding artefact of the grand total.

### Design counts

Verified **directly against the primary master**, not against any summary:

```bash
ssh euler 'cut -d, -f2 /cluster/work/wicky/lvalsecchi/rna_binder/analyses/cross_campaign/csy4_master_foot.csv | sort | uniq -c'
```

Returned RFD3+LiMPNN 1,524,796 / BG+LiMPNN 1,519,958 / BG+BIF 1,480,140 /
BD1+LiMPNN 1,444,614 / HD+LiMPNN 1,023,007, plus one header line, summing to
**6,992,515**. These match `inputs/Fig3e_gpu_efficiency.csv`'s `generated`
column, which is what the generator reads.

### Gate pass rate - and which BoltzIF variant

The applied gate produced **55,706** survivors from 6,992,515. Confirmed in the
header of the primary export, not inferred:

```
# binder_filter export  n=55706  store_n=6992515
```
(line 1 of `euler:.../analyses/foldseek_map/passing_v6_final.csv`; the filter
JSON is line 2, and a copy is at
`code/Fig3e_gpu_efficiency/methods/sources/gate/campaign_gate_55706.json`).

`cascade_by_stack.csv` holds **two variants** and they must not be mixed. T3
uses **`old_bif`**, for two independent reasons:

1. It is the only variant summing to exactly 55,706
   (27,881 + 14,752 + 9,042 + 3,229 + 802). `new_bif` - the true-backbone
   rescore - gives BoltzGen+BoltzIF **9,515** instead of 9,042 and totals
   **56,179**. The generator asserts the per-route sum equals the file's own
   `ALL` row for the chosen variant, so swapping variants cannot pass silently.
2. It is what the **library selection actually consumed**: the `eligible` column
   of `LIBRARY_SELECTION_REPORT.md` is 27,881 / 14,752 / 9,042 / 3,229 / 802,
   i.e. the `old_bif` numbers.

`old_bif` is campaign history. If a future table wants the corrected
BoltzGen+BoltzIF rescore, switch `CASCADE_VARIANT` in the generator and say so in
the caption - do not mix rows from the two.

### Panel share

**Counted from the delivered library itself**, not transcribed from the audit:
`inputs/csy4_binder_library_12k_v6_k30_CORRECTED.csv`, 12,000 rows, `stack`
column -> RFD3+LiMPNN 3,132 / BG+LiMPNN 3,266 / BG+BIF 2,400 / BD1+LiMPNN 2,400 /
HD+LiMPNN 802 = **12,000**. This independently reproduces the quotas in
`LIBRARY_SELECTION_REPORT.md` exactly.

### Footnote T3 needs: HalluDesign's panel share is a ceiling, not a choice

**HalluDesign+LigandMPNN fell below its own protected floor.** The allocation
rule is a protected floor of 1,200 per route, topped to an equal share of 2,400,
remainder distributed by ipTM under a cap of 1.5x base (3,600). Only **802**
HalluDesign designs survived the gate against that 1,200 floor, so the floor
could not be honoured and the **entire** HalluDesign pool was taken. Its
**802/12,000 (6.7%)** panel share is therefore a ceiling, not an allocation
decision, and the mechanism degrades silently to "take everything" with no
warning raised. Audit: `/Users/leo/LIBRARY_SELECTION_REPORT.md`
(2026-09-04 11:24), section "What the audit found worth knowing".

### One correction to the brief this table was built from

The brief described the ledger and the benchmark as differing by **"~3.7x"**.
Per route the ratio of benchmark to ledger designs-per-GPU-hour is **not
constant**: 3.74x (RFdiffusion3+LigandMPNN), 3.78x (BoltzGen+LigandMPNN), 2.95x
(BoltzGen+BoltzIF), **2.40x** (BoltzDesign1+LigandMPNN), 3.64x
(HalluDesign+LigandMPNN). A single factor would misdescribe the table, and the
spread is itself consistent with METHODS section 1: the ledger's waste and
GPU-generation confound are distributed unevenly across routes. The generator
asserts the range rather than a single factor.

---

## T4 - container and version inventory

### VERIFY: the count is nine - confirmed, and every digest re-computed

`container_digests.txt` holds **exactly nine** sha256 lines plus a `DONE`
trailer. The generator asserts `len(entries) == 9` and that the image set matches
its version table, so the count cannot drift unnoticed.

**No digest was retyped.** The generator reads them out of
`sources/container_digests.txt` with a regex, so a transcription slip is
impossible by construction.

Stronger than transcription: all nine were **independently re-computed** on
Euler with `sha256sum` over the live 63 GB of `.sif` files, and **all nine
reproduce exactly** (zero-line `diff` against the recorded file).

These are `sha256sum` digests **of the `.sif` files on disk**, not OCI image
digests. Do not describe them as image digests.

### Versions - read out of the images, not from notes

Every version in T4 was obtained by running a probe **inside the image**
(`ssh euler 'apptainer exec <image> bash -lc "<probe>"'`); the exact probe is in
T4's `version_probe` column. Verified this way:

| image | tool | version |
|---|---|---|
| `foundry_0.1.9.sif` | rc-foundry (carries RFdiffusion3, IPD) | 0.1.9 (torch 2.7.1) |
| `boltzgen.sif` | BoltzGen | 0.3.2, commit `31d9d9b9c72245b4ed6fe8742d6fbf4e1a3552a0` (torch 2.12.0) |
| `boltz_design1.sif` | BoltzDesign1 | commit `31aeedd7f187ba88dd9e92514070aee13381cb93` only |
| `LigandMPNN.sif` | LigandMPNN | commit `091ab1ff5fb4d13854cf6a7c41ec531e1d9d3e67` only (torch 2.2.1) |
| `alphafold3.sif` | AlphaFold3 | 3.0.0 (jax 0.4.34) |
| `boltz2.sif` | Boltz-2 | 2.2.1 (torch 2.8.0) |
| `esmc.sif` | ESMFold2 (ESMC-6B backbone) | esm 3.3.0 (transformers 4.57.6, torch 2.6.0+cu124) |
| `pyrosetta.sif` | PyRosetta | 2025.37+release.df75a9c48e |
| `af3_protenix_esm.sif` | HalluDesign predictor image | alphafold3 3.0.0, fair-esm 2.0.0 (torch 2.6.0+cu126) |

The AF3 3.0.0 / Boltz-2 2.2.1 figures quoted in
`code/Fig3e_gpu_efficiency/README.md` had no machine-readable source behind them
(prose only); they are now confirmed from the images themselves.

### Which image each route used - from the configs, not assumed

`grep -rn 'sif' euler:.../rna_binder/code/configs` gives an unambiguous mapping:

- `rfd3/default.yaml`, `rfd3__multif/default.yaml` -> `foundry.sif` + `LigandMPNN.sif`
- `boltzgen/default.yaml`, `boltzgen__multif/default.yaml` -> `boltzgen.sif` + `LigandMPNN.sif`
- `boltzdesign/default.yaml`, `boltzdesign__multif/default.yaml` -> `boltz_design1.sif` + `LigandMPNN.sif`
- `halludesign/default.yaml`, `halludesign/random_init.yaml` -> `af3_protenix_esm.sif` **only**
- `refold__af3` -> `alphafold3.sif`; `refold__b2` -> `boltz2.sif`; `refold__esm` -> `esmc.sif`

`foundry.sif` is a **symlink** -> `foundry_0.1.9.sif`
(`readlink -f` confirms), which is why the configs name `foundry.sif` and the
digest file names the resolved target. They are the same image.

`pyrosetta.sif` appears in **no** config. Searched for separately, it turns up
only in the metrics/aggregation layer - `ops/submit/{fullmetrics_chunk,hbbb_array,aant_chunk,iface_time,seed_metrics_array}.sbatch`
and `ops/agg/{hbond_map_job,hbond_map_af3_job,annotate_dockq_hbonds}.slurm`. It
is an analysis image (interface H-bond and packing metrics), not part of the
design or refold pipeline, and T4 says so rather than implying it ran in the
campaign.

### Two coverage gaps this inventory does NOT close

1. **HalluDesign's own source is not in any of the nine images.** The
   HalluDesign configs mount only `af3_protenix_esm.sif`, which is the
   predictor. HalluDesign's code - and its **vendored** LigandMPNN, which is
   what inverse-folds that route rather than `LigandMPNN.sif` - lives in an
   on-disk checkout at `euler:/cluster/work/wicky/lvalsecchi/halludesign_upgrade`
   (branch `upgrade-1b08312`, HEAD `42ccb00d0c786731e23e40e25f5ed2abfcbc1eca`).
   Cite that commit alongside the digest or the HalluDesign route is
   under-specified.
2. **Model weights mounted from outside the images are not covered by any
   digest**: the LigandMPNN checkpoint
   `/cluster/work/wicky/params/foundry/ligandmpnn_v_32_010_25.pt`, the AF3
   weights blob `/cluster/work/wicky/params/af3_params/af3.bin` (carries no
   version string), and the ESMFold2 Hugging Face cache
   `rna_binder/params/esmfold2_hf_cache` (checkpoint revision not recorded). A
   digest of the image does not pin the weights it loads.

### Not recorded inside the images - the digest is the identifier

These are absences, not gaps to fill, and the table says so in final wording
rather than carrying a drafting marker:

- `boltz_design1.sif` - no version string anywhere in the image; the git commit
  is the only identifier it carries.
- `LigandMPNN.sif` - same; commit only.
- `foundry_0.1.9.sif` upstream **commit** - the image has no git metadata and no
  `git` binary, so the commit cannot be read out of it. The version 0.1.9 is
  solid (dist-info `METADATA`, and the filename agrees).
- Protenix version/checkpoint inside `af3_protenix_esm.sif` - no pip package, no
  version string; weights ship inside the image, so the digest is the only
  identifier.
- `esmc.sif` build recipe - built from a `docker-archive` tar that no longer
  exists.

### Flagged: the version-control claim about the Euler code

**Euler's `rna_binder/code` is not a git repository, and T4 is the only way the
computational environment can be specified.** (This is about the Euler pipeline
tree, not about `LV_ThesisFigures/code/`, which is in this repo's git.) No
history was recovered with `git log`; none exists to recover.

Checked directly, and the picture is more precise than "Euler was never in git":

| path on Euler | git? |
|---|---|
| `/cluster/work/wicky/lvalsecchi/rna_binder` | **not a git repo** |
| `/cluster/work/wicky/lvalsecchi/rna_binder/code` | **not a git repo** |
| `/cluster/work/wicky/lvalsecchi/rna_binder/ops` | **not a git repo** |
| `euler_work` | **not on Euler at all** (checked 2026-11-19): it is the NAS directory `/Volumes/wicky/users/lvalsecchi/euler_work`, and it has no `.git` either |
| `/cluster/work/wicky/lvalsecchi/halludesign` | is a git repo |
| `/cluster/work/wicky/lvalsecchi/halludesign_upgrade` | is a git repo (branch `upgrade-1b08312`) |

So: **any Methods sentence claiming version control for the project's own Euler
pipeline code is false** - `rna_binder/code` and `rna_binder/ops` on Euler, and
the `euler_work` directory on the NAS, have no VCS at all. The two **vendored
upstream HalluDesign checkouts** are genuine git repos with recoverable commits,
so a blanket "Euler was never in
git" is also slightly too strong, and the HalluDesign commit above is citable
because of it. Write the narrow version of the claim.

---

## T1 - construct inventory

**87 plasmids**, measured from the DNA records rather than from journal labels,
which are unreliable here. Three scope calls, settled by the author: ordered
synthesis fragments are out (ordered DNA, not constructs), `pART255` is in although
it was built and then dropped from the panel, and the table is scoped to this work
rather than to the group part registry.

`size_bp` is the full plasmid length throughout. The one derived value is
`pART262`, 2,960 bp, whose only as-built record is a 795 bp insert window.

**Eight rows carry `Redacted for privacy purposes`.** They are vectors supplied by a
collaborator; the designs are theirs, not ours to publish. The identifiers stay, so
that a reader tracing one of our constructs back to its backbone finds the row and
knows why it stops there. The insert description, the plasmid length and the four
type IIS counts are withheld.

One check is still open and owned elsewhere: the ordered-scaffold list, the only
thing that can confirm each `pART243`-`pART254` file name is the number that
construct was built under. It does not block the table, since every row is measured
from sequence.

## T2 - primers

**9 rows, 9 sequences.** Complete as of 2026-11-23. Two things are worth carrying
forward from how it got there.

**The sequences were never in this project's tree.** They are in the group primer
register on Matteo's side of the NAS,
`mlampis/Matteo_Lampis/Matteo_Lampis_Cloning/Matteo_Lampis_PrimerList.txt`, with a
second copy in `2022_07_19_Labprimers.csv`. An earlier audit concluded the list was
"not on the NAS" after searching only `Leo_Valsecchi/`. `MicroSynthStandardPrimers.csv`
sits in the same directory and holds M13; the `PR####` primers are lab stock, not
provider primers.

**A reverse primer's row gives the oligonucleotide, not the site it reads.** PR249
and PR250 were first entered as the plus-strand sites annotated in `pLV005.gb`; the
register has the primers themselves, which are the reverse complements. The renderer
asserts that relation so the annotation and the register check each other.

M13 is listed once, forward only: every read in `LV_Sequencing` is
`LV1.NN.MM_<plate>_M13.ab1`, 56 of them, and no `M13r` read exists anywhere in the
tree. Every search, including the ones that came back empty, is in `T2_findings.md`.

## T5 - the 32 candidates

### ⚠️ It was listed as blocked on screen data. It is not, and nothing in it is `NEED:`

The build plan and `SM_BUILD_STATUS_2026-11-19.md` both mark T5 `PENDING` on screen
data. The pooled screen data has been on the NAS since 2026-10-19 at
`LV_BinderDesign/LV2.1.1_RNAbinders_PoolCraft_V1/11_alignment/counts__*.csv`. **Every
column the spec asks for is filled from a primary file**, so T5 ships with no `NEED:`
column at all. The two caveats it does carry are convention and provenance caveats,
below, not gaps.

### The gate is reproduced from the count tables, not transcribed

`reads >= 75` in POS **and** matched NEG, `splicing ratio >= 0.07`,
`fold change >= 10`, with the screen report's own pseudocount (NEG splicing =
1/NEG reads when no NEG read was spliced). Off the eight count tables this returns
**7 (context A) + 25 (context B) = 32**, matching `gated_Q1_75Reads_OrderedHits.xlsx`
row for row. `screen()` asserts it. The screen's *default* gate (`min_reads` 100,
`min_splicing` 0.15, `fc` 5) returns 6, from
`12_analysis/report_render_provenance.json`.

The count tables key on the alignment `reference_id`, not on `design_id`;
`11_alignment/reference_manifest.csv` is the only file tying the two together, and the
tie is **per context**. Every context carries its own 12,000 references with
`context_equivalence_size == 1` on all 48,000 rows, so the mapping is one-to-one
inside a context. The generator asserts that too.

**The 32 are a selection on ranked splicing ratio and fold change, not a hit rate.**

### Two conventions kept in separate columns, never merged

1. **ligand_RMSD and sc_DockQ.** The `stored` columns are what the campaign gate actually
   saw and are filled for all 32. For BoltzGen+BoltzIF those were scored against
   BoltzIF's own Boltz re-prediction rather than against the generated backbone, so
   the `corrected` columns carry the true-backbone rescore
   (`inputs/large/bif_true_per_backend.csv.gz`, 2026-11-17 15:18) and are **blank for
   every other route**, which has no defect to correct. The generator asserts that all
   nine BoltzGen+BoltzIF rows have a rescore.
2. **Fold change.** `fold_change_uncapped` is the raw ratio; `fold_change_capped50` is
   the screen report's export cap at 50, which is what the xlsx exports. **10 of the
   32 sit at the cap**, so the columns are not interchangeable and a median or a
   ranking differs between them.

The other selection metrics follow the gate's own combination rule across the two gate
backends: `sc_rmsd` is the **min** and `dockq` the **max** over AF3 and ESMFold2
(both were OR terms), `footnb_010` is the **max** (an AND term), `iptm_max` and
`seq_run` are single-valued.

### Sources, with mtimes

| file | mtime | what is taken |
|---|---|---|
| `.../11_alignment/counts__{POS,NEG}_{A,B}_rep{1,2}.csv` | 2026-10-19 10:50 - 16:37 | reads and splicing per member, 8 x 12,000 rows |
| `.../11_alignment/reference_manifest.csv` | 2026-09-30 16:58 | `reference_id` -> `binder_design` per context |
| `inputs/candidates32_mapped.csv` | 2026-11-09 18:56 | `twist_id` (= `ID_Number` = the Q1 number), generator, inverse folder, flow well |
| `inputs/csy4_binder_library_12k_v6_k30_CORRECTED.csv` | 2026-11-17 22:05 | length, fold cluster, epitope cluster |
| `inputs/large/gated_pool_v6_gate_metrics.csv.gz` | 2026-11-19 14:42 | the five gate terms, including `footnb_010` which no other local file carries |
| `inputs/large/bif_true_per_backend.csv.gz` | 2026-11-17 15:18 | the corrected BoltzGen+BoltzIF ligand_RMSD and sc_DockQ |
| `inputs/lv2_34_1_flow_on_off.csv` | 2026-11-09 18:47 | arrayed on/off, derived from `expdata/flow/lvalsecchi/LV2_34_1/LV2_34_1_Analysis.xlsx` |

`gated_pool_v6_gate_metrics.csv.gz` was extracted for `code/supp/metric_significance`
from `euler:.../analyses/foldseek_map/passing_v6_final.csv` (header line
`# binder_filter export  n=55706  store_n=6992515`); the exact `awk` command is in that
figure's README. All 12,000 ordered designs are in the 55,706, so the join is total.

### What the table says

Candidates per design route, which independently reproduces the per-route counts in
`MAINTEXT_DATA_REQUESTS` A4 and is the denominator `bif_quasicontrol` was waiting on:

| design route | candidates | panel members |
|---|---|---|
| RFdiffusion3 + LigandMPNN | 4/3,132 (0.13%) | 3,132 |
| BoltzGen + LigandMPNN | 11/3,266 (0.34%) | 3,266 |
| BoltzGen + BoltzIF | 9/2,400 (0.38%) | 2,400 |
| BoltzDesign1 + LigandMPNN | 7/2,400 (0.29%) | 2,400 |
| HalluDesign + LigandMPNN | 1/802 (0.12%) | 802 |

- Screen: POS reads 82 to 664, splicing ratio 0.072 to 0.505, 10/32 at the FC cap.
- Arrayed flow: on/off 0.478 to 1.563; **15/32 (46.9%) are proposed binders**
  (on/off >= 1.2), the other 17/32 (53.1%) inactive designs. Never a bare
  "binder"/"non-binder".
- BoltzGen+BoltzIF under the corrected convention: median sc_DockQ 0.740 -> **0.059**,
  median ligand_RMSD 2.08 -> **19.30 A**, and **1/9** would still clear the gate's
  ligand_RMSD and sc_DockQ terms on true-backbone values. Those nine entered the panel under
  a materially weaker filter, which is the `bif_quasicontrol` argument - but they
  still passed the other three terms and were still ipTM-ranked within cluster, so it
  is a weakened filter, not no filter.

### Caveats that belong in the caption

- **The flow readout is a reporter proxy for receptor dimerisation**, and for avidity
  on a multivalent scaffold. It is not affinity: no Kd was measured anywhere in this
  work.
- **One well per condition, no replicates**, and binder plasmid went in at 84 ng
  against a planned 100 ng, uniform across variant wells.
- ✅ **SETTLED 2026-11-20: the pLV15 / pLV16 question was not a substitution and
  not an accident.** pLV15 exists (`pLV015/pLV015.gb`, 5,879 bp, 8 NegCtrlLoops,
  from pART256) and is posB's matched negative - `5,951 - 5,879 = 72 = 8 x 9`,
  the same per-site arithmetic as pART252 -> pART256. But it was **never planned
  and never transfected**: 0 occurrences in the LV2.33 plan and in the
  as-performed record, against **41** for pLV16 as the planned MM3. So the screen
  used **one scaffold negative for both positives, by design**. Nothing was read
  against the wrong negative and no ON value inherits an error; the earlier
  "by accident" reading is withdrawn. T5 is not blocked on this.
- **Context A is sequenced to 15.4 M reads against context B's 40.8 M**, which is
  part of why the split is 7 A against 25 B. See `code/supp/pooled_readout/`.

## T6 - per-stage compute requests

Seven rows: four generators (the two BoltzGen routes share one generation job) and
three re-prediction backends. Built by `make_T6_compute_requests.py`.

### Where every value comes from

Nothing is retyped. The generator reads the `slurm:` block of the archived
production config for each stage out of `~/rna-binder-csy4` at tag
`v1.0-thesis-freeze`, and `CONFIGS` in the script names the file per row:

| Stage | Config read |
|---|---|
| RFdiffusion3 | `configs/rfd3__multif/default.yaml` |
| BoltzGen (both routes) | `configs/boltzgen__multif/default.yaml` |
| BoltzDesign1 | `configs/boltzdesign__multif/default.yaml` |
| HalluDesign, de novo | `configs/halludesign/random_init.yaml` |
| AlphaFold3 | `configs/refold__af3/default.yaml` |
| Boltz-2 | `configs/refold__b2/default.yaml` |
| ESMFold2 | `configs/refold__esm/default.yaml` |

The `__multif` variants are the production configs: `MIDRUN_V6_AS_RUN_CONDITIONS`
names them per stack. HalluDesign's production config is `random_init.yaml`, the
de-novo path; `halludesign/default.yaml` is the poly-ALA-helix fallback and is not
the archived production route.

### The two judgement calls

**Configured, not as-run.** Campaign submission overrode the scheduler fields, so
partition, GPU request and wall time in the table are what the config asked for. The
`as_run_deviation` column carries the routing that actually ran, each entry sourced to
`MIDRUN_V4_AS_RUN_CONDITIONS_2026-07-05.md` or `MIDRUN_V6_AS_RUN_CONDITIONS_2026-07-08.md`
in the archive. The largest deviations: the Turing-capable generators moved to the
public `gpu.*` pool with `gpu:1`, and all three re-predictions ran on `.bulk` tiers at
`06:00:00`.

**Local scratch.** No generation config requests any; the re-prediction stages ran
through node-local scratch workers at `--tmp=30G`, which is a property of the as-run
worker and not of the config, and the cell says so.

### Checks

- The generator asserts exactly 7 rows, so a config going missing fails loudly rather
  than silently shortening the table.
- Processor counts record which field the config used (`cpus_per_task`, `ntasks` or
  `ntasks_per_node`) rather than flattening them, because HalluDesign is the one stage
  that tolerates `--ntasks` and the others hard-abort on it.
- Memory per job is computed, not stated: processors x `mem_per_cpu`.

## Naming, settled 2026-11-20

**`sc_DockQ` and `ligand_RMSD (Å)`**, matching the main text (Fig 2c), wherever these two
quantities are an axis label or appear in prose. The bare `DockQ` and the hyphenated
`ligand-RMSD` are not used as display forms anywhere in `code/supp`.

Renamed **`sc_RMSD` -> `ligand_RMSD` on 2026-11-21**: the gate's first term is the
RNA-aligned *ligand* RMSD (`dockq_LRMSD_af3` / `_esm`), a pose metric, not the Ca fold
RMSD it was once written up as. Same measurement throughout, new name.

⚠️ **CSV column names are untouched and stay untouched** — `sc_dockq`, `sc_rmsd`,
`dockq_af3`, `dockq_esm`, `rmsd_af3`, `rmsd_esm`, `dockq_self`, `dockq_LRMSD`,
`dockq_fnat`, `dockq_iRMSD` and the `_true` variants all keep their lower-case spelling,
because they are what is on disk. The sweeps only ever matched capitalised display forms,
never a column name.

## Vocabulary, settled 2026-11-19, swept 2026-11-20

**"re-predictor"** for AF3 / Boltz-2 / ESMFold2, never "refolder", "predictor" or
"repredictor". **"refolding seeds"** for the five seeds. **"design route"**, never
"stack". In all prose, labels and legends.

⚠️ **Identifiers are deliberately exempt. Do not "fix" them — each one breaks a read or
an import:** the `refolder` and `stack` columns in the source CSVs, `palette.REFOLDERS`,
`palette.REFOLDER_ORDER`, `palette.STACK_ORDER`, `palette.SHORT`, the source file names
that contain `_by_stack_refolder`, and the directory name `Fig5fg_h1null_refolder`.
