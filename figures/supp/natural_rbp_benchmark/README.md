# natural_rbp_benchmark

*Was `F8` in the 2026-10-02 set. Built under its slug: supplementary numbers are assigned
in one pass once the main text stops moving, so no S-number appears in any file name.*

Designed interfaces against a natural RBP–RNA benchmark. **Four** bare panels.
Rebuilt 2026-11-18 from `~/rna_binder_plots/thesis_figures/fig_F8.py` (2026-10-02).
Revised 2026-11-19 on Leo's comment — the panels are normalised, see below — and the
H-bond distribution panel was **removed** the same day on his call, see
*What the removed panel carried* below. A panel describing the benchmark set itself was
**added 2026-11-20**, see *The benchmark set, drawn as (d)* below.

```bash
~/.venvs/binder_map/bin/python3 code/supp/natural_rbp_benchmark/make_inputs.py   # once
~/.venvs/binder_map/bin/python3 code/supp/natural_rbp_benchmark/make_figure.py
```

Emits `output/natural_rbp_benchmark/panel_{a,b,d,e}_*.{pdf,png,svg}`, 85 mm wide.
`make_inputs.py` needs re-running if the panel (a) size table **or** the panel (d)
selection table changes; it writes both.

⚠️ **Letter ≠ file stem, for two of the four panels.** `code/supp/make_supp_book.py`
letters panels by their position in **filename order** at assembly
(`sorted(d.glob("*.pdf"))`), and the stems still carry the letters they had when a
`panel_c` existed. The stems are deliberately *not* renamed: panel files are linked by
name from an Illustrator assembly, and renaming one breaks the link silently. The new
panel is `panel_e_*` for the same reason — naming it `panel_c` or `panel_a0` would have
re-lettered every existing panel and every caption clause with it. Read every "(c)" and
"(d)" below as the file, not the stem.

| drawn as | file | why the letters differ |
|---|---|---|
| (a) | `panel_a_size_confound.*` | — |
| (b) | `panel_b_normalisation.*` | — |
| (c) | `panel_d_base_reading_residues.*` | the old `panel_c` was deleted 2026-11-19 |
| (d) | `panel_e_dataset_overview.*` | added last so (a)–(c) keep their letters |

**Page layout, measured from the emitted SVGs** (the book packs greedily into
A4_W 210 − 2 × MARGIN 14 = **182 mm** rows with a **6 mm** gutter, in filename order):
(a) 84.3 + 6 + (b) 82.0 = 172.3 mm share the top row; (c) 79.2 + 6 + (d) 84.2 =
169.4 mm share the second. Two rows of two, 140 mm of the 283 mm page body. Before the
new panel the page was 84.3 + 82.0 on one row and (c) alone beneath it, so **(d) fills
a slot that was empty** and nothing reflowed.

## What the removed panel carried

> "I think I want to remove C." — Leo, 2026-11-19

The deleted panel drew the interface-H-bond **distribution** twice over — raw per
complex on a short strip, per interface residue on a taller one below it — with the
natural and design medians marked on each. Code (`panel_c`, `_hist`) and the emitted
`panel_c_hbond_distribution.{pdf,png,svg}` are gone.

🚨 **`git checkout` will not bring it back — these scripts are untracked** (`git ls-files
--error-unmatch` on this file: *"did not match any file(s) known to git"*). The deleted
function is kept **verbatim in `code/RUNLOG.md`**, under the 2026-11-19/20 entry for this
figure, with the one call site to restore. That is the only copy.

**Nothing it measured is lost, but two of its numbers are now text only.** These were
drawn in that panel and in no other:

| number | where it lives now |
|---|---|
| 43,378/99,268 (43.7%) of designs make **no** interface H-bond at all | `report()` stdout, and *Numbers* below |
| 5,372/99,268 (5.4%) reach the natural 5th percentile (0.151/residue) | `report()` stdout, *Numbers* — pulled from the caption, which now describes only drawn panels |
| design maximum 7 raw H-bonds against a natural maximum of 53 | `report()` stdout, *Numbers* |
| the shape of the two distributions and where they overlap | **not plotted anywhere** |

The **result** of the panel — that H-bonding is the one deficit surviving the size
correction — is unaffected: it is panel (b)'s 0.11× row, which is where it was always
argued from. `report()` still recomputes the per-residue medians on every run precisely
because reproducing that 0.11× is the check that this code path and the stored
`comparison_normalized.csv` agree; dropping the panel did not drop the check.

If the overlap statistic should be *seen* rather than quoted, it wants a home in another
figure — it does not need this one back.

🚨 **`NEED:` a citing sentence.** This figure is currently uncited and supports a claim the
main text does not presently make at all, so it needs one or it is an orphan.

## 2026-11-19 revision

> "here you need to normalize for the interface residues. but also we are analysing
> hbonds and very specific contacts for things we said befroe we are not super sure
> about" — Leo, on the panels

Both halves acted on.

**Normalisation.** The first half caught an internal contradiction. Panel (a)'s whole
argument is that an extensive count cannot be compared when the natural protein interface
is 2.18× larger, and panel (b) applies that correction — but the H-bond panel (since
removed) and the base-contact panel were still **raw counts**, committing exactly the
error (a) raises. Both were **divided by `n_contact_res_prot` per complex, before the
median or mean**, the same estimator (b) uses and `normalize_metrics.py` upstream, and the
surviving panel says so on its axis.

This changes what the base-contact panel says, and it is not a cosmetic change — see
*What normalisation did to the base-contact panel* below.

**Over-claiming.** The second half is a caution, and the numbers back it: H-bonds and
π-stacking are the two **least reproducible** metrics measured anywhere in this work.
Across five refolding seeds their ICC is **0.163–0.514** and **0.151–0.519**, against
**0.634–0.928** for burial (`code/supp/metric_reliability/`). So a per-complex H-bond count
is a noisy measurement and the figure must not read as though it were precise:

- the caption carries that caveat in one sentence, with the ICC ranges;
- the H-bond framing moved from *"most designs make none or one"* — a claim about
  individual complexes — to a statement about the two **distributions** and their
  overlap, and then off the page entirely when the panel was removed;
- the result itself stays. It is the one metric that survives the size correction, which
  is the figure's point, and it is carried by panel (b)'s single corrected ratio. It is
  simply no longer the firmest-looking thing on the page.

## The argument, which the panels must keep in this order

Panel (a) is **not a result — it is the objection**, raised before any of the data is used.
Natural RBPs are large multidomain proteins, so a raw count comparison is rigged. Because
the RNA-side breadth is nearly equal (1.10×) while the protein side differs 2.18×,
**protein interface residues** is the legitimate denominator.

Panel (b) applies it and **the correction destroys most of the apparent deficit** — base
contacts go 0.39× → 0.86×, and backbone and phosphate invert to designs doing *more*. Only
the H-bond row survives, 0.07× → 0.11×, which is why it alone is emphasised, and why the
third panel then looks at the residues that read bases.

**Do not flatten this into "designs make fewer H-bonds".** The claim is that one deficit
survives a correction that kills the rest. The normalisation is applied **per row before
the median**, never as a ratio of medians — in every panel, not just (b).

Panel (d) stands **outside** that chain — it describes the comparison population rather
than arguing anything with it — which is why it sits last rather than first. Logically it
is scene-setting and would read well in front of (a), but the book letters panels by
filename order, so moving it there would re-letter (a), (b) and (c) and every caption
clause that names them, for no gain in the argument. a → b → c is unchanged.

A second confound — `sc_value`, `packstat`, `worst_clash`, `n_interface_clashes`, all AF3
model vs crystal — is **excluded from the figure rather than normalised**, because
normalisation does not fix it. Only a staged natural refold would.

## What normalisation did to the base-contact panel, drawn as (c)

The raw panel read "designs retain arginine but lose the aromatics and base-edge amides".
Per interface residue that is **wrong on both ends**, and the corrected panel says
something sharper:

| | raw (nat ÷ des) | per interface residue |
|---|---|---|
| ARG | 1.2× | **0.6×** — designs make *more* per residue |
| LYS | 1.6× | **0.7×** — same inversion |
| PHE | 13.7× | **7.4×** — survives, still the largest deficit |
| GLN | 7.0× | **3.2×** |
| HIS | 6.0× | **3.0×** |
| TYR | 6.6× | **2.8×** |
| ASN | 3.7× | **1.4×** — mostly size |
| GLY, LEU | 1.4×, 2.9× | 1.0×, 1.1× — parity |

⚠️ **The raw column is 13.7× / 7.0× / 6.6×, not the 13.8× / 7.1× / 6.7× that
`aa_base.csv` and the earlier caption carry.** Both are right for their row set: the
stored values are over all 100,000 designs, these are over the 99,268 that have a
denominator. Dropping 732 zero-interface designs removes 732 zero base-contact rows and so
*raises* the design mean (ARG 23.38 → 23.56), which lowers the factor. Quoted on the
subset here so that the raw and corrected columns compare the same complexes. `report()`
prints both columns on every run.

So: the **aromatic** deficit and **Gln/His** survive the correction; the amide picture
splits (Gln yes, Asn largely not); and arginine, the single biggest base-contacting residue
in both populations, **inverts** — designs over-use it per residue rather than merely
matching naturals. Every row is annotated with its factor now, in both directions: a rule
that labelled only the deficits would have hidden exactly the inversions.

The ranked set also changes. Top twelve by **natural per-residue** mean is
`ARG PHE TYR SER LYS HIS ASN GLY PRO GLN GLU LEU` — ASP drops out and GLY enters, against
the raw ranking's `ARG PHE TYR ASN LYS SER HIS LEU GLN GLU ASP PRO`.

The estimator stays the **mean**, not the median: `aa_base_X` is zero for most complexes,
so a per-residue median would be zero for every row but arginine. The axis says `(mean)`
because every other panel here reports a median.

## Sources

The base-contact panel, and `report()`'s H-bond check, read the **primary per-complex
tables** —
`~/rbp_rna_db/subset59/metrics_natural_59.csv` and
`~/Documents/metric_explorer/data/metrics_100k.csv` — not the pre-aggregated
`nat_vs_design_findings/distributions.json` and `aa_base.csv`, which are **raw counts only**
and cannot be normalised after the fact. Those two files are no longer read by this figure.
Panel (b) still reads `nat_vs_design_findings/comparison_normalized.csv`; panel (a) reads
`inputs/natural_rbp_benchmark_size_confound.csv` and panel (d)
`inputs/natural_rbp_benchmark_selection.csv`, both from `make_inputs.py`.

Panel (d)'s table is built from the benchmark set's **own build artefacts** —
`~/rbp_rna_db/{cif,cif_dropped}/`, `manifest.csv`, `chain_map.tsv`,
`subset59/{manifest_59.csv,metrics_natural_59.csv}` — because the selection had no
machine-readable source at all: it existed as prose in `~/rbp_rna_db/README.md` and as the
size of two directories. Each count is recomputed with the rule its own build script used
and asserted against what that script wrote, so no count is copied out of the README.

The recomputed raw medians reproduce the stored ones exactly (natural 14, designs 1), and
the recomputed per-residue medians reproduce panel (b)'s 0.11× ratio, which is the check
that the new code path and the stored table agree.

## Colour

No route appears here, so no route colour is used. Natural is the reference distribution →
`NEUTRAL` (COLORS.md rule 3); designs → `MAIN` (rule 1). The one surviving metric is
emphasised by **weight and marker size**, never a fifth hue (STYLE_NOTES §4).

**No hatch, and there never was any.** The dot/stripe encoding withdrawn in palette v0.6
(COLORS.md rule 7, 2026-11-19) keyed *re-predictors*, and no re-predictor appears in this
figure — the design side is the 100k table, not a per-re-predictor split. Nothing to convert.

Panel (d) is entirely about the natural set, so every bar is a neutral: the funnel steps
take `NEUTRAL["300"]`, and the final step — the 58 this figure actually uses — takes
`NEUTRAL["500"]`, which is the exact colour naturals wear in (a) and (c). The ΔSASA
subset is an open bar with that same edge. No fifth hue, and the step that matters is
emphasised by **weight** (bold label) as in (b).

**Palette v0.9 (2026-11-20) does not touch this figure.** Its only change is
`scaffolds.mod1`–`mod4` moving to theme 600 `#0678a6`; this figure reads `MAIN`,
`NEUTRAL` and `SURFACE` and never touches `palette.SCAFFOLDS`. Re-rendered against v0.9:
(a), (b) and (c) are unchanged, down to their emitted sizes.

## The benchmark set, drawn as (d)

> "Can you add a general overview of the dataset we used? Either the way we selected
> them or just their size distribution or something like that." — Leo, 2026-11-20

**The selection, not the size distribution** — because the *size* of this set is already
the subject of panel (a), where it is the whole objection the figure is built on, while
nothing on the page said what the set is or where 58 complexes came from. Drawing the
size spread again would have restated (a); the funnel is the thing the figure used
without ever stating. What the set *contains* goes in the caption, below.

**The selection steps are on disk**, so nothing is invented. `~/rbp_rna_db/` holds the
three build scripts and every intermediate they wrote, and `make_inputs.py` re-derives
each count with that step's own rule rather than reading `~/rbp_rna_db/README.md`:

| step | n | re-derived from | what the step does |
|---|--:|---|---|
| RCSB search | 329 | `cif/` + `cif_dropped/` — `refine_rna_len.py` **moved** its rejects aside rather than deleting them, so the two directories together still *are* its input | ≥ 1 protein and ≥ 1 RNA entity, resolution ≤ 2.5 Å, ≤ 1,500 deposited polymer residues (which is what excludes ribosomes and spliceosomes), one representative per RCSB **30% protein-identity** cluster |
| longest RNA ≥ 10 nt | 180/329 (54.7%) | `manifest.csv`, which `refine_rna_len.py` rewrote to the kept set | `MIN_NT = 10`; drops 149 short-RNA edge cases |
| one protein + one RNA chain | 59/329 (17.9%) | `chain_map.tsv`, re-filtered with `build_subset59.py`'s own rule and asserted against the `manifest_59.csv` that script wrote | no DNA, no second copy of either polymer, so every complex has exactly one interface and a uniform A = protein / B = RNA layout |
| interface metrics computed | **58/329 (17.6%)** | `metrics_natural_59.csv` minus the rows with a `dsasa_biotite_error` | 8JOZ drops: its biotite SASA failed on a modified residue (`CM0` has no `OP3`) |

The 121 complexes the chain rule removes, in that script's own order: 23 carry DNA or
another polymer, 64 have multiple copies of **both** protein and RNA, 21 multiple protein
chains only, 13 multiple RNA chains only.

`assert list(out.n) == [329, 180, 59, 58]` runs on every `make_inputs.py` run, as do the
two reproduction checks above, so the funnel fails loudly rather than drifting.

**What the panel does not draw, and the caption carries instead.** All 58 are **X-ray**
(no cryo-EM survives; 11 were in the search), resolution 1.5–2.5 Å (median 2.0). The set
is **heterogeneous by design**: RNA 10–97 nt (median 29), deposited polymer 104–1,492
residues (median 347), protein interface residues 7–147 (median 37), RNA interface
residues 3–41 (median 11), ΔSASA 634–11,426 Å² (median 2,458, n = 53). That spread is
precisely why panel (a)'s correction is needed and is the one content fact worth stating.

## n, which is not uniform

**58 non-redundant natural RBP–RNA complexes** (X-ray, 1.5–2.5 Å). Four numbers circulate
and mean different things:

| n | what it is |
|---|---|
| 329 | RCSB search representatives before the RNA-length filter. Drawn as the first bar of (d); never a comparison population. |
| 180 | CIF **files** fetched, the whole set before subsetting. Never used as a population here. |
| 59 | the curated non-redundant subset (one protein + one RNA chain, no DNA). |
| **58** | what the figure uses — **8JOZ is dropped**, its biotite SASA failed. |

Write "58 non-redundant natural RBP–RNA complexes". Not "58 CIFs" (there are 59 on disk),
not "59". "Complexes" rather than "RBPs": 58 distinct entries is established, 58 distinct
*proteins* only to the 30%-identity clustering.

### ⚠️ The two n's on this page, and which applies where

They are **not** a disagreement — they are a set and a per-metric subset of it, and since
2026-11-20 panel (d) draws both so the page answers the question itself.

| where | n | why |
|---|--:|---|
| the benchmark set, everywhere except one bar | **58** | the 59-complex subset minus 8JOZ |
| panel (a)'s **ΔSASA bar only** | **53/58 (91.4%)** | `dsasa_total` (Rosetta InterfaceAnalyzer) is empty for **2VNU, 2XNR, 3ZJT, 4WAL, 6JVX** — a measurement that is missing, not a selection step, which is why (d) draws it as a subordinate open bar under the final step rather than as a fifth step |
| panel (a)'s other two bars | 58 | marked only where it differs; the panel labels the ΔSASA bar `n = 53` |

A blanket "n = 58" across panel (a) is wrong. The **99,268** that appears beside 58 is the
*design* n — 99,268/100,000 (99.3%) designs with at least one protein interface residue —
and has nothing to do with the natural count; see *The design n* below.

### The design n

⚠️ **The design n is no longer 100,000 in the normalised panel.** Normalising needs a non-zero
denominator, and **732/100,000 (0.7%)** designs have `n_contact_res_prot = 0` — no protein
interface residue at all, so no interface to normalise by. They are dropped by the same
`.where(den > 0)` rule `normalize_metrics.py` uses, leaving **99,268/100,000 (99.3%)**.
`report()`'s H-bond check uses the same subset. Panel (a) is untouched and still quotes
all 100,000. Dropping them does not move the raw medians (designs 1 either way).

All 58 natural complexes have ≥ 1 protein interface residue (minimum 7), asserted in code.

## Numbers, all re-derived on every run

Panel (d)'s funnel — 329 → 180 → 59 → 58, and 53/58 (91.4%) with ΔSASA — is printed by
`report()` with each step's reason, and the counts themselves are asserted in
`make_inputs.py`. Tabulated under *The benchmark set, drawn as (d)* above rather than
repeated here.

Panel (b), raw → per-residue (design ÷ natural): H-bonds **0.07 → 0.11**, base 0.39 → 0.86,
sugar 0.56 → 1.11, heavy-4 Å 0.61 → 1.20, backbone 0.97 → 1.59, phosphate 1.16 → 2.56.

H-bond check, no longer a panel: raw median 14 natural vs 1 design; per interface residue **0.392 vs 0.043**
(ratio 0.11, reproducing b). Design maximum is **7** raw H-bonds against a natural maximum
of 53. **43,378/99,268 (43.7%)** of designs make no interface H-bond at all, and
**5,372/99,268 (5.4%)** reach the natural 5th percentile (0.151 per interface residue) —
the 5th percentile rather than the natural minimum, which is unsteady at n = 58.

Base-contact panel, drawn as (c), per-interface-residue factors (natural ÷ design): PHE **7.4×**, GLN **3.2×**,
HIS **3.0×**, TYR **2.8×**; ARG **0.6×** and LYS **0.7×** inverted. `report()` prints the
raw factor beside each one — PHE 13.7×, GLN 7.0×, TYR 6.6×, HIS 6.0× on the same 99,268
designs, against the stored 13.8× / 7.1× / 6.7× / 6.0× over all 100,000 — so the raw column
is re-derived rather than quoted. Neither is what the panel draws.

## Statistics

scipy is unusable in `~/.venvs/binder_map` (`_propack` fails to dlopen), so everything is
hand-rolled on numpy/pandas — medians, means, quantiles and counts. **No test is run, and
none should be**: with ICC 0.163–0.514 on the H-bond metric, a p-value on a per-complex
H-bond count would claim a precision the measurement does not have. Descriptive only, as
`metric_reliability` concluded for these metrics.

## Caption draft

*Style: lead assertion, then panel-by-panel in order, then the caveat and the n. 233 words.*

**Interface comparison against natural RNA-binding proteins.** One deficit survives a size
correction; the rest were interface size. (**a**) Natural complexes bury a 2.18× larger protein
interface than designs at near-equal RNA breadth (1.10×), so protein interface residues is the
denominator throughout. (**b**) Design ÷ natural median ratio for six extensive metrics, raw
(open) and per protein interface residue (filled), normalised per complex before the median:
the correction removes the base-contact deficit and inverts backbone and phosphate, leaving
hydrogen bonding the only metric outside the natural range. Shape complementarity, packstat and
clash counts are excluded, not normalised: they are model-versus-crystal confounded. (**c**)
Mean base contacts per interface residue, the twelve top natural contributors: Phe (7.4×),
Gln (3.2×), His (3.0×) and Tyr (2.8×) survive the correction, while Arg and Lys invert.
Hydrogen bonds are among the least reproducible metrics here (ICC 0.163–0.514 over five
refolding seeds, against 0.634–0.928 for burial), so (**b**)'s hydrogen-bond row is read as
an order-of-magnitude deficit, not an exact ratio. (**d**) How the benchmark set was
selected: RCSB entries with at least one protein and one RNA entity at ≤ 2.5 Å and
≤ 1,500 deposited polymer residues, one per 30% protein-identity cluster, give 329;
requiring a longest RNA chain of ≥ 10 nt leaves 180/329 (54.7%), a single protein and a
single RNA chain with no DNA leaves 59/329 (17.9%), and interface metrics were computed
for 58/329 (17.6%), with ΔSASA measurable for 53/58 (91.4%) of those. All 58 are X-ray
structures at 1.5–2.5 Å spanning 10–97 nt of RNA and 104–1,492 deposited residues, which
is the heterogeneity (**a**) corrects for. n = 58 natural complexes and
99,268/100,000 (99.3%) designs with at least one protein interface residue; (**a**) uses all
100,000, its ΔSASA bar n = 53.

## Open

- SI numbering, above.
- Check `code/Fig5c_interface_chemistry/` and `code/Fig5def_asn_interface/` for overlap, so
  the SI does not restate a main-text panel. The base-contact panel here (drawn as (c))
  and the Fig 5 Asn work are
  adjacent: this one is the natural-complex comparison, Fig 5's is within the campaign.
  Asn now reads **1.4×** per interface residue here, not 3.7× raw — worth checking against
  what Fig 5 claims before both go in.
- `sum(aa_base_*)` equals `n_contacts_base` for every design but for only **55/58 (94.8%)**
  natural complexes (maximum shortfall 10 contacts), presumably non-standard residues
  outside the 20 columns. It does not affect the ranking; noted so nobody rederives it and
  thinks something is broken.
