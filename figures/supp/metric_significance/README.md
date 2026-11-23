# metric_significance

*Built under its slug. Supplementary numbers are assigned in one pass once the main
text stops moving, so no S-number appears in any file name.*

Whether the metrics that drove selection say anything about which designs the two
screens picked out. **Both panels are built.** Panel (a) was added 2026-11-19 as
precision-recall, on Leo's call (*"can you do precision recall curves? or else auc is
fine."*), then **rebuilt twice on 2026-11-20** — first on new populations, which cost it
the PR curve, then on Leo's explicit metric list and per-design backend rule, which
replaced the metric set entirely. See *Panel (a)* below.

```bash
# once, to build the two per-re-predictor input tables (~6 min the first time, then seconds)
~/.venvs/binder_map/bin/python3 code/supp/metric_significance/export_data.py
# the figure
~/.venvs/binder_map/bin/python3 code/supp/metric_significance/make_figure.py
```

Emits three panels into `output/metric_significance/`:

| file | size | what it is |
|---|---|---|
| `panel_a_auc_arrayed_vs_nonhits.{pdf,png,svg}` | 178 × 96 mm | the forest: rank AUC per metric, 17 rows, two contrasts |
| `panel_b_truncation_pool_vs_panel.{pdf,png,svg}` | 178 × 42 mm | the truncation ECDFs, unchanged since 2026-11-19 |
| `panel_c_precision_recall.{pdf,png,svg}` | 178 × 48 mm | four PR curves on contrast A, restored 2026-11-20 |

Panel letters are assigned at assembly (STYLE_NOTES §11); the stems carry the intended
reading order. ~30 s of compute, but **several minutes of wall time writing nine files to
the NAS mount** — 7% CPU on the last run, all of it I/O.

⚠️ **`panel_a_precision_recall.*` was deleted, not left beside the new file.** Its
populations no longer exist in this folder, and an assembly that picked it up would place
a curve built on 32-against-11,968 under a caption written for 15-against-8,020.

## Panel (b) — the truncation evidence

Each of the five campaign gate terms, distributed in the **gated pool** (55,706
designs, every survivor of the gate) against the **ordered panel** (the 12,000 taken
forward), as empirical cumulative distributions with the gate threshold drawn on the
axis. An ECDF is the form that shows truncation directly: every curve runs into its
own threshold and stops there.

This is the panel that makes the figure publishable rather than a failed analysis,
and it needs no statistical decision — only the two populations.

**The five terms are read out of the gate JSON**
(`code/Fig3e_gpu_efficiency/methods/sources/gate/campaign_gate_55706.json`, the
`"on": true` atoms), with the gate's own combination rule across the two gate
backends: `rmsd` and `sc_DockQ` were **OR** over AF3 and ESMFold2, so the better backend
decides; `footnb_010` was **AND**, so the worse one does. `seq_run` and `iptm_max`
are single-valued. Taking the other end of either consensus would not be the applied
gate, so the script fixes the rule and asserts that no design violates its own term.

| term | gate | pool median | panel median | within 5% of the threshold, pool → panel |
|---|---|---|---|---|
| ipTM (max) | ≥ 0.60 | 0.727 | 0.770 | 4.4% → 2.3% |
| ligand_RMSD (min, Å) | ≤ 3.0 | 2.241 | 2.279 | 6.7% → 7.7% |
| sc_DockQ (max) | ≥ 0.70 | 0.765 | 0.763 | 7.9% → 8.2% |
| seq_run | ≤ 5 | 3.0 | 3.0 | 10.2% → 6.4% |
| footnb_010 (max) | ≤ 5 | 1.0 | 1.0 | 10.3% → 9.4% |

**ipTM is the only term that moves**: median 0.727 → 0.770, and ≥ 0.8 in
8,725/55,706 (15.7%) of the pool against 4,334/12,000 (36.1%) of the panel. That is
the ipTM-ranked top-up inside each cluster showing through. The other four are
visually indistinguishable between the two populations — selection into the panel
added nothing on them.

## ⚠️ Two of the brief's supporting facts do NOT reproduce. Files win.

The brief for this figure carried five established facts to confirm. Three reproduce
exactly; two do not, and they are the two that would have been used to argue a metric
is independent of the gate.

*Every ρ and AUC in this section is measured over the **12,000 ordered panel**, which was
panel (a)'s population until 2026-11-20 and is no longer. The section is kept as written
because it is a record of what did and did not reproduce on that population. The same
correlations recomputed on the new 8,020 negatives are in **How the seven were chosen**,
below: they move by at most 0.030, so nothing here is overturned by the change.*

**Reproduce exactly**

* All five gate terms are 100% satisfied in the panel. (They are also 100% satisfied
  in the pool, necessarily — the pool *is* the gate's output. State it that way:
  saturation is by construction, not a finding.)
* ipTM 0.727 → 0.770; ≥ 0.8 15.7% → 36.1%.
* `dsasa_af3` ρ = 0.401 and `dsasa_polar_af3` ρ = 0.421 against `footnb_010`;
  `footnb_075_af3` ρ = 0.673 and `footnb_100_af3` ρ = 0.651 against the same. Those
  four are **not** independent of the gate and must not be presented as such.

**Do not reproduce**

* **`sc_value_esm` is not cleanly independent.** The brief has |ρ| ≤ 0.084 against
  every gate term. Measured over the panel it is ρ = **0.195** against ipTM
  (−0.025 rmsd, 0.065 sc_DockQ, −0.093 footnb_010, −0.008 seq_run). Still the least
  entangled of the non-gate metrics looked at, but not at 0.084.
* **`sc_value_esm` AUC is not 0.664.** Candidates against the rest of the panel it is
  **0.509**. Four constructions were tried and none gives 0.664: the gate's own
  OR/AND consensus (0.509), the worst-backend consensus (same column, unchanged),
  the 15 proposed binders instead of the 32 (0.573), `sc_value_af3` instead of ESM
  (0.474), and panel membership inside the pool rather than candidacy inside the
  panel (0.535).
* **`seq_max_aa_frac` is ρ = 0.396 with `seq_run`, not 0.608**, and its AUC on the raw
  value is **0.304**, not 0.50. ⚠️ **0.304 and the 0.696 in panel (a)'s table are the same
  measurement**: AUC on the reversed direction is 1 − AUC, and panel (a) scores this
  metric low-is-binder-like, which is the a priori direction. So this is not a metric
  with no signal — on the orientation its own meaning implies it is the only one of the
  seven that separates, and the caveats on it are in panel (a)'s results above. ⚠️ On the
  2026-11-20 populations it is 0.739, but the decomposition there shows that is **the same
  pooled-screen effect**, not a stronger one.

`VERIFY:` whoever produced 0.084 / 0.664 / 0.608 should say which population,
consensus rule and candidate set they used. Until then those three numbers must not
be quoted. Everything above is reproducible with
`inputs/large/gated_pool_v6_gate_metrics.csv.gz` and
`inputs/candidates32_mapped.csv`; the correlations are Spearman ρ hand-rolled on
ranks, and the AUC is the Mann–Whitney statistic over panel ranks — **scipy is
unavailable in `~/.venvs/binder_map`** (its `_propack` extension fails to dlopen).

## Sources, with mtimes

| file | mtime | what is taken |
|---|---|---|
| `inputs/large/gated_pool_v6_gate_metrics.csv.gz` | 2026-11-19 14:42 | 55,706 rows, 16 columns; the gated pool with all five gate terms |
| `inputs/csy4_binder_library_12k_v6_k30_CORRECTED.csv` | 2026-11-17 22:05 | the 12,000 ordered `design_id`s |
| `inputs/candidates32_mapped.csv` | 2026-11-19 | the 32 ordered hits; since 2026-11-20 they are **not** panel (a)'s positives — the file is kept as the assertion target that the parsed Fig 4c gate reproduces exactly |
| `inputs/screen_points.csv` | 2026-11-02 | 21,514 (design, comparison) rows: `reads`, `neg_reads`, `splicing_ratio`, `fold_change_raw`. Panel (a)'s negatives |
| `inputs/lv2_34_1_flow_on_off.csv` | — | the arrayed flow panel: 34 rows = 32 designs + 2 Csy4-H29A controls. `ONOFF_avg` is panel (a)'s positive class at ≥ 1.2 |
| `code/Fig4c_fc_vs_splicing/make_figure.py` | — | `R_GATE, S_GATE, FC_GATE` **parsed from the source**, so the pooled gate cannot drift between the two figures |
| `code/Fig3e_gpu_efficiency/methods/sources/gate/campaign_gate_55706.json` | — | the five `"on": true` gate atoms and their thresholds |
| `inputs/if_interface_aa_12k.csv.gz` | 2026-11-20 | **built here** by `export_data.py`: `if_N`, `if_arom`, `if_neg`, `L` and `n_if` per re-predictor for all 12,000, computed from the local AF3/ESMFold2 CIFs. Validated against `hits32_geom_12k.csv.gz` (`if_arom` to max \|d\| 5e-5, `n_if` exactly) |
| `inputs/metric_significance_backends_12k.csv.gz` | 2026-11-20 | **built here** by `export_data.py`: every other scored metric per re-predictor, from the local 12k master `m32.pkl`, with the BG+BIF true-backbone rescore substituted in |
| `/Users/leo/Downloads/design29_742db4d3/an32/` | — | ⚠️ **outside the repo**: `m32.pkl` (12,000 × 281) and `cif/` (24,032 structures), the same tree Fig 5d-f and Fig 5f-g read. `export_data.py` is the only script here that touches it, and only to write the two `inputs/` tables above |

`gated_pool_v6_gate_metrics.csv.gz` is **new, extracted for this figure**: columns
1, 4, 10, 11, 12, 16, 47, 48, 63, 65, 71, 79, 80, 84, 85, 86, 87 of
`euler:/cluster/work/wicky/lvalsecchi/rna_binder/analyses/foldseek_map/passing_v6_final.csv`
(50,007,971 bytes, 2026-09-01 17:12), whose own header line reads
`# binder_filter export  n=55706  store_n=6992515`. One `awk -F,` pass, skipping the
two comment lines. **This is the only local copy of `footnb_010`**, which the campaign
gate used and which no other file in `inputs/` carries —
`inputs/large/metric_panel_dsasa_foot.csv.gz` holds `foot_*`, a transplant *distance*,
and the two are not substitutable. Re-fetch with:

```bash
ssh euler 'awk -F, "NR>2{print \$1\",\"\$4\",\"\$10\",\"\$11\",\"\$12\",\"\$16\",\"\$47\",\"\$48\",\"\$63\",\"\$65\",\"\$71\",\"\$79\",\"\$80\",\"\$84\",\"\$85\",\"\$86\",\"\$87}" \
  /cluster/work/wicky/lvalsecchi/rna_binder/analyses/foldseek_map/passing_v6_final.csv | gzip -c' \
  > inputs/large/gated_pool_v6_gate_metrics.csv.gz
```

All 12,000 ordered designs are present in the 55,706 — the script asserts it — so the
panel really is a subset of the pool and the two ECDFs are comparable.

## Panel (a) — rebuilt 2026-11-20 on Leo's populations

Leo: *"We need to use the subset of the 12k that have >75 reads and dont have fc >10 and
dont have SR larger than 0.07 vs the proposed candidates in our one by one screen."* Both
sides of panel (a) are redefined by that sentence.

### One gate, shared with Fig 4c, parsed not retyped

`75 / 0.07 / 10` are not new numbers: they are the pooled screen's own hit gate, which
`code/Fig4c_fc_vs_splicing/make_figure.py` carries as `R_GATE, S_GATE, FC_GATE = 75, 0.07,
10.0`. Leo is asking for its **complement** as the background. So `make_figure.py` **parses
those three values out of Fig 4c's source** rather than restating them, and then asserts
that applying them here returns `inputs/candidates32_mapped.csv` **exactly**. The two
figures therefore cannot drift: change the gate in Fig 4c and either this figure follows or
it stops with an assertion.

The read floor applies to **both arms** of a comparison — `reads` and `neg_reads` — because
that NEG cut is what made the ordered hit set 32 rather than 34 (Fig 4c's own comment: two
POS_B hits had 64 and 71 NEG reads).

### `screen_points.csv` has two rows per design, so a design-level rule is stated

21,514 rows, 11,270 distinct designs: 10,244 appear in both comparisons and 1,026 in one.
`POS_A/NEG_A` and `POS_B/NEG_B` are **two backbone contexts carrying the same designs, not
two replicates of a design** (see *The independent unit is the design*, below). The rule,
written down rather than left to whichever row sorts first:

* **Measured** — at least one comparison with `reads ≥ 75` in both its POS arm and its
  matched NEG control. A comparison below the floor is not evidence either way, so it is
  dropped from the design entirely; it does **not** count as a non-hit. 8,895/12,000
  (74.1%) of the ordered panel is measured.
* **Pooled hit** — at least one *measured* comparison with `splicing_ratio ≥ 0.07` **and**
  fold change `≥ 10` in the **same** comparison. That is Fig 4c's per-point gate, and a
  design clearing it in one context is a hit, which is how the 32 were picked.
* **Negative** — measured, and **every** measured comparison has FC ≤ 10 **and** SR ≤ 0.07.

`fold_change_raw` is used, because that is the column Fig 4c gates on. `matching_fold_change`
is the same quantity capped at 50 (65 of 21,514 rows differ) and gives identical membership
at a threshold of 10; the script asserts it, so the choice is not load-bearing.

⚠️ **Leo's wording is stricter than "not a pooled hit", and the difference is 843 designs.**
The literal reading — FC ≤ 10 *and* SR ≤ 0.07 — also drops designs that are high on one axis
and low on the other. The looser complement (simply: not a hit) is **8,863**; the literal one
is **8,020**. The 843 between them are the designs sitting nearest the hit gate, so the
literal set is "comfortably non-hit" rather than merely "non-hit", and removing the
near-misses can only make a metric look better at separating. **The literal reading is what
is built**, because it is what Leo asked for; the looser count is printed on every run so the
cost is visible.

### The two sides, and the n

| | n | from |
|---|---|---|
| **Positives** — proposed binders, arrayed screen | **15** | `lv2_34_1_flow_on_off.csv`, `ONOFF_avg ≥ 1.2` |
| inactive designs (the rest of the arrayed panel) | 17 | same file, `ONOFF_avg < 1.2` |
| the two Csy4-H29A controls | 2 | excluded — they carry no `design_id` and are not designs |
| **Negatives** — measured, not pooled-screen hits | **8,020** | `screen_points.csv`, literal rule above |
| prevalence | **15/8,035 (0.187%)** | was 32/12,000 (0.267%) |

The 1.2 cut is not marginal: the 32 `ONOFF_avg` values run 0.478 … 1.188, **1.223** … 1.563,
so nothing sits within 0.023 of the threshold on either side and the 15/17 split is stable
to a hair's movement of it.

### ⚠️ The overlap, and what it does to the contrast

**Zero overlap, and not by luck.** All 15 proposed binders are pooled hits — they have to be,
the arrayed panel *is* the 32 pooled hits — so FC > 10 excludes every one of them from the
negative set by construction. The script prints this on every run.

The consequence is the thing to say out loud: **the primary contrast stacks two screens.**
It is not "proposed binder vs not"; it is "pooled hit *and* arrayed-confirmed" against
"measured and not a pooled hit". A separation there could come from either screen, and the
panel cannot attribute it. So a **second contrast is drawn beside it**, inside the 32, where
the pooled screen is held constant and only the arrayed readout varies:

* **A** — 15 proposed binders vs 8,020 measured non-hits (pooled **and** arrayed)
* **B** — 15 proposed binders vs 17 inactive designs (arrayed only)

### Precision-recall: not the main panel, but drawn — panel (c)

The first rebuild dropped PR entirely. Leo asked for it back (*"can you do the precision
recall too?"*, then *"so just do precision recall for a then."*), so it is **panel (c)**,
beside the forest rather than instead of it. The reasoning for not making it the main panel
is unchanged and is now what the panel itself demonstrates.

**Why the forest is still the panel that carries the argument.** At 15 positives:

* Recall advances in steps of **1/15 = 6.7%**, so the left-hand end of every curve rests on
  one or two designs.
* The 95% stratified-bootstrap band, measured as upper ÷ lower:

  | metric | at recall 0.1 | at 0.3 | AP | AP 95% |
  |---|---|---|---|---|
  | `seq_max_aa_frac` | **41×** | 27× | 1.14% | 0.36–5.51% |
  | `if_N` | 11× | 8× | 0.26% | 0.18–0.53% |
  | `sc_value` | 11× | 5× | 0.29% | 0.18–0.75% |
  | ipTM | 7× | 4× | 0.16% | 0.14–0.27% |

**The band is drawn at full width and nothing is clipped or smoothed** — the script asserts
the axis contains it. That width is the finding, not a defect: it is why an effect size with
an interval is the honest summary and a curve is not.

**Contrast A only**, on Leo's call. PR exists to show what imbalance does to attainable
precision, and contrast A's prevalence is 0.187%; contrast B is 15 against 17 (46.9%), where
there is no imbalance to correct for and the no-skill line would sit near 0.469. Drawing both
would also repeat, on the no-skill line, exactly the per-contrast misread that the forest's
shaded bounds were fixed for. **Contrast B stays in the forest**, which is where it carries
the `if_N` 0.588 → 0.784 comparison.

**Four metrics, chosen from the forest's verdicts before any curve was drawn** — never from
which curve looks best. The rule, written into `PR_METRICS`:

| metric | why it is in |
|---|---|
| `seq_max_aa_frac` | the only metric that clears contrast A's resolvable bound |
| `if_N` | the metric that leads contrast B, and this project's standing result |
| `sc_value` | the highest-AUC purely structural non-gate metric in contrast A |
| ipTM | the gate term a reader most expects to work, and which does not |

Sixteen curves would be unreadable and most of them sit on the no-skill line anyway.

**The no-skill line is drawn in every cell**, at that cell's prevalence — one prevalence
only, since only contrast A is drawn. A PR curve cannot be read without it.

### The metric set, settled 2026-11-20 (second pass)

The first rebuild scored "the seven that were not gate terms". Leo replaced that with an
explicit list **and a backend rule**, so the mechanical rule is gone and both the choice of
metrics and the choice of model are now decisions that can be written down.

#### The backend rule

> Only score a metric if it exists from **both** AF3 and ESMFold2. Per design, take the
> model with the lower **ligand_RMSD**, and read **every** metric from that one model.

This is a per-**design** model choice, not a per-**metric** best-of-two, and the difference
is a scar in this repo: `Fig2b_map_callouts` picked on best-of-two and drew an AF3 pose 29 Å
off its own call-out. The fix there was an explicit `model` column, and that is what
`metric_table()` carries here. Over contrast A's 8,035 designs the split is **AF3 3,814
(47.5%) / ESMFold2 4,221 (52.5%)** — close enough to even that an AF3-always or ESM-always
panel would be a different panel.

**Boltz-2 is dropped throughout.** It was never a gate backend and it fails the
both-backends requirement for several of these metrics.

Every metric in the list below is present from both backends on all 12,000 ordered designs
(0 missing, asserted), so nothing is dropped for coverage.

⚠️ **`seq_run` and `seq_max_aa_frac` have no backend.** They are sequence metrics; the model
choice simply does not apply to them and they are the same number whichever model is picked.
Flagged `noback` in the code and in the tables.

#### ligand_RMSD stays a scored metric — but read from a FIXED re-predictor

Leo: *"[ligand_RMSD] should still be a valid metric since the selection took place before this set
was created."* He is right about the campaign gate: the 12,000 were selected long before this
analysis, so scoring a gate term here is not circular in that sense.

The circularity measured here is a different one, and it comes from **this figure's own
backend rule**. The model is chosen by minimising ligand_RMSD, so on the selected model the
scored ligand_RMSD is `min(af3, esm)` **by construction**:

| ligand_RMSD read from | contrast A AUC |
|---|---|
| fixed AF3 | **0.458** |
| fixed ESMFold2 | **0.553** |
| the selected model | 0.602 |

**The selected value is higher than either backend alone** — that is the selection
manufacturing signal, not the metric carrying one.

**So ligand_RMSD is kept and scored, but from a fixed backend**, and because the two fixed
values are close and both sit at no effect, **both get an ordinary row** (`ligand_RMSD (Å), AF3`
and `ligand_RMSD (Å), ESMFold2`) rather than a special marker. Their agreement is itself the
answer. Everything else on the panel still comes from the per-design selected model, exactly
as Leo specified; ligand_RMSD is the one exception and the panel says so in a note.

The 0.602 is **kept in the README and in the run log as the record of what the choice
does** — a useful demonstration of the effect, not a result, and not drawn.

⚠️ **This is the established handling in this project, not a one-off.** The recorded note is
*"ligand_RMSD is NOT a predictor — I tested the selected min = circular"*, and
`code/Fig5fg_h1null_refolder/export_data.py` already says the same thing in its own header:
*"Panel a uses a FIXED backend per mark — never the selected refold. Testing the better of
two after selecting on it is circular and is what produced the spurious ligand_RMSD result in
the first version."*

Counting: **17 rows scoring 16 distinct metrics** — `dockq`/`dockq_self` are one column (see
below), and ligand_RMSD is one metric measured on each backend.

#### ⚠️ `dockq` and `dockq_self` are the same column

Leo's list names both. `export_data.py` asserts that `dockq_self_af3 == dockq_af3` and
`dockq_self_esm == dockq_esm` to max |d| = 0.0000 over all 12,000. **17 rows score 16
distinct metrics** — one duplicate removed, and ligand_RMSD counted once but measured on each
backend.

#### Directions, fixed a priori

Still never chosen from whichever way gives the larger AUC; reversing one maps AUC to
1 − AUC, so the direction is part of the result and is drawn as an arrow on every row.

⚠️ **`footnb_010` is drawn with a `?`, not an arrow.** The foot-clash count has **inverted**
as a filter in this project — it has behaved as the opposite of what its meaning implies —
so its direction is **not** assumed. It is scored low-is-binder-like, and the sign of its
result must not be read as confirming that direction.

#### Size confounding, and how it is handled

**Length is not a confound in this population**: every ordered binder is exactly 80 aa, and
`metric_table()` asserts it. **Interface size is.** ΔSASA and an interface H-bond count both
grow with the number of interface residues, so both are divided by the chosen model's own
`n_if` — the repo's existing idiom (`red_total_per_ifaceres_*`). What that costs:

| metric | ρ with `n_if` | raw AUC | per-interface-residue AUC |
|---|---|---|---|
| `dsasa_total` | +0.823 | 0.587 | **0.552** |
| `n_h_bonds_interface` | +0.284 | 0.492 | **0.492** |

ΔSASA is 0.82 correlated with interface size, so most of its raw separation was interface
size; normalising removes it. The H-bond count barely moves. `dsasa_polar / dsasa_total` is
Leo's ratio and is size-free already, so it is scored as given.

#### 🚨 The gate terms are scored now — but they are STILL range-restricted

Leo wants the five gate terms scored rather than reported as "not evaluable". They are. But
the premise for the change — that under the new populations the negatives are no longer the
gate's own output, so the terms are no longer truncated — **does not hold, and the script
measures it rather than assuming either way**:

| term | gate | passes in BOTH arms | gated-pool range | contrast range | range kept |
|---|---|---|---|---|---|
| ipTM | ≥ 0.60 | 8,052/8,052 (100.0%) | 0.600–0.922 | 0.600–0.922 | 100.0% |
| ligand_RMSD (Å) | ≤ 3.0 | 8,052/8,052 (100.0%) | 0.548–3.000 | 0.568–3.000 | 99.2% |
| sc_DockQ | ≥ 0.70 | 8,052/8,052 (100.0%) | 0.700–0.938 | 0.700–0.935 | 98.7% |
| seq_run | ≤ 5 | 8,052/8,052 (100.0%) | 1.000–5.000 | 1.000–5.000 | 100.0% |
| footnb_010 | ≤ 5 | 8,052/8,052 (100.0%) | 0.000–5.000 | 0.000–5.000 | 100.0% |

**Both arms are drawn from the 12,000 ordered designs, which are drawn from the 55,706 gate
survivors.** Redefining the negatives by the screen does not take them outside the gate — it
re-cuts the survivors. So all five terms are 100% satisfied on both sides and each keeps
98.7–100% of the gated pool's already-truncated range.

What changed is therefore the **label, not the statistics**: they are now scored and drawn
in a block of their own, headed *"Selection criteria — range-restricted: every design on
both sides is a gate survivor"*, instead of being excluded as *not evaluable*. The effect
each one carries is the effect **inside the surviving range**, which is not the effect of
the metric. **Panel (b) is unchanged and is still the picture of that restriction.**

#### ⚠️ What this does to "how were the metrics chosen?"

The first rebuild answered this with a correlation table and found the independence problem
intact. Leo's list changes the answer in one direction only:

* **The choice is no longer mechanical.** It is a stated list, which is what Methods needs.
* **Entanglement with the gate is much lower than the old seven's.** On the new negatives
  the largest |ρ| against any gate term among the eleven non-gate metrics is **0.421**
  (`seq_max_aa_frac` vs `seq_run`); the next is 0.180 (`dsasa_per_ifaceres` vs `seq_run`)
  and everything else is ≤ 0.175. The old set's `footnb_075_af3` (0.664) and
  `footnb_100_af3` (0.640) are **not on the new list at all**, and the new structural
  metrics are the cleanest in the folder so far: `pi_stack_score` 0.033, `buns_vsasa_opth`
  0.041, `packstat` 0.078, `if_neg` 0.079, `if_arom` 0.110, `if_N` 0.125.
* **But the five gate terms are now scored beside them**, and those are not independent of
  the gate in any sense — they *are* the gate. The independence question has moved from
  "are these metrics contaminated?" to "the top block is the gate, the bottom block is not",
  which is exactly what the panel's two blocks say.

Full table, recomputed on both populations, in the run log; the maxima are above.

### Results, re-derived on every run

**Contrast A — 15 proposed binders against 8,020 measured non-hits.** No-skill AUC 0.5;
smallest resolvable effect 0.709.

| metric | dir | AUC | AUC 95% boot | AP ÷ no-skill | |
|---|---|---|---|---|---|
| ipTM | ↑ | 0.418 | 0.287–0.544 | 0.87× | gate |
| ligand_RMSD (Å), AF3 | ↓ | 0.458 | 0.308–0.621 | 1.79× | gate, fixed backend |
| ligand_RMSD (Å), ESMFold2 | ↓ | 0.553 | 0.384–0.723 | 1.65× | gate, fixed backend |
| sc_DockQ | ↑ | 0.443 | 0.287–0.609 | 1.14× | gate |
| seq_run | ↓ | 0.499 | 0.350–0.651 | 1.03× | gate |
| footnb_010 | ? | 0.450 | 0.303–0.591 | 0.95× | gate, direction uncertain |
| packstat | ↑ | 0.418 | 0.299–0.531 | 0.84× | |
| sc_value | ↑ | 0.563 | 0.406–0.711 | 1.55× | |
| Buried unsat. (optH) | ↓ | 0.488 | 0.362–0.616 | 1.29× | |
| H-bonds / iface res. | ↑ | 0.492 | 0.350–0.634 | 2.78× | |
| π-stacking | ↑ | 0.481 | 0.443–0.552 | 1.08× | |
| ΔSASA / iface res. | ↑ | 0.552 | 0.374–0.722 | 2.43× | |
| ΔSASA polar fraction | ↑ | 0.396 | 0.285–0.514 | 0.81× | |
| `seq_max_aa_frac` | ↓ | **0.739** | 0.626–0.854 | 6.11× | **clears 0.709** |
| Asn at interface | ↑ | 0.588 | 0.459–0.707 | 1.41× | |
| Aromatic at interface | ↑ | 0.579 | 0.446–0.711 | 1.62× | |
| Acidic at interface | ↓ | 0.432 | 0.303–0.547 | 0.85× | |

**Contrast B — the same 15 against the 17 inactive designs**, where the pooled screen is held
constant and only the arrayed readout varies. Smallest resolvable effect 0.791.

| metric | dir | AUC | AUC 95% boot |
|---|---|---|---|
| ipTM | ↑ | 0.518 | 0.322–0.726 |
| ligand_RMSD (Å), AF3 | ↓ | 0.522 | 0.321–0.737 |
| ligand_RMSD (Å), ESMFold2 | ↓ | 0.620 | 0.408–0.827 |
| sc_DockQ | ↑ | 0.588 | 0.388–0.796 |
| seq_run | ↓ | 0.422 | 0.255–0.616 |
| footnb_010 | ? | 0.459 | 0.263–0.647 |
| packstat | ↑ | 0.437 | 0.247–0.643 |
| sc_value | ↑ | 0.690 | 0.506–0.871 |
| Buried unsat. (optH) | ↓ | 0.398 | 0.186–0.602 |
| H-bonds / iface res. | ↑ | 0.624 | 0.427–0.806 |
| π-stacking | ↑ | 0.478 | 0.382–0.576 |
| ΔSASA / iface res. | ↑ | 0.541 | 0.318–0.765 |
| ΔSASA polar fraction | ↑ | 0.325 | 0.149–0.510 |
| `seq_max_aa_frac` | ↓ | 0.631 | 0.400–0.804 |
| **Asn at interface** | ↑ | **0.784** | 0.621–0.922 |
| Aromatic at interface | ↑ | 0.592 | 0.396–0.784 |
| Acidic at interface | ↓ | 0.573 | 0.375–0.765 |

**Nothing in contrast B clears 0.791.** `if_N` comes closest, at 0.784.

#### ⭐ Asn at the interface is the largest honest effect in the arrayed contrast

`if_N` is 0.588 against the 8,020 measured non-hits and **0.784 against the 17 inactive
designs** — the opposite pattern to `seq_max_aa_frac`, and the one that an arrayed-screen
effect should have: it is largest where the pooled screen is held constant. It is also the
cleanest metric on the list with respect to the gate (|ρ| ≤ 0.125 against every term), it
does not depend on the model-selection variable, and it reproduces, from a different
population and a different contrast, **the one result this project already has** — interface
Asn, recorded in the Fig 5d-f write-up as surviving every stress test on both backends.

Three things still travel with it:

1. **Its interval is 0.621–0.922 and the bound is 0.791**, so it does not clear what 15
   against 17 designs can resolve. It is the best candidate on the panel, not a result.
2. **It is one of sixteen**, and nothing here corrects for that.
3. **The mechanism is not established** — the existing Fig 5d-f note is explicit that the
   presence of interface Asn holds while the base-reading mechanism does not. Nothing in
   this panel adds to the mechanism.

#### ⚠️ `seq_max_aa_frac` is still the pooled screen, not the arrayed one

Unchanged by the new metric set, because it has no backend and is the same number as before:
AUC 0.696 for the 32 ordered hits against 11,968, **0.700** for the same 32 against these
8,020, 0.739 for the 15 against 8,020, and **0.631** in contrast B. Changing the negatives
moves it by 0.004. It is the only metric to clear contrast A's bound and it is the one whose
effect **belongs to the pooled screen**. It is also the most gate-entangled thing on the new
list, ρ = 0.421 with `seq_run`.

#### Metrics the first rebuild scored that are now gone

`sc_value_esm` (replaced by `sc_value` on the chosen model), `dsasa_af3` and
`dsasa_polar_af3` (replaced by the per-interface-residue and ratio forms on the chosen
model), `footnb_075_af3`, `footnb_100_af3` and `seq_kmer` (not on Leo's list).

⚠️ **`seq_kmer`'s reversal is therefore no longer on the panel.** It was the first rebuild's
other finding — AUC 0.333 in contrast A and 0.210 in contrast B, both against an a priori
direction of *less* repetitive. It is not in Leo's metric set, so it is not drawn; the
numbers are kept here and in the run log so the lead is not lost.

### The independent unit is the design

One designed sequence = one library member = one unit. Specifically:

* **Reads are not replicates.** A member's reads are draws from one transfected pool of
  one library; a deeper-read member is measured more precisely, it is not measured more
  times. Treating reads as n would inflate the sample size by three to four orders of
  magnitude.
* **`rep1` and `rep2` are sequencing replicates of the same library**, summed by the
  screen's own report before any ratio is formed. They are not biological replicates.
* **Contexts A and B are two backbones carrying the same designs**, not two biological
  replicates of a design. This is why the design-level rule above collapses a design's two
  comparison rows into one unit instead of scoring them separately.
* **There is one biological replicate** of each screen.

### The two options not taken

**Cliff's delta** is the same number as AUC on a −1…+1 scale — δ = 2·AUC − 1 — so showing it
beside AUC would be the same quantity twice. **Permutation with an FDR correction** yields a
p-value rather than an effect size, and at n₁ = 15 the attainable resolution is coarser still
than it was at 32; it was not run, and no p-value appears in this folder. If it is wanted
later it should be *added to* the effect sizes, never shipped alone.

### Things to watch in the artefact

* **The shaded bound is PER CONTRAST, in half-row strips.** Contrast A resolves at AUC 0.709
  and contrast B at 0.791, so one band drawn at 0.709 would have put `if_N`'s light marker
  (0.784) visibly *outside* the shading while it is *inside its own* contrast's bound — a
  reader would take it as the one metric that cleared, which is the opposite of the result.
  Each marker now sits on a strip drawn at its own bound, and the step in the shaded edge
  always falls between the two markers of a row. ⚠️ **Do not "tidy" this into one rectangle.**
  The honest reading of the result: `seq_max_aa_frac` (dark) is the only mark outside its own
  strip; `if_N` (light) is the closest to leaving its own and does not.
* **The two contrasts are separated by lightness alone** (theme 700 / theme 400), which never
  passes a categorical check on its own (STYLE_NOTES §4, §9). What carries the identity is
  **position** — contrast A is always the upper interval of a row, B always the lower — plus
  the n/N in the key. At 178 mm they do read apart. ⚠️ **Known constraint if the panel is ever
  reused at a smaller size, or in a layout where that row order is not obvious: give B a
  different marker, not a different blue.** The stepped shading helps here too, since the
  wider strip is always the lower one.
* **ligand_RMSD occupies two ordinary rows, not one special one.** It is the only metric read
  from a fixed backend, and the note at the top right says so. There is no hollow marker any
  more; an earlier version used one, with an inline ⚠ that rendered as **tofu** — Arial has no
  warning sign, so no glyph outside Latin-1 goes on these panels.
* **The widest contrast-B intervals span more than 0.4 of the axis** (`Buried unsat.`
  0.186–0.602, `sc_value` 0.506–0.871). That is what 15 against 17 designs buys. It is drawn
  rather than implied, and it is the reason contrast B carries no conclusion of its own.
* **Panel (c)'s y axis stops just below the lowest point any band reaches**, not at the
  decade below it. Flooring to a decade added an empty decade under the data and squashed
  the part that matters; the script asserts no band is clipped.

## Decisions taken on Leo's behalf, 2026-11-20

⚠️ **These three were left open across several turns and were then defaulted. They are NOT
Leo's calls — they were made for him so the figure could be finished, and they are recorded
here as such so any of them can be reversed cleanly if he disagrees.** Each says what would
have to change.

**1. Contrast B stays on the forest.** It is the only contrast that holds the pooled screen
constant, and that is the question this figure exists to answer: without it the panel cannot
tell an arrayed-screen effect from a selection artefact. It is what exposed
`seq_max_aa_frac` as the pooled screen's (0.739 in A, 0.631 in B, moving 0.004 when the
negatives change) and what gives `if_N` its 0.588 → 0.784. ⚠️ Its intervals do span more
than 0.4 of the AUC axis — `Buried unsat.` 0.186–0.602, `sc_value` 0.506–0.871 — and the
caption says so rather than leaving a reader to find out. *To reverse:* drop `B` from
`panel_a()` and the key; the shaded strips go back to one rectangle at 0.709, and the
`seq_max_aa_frac` and `if_N` decompositions lose their evidence and must come out of the
caption too.

**2. The negative rule stays at 8,020** — the literal reading, FC ≤ 10 **and** SR ≤ 0.07 in
every measured comparison. It is what Leo actually specified, and the 843 designs between it
and the looser 8,863 are high on one axis and low on the other, so they are neither hits nor
clean non-hits. **Sensitivity check, now run on every execution and asserted:** rebuilding
contrast A on the looser 8,863 moves no metric by more than **|ΔAUC| = 0.003** (the largest
is `dsasa_polar_frac`), leaves the resolvable bound at **0.709**, and leaves
**`seq_max_aa_frac` as the only metric clearing it**. The rule is a presentation choice, not
a result. *To reverse:* pass `sets["not_hit"]` instead of `sets["neg"]` in `main()`; no
conclusion in this README changes.

**3. `seq_kmer` stays off the panel.** It was considered and excluded, and this is a
decision, not an oversight. It is not on Leo's metric list; it is gate-entangled at ρ = 0.305
with `seq_run`; and both its values run **against** its a priori direction (0.333 in contrast
A, 0.210 in B). Adding a metric back *after* seeing it score is exactly what makes a result
unpublishable, and the direction-fixed-a-priori rule is load-bearing for this whole panel —
it is what lets an AUC below 0.5 be read as a reversal rather than quietly flipped. The
numbers are kept in *Metrics the first rebuild scored that are now gone*, above, so the lead
is not lost. *To reverse:* it would have to be declared, with its direction, before being
scored again — not restored on the strength of the numbers above.

## Caption draft

*Style: lead assertion, then panel-by-panel, then the n. 290 words. Rewritten 2026-11-20
for the new populations and again for Leo's metric set; it is long because the backend rule,
the fixed-backend ligand_RMSD rows and the range restriction all have to be stated.*

**Neither the criteria that drove selection nor the structural metrics that did not
separate the proposed binders from the designs the pooled screen left behind.** (**a**)
Rank AUC per metric, with 95% intervals over 1,000 stratified bootstrap resamples of
designs; every metric is read from one re-predictor per design, the one with the lower
ligand_RMSD (AF3 for 3,814/8,035, 47.5%). Dark, the 15 proposed binders (arrayed ON/OFF ≥ 1.2)
against the 8,020 designs measured in the pooled screen that are not its hits (reads ≥ 75
in both arms, fold change ≤ 10, splicing ratio ≤ 0.07 — the complement of the gate in
**Fig 4c**); light, the same 15 against the 17 inactive designs of the same arrayed panel,
whose intervals span more than 0.4 of the axis at 15 against 17 designs; dashed, no skill (AUC 0.5, prevalence 15/8,035 = 0.187%); arrow, the direction scored as
binder-like, fixed a priori, and `?` where it is not assumed; shading, effects below what
each contrast's own group sizes resolve (AUC 0.709 dark, 0.791 light, the step in the shaded
edge falling between the two marks of a row). ligand_RMSD is the one metric read from a fixed
re-predictor rather than the selected one, on each backend, because it is the variable the
choice is made on. Above the rule, the five campaign gate terms: these are scored but range-restricted, because both arms are gate survivors and each
term keeps 98.7–100% of the gated pool's own truncated range. Precision-recall is not drawn:
at 15 positives its bootstrap band is 41× wide at recall 0.1. Only `seq_max_aa_frac` clears
a resolvable effect (AUC 0.739), and it is 0.700 for the 32 pooled hits against the same
negatives, so it is the pooled screen's effect and not the arrayed one's; interface Asn is
the largest effect where the pooled screen is held constant (0.784, interval 0.621–0.922
against a bound of 0.791). (**c**) Precision against recall for four of those metrics on the
first contrast, chosen from (**a**)'s verdicts: step, the empirical curve; band, 95% over
1,000 stratified bootstrap resamples of designs; dashed, the no-skill rate 15/8,035
(0.187%). With 15 positives recall advances in steps of 1/15 = 6.7%, so the left of every
curve rests on one or two designs and the band is 41× wide at recall 0.1 for
`seq_max_aa_frac` and 7–11× for the others — which is why (**a**) and not (**c**) carries
the argument. (**b**) Why the gate terms are range-restricted: cumulative
distributions of each in the gated pool (grey) and the ordered panel (blue), threshold
dashed — every curve runs into its own threshold and stops. Only ipTM shifts, median 0.727
to 0.770 and 8,725/55,706 (15.7%) to 4,334/12,000 (36.1%) at or above 0.8. n = 55,706 gated
designs, 12,000 ordered from them, 8,895 measured, 32 pooled hits and 15 proposed binders;
one biological replicate of each screen, so no interval covers experimental
reproducibility.

## Open

* ~~`NEED:` the statistical approach for panel (a)~~ — **chosen and built 2026-11-19**:
  precision-recall with average precision and a stratified bootstrap band, rank AUC
  carried beside it, no test run. **Superseded 2026-11-20**: at 15 positives PR does not
  fit and the panel is the rank AUC with its interval, two contrasts, still no test.
* ~~`VERIFY:` the rule that selected the 32~~ — **moot for panel (a) since 2026-11-20**,
  which no longer uses the 32 as its positive class. It still matters for **Fig 4ef** and
  for any sentence that calls the 32 "the ordered hits", so the rule should still be in
  Methods; panel (a)'s own classes are now both defined by thresholds that are written down
  (Fig 4c's gate, and ON/OFF ≥ 1.2).
* `VERIFY:` the 0.084 / 0.664 / 0.608 figures above, which do not reproduce.
* ~~🚨 `DECIDE:` ligand_RMSD cannot be both the model selector and a scored metric~~ —
  **settled 2026-11-20.** Leo: *"[ligand_RMSD] should still be a valid metric since the selection
  took place before this set was created."* It is kept and scored, but **fixed-backend**, on
  each backend as its own row; the selected-model 0.602 stays in the record as the
  demonstration of what the choice does. See *ligand_RMSD stays a scored metric* above.
* `DECIDE:` the five gate terms are drawn as range-restricted rather than excluded, which
  is what Leo asked for — but they are still 100% satisfied on both sides. If the figure is
  meant to say "the selection criteria carry no information", that sentence is not
  supportable from this panel; what is supportable is "not within the surviving range".
* ~~`DECIDE:` whether contrast B stays~~ — **defaulted 2026-11-20, not Leo's call.** It
  stays; see *Decisions taken on Leo's behalf*.
* `NEED:` a citing sentence in Section 5, at the retrospective metric analysis.
* ~~`NEED:` whether precision-recall comes back~~ — **settled 2026-11-20.** Leo: *"can you
  do the precision recall too?"* and *"so just do precision recall for a then."* Drawn as
  panel (c), contrast A only, four metrics, full-width band.
* ~~`NEED:` whether `seq_kmer` should be back~~ — **defaulted 2026-11-20, not Leo's call.**
  It stays off; see *Decisions taken on Leo's behalf*.
* ~~`NEED:` Leo's call on the negative rule~~ — **defaulted 2026-11-20, not Leo's call.**
  8,020, the literal reading; the sensitivity check against 8,863 runs on every execution and
  nothing moves (max |ΔAUC| 0.003). See *Decisions taken on Leo's behalf*.

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
