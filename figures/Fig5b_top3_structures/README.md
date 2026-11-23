# Fig5b_top3_structures

The three LV2.33 candidates with the highest reporter fold change, rendered in
the house ChimeraX style. 15 PNGs: three overview angles and two contact
close-ups per binder.

## READ THIS BEFORE USING THE PANEL

**None of these three is distinguishable from noise.** This is a structure
panel of the three top-ranked candidates, not evidence that they bind.

| | |
|---|---|
| Best fold change in the panel | **1.38x** (V15) |
| Csy4 positive control, same readout | **14.0x** (`l2fc_pct` 3.80) |
| Benjamini-Hochberg q, every candidate | **1.00** |
| Best empirical p | 0.26 (V15, z = 2.12) |
| Candidates with flow data | **20 of 32** |

`lv2_33_stats.py` builds the empirical null from log2(Neg_noScaf / NegA), the
same statistic on a contrast where no biology is expected, and nothing in the
panel clears it. `lv2_33_screen.py` itself asserts
`res.loc[cand, "l2fc_pct"].max() < 1.0`, i.e. that no candidate reaches a
doubling. There is one well per condition, so there is no replicate-based test
available at all.

Per `LV_writing/style/claim-boundaries.md`: the SMS receptor reports binding
indirectly through dimerisation-driven splicing, so the defensible wording is
"produces a reporter response", never "binds" and never "affinity".

**V1 to V12 have no flow data** (not yet re-exported), so "top three" means top
three of the 20 that were scored, not of the 32 ordered.

## Which three, and why

Ranked on what was asked for: mean %ON over the two hairpin targets (posA,
posB) divided by the mean over the two negatives (NegA, Neg_noScaf).

| Rank | Variant | Twist | Design | Generator | FC |
|---|---|---|---|---|---|
| 1 | V15 | 15 | `csy4free_BG_BIF_len80_25a69495_design18_seq0` | BG+BIF | 1.380 |
| 2 | V16 | 16 | `csy4free_RFD3_LiMPNN_len80_73d23e78_design49_seq4` | RFD3+LiMPNN | 1.352 |
| 3 | V26 | 26 | `csy4free_BD1_LiMPNN_len80_d4d5ef6f_design56_seq1` | BD1+LiMPNN | 1.335 |

One from each of three different generators, which is convenient but is not a
result at these effect sizes.

### Three numbering traps, all of which bite
1. **The project's own metric uses NegA alone as the denominator**, not the mean
   of both negatives: `l2fc_pct = log2(pct_hairpin / pct_NegA)`. On that metric
   the order is V15, V26, V16. The top three are the same set either way, only
   2nd and 3rd swap. NegA is the valency-matched negative, which is why the
   pipeline uses it; if the panel is captioned, say which denominator.
2. **CORRECTED 2026-11-13: `index.tsv`'s `flow_FC` is not untraceable.** This
   README previously said its provenance was unknown. It is
   `inputs/lv2_34_1_flow_on_off.csv` `ONOFF_avg` = mean(PosA, PosB) /
   mean(FullNeg, NegA) from **LV2.34.1**, a later and more complete experiment,
   verified identical on all 32 to 4 decimal places. LV2.33 scored only 20 of
   32; LV2.34.1 has all 32, and the handover treats it as ground truth.
   **On LV2.34.1 the top three are twist 23, 15, 16**, so the V26 panel here
   should be twist 23 if this figure is kept. See
   `code/Fig5c_interface_chemistry/README.md`.
3. **`working/structures/hits32/index.tsv` has a column called `ON_OFF_ratio`
   that is the HTS sequencing ratio, not flow.** V15 reads 50.6 there.

## Structures
`LV_BinderDesign/LV_32Candidates/Q1_*.pdb`, three models each: 1 design,
2 AF3 refold, 3 ESMFold2 refold, already superposed on the RNA.

**`working/structures/hits32/*.pdb` carries the same REMARK header claiming
three models but holds only model 1.** A refold read from there is silently the
design. Checked 2026-11-13; use `LV_32Candidates`.

| Panel kind | Model | Why |
|---|---|---|
| overview | design | what was designed |
| contact | AF3 refold | V16 and V26 designs are backbone only, 320 atoms, no CB, so their sidechains exist only in a refold |

V15 is BG+BIF and its design does carry sidechains (642 atoms), but the contact
panels use the AF3 refold for all three so the three are comparable.

## What the contact panels show
Measured from the AF3 model, not chosen by eye: polar contacts at or below
3.5 A, and aromatic / base ring-centroid pairs at or below 5.0 A.

| Panel | Residues | Reads |
|---|---|---|
| V15 hub | ARG67, ARG64, SER63 | G7, U10, G12, C5, A9 |
| V15 stack | HIS60, SER63 | stacks on G7 at 4.74 A |
| V16 hub | ARG47, GLN46, ARG50 | ARG47 makes 11 contacts to A1; GLN46 and ARG50 read C14, A15, G16 |
| V16 loop | HIS27, GLU34 | both reach U10 |
| V26 stack | TRP37, TYR32 | TRP37 on U10 at 4.33 A, TYR32 on U8 at 4.53 A |
| V26 hub | ARG20, ARG18, ARG50 | C6, G7, A9, U10, A11 |

`VERIFY:` several predicted polar contacts are shorter than real hydrogen bonds
(V16 ARG47 to A1 at 1.2 A, V15 HIS60 at 1.36 A). Those are model clashes, not
bonds. The drawn pseudobonds come from ChimeraX `hbonds` with its own geometry
test, which is why fewer dashes appear than the contact counts suggest. Do not
quote the raw contact counts as hydrogen bonds.

## House style
`code/common/CHIMERAX_RENDERING.md`, followed: named colours from
`config/thesis_colors.cxc` (no hand-picked hexes), binder in its stack colour,
RNA a neutral, soft lighting, neutral-700 silhouettes at width 1.5, white
background, transparent PNG, supersample 3, `view` called once per group.
All nine submodels are aligned independently onto one reference RNA (V15's
design) by C1', so the hairpin is in the identical pose in every panel.

The V26 design carries an 18-nt hairpin while its refolds carry the 16-nt core,
the documented gotcha; the range is mapped explicitly per model (design 2-17).

## Rebuild
```
python code/Fig5b_top3_structures/render.py output/Fig5b_top3_structures --width 3000
```
Writes 15 PNGs plus the generated `render_top3.cxc`. Width 1800 for screen.
