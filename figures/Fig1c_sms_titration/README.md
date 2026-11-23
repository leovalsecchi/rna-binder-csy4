# Fig1c_sms_titration

Csy4-H29A-SMS scaffold titration, three rows on one column grid: the scaffold
cartoon, the mScarlet micrograph, and the normalized mScarlet signal. Output:
`output/Fig1c_sms_titration/sms_titration.{png,pdf,svg}`.

## Files
| file | what it does |
|---|---|
| `make_panel.py` | draws the panel, reading only `inputs/` |
| `make_inputs.py` | rebuilds those inputs from Matteo's files on the share |

`make_panel.py` never touches his tree - that is the rule in `common/cartoon.py`,
his working copies are revised continuously. `make_inputs.py` is the bridge and
the record of where the inputs came from: run it with no arguments to check the
inputs still match the share, with `--write` to regenerate them. It asserts his
`FoldChange` column still reads 1 at the negative and 135.65 at x8, so if he
ever reorders the plate the script fails instead of silently relabelling.

**The figure number is a guess**, as in fig01b. It is one `paths.fig_output`
string in `make_panel.py`.

## Columns
`Csy4-H29A only` (no scaffold plasmid), `Neg Ctrl` (F30-2xNegCtrlLoop),
`1x`, `3x`, `8x` ((F30-2xCsy4Hairpin)xN). One F30 unit carries two Csy4 sites,
so 8x is 16 sites - the same construct fig04 calls the pre-existing F30 array.

The columns are on a regular pitch, each micrograph in its own axes, so a
cartoon and a bar always sit over their own image.

## Micrographs
`inputs/micrographs/sms_titration/`, five panels cut by `make_inputs.py` from
Matteo's own eight-panel montage
`mScarletZoomIn_xy09_10_11_12_13_14_15_16.png`
(`Matteo_Lampis_Microscopy/MLS11.104.1/MLS11_104_1_Plate3_HEK293H_Csy4MultipleHairpins_Part3/Montages/`).
xy09-16 are the 10 ng wells B3-C4 in the same order as the flow table, so
xy09 / xy10 / xy11 / xy13 / xy16 are the five columns here.

Taking them from his montage rather than from the raw TIFFs is deliberate: the
montage already carries **one display range across all eight panels**, which is
the thing that matters. Scaling each panel to its own maximum would put the
titration in the contrast instead of in the cells. The raw 16-bit stitched
wells are beside it in `mScarlet_10x/` (93 MP each) if the panels ever need
re-exporting; the same shared range has to be applied by hand there.

## Cartoons
`common.scaffold_glyphs.array`, the glyph fig04 uses for the F30 variation
("n Y units on one strand, two hairpins each"), with n = 1 / 1 / 3 / 8. No new
art was needed. Two choices worth keeping:

- **One scale for every column.** The 1x scaffold is drawn small because it is
  short. Drawing each column to fill its slot would make eight units look the
  same size as one, which is the thing the panel is about.
- **The control slot is empty on purpose.** `Csy4-H29A-Only` has no scaffold
  plasmid, so there is nothing to draw. The absence is the information.

Colour (COLORS.md): the scaffold is RNA, so `synbio.rna` rather than the
neutral fig04 gives F30 - there F30 is a reference against designed geometries,
here it is the subject. Hairpin tips keep fig04's convention, `DIVERGING[-1]`
for a Csy4 site and `NEUTRAL["500"]` for a neutral loop. The bars are
`synbio.reporter_alt`, the mScarlet colour, so they read as the channel imaged.

## Numbers
`inputs/sms_titration_flow.csv`, the 10 ng group of Matteo's **E.MLS11.104.1**
(`Csy4mRNAsensingMultipleHairpins`), wells B3-C4, from
`Matteo_Lampis_FlowCytometry/E.MLS11.104.1/E_MLS11_104_1_Analysis.xlsx`.
`norm_mscarlet` is his own `Norm mScarlet` column, copied unchanged by
`make_inputs.py`.

| column | well | norm. mScarlet | his fold change |
|---|---|---|---|
| Csy4-H29A only | B3 | 0.001646 | - |
| Neg Ctrl | B4 | 0.001770 | 1 |
| 1x | B5 | 0.003714 | 2.098 |
| 3x | C1 | 0.222802 | 125.868 |
| 8x | C4 | 0.240124 | 135.654 |

The CSV keeps all eight constructs at both doses. Construct names are not in his
workbook - only well ids - so they come from the x axis of his chart, in order.
His `FoldChange` column confirms the mapping: it is 1 at `F30-2xNegCtrlLoop` and
2.098 / 74.97 / 135.65 at the hairpin, x2 and x8, which are exactly the 2x, 75x
and 135x brackets he drew.

## Open points
- **No scale bar, by decision** (Leo, 2026-11-11). The montage is a crop of a
  10x field at a zoom that is not recorded in the file, so the micron-per-pixel
  would have to come from the acquisition settings. Do not add one from a guess;
  if it is wanted later, get the value from the acquisition metadata first.
- The montage panels are 376 x 348 px, which is close to 1:1 at 28 mm and
  300 dpi. They will not take much more enlargement.
- `Csy4-H29A only` is his first bar in the group. On the slide that tile is
  labelled `10ng`, which is also the dose for the whole strip.
- x2, x4 and x5 were measured but not imaged, so they are in the CSV and not in
  the panel. x4 (0.204) sits below x3 (0.223) at this dose.
- At 100 ng the whole series compresses to a 21x window; 10 ng is the dose that
  separates the constructs.
