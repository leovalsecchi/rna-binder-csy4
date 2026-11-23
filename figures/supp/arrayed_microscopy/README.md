# Arrayed microscopy — LV2.34, the 32 ordered candidates at 48 h

Imaging counterpart to the LV2.34 flow run. Same plates, same wells, same
48 h timepoint; this figure draws the images, `supp/arrayed_flow` owns the
quantitative panels. Nothing here is quantified — no intensity is measured, no
well is ranked. The fold-change numbers that appear on panel (c) labels are read
straight out of `inputs/lv2_34_1_flow_on_off.csv` to say which well is which.

## Run

```bash
PY=~/.venvs/binder_map/bin/python3
$PY code/supp/arrayed_microscopy/make_inputs.py     # ~17 min, 3 GB off the NAS
$PY code/supp/arrayed_microscopy/make_figure.py     # seconds, repo-only
```

`make_inputs.py` skips itself when its cache is newer than the six plates, so
`run_all.sh` does not re-read the NAS every time; pass `--force` to rebuild.
The whole-well maps are an *input* to field selection, not an output of it,
and they cost the one strided read per well that this script does — the
slowest thing in it, and unusable when the NAS is busy (it fell from 11 MB/s
to 2.9 MB/s mid-build). `--force` reuses them when they are newer than the
plates; `--rebuild-maps` forces them too. Field placement always reads the
*stored* map, so a reusing rebuild places fields identically to a cold one.
`nd2lite.py` is a reader, not a figure — `python3 nd2lite.py <plate>.nd2` runs
its self-check and prints what it found.

## Sources

| file | mtime | what |
|---|---|---|
| `/Volumes/wicky/expdata/Microscopy/LV2_34_{1..6}_Plate{1..6}_HEKwt_32RNAcandidateBinders_48h*.nd2` | 2026-11-06 11:37 – 14:38 | 45 GB, six plates, read-only |
| `/Volumes/wicky/.../LV_Transfection/LV2.33_transfection_AS-PERFORMED_2026-11-04.html` | 2026-11-05 10:43 | plate maps, doses, deviations |
| `/Volumes/wicky/expdata/flow/lvalsecchi/LV2_34_1/LV2_34_1_Analysis.txt` | 2026-11-09 14:10 | FlowJo export; the per-well evidence for plate 6 |
| `inputs/lv2_34_1_flow_on_off.csv` | 2026-11-09 18:47 | fold change per candidate |

Caches written to `inputs/micrographs/arrayed/`: `tiles.npz` (the 160 px field
per well), `wholewell.npz` (153 px whole-well maps), `closeups.npz`,
`wells.csv`, `acquisition.txt`.

### Why there is a hand-written ND2 reader in this folder

No ND2 reader is installed on this Mac — not `nd2`, `nd2reader`, `pims`,
`aicsimageio`, not Fiji, not bioformats — and installing one was not available
during the build. The plates are plain uncompressed ND2, so `nd2lite.py` walks
the container directly. Its traps and its self-check are documented in its
docstring; the one that matters is that chunk-map names carry no terminator.

## Verified

Everything below was read out of the files, not assumed.

* **Acquisition.** Nikon Ti2 widefield, Plan Apo λD 10× OFN25 NA 0.45,
  Photometrics Kinetix, 1×1 binning, 16-bit. **0.649403 µm/px**, recorded in
  `ImageCalibrationLV`, so the scale bars are measured, not estimated.
  Each well is one stitched **7344 × 7344 px = 4.77 × 4.77 mm** scan.
* **Channels**, in component order: `575_mCherry_500ms` (mScarlet),
  `390_mTagBFP` (500 ms), `475_mEGFP_300ms`. Named in `ImageTextInfoLV` and in
  the per-channel exposure chunks; both agree.
* **140 wells**, not the 136 the plan budgets: 24 on each of plates 1–5 and
  **20 on plate 6**.
* **Frame → well** is matched on recorded stage coordinates, agreeing with the
  point list to under 0.5 µm. Not on list order: plate 6 stops after 20 of the
  24 listed points, and index order would have mislabelled 16 wells.
* **Pixel decode** is exact: every `ImageDataSeq` chunk is
  `7344 × 44064 + 8` bytes, i.e. an 8-byte lead and no trailer, with the three
  components interleaved per pixel. Asserted for all 140 frames.
* **All three channels are drawn**, as red / green / blue: mScarlet
  `synbio.reporter_alt` `#e03127`, mNeonGreen `synbio.reporter` `#3fbf4f` —
  both matched to the real emission so these panels and the construct maps
  agree (COLORS.md rule 8) — and mTagBFP on theme step 400 `#65aed8`, because
  the palette has no mTagBFP entry and that is the lightest theme blue that
  still reads as blue against black. No hue is written in the figure script.
* **Every tile is the densest 480 × 480 px field** of its well —
  **311.7 × 311.7 µm** exactly, at the recorded 0.649403 µm/px. Tightened from
  640 px (415.6 µm) so single cells are resolvable at tile size: a 20 µm cell
  is 31 px of 480, which is 0.29 mm across a 4.55 mm array tile and 1.7 mm
  across a 26 mm close-up.
* The field is the window holding the most mTagBFP **edge**, found with a
  summed-area table over the gradient of the whole-well map. The same
  deterministic rule in all 140 wells, and scored on the **transfection
  marker, never on either reporter** — mNeonGreen is a reporter too, so
  scoring on it would let placement reach for a field because of the readout.
  Every well's chosen top-left pixel is in `wells.csv` (`field_y`, `field_x`).
* **A well with no cells gets a fixed field, not an argmax.** The
  untransfected control has nothing to find, so an edge scorer turned loose on
  it picks the brightest flare near the well wall — which is exactly how it
  first came out looking brightly transfected, and twice more after that as
  each guard pushed the field somewhere else. A well is taken to hold no cells
  when fewer than **5%** of its whole-well pixels stand more than 20% above
  their own 5-pixel neighbourhood. Measured: the four untransfected wells are
  the four lowest of all 140 at 0.0252, 0.0307, 0.0309 and 0.0315, and the
  lowest transfected well is 0.0730 — more than double — so 0.05 sits in the
  middle of a wide gap. One map pixel is ~31 µm, about one cell, so a cell is
  point-like at this scale and a flare is not; brightness alone cannot tell
  them apart. Such a well then takes a stated constant position, 1.2 mm due
  east of the well centre, because nothing in it is more representative than
  anything else. All four then sit at 1.08–1.12× their own well median, i.e.
  background.
* **No field may contain unscanned pixels.** A window holding zeros puts a
  black bar in the panel and drags the whole set's display low to 0, which
  washes out every other tile in it — seen once, guarded against since.
* **A field dominated by a bright artefact is rejected**: its mean mTagBFP
  may not exceed **1.4×** the well's own median. Also measured. With the
  exclusion disc but no cap, the four untransfected wells take the four
  highest field/well median ratios of all 140 — 1.78, 2.57, 2.80, 2.88 —
  against a median of 1.22 across the set and a maximum of 1.55 among the 136
  transfected wells. They contain no cells, so the only thing an edge scorer
  can find in them is a flare, and the negative control came out looking
  brightly transfected. At 1.4 every untransfected field lands at or below its
  own well median (0.91, 0.94, 1.00, 1.15); the six transfected wells the cap
  constrains simply take their next-densest field, and no well is left without
  a candidate.
* **The middle of the well is excluded outright**, not dodged by luck: no
  field may be centred within **1.2 mm** of the well centre. The radius is
  measured, not guessed. Over the 136 transfected wells, the probability that
  a point is cell-free (mTagBFP below a quarter of its own well's 90th
  percentile) runs 0.10 at the dead centre, peaks at **0.60 at r = 0.75 mm**,
  and only falls back to its outer-field floor of ~0.29 beyond r = 1.5 mm.
  1.2 mm is where the risk is back within a few points of that floor while
  still leaving the whole 1.2–2.2 mm annulus of a 4.77 mm scan to choose from.
  Checked by eye on eight wells: the excluded disc covers the hole in all
  eight and every chosen field lands in cell texture outside it.
* Scored on edge rather than brightness because several wells carry a
  **saturated artefact** a few hundred microns across — a bubble or a
  fluorescent speck — and a solid bright disc beats any field of cells on total
  brightness. It won in plate 2 B6 and filled that tile with a white blob. A
  disc is flat, so it has edge only at its rim; a field of cells has edge at
  every cell boundary. Checked by eye on six wells: edge lands in cell texture
  in all six, brightness lands on the artefact in four.
* Panels (a) and (b) are the **same field of the same well**. (b) is mScarlet
  alone because the merge cannot carry the result on its own: **mNeonGreen
  marks every transfected cell and saturates**, so a three-colour tile reads
  green whatever mScarlet is doing. The transfection control is the blue
  channel of (a) and the third column of (c). (c) is six of those same fields
  at full resolution, each channel and the merge.
* **One display range per channel per set** (low percentile to 99.9th of the
  pooled pixels), per the rule in `inputs/micrographs/README.md`. Two sets,
  because a tile is 4× block-averaged and a range taken from it would blow out
  the full-resolution close-ups. Counts, 16-bit:

  | channel | tiles (a, b) | close-ups (c) |
  |---|---|---|
  | mScarlet | 14,317 – 65,535 | 13,415 – 65,535 |
  | mNeonGreen | 36,395 – 65,535 | 25,662 – 65,535 |
  | mTagBFP | 26,173 – 65,535 | 21,988 – 65,535 |

  Every high is 65,535, which is the point: **all three channels clip in the
  acquisition itself**, not in the display. mNeonGreen is the worst — it marks
  every transfected cell and saturates across whole fields, which is why the
  merge reads green and panel (b) carries the readout alone.

### Three things the images say that nothing else in the repo does

* **Many wells have a cell-free hole in the middle**, a millimetre or more
  across — plate 4 well C6 is the clearest. The first version of this figure
  took the centre 640 × 640 px of every well and fell straight into it: it
  scored 0.0000 bright-mTagBFP area in wells the flow run puts at 43%
  BFP-positive, and tile brightness correlated with flow BFP-positivity at only
  r = 0.33. **Any future panel that crops these wells has to deal with this**;
  a centred field is not representative. The whole-well maps that show the
  holes are cached in `wholewell.npz`, 153 px per well.
* **mTagBFP is saturated in the acquisition itself** — 8.3% to 25.5% of pixels
  at 65535 across the plate-6 wells measured, with well medians far below it,
  so it is bright cells clipping rather than a flooded field. mTagBFP is usable
  as present/absent and must not be read as an intensity.
* **The untransfected wells carry a few saturating specks in mTagBFP** —
  debris or autofluorescence, not cells: their mTagBFP median is background
  (≈23,000) while a handful of pixels hit 65535, and mScarlet is flat
  (99th percentile 15,800–17,300, i.e. no signal at all). The flow run agrees,
  at 0.0037–0.093% BFP-positive. Worth knowing before anyone reads those specks
  in (b) as transfection.

A whole-well tile was also built and discarded: at 153 px one pixel is 31 µm,
so a block mean averages the sparse mScarlet cells away and a block max turns
every tile into speckle. Neither resolves a cell. The density-chosen field is
the one rule that both misses the hole and still shows cells.

### Where the files contradict the plan

The as-performed record maps plate 6 as 16 wells — V31, V32, Csy4 100 ng,
Csy4 10 ng in columns 1–4. **The data say otherwise, and the data win.**

* The ND2 holds 20 frames at columns **1–5**, rows A–D (stage coordinates).
  The FlowJo export has the same 20 wells, `Plate6_A1`…`Plate6_D5`.
* **Column 3 is an untransfected control the plan never mentions.** In flow it
  is 0.0037–0.093% BFP-positive against 36.3–43.0% in every other plate-6 well, and
  0–0.0088% mScarlet-positive; in the imaging all four of its wells are flat
  background in all three channels. Two independent sources, same answer.
* So plate 6 is V31, V32, untransfected, Csy4, Csy4 across columns 1–5, and
  140 = 136 planned + 4 untransfected.
* `inputs/lv2_34_1_flow_on_off.csv` agrees: it carries `P6_C1`, `P6_C2`,
  `P6_C4`, `P6_C5` and no `P6_C3`.

## NEED

* **`NEED:` which Csy4 column is which dose.** The two doses are 100 ng and
  10 ng (as-performed record). Because the whole block shifted one column
  right, the plan's column numbers cannot be used, and nothing found settles
  the order — not the FlowJo export, not the workbook, not Notion. The panels
  therefore say `Csy4 1` (plate 6 col 4) and `Csy4 2` (col 5) and the caption
  does not assign doses. If order was preserved, col 4 is 100 ng and col 5 is
  10 ng; col 4 does carry the higher mScarlet-positive fraction (22.9% / 21.4%
  on posA / posB against 15.7% / 10.4%), which is consistent with it being the
  higher dose, but consistent is not recorded.
* **`NEED:` the LV2.33 / LV2.34 journal entry.** The Notion journal is
  "LV Wet-Lab Journal (LV01.01 – LV2.18)" and stops at LV2.18; searches for
  LV2.33 and LV2.34 return the thesis build-plan pages, not a bench record.
  The as-performed HTML is the only written account, and it is the one the
  plate-6 layout contradicts.
* `NEED:` the as-performed record marks its own common-mix volume
  "TO CONFIRM against the bench record".
* Not a `NEED:`, but worth one line: the FlowJo filenames say `85ng` where the
  as-performed record computes 84 ng. Same number rounded differently
  (100 ng × 100 µL / 119.6 µL = 83.6 ng), not a second dose.

## Caption draft

> **Arrayed microscopy of the 32 ordered candidates, 48 h after transfection.**
> Each tile is the 311.7 µm field of its well holding the most cells, placed on
> mTagBFP edge, never on either reporter, and never within 1.2 mm of the well
> centre, which is cell-free in many of these wells; one display range per
> channel across all 140 wells.
> **a**, all three channels merged across 140 wells — mScarlet red,
> mNeonGreen green, mTagBFP blue. Columns are candidates V1–V32 and then,
> after the gap, the untransfected and the two Csy4-dose controls; rows are the
> four targets — posA and posB carry the Csy4 scaffold, NegA replaces it, and
> the fourth carries no scaffold. **b**, the same wells in mTagBFP, the
> co-transfected marker, so a difference in **a** that is not also here is not a
> difference in transfection. mNeonGreen marks every transfected cell and
> saturates, so the merge reads green and **b** carries the readout.
> **c**, the same fields at full resolution, each
> channel and the merge, for both Csy4 doses, the
> untransfected control, the highest and lowest fold-change candidates on posA,
> and the highest on the no-scaffold target. Scale bars 100 µm.
> mScarlet reports splicing of the dimerisation reporter — a proxy, not
> affinity; no Kd was measured. One biological replicate, and binder plasmid
> reached the cells at 84 ng against a planned 100 ng, uniform across the
> variant wells. The 32 are a selection on ranked splicing ratio and fold
> change, not a hit rate; 15/32 (47%) reach fold change ≥ 1.2 in the paired
> flow run and are the proposed binders, the other 17/32 (53%) inactive designs.

## Panels

| file | what |
|---|---|
| `panel_a_merge_array` | 4 targets × 36 columns, three-colour merge |
| `panel_b_mscarlet_array` | the same fields, mScarlet alone |
| `panel_c_controls` | 6 close-ups × mScarlet / mNeonGreen / mTagBFP / merge |

Each is emitted as `.pdf`, `.png` and `.svg` (page assembly consumes the SVG)
at 600 dpi, not the house 300: a vector save resamples every imshow raster to
the figure dpi, and a 4.55 mm tile drawn from 160 px would be thrown away at
300 (54 px).

Merge colours come from the palette (`synbio.reporter_alt` for mScarlet,
`MAIN` for mTagBFP); the single-channel images are greyscale, which is what
keeps 140 tiles comparable under one display range and prints unchanged.
No panel letters are drawn — they go on at assembly.
