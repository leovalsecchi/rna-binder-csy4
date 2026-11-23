# library_pcr_cycle_gel — the oligo-pool PCR, and the two bands cut from it

One bare panel in `output/library_pcr_cycle_gel/`, **115 × 66 mm**.

```bash
~/.venvs/binder_map/bin/python3 code/supp/library_pcr_cycle_gel/make_figure.py
```

No derived table: the gel is an image and the numbers are the worksheet's and
the journal's.

**Split out of `library_build` on 2026-11-19** (Leo: *"keep B, but rename it as
like library PCR cycle gel, and then improve the labeling for the two bands,
like the extracted and backup"*). It was panel (b) there.

## Sources

| source | what it gives |
|---|---|
| `LV_Gels/gel_LV2.2_oligopool_2026-09-30_inverted.png` | the image |
| `LV_PCR/LV2.2_Oligo_RNA_Binders/LV2.2.xlsx` (2026-09-30) | the reaction, the cycle counts, the lane plan |
| journal LV2.2 (Notion) | the pool, the resuspension |
| journal LV2.4 (Notion) | the two extractions and their yields |

## What the two bands are, and where that is written

Leo asked for the two band labels to tell the **extracted** band from the
**backup** band. Both identities come from journal LV2.4, transcribed in the
retired `library_build` README, now at
`archive/library_build_assembly_scheme_2026-11-20/README_library_build.md`,
which is the only place on disk that carries them:

> **Extraction, journal LV2.4.** Zymoclean. Band 2, the 20-cycle lane,
> **31.3 ng/µL**, taken forward; band 3, the 24-cycle lane, **87.0 ng/µL**, kept
> as a backup; the 16-cycle lane was not extracted.

The **cycle counts themselves** are the worksheet's, not the journal's.
`LV_PCR/LV2.2_Oligo_RNA_Binders/LV2.2.xlsx` lists three reactions off one
master mix — `LV2.1.1 | MM1 | 16x`, `LV2.1.2 | MM1 | 20x`, `LV2.1.3 | MM1 |
24x` — and states the manoeuvre:

> *"Made a single MasterMix, aliquoted in 3 physically separated tubes. Then
> used 1 thermocycler machcine and removed the tube at the wanted cycle by
> pausing the PCR machine during the last 5second of elongation of that
> cycle"*

So the three lanes differ in cycle number and in nothing else. The same sheet
gives TWIST oligopool 0006 (2,000 × 300 nt) at 5 ng/µL, **2 ng of input per
reaction**, KAPA HiFi HotStart ReadyMix, PR21 (`ATTCAAGGGTTGGACGACTC`) and
PR33 (`AAGGCCCAGAAGGATACAAC`), **302 bp expected**, 1.8% gel at 130 V for
120 min.

`VERIFY:` that sheet's own lane plan does **not** match the capture. It lists
the three samples in lanes 1, 3 and 5 across *two* gels with the ladder in
lane 9; the image is one gel with the ladder leftmost and the three products
in ascending cycle order to its right, which is how the panel labels them.
The image is taken as the record of what was loaded. The lane x positions were
measured off the raw capture by eye.

## The labels, and what changed

Old panel (b) wrote `extracted, 31.3 ng/µL` and `backup, 87.0 ng/µL` on two
leader lines that ran out of the bottom-right corner of the gel. Neither said
which lane it came from or what "backup" was a backup *of*, and the second one
sat on the gel's own edge. Now:

| band | label |
|---|---|
| 20-cycle lane | **Excised, taken forward** · 20 cycles · 31.3 ng/µL |
| 24-cycle lane | **Excised, kept as backup** · 24 cycles · 87.0 ng/µL |
| 16-cycle lane | `not excised`, on the lane label itself |

Three things make them read:

* The labels sit in a **column to the right of the gel**, not under it, so no
  leader crosses the image lengthways and nothing is written over the gel.
* **Taken-forward is above backup**, which is the one order in which the two
  leaders cannot cross: both leave their band going up-and-right in lane order.
  With backup on top (the first version), they crossed in the middle of the
  gel.
* Each leader stops **beside** its label block, between the two lines. Run to
  the text anchor instead and it strikes through whichever line it reaches —
  which the first render did, through `20 cycles · 31.3 ng/µL`.

The third lane is named rather than left unexplained: without `not excised`,
the 16-cycle lane is the only one with no annotation and nothing says whether
that is an omission.

## Deliberate

* The gel is **flipped vertically** so the wells sit at the top and migration
  runs downward. That is the only change to the pixels; no gamma, no contrast
  work. Grey levels are the crop's own range.
* **No on-canvas note.** Panel (b) carried `TWIST pool 0006, 2 ng input, 302 bp
  expected; wells at the top, image flipped vertically` along the bottom. That
  is caption material (STYLE_NOTES §0) and is now in the caption.
* **No size calibration.** The ladder lane is shown but not read: the one claim
  the panel makes is which band was cut, and the expected 302 bp is in the
  caption. Interpolating sizes would be a second, unneeded claim.
* `VERIFY:` **the pale rectangles above the two strong bands are not cut
  marks.** They look like excision holes, and the first reading of this panel
  took them for some. They are not: the ladder lane and the un-excised
  16-cycle lane carry the same lighter-than-background band in the same rows
  (0.65 and 0.61 against 0.50 for the open gel, inverted image), so it is a
  lane artefact and nothing about it is drawn or claimed.

## Still open

* `NEED:` **colony counts.** Carried over from the retired `library_build`
  README, because `code/supp/pooled_readout/` cites it and it is library
  construction. Journal LV1.40 records the intervention (5 µL DNA, 100 µL
  cells, 100 µL KCM, 400 µL SOC, plate everything) and Leo's reason (*"for
  some plasmids like pLV004/008 etc. the colony count was super low at
  times"*), but no number. LV2.3 then blames the pLV008 failure on *"low
  colony count led to faulty cell being picked twice"*, so the count is the
  evidence behind that conclusion and it is missing. LV1.44 records the
  opposite failure: *"Colonies overgrew and plate had to be discarded"*.
* `NEED:` **achieved library complexity.** Methods 2.10.1 says "approximately
  12,000 members" while 2.12.4 and 2.12.6 compute from **24,682**. LV2.2
  records only the resuspension — 256 ng at 5 ng/µL into 51.2 µL — and no file
  on the NAS settles it. The panel states no complexity.

## Caption draft

**One master mix, three cycle counts, and the band that became the library.**
The TWIST oligo pool 0006 (2,000 × 300 nt) amplified from 2 ng of input with
KAPA HiFi HotStart and primers PR21 + PR33, 302 bp expected, on a 1.8% gel run
at 130 V for 120 min. A single master mix was split into three tubes and each
was lifted out of the cycler at its own cycle number, so the lanes differ in
cycle count and in nothing else: 16 cycles gives no visible product, 20 a
single clean band, 24 the same band over-amplified, with a heavier species
above it. The 20-cycle band was gel-extracted (Zymoclean) at 31.3 ng/µL and
carried into the BsmBI oligo-pool assembly; the 24-cycle band was excised as
well and kept as a backup at 87.0 ng/µL, and was not used. The 16-cycle lane
was not cut. Sizes are not read off the ladder, which is shown for orientation
only. The image is flipped vertically so the wells sit at the top and migration
runs downward; nothing else is changed. n = 3 reactions, one per cycle count,
from one master mix on one gel; 2 of 3 lanes (67%) were excised.
