# missing_gel_records — the two readout gels whose captures are not archived

Two bare panels in `output/missing_gel_records/`, **85 × 46 mm** each
(`csy4_retest`, `double_band`).

```bash
~/.venvs/binder_map/bin/python3 code/supp/missing_gel_records/make_figure.py
```

No image and no derived table. Each panel is a dashed empty frame carrying the
written record that survives the missing picture.

**Split out of `library_build` on 2026-11-19** (Leo: *"and then the others can
be separate figures as well"*). These were the two right-hand frames of panel
(c) there, where they said only *"NEED: image not on the NAS"*.

## Why this is a figure at all

Two gels were run and read, and their captures are not on the NAS.
`LV_Gels/` holds **4 captures in total** — the 2026-07-15 primer test, the
LV2.2 oligo pool, and the two hairpin-test gels — and neither of these is
among them. The gap is drawn rather than cropped out, because one of the two
is the gel carrying the **unexplained double band** that was nonetheless
extracted and carried forward, which makes its absence a fact about the
library and not only about the filing.

## Sources

| source | what it gives |
|---|---|
| `LV_PCR/LV1.27.xlsx` | panel (a): five reactions, the primer pair, the cycle number, the lanes |
| journal LV2.18, 13.10 (Notion) | panel (b): the gel conditions, the double band, the eight extractions |
| journal LV2.11, LV2.9 (Notion) | the yields quoted with it |
| `LV_Gels/` | the four captures that *do* exist, which is how the absence is established |

## Verified

**Panel (a), LV1.27.** `LV1.27.xlsx` is headed *"Splicing Assay integrated
Csy4 (H29A) Scaffold for Csy4 retest"* and lists five reactions
(`LV1.27.1_PCR` to `LV1.27.5_PCR`) off one master mix: KAPA HiFi, 50 µL,
`PR10148` + `PR10137`, `~850 bp` expected, `> step 2, 24x`, in lanes 1 to 4
and 10. The same amplicon was later sequenced — `LV_Sequencing/LV1.27.1_ONT/`
— so the reaction is on record twice over; only its gel is missing.

**Panel (b), LV2.18.** The journal entry (13.10) reads:

> *"Both cycle counts look good visually. Some kind of double band
> (splicing?). Extracted both bands for all 8 samples."*

with 1.25% agarose, 250 mL TAE, 25 µL SYBR Safe. The eight samples are pLV009
to pLV012 at 24 cycles and at 29. Both bands were cut for **8/8 samples
(100%)**. The 24-cycle extractions gave 87.3, 61.9, 84.0 and 67.5 ng/µL, each
annotated *"Concentration good but purity not optimal"*; the 29-cycle set was
not measured.

`VERIFY:` the cause of the double band. The journal's own guess, *splicing?*,
is the only one on record and nothing tested it.

## Deliberate

* **A dashed frame, at the size the image would have had.** An omitted panel
  reads as a figure that was never planned; an empty frame reads as a record
  that is missing, which is the true state.
* **The frames carry the surviving record, not just "not on the NAS".** That
  is what changed in this split: the two frames used to be empty but for the
  `NEED:`, which is one word of information per A4 panel. What is written down
  about each gel is drawn instead, so the panel is worth its page.
* **Muted type inside, ink type for the gel name.** The frame contents are a
  record, not a measurement; only the identification of which gel is missing
  is set in ink.
* No count, ratio or yield here is plotted — they are quoted as the journal
  wrote them, including the four 24-cycle yields, which have no 29-cycle
  counterpart to compare against.

## Still open

* `NEED:` **both captures.** If either turns up, it replaces its frame
  directly: the panel size is the image's.
* `NEED:` **what the second band was.** Sequencing of the LV2.18 extractions
  would settle it; nothing on the NAS does.

## Caption draft

**Two gels in the construction record have no surviving capture, and one of
them is the one with the unexplained band.** Each panel stands in for a
missing image at the size that image would have occupied, and carries what the
written record still gives. (a) LV1.27, the Csy4 retest readout PCR: five
reactions from one master mix, KAPA HiFi in 50 µL with primers PR10148 +
PR10137, ~850 bp expected, 24 cycles, loaded in lanes 1 to 4 and 10. The same
amplicon was sequenced afterwards, so only its gel is lost. (b) LV2.18, run on
a 1.25% gel in 250 mL TAE with 25 µL SYBR Safe: eight samples, pLV009 to
pLV012 at 24 and at 29 cycles, which the journal records as looking good at
both cycle counts but carrying *"some kind of double band (splicing?)"*. Both
bands were excised for 8 of 8 samples (8/8, 100%); the 24-cycle extractions
gave 87.3, 61.9, 84.0 and 67.5 ng/µL, all annotated as good in concentration
but not in purity, and the 29-cycle set was not quantified. `LV_Gels/` holds
four captures and neither of these two is among them. n = 0 images; 5 recorded
reactions in (a) and 8 in (b).
