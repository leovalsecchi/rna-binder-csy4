# readout_pcr_cycle_test — the readout PCR at 24 and at 29 cycles

Two bare panels in `output/readout_pcr_cycle_test/`, **85 × 92 mm** each
(`cycles24`, `cycles29`).

```bash
~/.venvs/binder_map/bin/python3 code/supp/readout_pcr_cycle_test/make_figure.py
```

No derived table: both panels are images, and the reaction behind each is a
worksheet.

**Split out of `library_build` on 2026-11-19** (Leo: *"separate like C, the
first two gels can be a separate figure"*). These were the two left-hand gels
of panel (c) there.

## Sources

| source | what it gives |
|---|---|
| `LV_Gels/gel_hairpinTests_24cycles_inverted.png` | panel (a) |
| `LV_Gels/gel_hairpinTests_29cycles_day2_inverted.png` | panel (b) |
| `LV_PCR/LV1.14.xlsx` | the 24-cycle run |
| `LV_PCR/LV1.15.xlsx` | the 29-cycle run, same ten samples |

## Verified

The two worksheets are the **same reaction twice**. Both list one master mix —
KAPA HiFi, 50 µL, `PR10148` + `PR10137` at 2.5 µL of 10 µM each, 45 µL per
tube plus 5 µL of template — and the same ten samples, `LV1.13.1_PCR` to
`LV1.13.10_PCR`, in lanes 1 to 10, `~850 bp` expected. The cycling blocks are
identical down to the second (95 °C 3 min; 98 °C 20 s / 65 °C 15 s / 72 °C
1 min; 72 °C 5 min; hold 4 °C) except for the one line that is the experiment:

* `LV1.14.xlsx`: `S4 Elongation | 72˚C | 01:00:00 | > step 2, 24x`
* `LV1.15.xlsx`: `S4 Elongation | 72˚C | 01:00:00 | > step 2, 29x`

So the comparison is cycle number against the 24-cycle run as the baseline,
with template, primers, polymerase and programme held fixed. 29 cycles was
adopted for the readout PCR on the strength of it.

The second capture's file name carries `Day2`: the 29-cycle run was done the
day after, on the same samples.

## What the panels show

At **24 cycles** the product row is barely above background — a faint band in
five lanes of the upper gel and effectively nothing in the lower one. At **29
cycles** the same five lanes carry a strong, single, sharp band, and the lower
gel's lanes now show theirs as well. Same samples, same lanes.

## Deliberate

* Both gels are **flipped vertically** so the wells sit at the top and
  migration runs downward. That is the only change to the pixels. Both
  captures are shown at the same crop and the same scale, which is what makes
  the two readable as one comparison.
* **No contrast work on the 24-cycle panel.** Stretching it to make its faint
  bands visible would destroy the only thing the pair says.
* `24 cycles` / `29 cycles` is drawn on each panel. It is the condition, not a
  title: it is the one quantity that differs between the two panels, so it
  belongs on the canvas (STYLE_NOTES §1) rather than only in the caption.
* **Lanes are not annotated.** `VERIFY:` both sheets list ten samples in lanes
  1 to 10, but each capture holds **two gels, one above the other**, which
  neither sheet describes. Assigning a sample to a lane would be a guess, so
  no lane is named.
* **No on-canvas note.** Panel (c) carried `flipped vertically so migration
  runs downward; each capture holds two gels…` along the bottom; that is
  caption material (STYLE_NOTES §0) and is now in the caption.

## Still open

* `NEED:` **which gel of each capture is which half of the sample set.** See
  above; it is the same gap in both captures.

## Caption draft

**Five cycles more is the difference between no band and a clean one.** The
splicing-readout PCR run twice on the same ten templates (LV1.13.1 to .10),
with one master mix per run — KAPA HiFi, 50 µL, primers PR10148 + PR10137,
~850 bp expected — and identical cycling but for the cycle number. (a) 24
cycles (LV1.14): the product row is barely above background. (b) 29 cycles
(LV1.15), run the next day on the same samples: a single sharp band in the
same lanes, and the second gel's lanes now resolve as well. Everything else
was held fixed against the 24-cycle run, so the difference between the panels
is cycle number alone; 29 cycles was adopted for the readout. Each capture
holds two gels, one above the other, and lanes are deliberately not labelled —
the worksheets list ten samples in lanes 1 to 10 but do not say how the two
gels divide them. Both images are flipped vertically so the wells sit at the
top and migration runs downward, shown at the same crop and scale, with no
contrast adjustment; nothing else is changed. n = 10 reactions per cycle
count, one gel capture each.
