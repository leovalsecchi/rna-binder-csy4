# Fig4d - screen hit rate per stack

    python3 code/Fig4d_stack_hit_rate/make_figure.py  -> output/Fig4d_stack_hit_rate/Fig4d_stack_hit_rate.{svg,pdf,png}

Hit rate = ordered hits (the 32, `inputs/candidates32_mapped.csv`) / designs of
that stack in the 12k library (`inputs/csy4_binder_library_12k_v6_k30.csv`).
Bars sorted by rate, best on top, on Leo's request: this panel deliberately
breaks the fixed stack order (STYLE_NOTES 9). Colour still follows the stack.

| Stack | Hits | Library | Rate |
| --- | --- | --- | --- |
| BG+BIF | 9 | 2,400 | 0.375% |
| BG+LiMPNN | 11 | 3,266 | 0.337% |
| BD1+LiMPNN | 7 | 2,400 | 0.292% |
| RFD3+LiMPNN | 4 | 3,132 | 0.128% |
| HD+LiMPNN | 1 | 802 | 0.125% |

Caption draft: *Screen hit rate by design stack: designs passing the gate
(>= 75 reads in POS and matched NEG, splicing >= 0.07, ON/OFF fold change >= 10)
as a fraction of each stack's designs in the 12,000-design library. Labels give
hits / library designs.*

Panel 4d (Leo, 2026-11-12). The ordered-32 pie and grid moved to 4e and 4f (`Fig4ef_ordered32`).

## Drawn twice, at two sizes (2026-11-13)
`make_figure.py` takes an optional `W_MM H_MM FONT_PT OUTDIR` argv so
`code/Fig4_assembly` can redraw it through `page_assembly.variant()`. It has
to: this panel sits beside the 4:1 epitope-columns strip, which leaves it
about 45 mm. Scaled into that slot it shrinks ~0.6x, and text re-set to
T = 5.5 pt then grows against the bars until "RFD3+LiMPNN" runs off the left
edge of the page. Drawn at the slot size instead, the labels fit.
Same failure and same fix as Fig 2's GPU panel (`_d_variant`).
