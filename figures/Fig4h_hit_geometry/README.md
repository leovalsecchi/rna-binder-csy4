# Fig4h - the geometry of the 32 hits against the 12k library

    PY=~/.venvs/binder_map/bin/python3
    $PY code/Fig4h_hit_geometry/make_shape.py       -> output/.../shape_space.{svg,pdf,png}
    $PY code/Fig4h_hit_geometry/make_percentile.py  -> output/.../percentile_strip.*
                                                       output/.../density_curves.*

Letter provisional, same note as `Fig4g_hit_epitope`.

## shape_space

**a** Elongation against radius of gyration: library as a log-scaled hexbin,
the 32 as dots in stack colour. The view is clipped to the library's 1st-99th
percentile on both axes and all 32 hits fall inside it, so the clip hides
nothing but library tail. Panel b's coloured stack names double as the key
for a.

**Neither axis separates the hits, and the panel must not be captioned as if
it did.** The raw medians differ a little - elongation 2.39 against 2.66, Rg
13.1 A against 13.6 A - but both distributions are long-tailed to the right,
so a median gap of that size is what a 32-design draw from the library looks
like. On the rank scale the hits land dead centre: median within-stack
percentile **51.7 for elongation and 51.4 for Rg**, against a null band of
33.6-66.5. They are the two most featureless metrics of the ten in
`percentile_strip`. Panel a's job is to show the hits sitting inside the
library's cloud, not to claim they sit anywhere particular in it.

**b** Share of designs with any beta-strand. The **pooled** comparison leads,
drawn at full weight (25.0% of the library against 31.2% of the hits); the
five per-stack rows beneath it are recessive, because they are the
explanation rather than the result. Bar is the library, dot is the hits.

## The strand result, in full

Pooled, the hits do look strand-rich: **25.0%** of the library has any strand
against **31.2%** of the hits (10/32), and mean strand content is 0.054
against 0.085. ESM agrees (0.053 against 0.083).

It does not survive looking at it per stack:

| Stack | Library with strand | Hits with strand |
| --- | --: | --: |
| RFD3+LiMPNN | 1.6% (49/3,128) | 0/4 |
| BG+LiMPNN | 35.5% (1,154/3,255) | 5/11 |
| BG+BIF | 39.7% (949/2,391) | 3/9 |
| BD1+LiMPNN | 26.3% (629/2,393) | 2/7 |
| HD+LiMPNN | 26.6% (213/801) | 0/1 |

**RFD3 almost never makes a strand, and RFD3 is the stack most depleted among
the hits** (4 of 32, against 26% of the library). Most of the pooled gap is
that mix shift. Ranking each design within its own stack and permuting the hit
labels inside each stack 20,000 times gives **p = 0.31** (AF3) and **p = 0.22**
(ESM). Panel b is drawn per stack so a reader sees the confound instead of
being told about it. Do not caption this as hits being beta-rich.

## percentile_strip and density_curves

Ten geometry metrics, each hit at its percentile **within its own stack**.
Pooling the ranks would manufacture a signal from the same mix shift. The grey
band is the 2.5-97.5 percentile of the median of 32 uniform draws (33.6-66.5,
40,000 simulations).

Eight of ten hit medians land inside the band. The two that do not are
marginal and are the same weak trends the original 220-feature analysis
already turned up and discarded:

| Metric | Hit median percentile | |
| --- | --: | --- |
| Interface atom pairs | 68.9 | just above the 66.5 band top |
| Helix fraction | 32.5 | just below the 33.6 band floor |

Ten metrics against a 95% band expect 0.5 outside by chance, neither clears BH
correction over the full feature set, and the helix one rides the same stack
mix as strand. Caption them as marginal, not as findings.

**Strand is deliberately not in these two.** 75% of the library has no strand
at all, so its ranks are one large tie and a smooth curve draws that point
mass as a feature - the dead end already recorded in `Fig3d`'s header. It gets
the binned panel in `make_shape.py` instead. For the same reason
`density_curves` switches to a step histogram for any integer metric with 30
or fewer levels; without that, nucleotides-touched combs into a spike per
integer.

scipy in `~/.venvs/binder_map` is broken (dlopen on `_propack`), so the KDE is
six lines of numpy with a Silverman bandwidth.

Caption draft: *Geometry of the 32 screen hits against the 12,000-design
library. (a) Elongation against radius of gyration, AF3 complexes; library as
a density, hits as points coloured by design stack. (b) Share of designs
containing any beta-strand: pooled across the library (top row), then per
design stack; bars are the library, points the hits, labelled hits/total.
Pooled, 31.2% of the hits contain a strand against 25.0% of the library, but
RFD3+LiMPNN produces a strand in only 1.6% of its designs and is the stack
most depleted among the hits, so most of that gap is the change in stack mix;
ranking designs within their own stack and permuting hit labels within stacks
gives p = 0.31. (c) Ten geometry metrics, each hit at its
percentile within its own stack; the band is where the median of 32 randomly
drawn designs falls 95% of the time. The hits are not a geometric subclass of
the library.*
