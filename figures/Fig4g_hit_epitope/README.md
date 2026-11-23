# Fig4g - the RNA epitope of the 32 hits against the 12k library

Two forms of one comparison. Build either with the venv python:

    PY=~/.venvs/binder_map/bin/python3
    $PY code/Fig4g_hit_epitope/make_hairpin.py   -> output/.../hairpin_footprint.{svg,pdf,png}
    $PY code/Fig4g_hit_epitope/make_columns.py   -> output/.../epitope_columns.{svg,pdf,png}

**Letter is provisional.** Figure 4 already runs a-f (`Fig4d_stack_hit_rate`,
`Fig4ef_ordered32`). These two and `Fig4h_hit_geometry` are the new
hits-against-library content; renumbering is a folder rename plus the deck.

## Campaign epitope moved to Fig2d

The campaign / gate / library version of this panel now lives in
`code/Fig2d_epitope_funnel/` (`make_columns.py` + `make_summary.py`), with its
own source and caveats. It is a different population AND a different
measurement - design-model coordinates, not an AF3 refold - which is why the
library reads U10 = 0.362 there and 0.231 here. Never plot them on one axis.

The shared hairpin topology those scripts and these use is
`code/common/csy4.py`.

## What is drawn

Per nucleotide, the mean share of a design's protein-RNA contacts landing on
it. From `geom.py`: heavy-atom pairs within 4.5 A, each nucleotide's count
divided by that design's total, so one design's 16 values sum to 1. A mean of
those shares is a footprint that does not reward designs for having a big
interface. `make_columns.py` splits each nucleotide into the base part
(`nt{i}_base`) and the rest, which is sugar + phosphate.

A three-way phosphate / sugar / base split **per nucleotide is not in the
data** - geom.py stores those three only summed over the whole hairpin. Those
whole-hairpin numbers are printed for the caption instead of drawn:

| | phosphate | sugar | base |
| --- | --: | --: | --: |
| Library | 0.403 | 0.232 | 0.365 |
| Hits | 0.382 | 0.203 | 0.415 |

## The topology is measured, not assumed

The 16-nt target is 4AL5 chain B residues 5-20, `ACUGCCGUAUAGGCAG`, confirmed
against the design CIFs, which number it 1-16. Base pairs are the five with an
N1/N3 distance under 3.3 A in 4AL5 itself (`C6-G20 U7-A19 G8-C18 C9-G17
C10-G16` in chain B numbering = `nt2-16 nt3-15 nt4-14 nt5-13 nt6-12`), so the
fold is a **5-bp stem, a 5-nt GUAUA loop at nt7-11, and an unpaired A at the
5' end**. Nothing in the drawing is a cartoon guess, and the stem/loop
brackets in the column form come from the same pairing.

## Numbers

One design of the 12,000 has no AF3 protein-RNA contact at all and is dropped,
so the library is n = 11,967 against 32 hits. Carried as NaN it silently
poisons every numpy mean; both scripts drop it explicitly.

- Hottest nucleotide in both sets is **U10**, the fourth loop position:
  0.231 of library contacts, 0.253 of hit contacts.
- The loop (nt7-11) takes **49.7%** of library contacts and **52.5%** of hit
  contacts.
- Largest hit-minus-library difference anywhere is **0.022, at U10**.
- Of 16 nucleotides, the library mean falls outside the hits' bootstrap 95%
  interval at **two** (G4 and U10), uncorrected, against 0.8 expected by
  chance at 16 tests.
- The 32 hits spread over **28 epitope clusters** (of 341 in the library).

So the hits read the same epitope as the library: the loop, U10 above all.

Caption draft: *RNA epitope of the 32 screen hits against the 12,000-design
library. Mean share of each design's protein-RNA heavy-atom contacts (4.5 A)
falling on each nucleotide of the Csy4 repeat hairpin; shares sum to 1 per
design. Secondary structure is the five base pairs resolved in 4AL5 (chain B
5-20). Error bars are bootstrap 95% intervals over the 32 hits; the library
(n = 11,967) is plotted without one. Both sets concentrate on the GUAUA loop
and on U10 in particular, and the hits are not reading a different epitope.*
