"""Csy4 repeat hairpin: sequence and secondary structure, measured not assumed.

The 16-nt target used throughout the campaign is 4AL5 chain B residues 5-20,
confirmed against the design CIFs, which number it 1-16.

Base pairs are the five with an N1/N3 donor-acceptor distance under 3.3 A in
4AL5 itself:

    C6-G20  U7-A19  G8-C18  C9-G17  C10-G16   (4AL5 chain B numbering)
    nt2-16  nt3-15  nt4-14  nt5-13  nt6-12    (target numbering)

so the fold is a 5-bp stem, a 5-nt GUAUA loop at nt7-11, and an unpaired A at
the 5' end. Anything drawn from these constants is the measured topology, not
a cartoon guess.

Shared by the epitope panels in Fig2d_epitope_funnel and Fig4g_hit_epitope.
It lives here rather than in one of them because a figure folder importing
from a sibling figure folder breaks as soon as either is renamed or moved.
"""

#: nt1..nt16, 5' to 3'.
SEQUENCE = "ACUGCCGUAUAGGCAG"

#: Watson-Crick pairs as (5' index, 3' index), 1-based in target numbering.
PAIRS = [(2, 16), (3, 15), (4, 14), (5, 13), (6, 12)]

#: The GUAUA loop.
LOOP = range(7, 12)

#: The unpaired 5' nucleotide.
DANGLE = 1

#: Spans for the stem / loop / stem brackets drawn under a per-nucleotide
#: axis, as 0-based bar positions. The three runs are contiguous, so they are
#: inset slightly or they read as one continuous rule.
BRACKETS = [(0.6, 5.4, "stem"), (5.6, 10.4, "loop"), (10.6, 15.4, "stem")]
