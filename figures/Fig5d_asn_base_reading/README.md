# Fig5d_asn_base_reading

Asparagine reading a base pair, under the 2026-11-14 framing that asparagine is
the residue giving the most signal.

Six PNGs: an overview with the asparagines highlighted, and a zoom on the
base contact with the measured distance, for three designs.

## The top three cannot carry this panel

This was the request, and the structures do not support it:

| | FC | Interface asparagines | Reading a base |
|---|---|---|---|
| twist 23 | 1.56 | **0** | 0 |
| twist 15 | 1.54 | **0** | 0 |
| twist 16 | 1.47 | 1 (ASN24) | 0, it reads backbone |

So the panel uses the three binders whose asparagine actually reads a base
**and** whose contact survives a second structure prediction:

| | FC | Model rendered | ligand_RMSD | Contact drawn | Also in |
|---|---|---|---|---|---|
| twist 19 | 1.34 | **design** | 0 by definition | ASN71 OD1 to C5 N4, **3.40 A**, pair C5-G13 | AF3 and ESM |
| twist 12 | 1.34 | **ESM** | **1.77 A** | ASN7 ND2 to G12 O6, **2.92 A**, pair C6-G12 | AF3, at 2.65 A |
| twist 07 | 1.41 | **ESM** | **2.39 A** | ASN54 ND2 to G7 O6, **3.39 A** | AF3, at 3.41 A |

twist 19 is the strongest case: the contact is in the designed model itself, not
only in a refold. twist 21 was considered and dropped: AF3 puts ASN14 on G13 but
ESM puts it on U8 at 1.82 A, a clash, so the two models disagree about what it
touches.

## Model choice: lowest self-consistency RMSD, not AF3

Each panel is rendered on the model that best reproduces the design, by
`dockq_LRMSD`, among those present in the file.

| twist | af3 | b2 | esm | rendered |
|---|---|---|---|---|
| 19 | 8.30 | 1.32 | 2.93 | design |
| 12 | 2.70 | 22.73 | **1.77** | ESM |
| 07 | 3.62 | 6.72 | **2.39** | ESM |

**CORRECTED 2026-11-14.** The first version rendered AF3 for all three. AF3 is
the worst available model for every one of them, and for twist 19 it sits
**8.30 A** from the design, far outside the 3 A the gate asks for. That is the
same mistake `Fig5c_interface_chemistry` documents for the acid-phosphate flag,
made again one figure later. Boltz-2 is closest for twist 19 but is not in
these PDB files; the designed model is used instead, which is better, because
twist 19 is BG+BIF and its design carries sidechains.

The contacts survive the change: all three are present in at least two models,
and the drawn distances move by less than 0.3 A.

## Does asparagine track the signal?

Computed over all 32 against LV2.34.1 flow FC, Spearman with a 20,000-draw
permutation test:

| | rho | p |
|---|---|---|
| `aa_N`, asparagine fraction of the sequence | **+0.356** | - |
| asparagines at the interface (count) | +0.312 | 0.081 |
| asparagines at the interface (fraction) | +0.319 | 0.076 |
| asparagines **reading a base** (count) | **+0.073** | **0.686** |

So the framing is supported at about the same strength as the Csy4-direction
score (+0.346), and `aa_N` is **not** one of the features in that score, so it
is not a restatement of it.

**But the base reading itself does not track the signal.** Asparagine being
present at the interface weakly does; asparagine reaching a base does not at
all. These panels show a mechanism that is plausible and measurable in the
models, not the thing that correlates with the readout. Say "these are the
asparagine-base contacts in the better binders", not "asparagine base reading
drives the signal".

n = 32 and p around 0.08 is not significant, and this is one of 408 features
the earlier scan found nothing in. Treat as a hypothesis for the next round.

## Definitions
- **Interface**: binder residue with any atom within 5.0 A of any RNA atom.
- **Reads a base**: ASN OD1 or ND2 within 3.5 A of a base atom, where base
  atoms exclude sugar and phosphate (so backbone binding does not count).
- Watson-Crick partners from the canonical frame `ACUGCCGUAUAGGCAG`,
  pairs (2,16) (3,15) (4,14) (5,13) (6,12) (7,11). G7-A11 is a sheared G-A,
  not a Watson-Crick pair; the G7 panel is labelled by the base, not the pair.

## What the blue dashes are, and what they are not

They are ChimeraX `distance` lines between two atoms **I named**, drawn in the
theme petrol so they read with the asparagine sticks, labelled with the
measured separation. They are **not** bonds that ChimeraX detected.

I switched from `hbonds` to `distance` because the H-bond test fired on only
one of the three panels. Running that test explicitly (2026-11-14, default and
relaxed `distSlop 0.5 angleSlop 25`, on the model each panel renders):

| Panel | Drawn | ChimeraX `hbonds` verdict |
|---|---|---|
| twist 12 | ASN7 ND2 to G12 O6, 2.92 A | **confirmed H-bond**, plus ND2 to G12 N7 at 3.29 A |
| twist 19 | ASN71 OD1 to C5 N4, 3.40 A | **not an H-bond.** The only one found is C**6** N4 to OD1 at 3.66 A, a different nucleotide |
| twist 07 | ASN54 ND2 to G7 O6, 3.39 A | **no H-bonds at all**, even relaxed |

These models carry no hydrogens, so the test is already the permissive
donor-acceptor form, and twist 07 still fails it.

**So only the twist 12 panel shows a hydrogen bond.** The other two show a
close polar contact that does not satisfy H-bond geometry. The donor/acceptor
chemistry is compatible in all three (ASN OD1 accepts, ND2 donates; cytosine N4
donates, guanine O6 accepts), so these are not nonsense, but the angles are
wrong.

`VERIFY:` before this figure is used, either relabel the dashes as "closest
polar contact" with the H-bond status stated per panel, or redraw twist 19 to
C6 (the contact that does pass) and reconsider whether twist 07 earns a panel.
The current labels invite the reader to see three hydrogen bonds.

## Camera

The zoom camera is chosen, not inherited from the shared frame. For each panel
the scene is turned so the measured contact lies along screen Y, so the drawn
distance is never foreshortened, and then swept about that same axis. Rotating
about the contact axis cannot change the distance, which is why the label reads
identically at every angle: it is a pure legibility choice about what occludes
the asparagine and the base.

Angles in use, picked from the sweep: twist 19 y60, twist 12 y30, twist 07 y120
(`CHOSEN` in the script). To see the options again:

```
python code/Fig5d_asn_base_reading/render.py OUT --sweep --width 1200
```

## Rebuild
```
python code/Fig5d_asn_base_reading/render.py output/Fig5d_asn_base_reading --width 3000
```
