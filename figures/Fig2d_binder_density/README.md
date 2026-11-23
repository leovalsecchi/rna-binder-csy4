# Fig 2d — binder density on the target, per backbone generator

Four renders of the Csy4 crRNA hairpin, each carrying the density of where one
backbone generator's **design backbones** sit on it. Built 2026-11-18.

    python code/Fig2d_binder_density/prep_inputs.py     # renders -> inputs/binder_density/
    python code/Fig2d_binder_density/contacts_table.py  # contact distributions -> inputs/
    python code/Fig2d_binder_density/make_figure.py     # -> output/Fig2d_binder_density/

## What a panel shows

A surface at the level where **30% of that generator's binders reach**. The
voxel value is the percentage of binders with a backbone atom within 2 Å of
that point, each binder counted once per voxel, on a 1 Å grid — so the contour
level *is* a percentage, not an arbitrary isovalue.

**Four panels, not five.** The inverse folder moves no atoms, so the density is
a property of the backbone generator: BoltzGen+BoltzIF shares BoltzGen's
backbones and is not a separate panel. Its stored design file is also a Boltz
re-prediction rather than the generated backbone.

**800 designs per generator**, equal N. HalluDesign has only 802 in the whole
55,706-design gate, which is what set the number. Equal N matters because the
value is a percentage of that generator's own binders: with 802 in one group
and 27,881 in another, the two scales would be quantised and smoothed
differently and part of the difference would be sample size.

Sampled from `passing_v6_final.csv`: RFD3 v7, BoltzGen v8, BoltzDesign v8,
HalluDesign pooled v3–v9 (it needs every campaign to reach 800).

## Volume enclosed, 1 Å³ voxels

| generator | peak | ≥20% | ≥30% | ≥50% | ≥70% |
|---|---|---|---|---|---|
| RFD3 | 91% | 6,142 | 3,278 | 834 | 294 |
| BoltzGen | 90% | 8,884 | 5,404 | 1,539 | 388 |
| BoltzDesign | 63% | 5,433 | 2,395 | 236 | **0** |
| HalluDesign | 88% | 7,944 | 4,149 | 789 | 210 |

30% is a deliberate mid-range choice: at 20% the shell is 5–9k Å³, the size of
the RNA itself, and hides it; at 50% it is the consensus core and BoltzDesign is
down to two specks. Above ~30% is also where the generators start to differ at
all — at ≥10% all four enclose 14,166–14,656 Å³, which is suspiciously equal and
probably set by binder size and target surface rather than by the method.

## Camera

Identical in all four, so the grey RNA is **pixel-identical across the row** and
only the density differs — a cheap whole-set check that the panels really do
share a camera. The hairpin stands on its own helical axis (PC1 of the 16 core
C1′ atoms), truncated ends at the bottom, loop tip at the top, leaning 25°
toward the viewer.

That lean is measured, not eyeballed: at 0° the hairpin already leans ~4 Å
*away* (tip depth −1.5 vs ends +2.6, because the loop sits off the PC1 axis), so
`turn x 15` merely straightens it and only +25 is a real forward lean.

## Pipeline and where it lives

**Every script is in this folder**, so the panel can be rebuilt from here; the
first three also have a copy on Euler, where they were actually run, in
`rna_binder/analyses/binder_density/` (that is also where the 3,200 harvested
CIFs, the four `.npz` and the four `.mrc` live — they are not in the repo).

| file | stage |
|---|---|
| `harvest_designs.py` | pulls the sampled design CIFs out of the campaign archives (Euler) |
| `extract_binders.py` | stage 1: binders into the target frame → `.npz`. `--self-check` runs without data |
| `render_methods.py` | the four ChimeraX scenes → PNG |
| `prep_inputs.py` | renders → the panel images in `inputs/binder_density/` |
| `make_figure.py` | panel images → `output/Fig2d_binder_density/binder_density_row.svg` |
| `chimerax/` | the four emitted `.cxc`, `thesis_colors.cxc`, and `upright.txt` (the computed upright transform) |

Not here: `binder_density_volume.py`, stage 2, which is Laetitia's and lives on
Euler and on the NAS under `users/llehman/`.

The original three stages on Euler:

| stage | script | what |
|---|---|---|
| harvest | `harvest_designs.py` | pulls the sampled design CIFs out of the nested campaign `tar.zst` archives |
| 1 | `extract_binders.py` | superposes each design's RNA onto the reference, applies the same transform to the binder, writes one `.npz` per generator |
| 2 | `binder_density_volume.py` | Laetitia Lehmann's, unmodified — turns the `.npz` into `.mrc` maps |
| render | `render_methods.py` | the four ChimeraX scenes |

Stage 2 is **not ours**: it is Laetitia's binder-density code from
`users/llehman/`, run unmodified (md5 `0bac577da7e3`). Her stage 1 was not
recoverable, so `extract_binders.py` is a rebuild against the `.npz` contract her
scripts read; it is verified by running *her* `campaign_parts`,
`candidate_references`, `frame_deviation` (0.000 Å against a 0.01 Å tolerance)
and `percentage_field` on its output.

A local copy of stage 1 and the renderer is in `~/binder_density/`.

## Gotchas

* **Crop all four renders to one common box**, not each to its own. They share a
  camera, so a common box keeps the hairpin in register; per-image trimming
  slides each panel by its own density extent. `prep_inputs.py` does this.
* **The renders are not grown to fill the block.** A render is taller than it is
  wide, so four at full width overflow the page: `_page` would shrink everything
  to fit ≈ 0.96 and the type would drop below T. They sit at a fixed `IMG_W`
  (95 pt — the largest that still gives fit = 1.0) with the slack as gaps.
* **Images are embedded as base64**, not linked: the panel body is inlined into
  the page, where a relative href would not resolve. Hence the 1200 px /
  128-colour inputs — four full-size RGBA panels would put 5.5 MB into the page.
* **The labels carry an unstroked `stroke-width`.** The assembly floors every
  stroke and takes `min` over the ones it finds; a panel with no `stroke-width`
  anywhere raises `ValueError` on an empty `min()`.

## The graph under each render — contacts per nucleotide

Mean protein→RNA heavy-atom pairs within 4.5 Å landing on each of the 16
nucleotides of the repeat hairpin, over the designs rendered above. 4.5 Å and
pair-counting to match the other epitope panels (`Fig2d_epitope_funnel`,
`Fig4g_hit_epitope`), not the campaign's 4.0 Å `n_heavy_contacts_4A`. The GUAUA
loop (nt7–11, from `common/csy4.py`) is a pale band behind the bars and all four
share one y scale.

The x axis is the **sequence itself** — `ACUGCCGUAUAGGCAG`, one base letter
under every bar, on every column (2026-11-19, Leo: *"id rather you put the nt
letters under the columns rather than 1-16 or the contacts per nt label"*).
Which nucleotide a bar belongs to is the whole read of this strip, so a scale
drawn once on RFD3 would make a reader count across four panels to use the
fourth. It replaces both the 1/16 ticks and the axis label, so what the bar
heights are gets its own line **under RFD3's sequence**: `mean contacts`, left
aligned with the plot.

It was briefly drawn *inside* that plot at top left, which does clear RFD3's
bars (2.55 pt, computed against the shared y max) but sits on the pale loop band
and reads as a clash — Leo, 2026-11-19. Under the letters it touches nothing.
It would not have worked on any other column anyway: BoltzGen's nt10 is 71 of
the shared max 73 and would run straight through the text.

The caption still has to carry the full definition — "mean protein–RNA
heavy-atom contacts (4.5 Å) per nucleotide, n = 614 / 674 / 729 / 800" — and
the y scale, which is shared across the four columns and tops out at 73
contacts; nothing on the panel gives the reader an absolute number.

### Every structure counted is full-atom, and that took work

Three of the five stacks write **backbone-only** design files (RFD3+LiMPNN,
BG+LiMPNN, BD1+LiMPNN ≈ 4.0 atoms/residue; HD+LiMPNN 7.2 and BG+BIF 7.8 are
full-atom). Counting contacts straight off the design files measures what each
generator wrote to disk, not what it designed. Measured on the union
fingerprint behind Fig 2d/4e, gated means are:

| stack | atoms on file | mean contacts (gate) |
|---|---|---|
| RFD3+LiMPNN | backbone only | 55.2 |
| BG+LiMPNN | backbone only | 52.9 |
| BD1+LiMPNN | backbone only | 82.0 |
| BG+BIF | full atom | 247.9 |
| HD+LiMPNN | full atom | 276.8 |

A 5× split falling exactly on the backbone-only / full-atom line. So the
backbone-only designs are read from the **LigandMPNN side-chain retro-pack**
(`analyses/library_selection/packed_designs/`) — sequence frozen, CA exact,
N/C/O transplanted back — and the rest from the design file itself.
`nt_contacts.py` decides per file by measuring atoms/residue and refuses
anything at or below 4.5, so no stack is hardcoded and a backbone-only
structure cannot slip through.

**Why packed designs and not a refold.** Side chains are a model's guess either
way. What the packed route keeps, and a refold does not, is the design's own
pose — which is exactly what the render above it shows. A refold may move the
whole binder (that is what ligand_RMSD measures), so its contacts would describe a
different placement from the density.

**BoltzGen's column is labelled `BoltzGen + LigandMPNN`**, the stack name folded
into the panel label rather than repeated under the graph: BoltzGen is the one
generator whose backbones feed two stacks, and a contact count belongs to a
sequence, so it has to name the stack it was measured on. The other three have
one stack each, so naming it there would be ink for nothing — but all four
groups ARE the +LigandMPNN stack, which belongs in the caption so the short
labels are not read as "these three are not".

### ⚠️ Coverage is not equal across the groups

The retro-pack was built for the `gate_lrfix` pool, a slightly earlier gate than
`passing_v6_final`, so some sampled designs have no packed structure and the per
group `n` is below 800 for the three backbone-only stacks. It cannot simply be
topped up: `pack_task.sh` resolves its inputs through per-uuid design archives
that the 2026-09-11 campaign repack folded into solid midrun tars, so every path
in `pack_manifest.csv` is now stale. Re-packing would mean rewriting the staging
step against `harvest_designs.py`'s archive reader. `n` is carried into the
caption for this reason.

The subset is **not biased**, which was checked rather than assumed: AF3-refold
contact totals for the designs that do have a packed structure sit within 0.4
contacts of each group's whole gated pool (RFD3 173.6 vs 173.2, BoltzGen 157.1
vs 157.4, BoltzDesign 149.9 vs 149.4, HalluDesign 152.5 vs 152.5). The missing
383 are a coverage accident, not a selection.

`inputs/binder_density/gate_contacts.csv` is NOT the source of this panel. It is
the independent AF3-refold cross-check - total `n_heavy_contacts_4A` per design
for all 55,706 gated designs, from one awk pass over
`cross_campaign/csy4_master_foot.csv` - and it agrees: RFD3's median is 171
against 144-154 for the other three, the same ordering the full-atom design
totals give (359 vs 268-278) from entirely different structures. Its backends
disagree about the other three, though: AF3 says BG > HD > BD1 and ESM says
HD > BD1 > BG, so no ranking below RFD3 is safe on refolds.

## Naming

Labels are the backbone generators: RFD3, BoltzGen, BoltzDesign, HalluDesign.
Panel c uses the stack shorthands (RFD3, BG+LiM, BG+BIF, BD1, HD) and panel b's
legend mixes the two, so the row does not match them. This panel is grouped by
generator rather than by stack, so the long names are the accurate ones here —
change the first element of each `PANELS` tuple to match c if that is preferred.
