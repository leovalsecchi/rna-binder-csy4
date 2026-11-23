# target_definition

What each design route was given, and where those coordinates come from. Two bare
panels, cited from the Introduction.

```bash
# (a) renders first - ChimeraX cannot render offscreen on macOS, so this steals focus
~/.venvs/binder_map/bin/python3 code/supp/target_definition/render_target.py --azimuth 60
~/.venvs/binder_map/bin/python3 code/supp/target_definition/make_figure.py
```

`render_target.py` with no `--azimuth` writes a 0/30/60/90 sheet to pick from; 60 is
the chosen frame and is pinned in `make_figure.AZIMUTH`. Emits
`output/target_definition/panel_a_complex_and_hairpin.{pdf,png}` (85 × 42 mm) and
`panel_b_route_inputs.{pdf,png}` (85 × 44 mm), plus `output/target_definition/renders/`.

> The supplementary NUMBER is not settled and appears nowhere in this folder or in the
> output path.

## Sources

| file | mtime | what is taken |
|---|---|---|
| `~/4al5.cif` (RCSB deposition) | 2026-04-29 | the structure, the citation block, and every geometric claim |
| `working/structures/target_definition/csy4_rna_only.pdb` | copy of the Euler target, 2026-07-03 09:57 | checked against 4AL5 chain B, atom for atom |
| Euler `rna_binder/code/configs/{rfd3__multif,boltzgen__multif,boltzdesign__multif,halludesign,refold__af3,refold__b2,refold__esm}/*.yaml` | read 2026-11-19 | what each route and each re-predictor was passed |
| Euler `rfdiffusion3/examples/rna/csy4_crrna.json` | read 2026-11-19 | the RFdiffusion3 contig |
| Euler `rna_binder/targets/csy4_crrna_free/template_{af3,protenix}.json` | read 2026-11-19 | the 18-nt sequence HalluDesign de novo was given |
| `code/common/csy4.py` | — | the measured base pairs, loop span and brackets |

The configs are not readable from this machine, so they are **transcribed** into
`ROUTES` in `make_figure.py` with the source path beside each line. The two things that
*are* local are re-derived on every run by `verify()`, which asserts rather than prints.

## Verified — all of it re-derived from `~/4al5.cif` on every run

- **Citation.** `_struct.title` = "Crystal structure of the Csy4-crRNA product complex";
  deposited **2012-03-02**; primary citation **Haurwitz, R. E.; Sternberg, S. H.;
  Doudna, J. A., *EMBO J* 31:2824 (2012)**, DOI 10.1038/EMBOJ.2012.107, PMID 22522703.
  The Haurwitz 2010 *Science* structure is **2XLK**, a different entry.
  ⚠️ **Two files in the project repeat the wrong year**:
  Euler `targets/_registry.csv` notes column says `"Haurwitz et al. 2010"` for both csy4
  targets, and `targets/csy4_crrna_free/notes.md` opens "18-nt crRNA from PDB 4AL5
  (Haurwitz et al. 2010)". Neither is a source for the figure; both should be fixed.
- **The campaign target file IS 4AL5 chain B.** 348 heavy atoms matched by
  (residue, name), **0 unique to either side, RMSD 0.000000 Å with no superposition.**
  So the coordinates every structure-conditioned route designed against are the
  **Csy4-bound product-complex conformation**. The registry row calling that target
  `csy4_crrna_free` — "refolded in isolation … starting from sequence" — does not
  describe this file, and `csy4_crrna_bound` was never instantiated.
- **16 vs 18 is not a base pair.** Chain B has 18 modelled residues, 4–21,
  `CACUGCCGUAUAGGCAGC`. Residue 4 is **one atom, O3′**; residue 21 is **four atoms,
  P / O5′ / OP1 / OP2**. Neither carries a base, so the two extra residues cannot be a
  preserved closing pair. The 16 complete nucleotides are residues 5–20,
  `ACUGCCGUAUAGGCAG` = `csy4.SEQUENCE`. (The deposited construct is 20 nt,
  `UUCACUGCCGUAUAGGCAGC`; residues 2 and 3 are unmodelled.)
- **Stem.** Five Watson-Crick pairs, N1/N3 2.64–3.02 Å: res 6–20, 7–19, 8–18, 9–17,
  10–16 = nt2–nt16 … nt6–nt12, exactly `csy4.PAIRS`.
- **No groove channel.** Those five pairs span **122.2° = 0.339 of a helical turn** over
  a 15.1 Å rise (C1′–C1′ vectors projected perpendicular to the axis through the pair
  midpoints). A groove reads as a channel only across a full turn, so major and minor
  must be shown by base-edge colouring, never by shape.
- **Three flipped bases, not five.** Ring-centroid gaps and interplanar angles give
  breaks at nt1→2 (9.20 Å, 45.8°), nt7→8 (8.22 Å), nt9→10 (8.55 Å, 70.4°) and nt10→11
  (10.41 Å, 77.7°): continuous stacks nt2–nt7 and nt11–nt16, so **U8, A9 and U10** leave
  the stack although the loop is 5 nt. The 5′ A1 is also unstacked.

### One correction to a repo value

`csy4-no-groove-channel` and the figures that cite it give **0.37 of a turn, 6 bp,
20.6 Å rise**, measured on AF3 complexes. On the crystal itself the five Watson-Crick
pairs give **0.339 of a turn, 15.1 Å rise**. The two are not in conflict: 4AL5 closes
the loop with a **sheared G·A** (res 11 G N2 – res 15 A N7 **2.95 Å**, N3–N6 3.57 Å),
which is the sixth pair the AF3 count includes. Either way the stem is far under one
turn and the conclusion is unchanged. Quote 0.34/5 bp when the source is 4AL5.

## What each route received (transcribed from the configs, printed on every run)

| route | form | residues | config |
|---|---|---|---|
| RFdiffusion3 + LigandMPNN | coordinates | 16 (res 5–20) | `rfd3__multif/default.yaml` → `csy4_crrna.json` contig `B5-20,/0,80-80` |
| BoltzGen + LigandMPNN | coordinates | 18 (res 4–21) | `boltzgen__multif/default.yaml`, `file: csy4_rna_only.pdb`, `chain: B` |
| BoltzGen + BoltzIF | coordinates | 18 (res 4–21) | same BoltzGen job, one backbone feeding both inverse folders |
| BoltzDesign1 + LigandMPNN | coordinates | 18 (res 4–21) | `boltzdesign__multif/default.yaml`, `pdb_target_ids: B` |
| HalluDesign + LigandMPNN, helix mode | coordinates | 18 (res 4–21) | `halludesign/default.yaml`, `input_file: csy4_helix.pdb`, `fix_chain_index: B` |
| HalluDesign + LigandMPNN, de novo | sequence | 18 nt `CACUGCCGUAUAGGCAGC` | `halludesign/random_init.yaml`, `random_init_chain_spec: "A:80"`, RNA from `template_{af3,protenix}.json` |
| AF3, Boltz-2, ESMFold2 re-prediction | sequence | 16 nt `ACUGCCGUAUAGGCAG` | `refold__{af3,b2,esm}/default.yaml`, `rna_sequence:`, all three `use_rna_msa: false` |

Two things worth keeping in the Methods rather than the legend: `csy4_helix.pdb` chain B
is the same 18 residues (verified residue by residue), its chain A being an 80-residue
CA-only poly-Ala placeholder; and **no generator received any site constraint** —
`targets/csy4_crrna_free/hotspots.json` is an empty stub (`"hotspots": []`) that nothing
reads, so where a binder lands is emergent.

## Panel notes

- **(a)** Both images are saved from **one camera** by `render_target.py`, so the second
  is the first with the protein hidden. They are **not drawn to a shared scale**: the
  hairpin is enlarged, which the caption states. Colours: **the hairpin takes the main
  figures' colours** — backbone theme 700 `#015d82`, base slabs theme 300 `#95c9e8` —
  with the protein neutral 300 (context, rule 3) and silhouettes `synbio.outline`. No
  hex is written by hand; every value comes from `config/palette.json`.
- **(a) the hairpin was `synbio.rna` until 2026-11-19**, changed on Leo's comment:
  *"the hairpin has a color which is weird. It should be the color that the hairpin has
  in the main figures."* It now takes the `BACKBONE`/`BASES` pair that
  `Fig2a_pipeline/render_hairpin.py` gives **this same molecule in this same canonical
  frame**, so the supplement and the main text show one hairpin, not two. Two faults
  were behind the "weird": `synbio.rna` is theme **400**, a pale sky blue that no
  main-figure structure render uses, and it was applied to the **cartoon and the slabs
  alike**, so backbone and bases were the same colour and the hairpin came out one flat
  tone. Against the grey protein the RNA went from ΔE 32.2 to **54.5**, so the subject
  now separates from its context.
- ⚠️ **(a) the two-tone needs `target c`, not `target a`** (fixed 2026-11-20, a day after
  the colour change). `nucleotides tube/slab` gives a slab its **whole residue's**
  colour, so `color {rna} {C_BASE} target a` painted the slabs **nothing** — the
  atom-subset idiom works on a molecular *surface*, not on a slab — and the preceding
  `color #1/B {C_RNA}` had already painted ribbon and slab alike. The working order is
  **residue, then revert the ribbon**:

  ```
  color #1/B:5-20 #95c9e8            # residue: ribbon AND slab both go pale
  color #1/B:5-20 #015d82 target c   # ribbon back to 700; the slab keeps the pale
  ```

  because in tube/slab the base **is** the slab and the sugar+phosphate **are** the
  ribbon. (The tube also needs the cartoon shown and the slabs need the atoms shown —
  line 115 does both.)
- **Verified from the rendered pixels, not from the hex values.** The two stages come off
  one camera, so the PNGs are pixel-aligned and the pixels that changed are exactly the
  slab faces — **36.9%** of the hairpin's solid pixels. Slab mean vs backbone mean in
  `hairpin_y060.png`: ΔE **5.9** with `target a` (i.e. no two-tone at all, just shading)
  → ΔE **34.1** with `target c`. The backbone itself moved ΔE **0.00** between the two,
  so the fix touched only the slabs. Measuring the palette hexes instead would have
  reported a healthy ΔE 42.1 for a render that had no two-tone in it.
- **This is a deliberate departure from the synbio molecular-identity rule**, which
  COLORS.md numbers **8** (the file has no literal rule 6 — the list runs 1–5, 7, 9, 8,
  and `render_target.py` had been citing the number that rule carried before the
  renumbering). That rule sets `rna` = theme 400 because *"the target is context"* — but
  the target is the **subject** of this figure, and the rule is written for schematics
  and imported cartoons. **No main-figure structure render follows it for the hairpin:**
  Fig 2a uses theme 700/300, Fig 2b theme 600, Fig 5c neutral 400. The convention this
  panel now follows is Fig 2a's, which is the closest analogue — same molecule, same
  frame, hairpin as subject.
- **The warm diverging end `#ae4e1a` is deliberately NOT used here.** It is the Csy4
  **site** colour on the **schematic** glyphs and construct maps (Fig 1c, 1f,
  `supp/scaffold_design`, `supp/negative_controls`), where it marks *where* a site sits
  on a linear construct. COLORS.md rule 4 confines diverging to signed data plus those
  glyphs, and `CHIMERAX_RENDERING.md` is explicit that the warm end *"is not in
  `thesis_colors.cxc` precisely because it has no categorical use"*, flagging its one
  appearance in a render as a rule-4 violation. No render script colours the hairpin
  warm, so warm would have matched the schematics while breaking from every structure
  panel in the thesis.
- **(a) shows the campaign file whole, and nothing else** (changed 2026-11-19 on Leo's
  comment: *"I want the hairpin target completely clean like its in the campaign, no
  subs and no other small molecules"*). The loaded structure is chain A protein plus
  chain B RNA only — the deposition's glycerol (`GOL` 1022) and its 128 waters are
  dropped by `prepare()`, and **eleven of those waters and the glycerol sit on auth
  chain B**, so "chain B" alone is not clean. What Leo saw as a small molecule was
  ours, not the deposition's: residues 4 and 21 were drawn as grey **ball atoms**, and
  residue 21 — a bare `P / O5' / OP1 / OP2` — reads exactly as a tetrahedral ligand
  bound under the stem. That block is gone.
- **Nothing was trimmed to achieve it.** The cartoon is now `#1/B`, the **whole 18
  residues** the campaign file contains, not the `5-20` slice. ChimeraX then draws what
  it can: `info residues #1/B attribute ribbon_display` reports **res 21 True** (its
  phosphate is carried as the short 3' ribbon tail) and **res 4 False** — one atom,
  `O3'`, gives no backbone path, so it cannot be ribboned at any setting. 17 of 18
  residues draw; the one that does not is a limit of the coordinates, not a choice.
  The stubs are panel (b)'s subject and are stated there explicitly, which is why
  nothing is added to (a) to mark them.
- **(b)** A filled cell is a residue whose coordinates were supplied; an open cell with
  a short internal rule is a residue supplied as **backbone atoms with no base** (the two
  stubs); a letter is a residue supplied as **sequence only**. In sequence form residues
  4 and 21 are ordinary C's and the model builds them whole, which is why those cells
  carry a letter — the stub is a property of the coordinates, not of the molecule.
  Stem/loop brackets come from `csy4.BRACKETS`, which is measured, not drawn by eye.

## VERIFY — not settled here

- **"No free-state experimental structure of this hairpin exists."** This is a universal
  negative over the PDB and cannot be established from a local file; nothing was queried
  externally. What *is* established is narrower and sufficient for the caption as
  written: the campaign's only structural target is the bound crystal, the registry's
  `csy4_crrna_bound` was defined and never run, and the two Csy4 entries on this machine
  (4AL5, 2XLK) are both protein-bound complexes. If the stronger sentence is wanted,
  someone has to run the PDB search and record it.

## Caption draft

**Every structure-conditioned route designed against the Csy4-bound conformation, and
the 16-vs-18-nucleotide difference between routes is two crystallographic stubs, not a
base pair.** (a) The Csy4–crRNA product complex (PDB 4AL5; Haurwitz, Sternberg & Doudna,
*EMBO J* 31:2824, 2012) and the repeat hairpin taken from it, from one camera, the
hairpin enlarged: the two share no scale. It is the campaign's target file whole — chain
B of this entry atom for atom (348 heavy atoms, RMSD 0.000 Å, no superposition), with
the deposition's glycerol and waters dropped; residue 4, a lone O3′, carries no ribbon.
Generation used this **protein-bound** product conformation while validation
re-predicted the hairpin free from sequence, so the two share no conformational
reference. `VERIFY:` if the PDB search is run and comes back empty, add "and no
free-state experimental structure of this hairpin exists" here. (b) What each design
route and each re-predictor received, over the 18 residues modelled in chain B. Filled =
coordinates; open with a rule = backbone atoms and no base; letter = sequence only.
Residues 4 and 21 are a lone O3′ and a bare 5′-phosphate, so the 16 complete nucleotides
(5–20, ACUGCCGUAUAGGCAG) are all that carries a base: the 16-vs-18 split is where a tool
draws the boundary, not a different RNA. Stem and loop are the five Watson-Crick pairs
and the GUAUA loop measured in 4AL5.
