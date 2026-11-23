# ChimeraX rendering recipe (design vs refold overlays)

How the structure panels for Fig03a were rendered. Everything here was run and
verified on 2026-11-05; the gotchas at the bottom are mistakes that actually
happened, not hypotheticals.

## Binary
```
/Applications/ChimeraX-1.11.1.app/Contents/MacOS/ChimeraX
```
Run a script with:
```
ChimeraX --exit --cmd "cd <dir>; open render.cxc"
```
**Do not add `--nogui`.** It parses and runs the script fine but silently saves
no images (no OpenGL context). The symptom is a clean log with the alignment
RMSD printed and zero PNGs on disk.

## Getting the structures: use the NAS mirror first
A **full mirror of the Euler tree** is on the NAS at
`/Volumes/wicky/users/lvalsecchi/euler_work/rna_binder/`, and much of it is
**unpacked**, e.g. `campaigns/csy4_crrna_free/<stack>/midrun_v8/` exists as a
real directory with `designs/`, `refold_af3/`, `refold_esm/` and `refold_b2/`
inside. Pulling one design from there means extracting a single small per-uuid
`designs.tar.zst`, not sweeping a 30 GB archive, and it needs no SSH. Check the
mirror before going to Euler.

The mirror also carries things the gate never used: **`refold_b2` (Boltz-2)**
alongside the two gate backends, and `_confidence.json` / `_ipsae.json` per
design (the confidence JSON has `iptm`, `ptm`, `complex_plddt`).

## Getting the structures off Euler (fallback)
Per-design CIFs live in **nested** `tar.zst` archives. The store's
`meta.json` has a `templates` list; `tpl_idx` for a design indexes into it, and
each entry has three paths (`af3_msaoff_cif_path`, `esm_cif_path`,
`designed_cif_path`) with `<U>` = uuid8 and `<ID>` = design_id.

```python
import json
m = json.load(open(".../filter_store_v2v9_hdtop/meta.json"))
for t in m["templates"][TPL_IDX]:
    print(t.replace("<U>", UUID8).replace("<ID>", DESIGN_ID))
```

The printed path looks like `campaigns/.../midrun_v8/designs/<uuid>.tar.zst!<ID>.cif`,
but that inner archive is itself inside `midrun_v8.tar.zst` (20-30 GB). Two-step
extraction, `--occurrence=1` so tar stops once it has the members:

```bash
D=.../campaigns/csy4_crrna_free/<stack_dir>
nice tar --zstd -xf $D/midrun_v8.tar.zst --occurrence=1 \
  midrun_v8/designs/<uuid>.tar.zst \
  midrun_v8/refold_af3/msaoff_<uuid>/designs.tar.zst \
  midrun_v8/refold_esm/msaoff_<uuid>/designs.tar.zst
tar --zstd -xf midrun_v8/designs/<uuid>.tar.zst <ID>.cif
tar --zstd -xf midrun_v8/refold_af3/msaoff_<uuid>/designs.tar.zst designs/<ID>_AF3_seed1_sample0.cif
```
Takes a few minutes per design; run it in the background. Stack directory names:
`boltzdesign__ligandmpnn`, `boltzgen__ligandmpnn`, `boltzgen__boltzif`,
`halludesign__ligandmpnn`, `rfdiffusion3__ligandmpnn`.

## Which frame to superpose in
This is the decision that makes the panel say what the filter measures.

| Panel | Superpose on | Command |
|---|---|---|
| ligand_RMSD | the **binder**, so fold error shows | `align #2/A@CA toAtoms #1/A@CA move structures` |
| DockQ | the **RNA**, so pose error shows | `align #2/B@C1' toAtoms #1/B@C1' move structures` |

`#1` is the design, `#2` the refold. `align` is a straight least-squares fit with
no outlier pruning, which is what you want for an honest picture. `matchmaker`
prunes and will flatter the result.

## Colours: follow COLORS.md, do not hand-pick hexes
`COLORS.md` governs ChimeraX renders as well as plots: **never hand-pick a hex in
a figure script.** ChimeraX has its own named-colour file, so load it first and
use the names:

```
open config/thesis_colors.cxc
color #1/A thesis_bg_limpnn
color #1/B thesis_grey_300
graphics silhouettes true width 1.5 color thesis_grey_700
```

Names available: `thesis_petrol_100`..`_900`, `thesis_grey_100`..`_900`, and one
per stack (`thesis_rfd3_limpnn`, `thesis_bg_limpnn`, `thesis_bg_bif`,
`thesis_boltzdesign_limpnn`, `thesis_hallu_limpnn`).

Rules that bite in a structure panel:
- **Stack colours mean stack identity and nothing else** (rule 2). Use a stack
  colour for a binder of that stack. Do not reuse one to mean "the other
  backend" or "the clashing part".
- **Diverging is only for signed data with a meaningful zero** (rule 4). It is
  not a source of contrast colours. The warm end `#ae4e1a` is not in
  `thesis_colors.cxc` precisely because it has no categorical use.
- **Several things of the same kind** (three refolds of one design) are a single
  series, so take theme steps: `thesis_petrol_800` / `_600` / `_400`. Monotone
  lightness also survives greyscale (rule 5).
- **Neutrals** for the RNA, reference geometry, axes and silhouettes (rule 3).

### NEED: synbio element colours are undefined
`COLORS.md` rule 6 says hairpin, binder, DNA and receptor colours are **"to be
defined; derive from theme and neutrals"**. Until Leo defines them, every
structure panel's colour choice is provisional:

| Element | Used here | Status |
|---|---|---|
| binder, single model | its stack colour | consistent with rule 2 |
| RNA hairpin | `thesis_grey_300` | neutral, consistent with rule 3 |
| three refolds of one design, labelled | `thesis_petrol_800/600/400` | theme steps, rule 1 |
| three refolds as one uncertainty cloud | `thesis_grey_400` at 55%, design `thesis_grey_800` | neutrals, rule 3 |
| **hairpin foot extension** | `#ae4e1a` (diverging warm) | **violates rule 4, NEED a defined colour** |
| **clashing residues** | `#ae4e1a` | same, and dropped from the final foot render |
| RNA *surface* in older renders | `#d7dcdf` | **hand-picked, not in the palette** |

The foot render `TR_m30_y330` deliberately keeps the warm foot because Leo asked
for orange/red, and the palette offers no sanctioned warm that is not a stack
colour or the diverging end. Resolve by adding a synbio entry to
`palette.json` + `thesis_colors.cxc` (rule 6), then re-render.

`width 1800` for screen, `width 3000` for print. `supersample 3` always.

### Turntable
Call `view` once, then `zoom 0.85` for margin, then `turn y 30` between saves
(24 views: two tilts, `turn x -30` and none, 12 rotations each). Do **not**
re-run `view` between frames: it refits the camera each time and the object
changes size from frame to frame. Earlier sets used `zoom 0.8` / `turn y 45`;
do not mix the two in one panel row, the object size differs.

## Gotchas

- **`set bgColor white` is still required with `transparentBackground true`.**
  The alpha is correct either way, but the RGB under transparent pixels is
  whatever the background was, so without it the PNG composites to black in any
  viewer that ignores alpha.
- **The PNGs really are transparent.** Corner alpha 0, ~86% of pixels fully
  transparent. Chat previews and Preview.app flatten onto white, which looks
  like it failed. Check with `sips -g hasAlpha`, or composite onto a colour.
- **The refold is drawn at 35% transparency, so it picks up the backdrop.** Fine
  on white, tinted on anything else. Make it opaque if the panel sits on colour.
- **Design and refold RNA are not always the same length.** A BoltzDesign design
  carried an 18-nt hairpin (`CACUGCCGUAUAGGCAGC`) while its refolds had only the
  16-nt core, which is design residues **2-17**. `align` then fails with
  "Unequal number of atoms to pair, 16 and 18". Map the range explicitly
  (`#1/B:2-17`) and surface the same range so panels stay comparable. Check the
  residue counts before assuming they pair.
- **Design CIFs are backbone only; sidechains live in the refolds.** A design
  binder has 4 atoms per residue (N, CA, C, O) and **no CB**, so no sidechain can
  be drawn from it. Render sidechains from an AF3 or ESM refold and label them as
  such. ChimeraX does not error on this: its `sidechain` selector falls back to
  CA when a residue has no sidechain atoms, so `show ... & sidechain atoms`
  reports a match (one atom per residue) and saves an image with no sidechains
  in it. Count atoms per residue before trusting it (design 320 vs refold 541
  for an 80-residue binder).
- **A `.cxc` file has no `;` comment syntax.** A `;` line is parsed as a command,
  fails, and aborts everything after it, which looks like the render silently
  producing nothing. Use no inline comments, or run once with `--nogui` first to
  surface command errors (it validates but cannot save images).
- **Silhouettes are GLOBAL, not per model.** `graphics silhouettes false models #2`
  and `graphics model #2 silhouette false` both error. Turning outlines off for
  one model means turning them off for the whole scene. This matters for faded
  overlays: with silhouettes on, a cartoon keeps its full outline at any
  transparency, so raising transparency past ~80% removes shading but leaves
  three traced ribbons. If something must genuinely dissolve, silhouettes have to
  go off everywhere. Mixing the two needs two renders on the same camera,
  composited.
- **`models` takes `#1-4`, not `#1,#2,#3,#4`.** The repeated-`#` form fails with
  `Invalid "models" argument: only initial part "#1" of atom specifier valid`.
- **`nucleotides` only acts on displayed atoms.** `nucleotides <spec> tube/slab`
  silently does nothing unless `show <spec> atoms` ran first, leaving a bare
  backbone ribbon with no bases and no error. Colour the bases with `target a`
  as well as the cartoon, or they keep the default colour. `tube/slab` is the
  house style for RNA; `ladder` is the rounder alternative.
- **`hide /B models` hides the WHOLE model**, because chain B belongs to it, so
  the image comes out blank with no error. To draw only some chains, use
  `hide cartoons` then `cartoon <spec>` for the ones you want.
- **Chain IDs are NOT always `A` = binder, `B` = RNA.** Corrected 2026-11-05
  against all 50 extracted triplets in `foldseek_map/browse50/`:

  | File | Chain A | Chain B | Count |
  |---|---|---|---|
  | `*_af3_*.cif`, `*_esm_*.cif` (refolds) | binder | RNA | 50 / 50 |
  | `*_design_*.cif`, every stack except RFD3 | binder | RNA | 39 / 39 |
  | `*_design_*.cif`, **RFD3 + LigandMPNN** | **RNA** | **binder** | 11 / 11 |

  So an RFD3 design file has the chains the other way round, deterministically,
  while its own refolds do not. Hardcoding `#1/A` = binder silently colours the
  RNA as the binder and aligns the wrong chains: it produces a picture, not an
  error. **Detect the chain by content, not by ID** (the RNA chain is the one
  whose residue names are all in `{A,U,G,C}`) before writing any `.cxc`.

  The `label_asym_id` column in these CIFs shows things like `Axp`, which is not
  what ChimeraX selects on; it uses `auth_asym_id`.

- **`transparency <res> 0 target ac` also makes that residue's CARTOON opaque.**
  Resetting a highlighted residue with `target ac` while the rest of the chain
  is faded leaves a solid, fully saturated ribbon wedge beside the sticks that
  reads as a nucleotide slab or a stray blob. It is the residue's own cartoon.
  Use `target a` to reset the atoms only. Found 2026-11-14 in Fig5d; Fig5c had
  the same artefact and was re-rendered.
- **`nucleotides <spec> fill` resets slab colour.** Call it before colouring,
  not after, or the per-residue colour is silently overwritten.

## Where the rendered examples live
**The registry of what was actually rendered, with per-design metrics, chosen
orientations and caption caveats, is `working/structures/REGISTRY.md`.**

`LV_ThesisFigures/working/structures/<name>/`, each with its CIFs, the exact
`.cxc` used, and the PNGs. Current ones:

| Folder | Filter it illustrates | Design |
|---|---|---|
| `_reference` | the transplant foot reference | `extended_csy4_model_0` (21 nt) |
| `rmsd_fail_1160e10b` | ligand_RMSD | `csy4free_HD_LiMPNN_len80_1160e10b_recycle13_sample2_pack1` |
| `dockq_fail_B_605cc505` | DockQ (chosen) | `csy4free_BD1_LiMPNN_len80_605cc505_design55_seq3` |
| `dockq_fail_A_7960a79f` | DockQ (alternative) | `csy4free_BG_LiMPNN_len80_7960a79f_design46_seq7` |
| `aarun_fail26_f5d0270f` | longest AA run (chosen) | `csy4free_RFD3_LiMPNN_len80_f5d0270f_design3_seq6` |
| `aarun_fail_886a362a` | longest AA run (alternative) | `csy4free_BG_LiMPNN_len80_886a362a_design49_seq1` |
| `foot_fail_d1557f4c` | foot clash | `csy4free_BG_LiMPNN_len80_d1557f4c_design49_seq1` |
| `iptm_fail_64e24616` | ipTM | `csy4free_BG_LiMPNN_len80_64e24616_design1_seq9` |

The foot-clash panel is **not** built from the recipe above: its structures come
from the metric's own code on Euler (`make_foot_fig.py`, which imports
`foot_clash/hairpin_foot_clash.py`), emitting chains A binder, B design RNA,
F foot, T full transplant. See the registry.

## Finding candidate designs
Query the store directly with numpy, no pandas needed. The gate that was
actually applied is in `code/fig07_filter_funnel/README.md` (gate provenance;
the funnel panel itself is `code/Fig3a_filter_funnel/`); to find a design
that fails exactly one term, require every other term and negate that one. Note
the ligand_RMSD and DockQ terms are **OR** across backends, so "fails DockQ" means
both `dockq_af3 < 0.7` **and** `dockq_esm < 0.7`.

For a globular-looking binder, filter on radius of gyration against the compact
globule expectation `2.2 * N^0.38`; `rog_af3 / that` under about 1.05 is compact.
A long single helix scores well above it and renders badly.
