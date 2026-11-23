# Figure 3b: one binder per generator

Four binders, one from each design stack, drawn in the **same shared RNA frame
and camera as the Fig 2b structure panels**, so the hairpin sits identically in
all of them and only the binder moves. Binder colour is the stack colour, so
colour means generator and nothing else.

The published panels are the **largest backbone cluster per generator**, all at
**y120**, 3000 px.

| File | Design | Cluster | H / E | Candidates |
| --- | --- | --- | --- | --- |
| `rfd3.png` | `8eaafbb8_design14_seq0` | RFD3 cl3 | 87.5 / 0.0, 3H/0E | 396 |
| `bg.png` | `81878aac_design15_seq0` (BG+LiMPNN) | BG cl9 | 36.2 / 26.2, 2H/3E | 14 |
| `bd1.png` | `20a7e2a8_design7_seq0` | BD1 cl0 | 15.0 / 45.0, 1H/4E | 20 |
| `hd.png` | `58be5280_recycle15_sample0_pack1` | HD cl1 | 61.2 / 0.0, 3H/0E | 44 |

Secondary structure was recomputed from each CIF and matches the table (HD
reads 61.3, which is 49/80 rounded the other way).

## Metrics
**All four pass every term of the `passing_v6_final` gate.** The HD pick is also
in the ordered 12k (`csy4_binder_library_12k_v6_k30.csv`); the other three are
not. Values are from the Euler store `filter_store_v2v9_hdtop`, AF3 unless noted.

| | BD1 | RFD3 | HD | BG |
| --- | --- | --- | --- | --- |
| ipTM (max) ≥ 0.6 | 0.685 | 0.67 | 0.82 | 0.76 |
| ligand_RMSD ≤ 3 Å | 2.62 | 1.52 | 1.91 | 1.88 |
| DockQ ≥ 0.7 | 0.765 | 0.838 | 0.796 | 0.837 |
| AA run ≤ 5 | 3 | 4 | 4 | 4 |
| Foot clash ≤ 5 | 1 | 5 | 0 | 0 |
| dSASA Å² | 1375 | 2600 | 1536 | 1166 |
| H-bonds af3 / esm | 2 / 1 | 6 / 7 | 1 / 4 | 1 / 1 |
| packstat | 0.675 | 0.568 | 0.673 | 0.581 |
| shape compl. | 0.529 | 0.389 | 0.313 | 0.445 |

## Sets that were rendered and dropped
`make_figure.py` keeps them so they are not picked up again by accident.

- `--set library`: four designs from the ordered 12k (`88b66226_design33_seq1`,
  `ef8ee7a5_design8_seq2`, `8bca6080_design31_seq0`,
  `1ef1a7dc_recycle11_sample1_pack1`), published at y060 until 2026-11-08.
  Replaced by the cluster set.
- `--set clusters`, `--set clusters2`: two earlier rounds of cluster picks.
  **None of those designs was in the 12k and almost none passed the gate**
  (ligand_RMSD 10 to 44 Å, DockQ 0.02 to 0.31), because they were the largest
  backbone clusters over all designs rather than over candidates.

## The BG row needs disambiguating every time
Both BG stacks contain a design called `81878aac_design15_seq0`, with the
identical fold, so the H/E in the cluster table cannot separate them. The gate
does: **BG+LiMPNN passes every term, BG+BIF fails ligand_RMSD (13.1 Å) and DockQ
(0.172)**, so a candidate pool can only mean BG+LiMPNN. This recurred in all
three rounds; in the first round the two differed in fold and the H/E settled
it, in the later two only the metrics could.

## Provenance
CIFs come from **several different midruns** (RFD3 v6, BG v6, BD1 v9, HD v3),
resolved through the store's `meta.json` templates and `tpl_idx`, then extracted
from the unpacked NAS mirror at
`euler_work/rna_binder/campaigns/csy4_crrna_free/`. The Euler tree still has
these midruns packed inside `midrun_vN.tar.zst`, so the mirror is the fast path.
Extracted copies live in `working/structures/fig02b_examples/{cifs,af3}/`.

## Rebuild
```
~/.venvs/binder_map/bin/python3 code/Fig2b_stack_binders/make_figure.py <out_dir> --width 3000 --angle 120
```
Drop `--angle` for the full 12 angle turntable to pick from, then re-run with
the chosen `--angle` at full width. `--set` selects a different panel set.

## Things that are deliberate
1. **Shared RNA frame with Fig 2b.** The reference hairpin is still panel a's
   design (`8b13398e_design50_seq0`), opened as model 1 and never drawn as a
   binder. Every other complex is aligned onto it by C1' and only the reference
   RNA surface is shown, so the surface in all panels is literally the same
   atoms.
2. **Cartoon only, no sidechains.** Design CIFs are backbone-only (N, CA, C, O),
   so no H-bond or stacking can be drawn from them. AF3 refolds for all four are
   already extracted next to the designs if contacts are ever wanted.
3. **Colour carries generator identity only**, per `config/COLORS.md`.

## Map callouts (`callout_*.png`, set `map8`, y120, 3000 px)

Eight binders, two per generator, to annotate the campaign UMAP. All eight pass
every term of `passing_v6_final` and all eight are in the ordered 12k.

| Panel | Design | Map x, y | Cluster | ipTM | DockQ | ligand_RMSD | dSASA |
| --- | --- | --- | --- | --- | --- | --- | --- |
| rfd3_a | `ccddfc0f_design11_seq6` | 7.0, -0.7 | 5 | 0.710 | 0.819 | 1.45 | 2764 |
| rfd3_b | `7c0a8297_design17_seq0` | 4.7, 5.4 | 4 | 0.830 | 0.838 | 0.96 | 2363 |
| bg_a | `aebd763e_design11_seq0` | 2.7, 4.0 | 1 | 0.918 | 0.896 | 1.64 | 1495 |
| bg_b | `e1e86ec3_design9_seq0` | 9.4, 9.6 | 14 | 0.854 | 0.873 | 1.11 | 1247 |
| bd1_a | `b6ccc04e_design28_seq0` | 4.1, 2.6 | 8 | 0.820 | 0.866 | 1.51 | 1997 |
| bd1_b | `731e2185_design57_seq0` | 4.1, 7.1 | 4 | 0.680 | 0.860 | 1.66 | 2603 |
| hd_a | `98a6ba23_recycle9_sample0_pack1` | 1.2, 1.8 | 7 | 0.859 | 0.829 | 1.96 | 2076 |
| hd_b | `13538ba3_recycle15_sample0_pack1` | 7.5, 3.8 | 3 | 0.768 | 0.862 | 1.80 | 2052 |

Selection: `pick8.py` on Euler (SLURM, 64 GB). Rebuilds the fold map, drops
`deg0` isolates, keeps the elected representative per backbone, applies the full
gate (4,396 survive: RFD3 2,575, BG 1,452, BD1 291, HD 78), scores on
DockQ, ipTM, ligand_RMSD, 12k membership and dSASA, then takes the top design per
generator plus the best one at least 4.0 UMAP units away (2.0 as fallback).

hd_a was originally `1ef1a7dc_recycle11_sample0_pack1`, the same backbone as the
published HD generator panel. Replaced with `98a6ba23` (job `hdalt.sbatch`, which
lists all 78 gate-passing HD representatives): fold cluster 7, so distinct from
the published HD panel (cluster 1), from hd_b (cluster 3) and from the dropped
pick (cluster 13).

Known limitation: the search sees only the elected representative per backbone,
so a backbone whose representative fails the gate is invisible even if a sibling
design passes.

Rebuild:

    ~/.venvs/binder_map/bin/python3 code/Fig2b_stack_binders/make_figure.py OUT \
        --set map8 --width 3000 --frames 12 --angle 120

Numbers 1-8 in `output/Fig2b_map_callouts/map_callouts.*` refer to these panels;
see `code/Fig2b_map_callouts/README.md` for the map positions and the placement rules.

`euler/pick8.py` and `euler/hdalt.py` are the selection scripts for the eight
map callouts, pulled off Euler into version control on 2026-11-09. Submit them
from `/cluster/work/wicky/lvalsecchi/figures/`.


## rfd3_a swapped for a pi-stacking design (2026-11-10)

`rfd3_a` was `4bcad7a8_design62_seq0`. Measured with
`code/Fig2b_map_callouts/analyse_examples.py`, **none of the original eight had a
single aromatic stack**, and in all eight no aromatic ring was even within 7.5 A
of a base. That matches the campaign rate: ~8% of designs stack in one backend,
0.09% in all three, against 71% for natural RBPs.

Replacement `ccddfc0f_design11_seq6` stacks **TYR24 on U10 at 4.21 A, interplanar
angle 13.2 deg**. U10 is a LOOP nucleotide (loop = nt 8,9,10) and is also this
design's most contacted base, 69 contacts. Same generator, same fold cluster 5,
0.65 UMAP units from the design it replaces, also in the ordered 12k, so the map
position and the composite's slot assignment are unchanged.

Costs: ipTM 0.839 -> 0.710 and H-bonds 6 -> 1. It trades hydrogen bonding for a
stack, which is coherent rather than alarming.

⚠️ **The stack is in the AF3 refold only.** Boltz-2 and ESMFold2 both report zero
stacks for this design, and backends barely agree on stacking anywhere. Caption it
as "TYR24 stacks on U10 in the AF3 model", not "this design stacks".

⚠️ **The search was shallow.** Only 180 of the ordered 12k appear in
`app_metrics_100k.csv`, the one table with pi_stack computed, so 98.5% of the
library was never tested. This is the best of 180, not the best of 12,000.
Note also the limitation already recorded above: pick8 sees only the elected
representative per backbone, and this design is `seq6` whose backbone's
representative is `seq0` -- which is exactly why it was invisible to the original
selection.

Backups: `README.md.bak_pre_pistack_swap`, `output/Fig2b_stack_binders/bak_pre_pistack_swap/`
(the previous eight renders), `inputs/map8_callouts.csv.bak_pre_pistack_swap`.
All eight renders were rebuilt, since the set shares one frame and one crop.
