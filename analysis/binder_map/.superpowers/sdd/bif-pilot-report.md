# BG+BIF pose-reference pilot (10 designs) — IN PROGRESS

Status: archives staging, DockQ runs not yet done. This file is written early and updated as
results land. If it still says IN PROGRESS, the run was interrupted before finishing.

## Question

For the BG+BIF stack the campaign stores BoltzGen's own re-prediction as `designed_cif_path`
(the stack is invoked `--steps inverse_folding folding` and harvests `refold_cif`). Every pose
metric for BG+BIF is therefore measured against a structure that already moved (RNA-aligned
protein RMSD median ~35 Å; binder centroid displaced 19–32 Å). The true generated backbone is the
BG+LiMPNN `designed_cif_path` at the same (uuid, design index) — LigandMPNN is pure inverse
folding and moves no atoms.

Pilot: recompute DockQ for 10 BG+BIF designs (all in the ordered 12k library) against the true
generated backbone, and count how many lose the passing_v4 / 12k gate.

Gate: `(sc_RMSD ≤ 3 Å on AF3 OR ESM) AND (DockQ ≥ 0.7 on AF3 OR ESM) AND seq_run ≤ 6`.

## Method

- `code/scripts/dockq_self.py` on the Euler host python3 (DockQ not reimplemented).
- Run 1 (control / reproduction): `--ref-col` = original BG+BIF `designed_cif_path`.
  Must reproduce stored `dockq_af3` / `rmsd_af3` (filter store) == `dockq_self_af3` /
  `dockq_LRMSD_af3` (master).
- Run 2 (corrected): `--ref-col` = BG+LiMPNN `designed_cif_path` at the same (uuid, design, seq).
- Inner archives restaged out of the post-2026-09-11 repack
  (`campaigns/csy4_crrna_free/boltzgen__{boltzif,ligandmpnn}/midrun_v8.tar.zst`) into a tree
  mirroring the pre-repack layout, one streaming pass per outer archive.

## Stored values (filter store == master, verified exact)

`sc_RMSD` == `dockq_LRMSD`; store `dockq_af3` == master `dockq_self_af3`. Confirmed for all 10
(max|diff| = 0).

| design_id | sc_RMSD af3 | DockQ af3 | sc_RMSD esm | DockQ esm | seq_run |
|---|---|---|---|---|---|
| 3ab1f5f0_design12_seq1 | 4.7765 | 0.5947 | 2.5668 | 0.7301 | 3 |
| 3ab1f5f0_design12_seq5 | 2.3701 | 0.7286 | 3.9111 | 0.5589 | 2 |
| 3ab1f5f0_design5_seq5 | 24.6243 | 0.0400 | 1.6687 | 0.7572 | 2 |
| 3ab1f5f0_design6_seq0 | 4.4508 | 0.6816 | 2.7656 | 0.7240 | 2 |
| 3ab1f5f0_design6_seq2 | 2.1619 | 0.7561 | 9.7629 | 0.4281 | 3 |
| 3ab1f5f0_design6_seq3 | 9.6965 | 0.3719 | 2.2804 | 0.7393 | 3 |
| 3ab1f5f0_design6_seq7 | 1.5221 | 0.7441 | 3.7652 | 0.6472 | 3 |
| 66dfd45f_design5_seq1 | 2.2564 | 0.7179 | 2.3999 | 0.6875 | 2 |
| 66dfd45f_design5_seq2 | 2.3537 | 0.7064 | 1.3987 | 0.7237 | 2 |
| 66dfd45f_design5_seq3 | 1.7011 | 0.7093 | 0.9942 | 0.7898 | 2 |

All 10 pass the gate on stored values, as expected (they are in the ordered library).

## Reproduction run

PENDING

## Corrected run

PENDING

## Verdict

PENDING
