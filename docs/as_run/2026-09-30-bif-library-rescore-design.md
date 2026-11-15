# BG+BIF library rescore — design

**Date:** 2026-09-30
**Question:** of the 2,400 BG+BIF designs in the delivered 12k library, how many
fail the campaign's own gate when scored against the **true BoltzGen backbone**
instead of against BIF's own re-prediction?
**Deliverable:** that rate, plus the per-design table behind it.

---

## 1. Why

`code/scripts/stages/boltzgen_multif.sh` runs one BoltzGen job that feeds both
inverse folders, so `BG+BIF` and `BG+LiMPNN` share a `uuid8` and the same
`spec_<D>.cif` backbones (header L2-5, backbones written L146-149).

- **BG+LiMPNN** stores `spec_<D>` + the MPNN sequence — LigandMPNN moves no atoms.
- **BG+BIF** stores a file harvested from `bg_output/.../refold_cif/` (L209-212,
  copied L228-246). `--steps inverse_folding folding` (L203) runs the *folding*
  step, and `boltzgen_patch/writer.py:92` writes that file from
  `pred_dict["coords"]` — the best of N samples by `0.8*ipTM + 0.2*pTM` (L72-74)
  — with pLDDT in the B-factor column (L85-89).

So for BG+BIF alone, `dockq_self`, `dockq_LRMSD` (= `sc_RMSD`) and
`ca_rmsd_vs_designed` were measured against a **prediction of the designed
sequence**, not against the generated backbone. Both are gate terms.

### Evidence from the 20-pair probe (`analyses/bif_backbone_probe/`)

Binder Ca vs `spec_<D>`: the BIF file has the right **fold** (0.82 A median after
fitting) but the wrong **pose** — 4.4 A median RNA-aligned, 5/20 beyond 10 A. So
it is not a rigid-body frame shift; the binder genuinely sits elsewhere on the
target. Rescored: sc_RMSD median 2.723 -> 5.655 A, DockQ 0.708 -> 0.545,
**14/20 flipped PASS -> FAIL** on the base gate.

Those 20 are not a sample to extrapolate from — but the 2,400 are not a sample at
all. They *are* the population the question is about, so conditioning on "passed
under the wrong reference" is correct here, not a bias.

---

## 2. Scope

**In:** the 2,400 `BG+BIF` rows of
`analyses/library_selection/csy4_binder_library_12k_v6_k30.csv`; AF3 and
ESMFold2 refolds; the `sc_RMSD <= 3` and `DockQ >= 0.7` gate terms.

**Out (explicit non-goals):**
- the other four stacks — only BG+BIF stores a re-prediction, so they cannot move
- the remaining ~1.478M BG+BIF designs, and the gated pool
- re-running the 12k selection, or producing a corrected library
- backfilling corrected columns into the master / filter store
- fixing `boltzgen_multif.sh` forward

Each is a separate decision that this result should inform first.

---

## 3. Established facts (measured 2026-09-30, not assumed)

| fact | value |
|---|---|
| BG+BIF rows in the 12k | **2,400** |
| distinct `uuid8` | 1,353 |
| campaign spread | v2 235 - v4 138 - v5 158 - v6 133 - v7 542 - v8 629 - v9 565 |
| rows with both backends already scored | 2,400 / 2,400 |
| `spec_<D>` backbone directly available | **2,292 / 2,400** |
| missing backbone | 108, **all v2** (the known v2 backbone cleanup) |
| of those, recoverable via a BG+LiMPNN design on the same backbone | **107 / 108** |
| genuinely unresolvable | **1** |

Archive structure (checked on `midrun_v6`):
- **one** zstd frame, 0.5 MB window — no random access, no intra-archive parallel decode
- tar order: `designs/ -> refold_af3/ -> _refold_arrays/ -> refold_b2/ -> refold_esm/ -> run_*/`

Two consequences: `refold_esm` sits *after* the `refold_b2` third we never use, so
we must decompress past it; and `designs/` comes first, so skipping the design
files buys nothing. The only available skip is the `run_*` tail.

---

## 4. Method

Reuse the campaign's own scorers unchanged — `code/scripts/dockq_self.py` (stock
DockQ) and `code/scripts/lrmsd_rna_aligned.py` (RNA-aligned LRMSD) — recombined
with the formula from `dockq_rna_recombine.py`:

```
LRMSD_rna = dockq_LRMSD_stock * (lrmsd_rna / lrmsd_prot)     # = sc_RMSD
DockQ_rna = [ fnat + 1/(1+(iRMSD/1.5)^2) + 1/(1+(LRMSD_rna/8.5)^2) ] / 3
```

The campaign-side numbers are taken from the library CSV rather than recomputed,
justified by the probe's **exact** reproduction: scoring the 20 pairs against the
BIF design reproduced the stored values to max |delta| **0.0005 A** (sc_RMSD) and
**0.0000** (DockQ) over 40 values.

### The reference

`spec_<D>.cif` cannot be handed to the scorers directly: it carries BoltzGen's
own co-designed sequence, only ~50% identical to the BIF sequence. That trips
`lrmsd_rna_aligned.py`'s `MIN_REGISTER = 0.9` residue-name gate (returns `None`)
and makes DockQ's internal sequence alignment unreliable.

**Fix — sequence transplant:** take the `spec_<D>` protein backbone atoms
(N, CA, C, O), set residue names positionally from the design's own sequence,
keep the RNA verbatim. This reconstructs the file BG+BIF *should* have written,
and has exactly the shape of a real BG+LiMPNN design file — which is already what
the campaign uses as ground truth for the three LiMPNN stacks.

**The recipe is validated, not assumed:** run with the LiMPNN sequence it
reproduces the real BG+LiMPNN design file atom-for-atom on 20/20 pairs
(0.0005 A over 320 backbone atoms, 320/320 residue names).

**Fallback for the 107 v2 designs with no `spec_<D>`:** transplant onto the
BG+LiMPNN design file at the same `(uuid8, design_idx)` instead. Exact, not
approximate — that file *is* the spec backbone, by the same validation.

### Why sidechains do not confound this

DockQ's LRMSD and iRMSD are backbone terms; only `fnat` reads sidechains. Measured
on the probe: stripping sidechains from the reference leaves sc_RMSD identical
(max |delta| 0.0015 A) and moves only fnat (+0.048). The reference is
backbone-only on both sides of the comparison that matters, so the DockQ change
is pose.

---

## 5. Architecture

Extraction is the expensive step and has been the recurring tax, so it is
separated from scoring and paid **once**.

| | bytes |
|---|---|
| must be **read** | ~174 GB (7 BIF archives + LiMPNN v2 for the fallback) |
| actually **wanted** | 2,400 x (af3 ~83 KB + esm ~104 KB) ~= **0.45 GB, 4,800 files** |

That 385:1 ratio is the design. The working set is small enough to keep on Lustre
permanently (0.45 GB, 4,800 inodes), so every later scoring run — bug fixes,
re-gating, a new metric — reads no archive at all.

Everything lands in a new `analyses/bif_library_rescore/`, leaving the 20-pair
probe untouched.

```
bif_library_rescore/
  manifest.tsv              2,400 rows: design_id, uuid8, design_idx, campaign,
                            sequence, ref_kind {spec|limpnn_fallback|unresolvable},
                            and the library's sc_rmsd/dockq per backend
  cifs/<design_id>_{af3,esm}.cif      the persistent working set (~0.45 GB)
  refs/<design_id>.pdb                transplanted references
  score_out/                          scorer output, per campaign
  rescore_2400.csv                    per-design result
  REPORT.md                           the rate + caveats
```

### Stages

1. **`build_manifest.py`** — the 2,400 rows, each tagged with its `ref_kind`.
   Pure metadata, no archive access.
2. **`extract_refolds.py`** — one SLURM task per outer archive (7 BIF + LiMPNN v2
   = 8, run in parallel, so wall time is the largest archive, not the sum).
   Streams once, pulls `refold_af3` + `refold_esm` members for its campaign's
   uuids, stops after the last needed `refold_esm` member. Idempotent: skips
   files already present, so a failed task is re-runnable alone.
   A `--with-design` flag additionally pulls `designs/<uuid>.tar.zst`; it is used
   **only for the T2 campaign**, to feed correctness gate 3. The design files sit
   earlier in the same stream, so the flag costs no extra read.
   Backbones come separately from the loose per-uuid `_backbones/` archives,
   which the 2026-09-11 repack did not touch.
3. **`build_refs.py`** — the transplant, extended with the fallback path.
4. **`score.py`** — drives `dockq_self.py` (4,800 runs) and
   `lrmsd_rna_aligned.py` (2,400 rows) off `cifs/` only, then recombines.
5. **`report.py`** — re-gate and emit the rate.

### Gate

A design flips PASS -> FAIL when, under the corrected reference, it stops
satisfying:

```
min(sc_RMSD_af3, sc_RMSD_esm) <= 3    AND    max(DockQ_af3, DockQ_esm) >= 0.7
```

`seq_run`, `foot` and `iptm_mean` are sequence- or refold-derived, do not touch
the reference, and are therefore held fixed.

---

## 6. Correctness gates

These are hard assertions, not warnings. Any breach stops the run.

1. **Sequence identity between reference and refold.** The transplant takes its
   sequence from the library CSV; the refold was folded from the design's own
   sequence. Assert they match for every design.
   This is not defensive padding: BG+LiMPNN sequences were mislabelled as BG+BIF
   in **v2** once before (`aggregate.py`), and v2 is 235 of these 2,400.
2. **Register.** Every design must return a finite `lrmsd_rna` and `lrmsd_prot`
   for both backends. A `None` means the name-agreement gate rejected the
   pairing and the transplant did not do its job; such designs are listed, never
   silently dropped.
3. **Reproduction spot-check.** On the T2 campaign only, also score against the
   BIF design (via `extract_refolds.py --with-design`) and confirm the stored
   campaign values reproduce, as the probe did. This is the single deliberate
   exception to approach A's "refolds only": it re-verifies the "trust the CSV"
   assumption on fresh data, and because `designs/` sits earlier in the same
   stream it costs no extra archive read. If v6 reproduces, the remaining
   campaigns take their campaign-side numbers from the CSV unchecked.
4. **Accounting.** `resolved + fallback + unresolvable + failed == 2,400`,
   printed. The headline rate is quoted over resolved designs with the
   unresolved count stated alongside, never folded in silently.
5. **Chains by composition, never by ID** — design CIF chain lettering is not
   consistent across stacks.

---

## 7. Test ladder

Small tests first; each stage gates the next.

| | what | measures / proves |
|---|---|---|
| **T1** | re-run the existing 20 probe pairs through the new scaled scripts | pure regression: must match the probe's numbers exactly |
| **T2** | **v6 end-to-end** — 133 designs, smallest archive (9.3 GB) | real full-stream time (the probe's 26-37 s is *not* extrapolatable — it broke early on low-hex uuids); validates extract -> ref -> score -> gate |
| **T2b** | the 107 LiMPNN-fallback designs (v2) | separate code path, tested separately |
| **T3** | remaining 5 archives | the run |

T2 is the point at which the extraction estimate becomes measured rather than
assumed. If it lands far off, re-plan T3 before committing.

---

## 8. Risks

| risk | handling |
|---|---|
| extraction slower than estimated | T2 measures before T3 commits; tasks are independent and idempotent |
| a scoring bug found after extraction | costs nothing — `cifs/` is persistent, re-scoring touches no archive |
| minimal ESM mmCIFs crash DockQ's parser | already handled by `dockq_self.py`'s gemmi fallback; proven on 40 structures |
| v2 sequence mislabelling | gate 1 catches it |
| RAM | one inner blob held at a time; probe peaked at 248 MB. Request 8 GB/task |

---

## 9. Success criteria

1. A flip rate over the 2,400, with the resolved / fallback / unresolvable
   accounting stated.
2. A per-design table (`rescore_2400.csv`) giving campaign and corrected
   sc_RMSD and DockQ per backend.
3. All five correctness gates passed, or every exception enumerated.
4. `REPORT.md` stating plainly what the number does **not** cover: designs
   outside the library, and the other four stacks.
