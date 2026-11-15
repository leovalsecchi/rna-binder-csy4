# Structure file locations — generated designs + refolds

> **OUTDATED for current runs.** This describes the OLD reference run
> `midrun_v2_f2653389`: gzipped `.cif.gz` designs + the deeply-nested AF3
> `outputs/<DESIGN>/<lc>/` layout. Since then the pipeline **canonicalizes**
> (flat `<DESIGN>.cif`, one AF3 `_AF3_top.cif` per design under
> `refold_af3/<leaf>/designs/`) and **packs to tar.zst** — multi-run designs are
> per-uuid slices `designs/<uuid>.tar.zst`, refold outputs are
> `<leaf>/designs.tar.zst`. See `~/rna_binder_archiving/phase2/AS_BUILT.md`.

Reference run: campaign `csy4_crrna_free`, stack `rfd3__ligandmpnn`, run
`midrun_v2_f2653389` (200 designs).

**Root** (`<RUN>`):
- Euler:  `/cluster/work/wicky/lvalsecchi/rna_binder/campaigns/<CAMPAIGN>/<STACK>/<RUN>`
- Mount:  `/Volumes/eulermount/lvalsecchi/rna_binder/campaigns/<CAMPAIGN>/<STACK>/<RUN>`

Placeholders: `<DESIGN>` = e.g. `csy4free_RFD3_LiMPNN_len80_f2653389_design0`;
`<lc>` = the **lowercased** design name (`csy4free_rfd3_limpnn_len80_f2653389_design0`).

---

## 1. Generated designs (the reference pose)
```
<RUN>/designs/<DESIGN>.cif.gz          # gzipped mmCIF  (+ <DESIGN>.json sidecar)
```
- 1 per design — **200** files. Original mixed case. Chains: **A = RNA, B = protein**.

## 2. AF3 refold — msa-off and msa-on
```
<RUN>/refold_af3/<MSA>/outputs/<DESIGN>/<lc>/<lc>_model.cif      # top-ranked model
<RUN>/refold_af3/<MSA>/outputs/<DESIGN>/<lc>/seed-1_sample-{0..4}/   # 5 seed samples
```
- `<MSA>` ∈ `{msaoff, msaon}`. 1 top-ranked `_model.cif` per design per MSA →
  **200** (msaoff) + **200** (msaon).
- Plus `ranking_scores.csv`, `*_confidences.json` in the `<lc>` dir.
- ⚠ The per-design dir name is the **lowercased** design name (`RFD3→rfd3`,
  `LiMPNN→limpnn`). Chains: **A = protein, B = RNA** (swapped vs designs/).

## 3. Boltz-2 (B2) refold — msa-off only
```
<RUN>/refold_b2/msaoff/designs/<DESIGN>_B2_sample{0..4}.cif      # 5 samples, flat dir
                                <DESIGN>_B2_sample{N}_confidence.json
```
- **5 samples per design**, flat directory, original mixed case →
  **1000** cif files (200 × 5). No `msaon`. No single representative `_model.cif`
  (pick "best" via the `*_confidence.json`). Chains: **A = protein, B = RNA**.

---

## Counts (this run)
| set | path | count |
|---|---|---|
| designs | `designs/*.cif.gz` | 200 |
| AF3 msaoff | `refold_af3/msaoff/outputs/*/*/*_model.cif` | 200 |
| AF3 msaon | `refold_af3/msaon/outputs/*/*/*_model.cif` | 200 |
| B2 msaoff | `refold_b2/msaoff/designs/*_B2_sample*.cif` | 1000 (200×5) |

## Naming gotchas (relevant to the conventions cleanup)
1. **Three different layouts**: designs = gzipped flat; AF3 = deeply nested + name
   used twice; B2 = flat with `_B2_sample{N}` suffix.
2. **Case inconsistency**: AF3 nested dirs/filenames are **lowercased**; `designs/`
   and B2 keep mixed case → any design↔refold join must be case-insensitive.
3. **Representative model differs**: AF3 has one `_model.cif` (best) + 5 seed dirs;
   B2 has 5 flat sample cifs and *no* single best file.
4. **MSA coverage differs**: AF3 has off+on; B2 has off only → 3 refold variants
   total (af3-on, af3-off, b2-off).
5. **Chain-letter swap**: designs use RNA=A/protein=B; both refolds use
   protein=A/RNA=B.
