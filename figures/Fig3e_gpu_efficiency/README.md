# Fig3e — GPU cost against design quality, by design route

`Fig3e_gpu_efficiency.{pdf,png}` — one point per design route. The folder name is
historical: this is **panel f** of figure 3, the letter `code/Fig3_assembly/`
(layout V2) and `captions/Fig3.md` both give it.

- **x** designs per GPU-hour = 1 / (GPU-h per design), **linear from zero**:
  0–300, ticks every 100. Log until 2026-11-17, changed at Leo's request.
- **y** pass rate = share of that route's own designs reaching the 55,706 set
- **no shading and no reference line.** The green/grey quadrant tints and the
  campaign-average pass-rate edge were removed from the panel on 2026-11-17; the
  five marks and their direct labels carry the comparison on their own, so the
  campaign-wide rate (0.797%) is no longer drawn anywhere on the canvas.

## Caption draft

> **GPU cost against design quality for the five design routes.** Throughput (x)
> is designs per GPU-hour measured on campaign-identical jobs pinned to NVIDIA
> RTX 4090s, covering generation and all three refold backends (AlphaFold3
> 3.0.0, Boltz-2 2.2.1, ESMFold2). Pass rate (y) is the fraction of each route's
> own designs that satisfy the campaign filter, which admitted 55,706 of
> 6,992,515 designs overall (0.80%). Neither that rate nor any quadrant shading
> is drawn on the panel. Refolding costs the same for every route
> (0.00364–0.00370 GPU-h per design, 1.6% spread), so the spread on the x axis
> is generation. The two BoltzGen routes share one allocation; each is charged
> the full shared backbone stage plus its own inverse folding. Combining both
> axes, one design passing the filter costs 0.22 GPU-h with RFD3 + LigandMPNN
> and 29 GPU-h with HalluDesign + LigandMPNN, a 133-fold range.

## Numbers

| design route | designs/GPU-h | pass rate | GPU-h per passing design |
|---|--:|--:|--:|
| RFD3 + LigandMPNN | 250 | 1.829% | 0.22 |
| BoltzGen + LigandMPNN | 233 | 0.971% | 0.44 |
| BoltzGen + BoltzIF | 179 | 0.611% | 0.92 |
| BoltzDesign + LigandMPNN | 101 | 0.224% | 4.41 |
| HalluDesign + LigandMPNN | 44 | 0.078% | 29.11 |

## Inputs

`inputs/Fig3e_gpu_efficiency.csv` — its header records the provenance of every
column. GPU-h from the `gpubench_v1` benchmark on Euler
(`rna_binder/bench_gpueff/`, 2026-11-12, see its README for the protocol and the
three measurement traps it had to avoid). Campaign counts are the `old_bif`
variant of `cascade_by_stack.csv`, i.e. BG+BIF scored on the stored
re-prediction values the campaign actually gated on — that variant is the one
summing to exactly 55,706.

## Two things worth knowing before reusing these numbers

1. **Refolding is measured on an equal 600-design pool per design route, not per uuid.**
   AF3 and ESM amortise one model load across `--input_dir`, and HalluDesign only
   makes 150 designs per uuid, so a per-uuid refold would have given it a
   quarter-size batch and inflated its cost for a reason unrelated to the route.
2. **BoltzGen is charged the full shared backbone stage to each inverse folder.**
   Measured stage split from a timestamped run: backbone 16.7%, BIF 72.6%,
   LigandMPNN 10.7% — BIF costs 6.8x LigandMPNN, because it runs a Boltz
   prediction per design. A 50/50 split of the whole job would put both routes at
   the same x, which is wrong by ~1.9x in opposite directions.
