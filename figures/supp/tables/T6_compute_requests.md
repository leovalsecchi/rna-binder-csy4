# T6 - per-stage compute requests

Configured requests for the seven production stages, read from the archived
production configs at tag `v1.0-thesis-freeze`. Campaign submission overrode only
the scheduler fields, so these are the configured defaults; the as-run routing,
where it differed, is in the CSV's `as_run_deviation` column and in the notes below.

| Stage | Processors | Memory | Local scratch | GPU request | Partition | Wall time |
|---|---|---|---|---|---|---|
| Backbone generation, RFdiffusion3 + LigandMPNN | 4 (cpus_per_task) | 20 GB | none requested; node-local $TMPDIR only | rtx_4090:1 | gpuhe.4h | 04:00:00 |
| Backbone generation, BoltzGen (shared by both BoltzGen routes) | 4 (cpus_per_task) | 64 GB | none requested; node-local $TMPDIR only | rtx_4090:1 | gpuhe.24h | 23:30:00 |
| Backbone generation, BoltzDesign1 + LigandMPNN | 4 (cpus_per_task) | 20 GB | none requested; node-local $TMPDIR only | rtx_4090:1 | gpuhe.24h | 23:30:00 |
| Backbone generation, HalluDesign + LigandMPNN (de novo) | 4 (ntasks) | 28 GB | none requested; node-local $TMPDIR only | a100_80gb:1 | gpuhe.24h | 23:30:00 |
| Re-prediction, AlphaFold3 | 4 (cpus_per_task) | 32 GB | 30 GB (--tmp=30G), requested by the as-run scratch workers and not by the config | 1 GPU + gpumem:24g | gpuhe.4h | 04:00:00 |
| Re-prediction, Boltz-2 | 4 (ntasks_per_node) | 32 GB | 30 GB (--tmp=30G), requested by the as-run scratch workers and not by the config | 1 GPU + gpumem:8g | gpu.4h | 04:00:00 |
| Re-prediction, ESMFold2 | 4 (cpus_per_task) | 32 GB | 30 GB (--tmp=30G), requested by the as-run scratch workers and not by the config | 1 GPU + gpumem:24g | bulk | 12:00:00 |

## Configured against as-run

- **Backbone generation, RFdiffusion3 + LigandMPNN** (`configs/rfd3__multif/default.yaml`) - moved to the public gpu.4h tier with gpu:1 (the gpumem request dropped, eu-g3 nodelist pinned) from midrun_v6 (`docs/as_run/MIDRUN_V6_AS_RUN_CONDITIONS_2026-07-08.md`).
- **Backbone generation, BoltzGen (shared by both BoltzGen routes)** (`configs/boltzgen__multif/default.yaml`) - gpu.4h with gpu:1 in midrun_v6; gpu.24h and gpu.120h in midrun_v4 (`docs/as_run/MIDRUN_V6_AS_RUN_CONDITIONS_2026-07-08.md`).
- **Backbone generation, BoltzDesign1 + LigandMPNN** (`configs/boltzdesign__multif/default.yaml`) - gpu.24h and gpu.120h with gpu:1; too long for the 4 h tier at a 6.4 h median (`docs/as_run/MIDRUN_V6_AS_RUN_CONDITIONS_2026-07-08.md`).
- **Backbone generation, HalluDesign + LigandMPNN (de novo)** (`configs/halludesign/random_init.yaml`) - gpuhe.24h with gpumem:24g and gpupr.24h with gpumem:40g; sm_80+ only, so never on the public Turing pool (`docs/as_run/MIDRUN_V4_AS_RUN_CONDITIONS_2026-07-05.md`).
- **Re-prediction, AlphaFold3** (`configs/refold__af3/default.yaml`) - gpuhe.bulk at 06:00:00 via the node-local scratch worker (`docs/as_run/MIDRUN_V6_AS_RUN_CONDITIONS_2026-07-08.md`).
- **Re-prediction, Boltz-2** (`configs/refold__b2/default.yaml`) - gpu.bulk at 06:00:00 via the node-local scratch worker (`docs/as_run/MIDRUN_V6_AS_RUN_CONDITIONS_2026-07-08.md`).
- **Re-prediction, ESMFold2** (`configs/refold__esm/default.yaml`) - gpuhe.bulk at 06:00:00 via the node-local scratch worker (`docs/as_run/MIDRUN_V6_AS_RUN_CONDITIONS_2026-07-08.md`).

Processors are per job: the generation and ESMFold2 configs set `cpus_per_task`, the HalluDesign config sets `ntasks` and the Boltz-2 config `ntasks_per_node`, and the bucketers emit `--cpus-per-task=4` for every stage. No generation stage requested local scratch; the re-prediction workers did, which is what kept the project inode count flat.
