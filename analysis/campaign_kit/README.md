# campaign_kit — YAML-driven, probe-calibrated cluster campaigns (SCAFFOLD)

Hand this to a labmate: they bring their own models + target, and drive the cluster at full
power without inheriting anyone's head-knowledge. Full design: `LAUNCHER_LOOPER_ARCHITECTURE.md`.

## Layers
- `cluster/euler.yaml`  ① static cluster inventory (pools, arch, VRAM, QOS caps, tiers) — set once.
- `models/*.model.yaml` ② one file per generator/refolder — container, entrypoint, params, I/O.
- `campaigns/*.campaign.yaml` ③ per campaign — target, stacks, counts, seeds, refold backends.
- `compat/euler.compat.yaml`  ← **probe output**: measured VRAM/runtime/where-each-model-runs.
  Everything downstream (planner, bucketer, looper, gpu_watcher) reads this one file.

## First-time setup (once per user)
```bash
cp site.yaml.example site.yaml       # then edit: user / rb_root / shared  (your cluster paths)
```
`site.yaml` is the ONLY account-specific file; every `${...}` in cluster/*.yaml and models/*.yaml
resolves from it, so nobody inherits anyone else's absolute paths. It is gitignored.

## Handoff workflow
```bash
export PATH="$PWD/bin:$PATH"          # or alias: campaign=bin/campaign
campaign models new mymodel          # author models/mymodel.model.yaml  (or reuse an example)
campaign probe --model mymodel       # measures it on every compatible pool -> compat.yaml
campaign new myrun                   # scaffold campaigns/myrun.campaign.yaml, edit target/stacks
campaign plan  myrun --dry-run       # see the exact submission plan
campaign doctor myrun                # guardrail preflight (must pass)
campaign launch myrun                # submit
campaign loop  myrun start           # auto-generate + start the refold/HD loops (via loopctl)
campaign status myrun                # tally / array_summary / campaign_status
```

## Status (this scaffold)
- **Working:** `campaign help|new|models`, the 3 YAML schemas + templates, filled `cluster/euler.yaml`,
  `compat/euler.compat.yaml` schema, `GUARDRAILS.md`, and `lib/ck/probe.py:read_peak_vram_util()`
  (the core sacct TRESUsageInMax reader).
- **Stubbed (print intent, don't execute):** `probe|plan|doctor|launch|loop|status` and the
  `lib/ck/{planner,doctor}.py` internals. Build order: probe → planner → cutover → monitors (see the
  design doc §9).

## Reuses (does not reinvent)
`code/scripts/launch.py` builders (executor) · `ops/loops` + `loopctl` (looper) ·
`campaign_status.py`/`array_summary.sh`/`recent_jobs.sh` (monitors) · the `tar.zst` archive lib ·
`ops/watchers/gpu_watcher.py` `parse_nodes()`/`class_runtimes()` (probe primitives).
