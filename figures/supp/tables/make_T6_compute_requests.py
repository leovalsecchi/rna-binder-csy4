#!/usr/bin/env python3
"""T6 - per-stage compute requests. Writes T6_compute_requests.{csv,md}.

Run from the repo root:
    ~/.venvs/binder_map/bin/python3 code/supp/tables/make_T6_compute_requests.py

Every scheduler value is read out of the archived production config in the code
archive (~/rna-binder-csy4, tag v1.0-thesis-freeze), never retyped: CONFIGS below
names the file and the parser pulls its `slurm:` block. The configured values are
what the table reports, because campaign submission overrode only the scheduler
fields; where the as-run routing differed, the `as_run` column says so and names
the archived as-run record it came from.

Seven production stages: four generators (the two BoltzGen routes share one
generation job) and three re-prediction backends.
"""
import csv
import re
import sys
from pathlib import Path

import yaml

ARCHIVE = Path.home() / "rna-binder-csy4"
OUT = Path(__file__).resolve().parent

# as-run routing, read from docs/as_run/MIDRUN_V{4,6}_AS_RUN_CONDITIONS_*.md in
# the same archive; the file is named in `as_run_source` so nothing here is bare.
V6 = "docs/as_run/MIDRUN_V6_AS_RUN_CONDITIONS_2026-07-08.md"
V4 = "docs/as_run/MIDRUN_V4_AS_RUN_CONDITIONS_2026-07-05.md"

NO_SCRATCH = "none requested; node-local $TMPDIR only"
REFOLD_SCRATCH = "30 GB (--tmp=30G), requested by the as-run scratch workers and not by the config"

CONFIGS = [
    dict(stage="Backbone generation, RFdiffusion3 + LigandMPNN",
         config="configs/rfd3__multif/default.yaml", scratch=NO_SCRATCH,
         as_run="moved to the public gpu.4h tier with gpu:1 (the gpumem request dropped, eu-g3 nodelist pinned) from midrun_v6",
         as_run_source=V6),
    dict(stage="Backbone generation, BoltzGen (shared by both BoltzGen routes)",
         config="configs/boltzgen__multif/default.yaml", scratch=NO_SCRATCH,
         as_run="gpu.4h with gpu:1 in midrun_v6; gpu.24h and gpu.120h in midrun_v4",
         as_run_source=V6),
    dict(stage="Backbone generation, BoltzDesign1 + LigandMPNN",
         config="configs/boltzdesign__multif/default.yaml", scratch=NO_SCRATCH,
         as_run="gpu.24h and gpu.120h with gpu:1; too long for the 4 h tier at a 6.4 h median",
         as_run_source=V6),
    dict(stage="Backbone generation, HalluDesign + LigandMPNN (de novo)",
         config="configs/halludesign/random_init.yaml", scratch=NO_SCRATCH,
         as_run="gpuhe.24h with gpumem:24g and gpupr.24h with gpumem:40g; sm_80+ only, so never on the public Turing pool",
         as_run_source=V4),
    dict(stage="Re-prediction, AlphaFold3",
         config="configs/refold__af3/default.yaml", scratch=REFOLD_SCRATCH,
         as_run="gpuhe.bulk at 06:00:00 via the node-local scratch worker",
         as_run_source=V6),
    dict(stage="Re-prediction, Boltz-2",
         config="configs/refold__b2/default.yaml", scratch=REFOLD_SCRATCH,
         as_run="gpu.bulk at 06:00:00 via the node-local scratch worker",
         as_run_source=V6),
    dict(stage="Re-prediction, ESMFold2",
         config="configs/refold__esm/default.yaml", scratch=REFOLD_SCRATCH,
         as_run="gpuhe.bulk at 06:00:00 via the node-local scratch worker",
         as_run_source=V6),
]

MEM = re.compile(r"^(\d+(?:\.\d+)?)\s*([GM])$")


def processors(slurm, path):
    """Processors per job, with the field the config actually used."""
    if "cpus_per_task" in slurm:
        return int(slurm["cpus_per_task"]), "cpus_per_task"
    if "ntasks" in slurm:
        return int(slurm["ntasks"]), "ntasks"
    if "ntasks_per_node" in slurm:
        return int(slurm["ntasks_per_node"]), "ntasks_per_node"
    sys.exit(f"no processor field in {path}")


def total_mem(per_cpu, n):
    m = MEM.match(str(per_cpu).strip())
    if not m:
        sys.exit(f"unparsed mem_per_cpu {per_cpu!r}")
    val, unit = float(m.group(1)), m.group(2)
    return f"{val * n:g} {unit}B"


def gpu_request(slurm):
    """The GPU class as requested: a pinned model, or a count plus a memory class."""
    g = slurm.get("gpus")
    gres = slurm.get("gres")
    if isinstance(g, str) and ":" in g:
        out = g
    else:
        out = f"{g} GPU" if g else "NEED: no GPU field"
    if gres:
        out += f" + {gres}"
    return out


rows = []
for spec in CONFIGS:
    path = ARCHIVE / spec["config"]
    slurm = yaml.safe_load(path.read_text())["slurm"]
    n, field = processors(slurm, path)
    rows.append({
        "stage": spec["stage"],
        "config_file": spec["config"],
        "processors": f"{n} ({field})",
        "memory_per_processor": str(slurm["mem_per_cpu"]).strip(),
        "memory_per_job": total_mem(slurm["mem_per_cpu"], n),
        "local_scratch": spec["scratch"],
        "gpu_request": gpu_request(slurm),
        "partition": slurm["partition"],
        "wall_time": str(slurm["time"]).strip(),
        "as_run_deviation": spec["as_run"],
        "as_run_source": spec["as_run_source"],
    })

assert len(rows) == 7, f"expected 7 production stages, built {len(rows)}"

cols = list(rows[0])
with (OUT / "T6_compute_requests.csv").open("w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=cols)
    w.writeheader()
    w.writerows(rows)

show = ["stage", "processors", "memory_per_job", "local_scratch", "gpu_request",
        "partition", "wall_time"]
head = ["Stage", "Processors", "Memory", "Local scratch", "GPU request", "Partition", "Wall time"]
md = ["# T6 - per-stage compute requests", "",
      "Configured requests for the seven production stages, read from the archived",
      "production configs at tag `v1.0-thesis-freeze`. Campaign submission overrode only",
      "the scheduler fields, so these are the configured defaults; the as-run routing,",
      "where it differed, is in the CSV's `as_run_deviation` column and in the notes below.",
      "", "| " + " | ".join(head) + " |", "|" + "|".join(["---"] * len(head)) + "|"]
for r in rows:
    md.append("| " + " | ".join(str(r[c]) for c in show) + " |")
md += ["", "## Configured against as-run", ""]
for r in rows:
    md.append(f"- **{r['stage']}** (`{r['config_file']}`) - {r['as_run_deviation']} "
              f"(`{r['as_run_source']}`).")
md += ["", "Processors are per job: the generation and ESMFold2 configs set `cpus_per_task`, "
       "the HalluDesign config sets `ntasks` and the Boltz-2 config `ntasks_per_node`, "
       "and the bucketers emit `--cpus-per-task=4` for every stage. No generation stage "
       "requested local scratch; the re-prediction workers did, which is what kept the "
       "project inode count flat.", ""]
(OUT / "T6_compute_requests.md").write_text("\n".join(md))
print(f"wrote T6_compute_requests.csv and .md: {len(rows)} rows, {len(cols)} columns")
for r in rows:
    print(f"  {r['stage'][:52]:54} {r['processors']:18} {r['memory_per_job']:7} "
          f"{r['gpu_request']:28} {r['partition']:12} {r['wall_time']}")
