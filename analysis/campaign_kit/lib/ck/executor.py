"""ck.executor — turn a plan into bucketed array sbatch scripts (and optionally submit).

Generalizes bucket_midrun_v<N>.py. The SBATCH HEADER is generated entirely from the plan
(the routing the probe measured): correct --cpus-per-task from the start (no post-hoc
--ntasks rebucket), --requeue/--open-mode=append, blackhole --exclude, node-local $TMPDIR.

The per-task BODY is now real and generic:
  * standard env for the unit (CK_UUID, CK_SEED, CK_INDEX, CK_OUT, CK_TARGET_PDB, CK_RNA);
  * the model's RUN CONTRACT — `model.run_contract` if given, else a default apptainer
    invocation built from the model spec (container + entrypoint + params + seed + --out);
  * generic PACK + VERIFY + PUBLISH: tar.zst -19 -T0 of the node-local output, `zstd -t`
    verified, then moved to the campaign pool as ONE file (the inode + integrity guardrails,
    now code not lore).
A new model usually only needs to set `run_contract` (one line) — or nothing, if its
entrypoint takes `--seed`/`--out`. Safe by default: nothing is submitted unless submit=True.
"""
from __future__ import annotations
import os
import subprocess

TIER_TIME = {"4h": "04:00:00", "24h": "23:30:00", "120h": "48:00:00", "bulk": "06:00:00"}
ARRAY_LIMIT = 200          # %concurrency cap; the QOS GPU cap is the real throttle


def _slug(s):
    return s.replace(":", "_").replace("/", "_")


def _flags(model):
    fs = model.get("flag_style") or {}
    bools = set(fs.get("boolean_flags") or [])
    out = []
    for k, v in (model.get("params") or {}).items():
        out.append(f"--{k}" if k in bools else f"--{k} {v}")
    return " ".join(out)


def _default_run_contract(model):
    c = model.get("container", {}) or {}
    binds = " ".join(f"--bind {b}" for b in (c.get("binds") or []))
    envs = " ".join(f'--env {k}="{v}"' for k, v in (c.get("envs") or {}).items())
    app = (f"apptainer {c.get('apptainer_cmd', 'run')} {'--nv ' if c.get('nv') else ''}"
           f"{binds} {envs} {c.get('sif', '<sif>')} {model.get('entrypoint', '<entrypoint>')}")
    seed = f'--{model.get("seed_flag", "seed")} "$CK_SEED"' if model.get("per_unit_seed") else ""
    wd = f'cd {model["workdir"]} && ' if model.get("workdir") else ""
    return f'{wd}{app} {_flags(model)} {seed} --out "$CK_OUT"'


def render_sbatch(w, campaign, model, inventory, rb_root, logdir):
    name = f"{campaign['name']}_{_slug(w['label'])}"
    tgt = (campaign.get("target") or {})
    label = _slug(w["label"])
    publish = f"{rb_root}/campaigns/{tgt.get('id')}/{label}/{campaign['name']}"
    seed_base = (campaign.get("stacks") or [{}])[0].get("seed_base", 1)
    excl = ",".join(inventory.get("blackhole_excludes") or [])
    njobs = w.get("n_jobs") or 1
    run_contract = model.get("run_contract") or _default_run_contract(model)

    header = [x for x in [
        "#!/usr/bin/env bash",
        f"#SBATCH --job-name={name}",
        f"#SBATCH --partition={w['partition']}",
        f"#SBATCH --qos={w['qos']}",
        f"#SBATCH --time={TIER_TIME.get(w['walltime_tier'], '04:00:00')}",
        f"#SBATCH --cpus-per-task={w['cpus_per_task']}",
        f"#SBATCH --mem-per-cpu={w['mem_per_cpu']}",
        "#SBATCH --gpus=1",
        f"#SBATCH --gres={w['gres']}" if w.get("gres") else "",
        f"#SBATCH --array=0-{njobs-1}%{ARRAY_LIMIT}",
        f"#SBATCH --exclude={excl}" if excl else "",
        "#SBATCH --requeue",
        "#SBATCH --open-mode=append",
        f"#SBATCH --output={logdir}/{name}_%A_%a.log",
    ] if x != ""]

    body = f'''
set -euo pipefail
module load stack/2024-06 python/3.12.8 2>/dev/null || true

# node-local scratch — extraction/outputs NEVER touch /cluster (inode-flood guardrail)
export TMPDIR="${{TMPDIR:-/scratch/$USER/${{SLURM_JOB_ID}}}}"; mkdir -p "$TMPDIR"

# --- per-unit identity + inputs ---
export CK_INDEX="$SLURM_ARRAY_TASK_ID"
export CK_UUID="$(printf '{campaign['name']}-{label}-%06d' "$CK_INDEX")"
export CK_SEED="$(( {seed_base} + CK_INDEX ))"
export CK_TARGET_PDB="{tgt.get('pdb', '')}"
export CK_RNA="{tgt.get('rna_sequence', '')}"
export CK_OUT="$TMPDIR/$CK_UUID"; mkdir -p "$CK_OUT"
PUBLISH="{publish}"; mkdir -p "$PUBLISH"
[ -f "$PUBLISH/$CK_UUID.tar.zst" ] && {{ echo "already done: $CK_UUID"; exit 0; }}   # idempotent/resumable

# --- run ONE unit (model run contract; must write into $CK_OUT) ---
{run_contract}

# --- pack + verify + publish (generic): node-local -> Lustre as ONE verified file ---
if [ -n "$(ls -A "$CK_OUT" 2>/dev/null)" ]; then
  ( cd "$CK_OUT" && tar -cf - . | zstd -q -19 -T0 -o "$TMPDIR/$CK_UUID.tar.zst" )
  zstd -tq "$TMPDIR/$CK_UUID.tar.zst"
  mv -f "$TMPDIR/$CK_UUID.tar.zst" "$PUBLISH/$CK_UUID.tar.zst"
  echo "published $PUBLISH/$CK_UUID.tar.zst"
else
  echo "ERROR: $CK_OUT is empty — run_contract produced no output" >&2; exit 1
fi
'''
    return "\n".join(header) + "\n" + body


def calibration_sbatch(model, pool, tier, units, logdir, target_pdb="", rna="", gres_override=None):
    """A tiny pool-pinned single job that runs `units` of the model, for the cold-probe.
    Pinned to the pool via its partition; requests the pool's FULL card VRAM by default so
    the probe can measure a memory-elastic model's true peak — but `gres_override` (e.g.
    gpumem:8g) lets a small model schedule on more (incl. lower-VRAM) idle nodes. Params
    come from probe.calib_params (a small override) so the calibration is quick."""
    card = int(pool.get("vram_g") or 24)
    gres = gres_override or f"gpumem:{card}g"
    tiern = "bulk" if tier == "bulk" else tier.rstrip("h")
    name = f"calib_{model['name']}_{pool['name']}"
    m = dict(model)
    m["params"] = {**(model.get("params") or {}),
                   **((model.get("probe") or {}).get("calib_params") or {})}
    rc = m.get("run_contract") or _default_run_contract(m)
    header = [x for x in [
        "#!/usr/bin/env bash",
        f"#SBATCH --job-name={name}",
        f"#SBATCH --partition={pool['name']}.{tier}",
        f"#SBATCH --qos={pool.get('qos_base', 'es_wicky/' + pool['name'])}/{tiern}",
        f"#SBATCH --time={TIER_TIME.get(tier, '04:00:00')}",
        "#SBATCH --cpus-per-task=2",       # calibration is tiny — keep the QOS-memory footprint small
        "#SBATCH --mem-per-cpu=4G",
        "#SBATCH --gpus=1",
        f"#SBATCH --gres={gres}",
        "#SBATCH --requeue",
        "#SBATCH --open-mode=append",
        f"#SBATCH --output={logdir}/{name}_%j.log",
    ] if x != ""]
    body = f'''
set -euo pipefail
module load stack/2024-06 python/3.12.8 2>/dev/null || true
export TMPDIR="${{TMPDIR:-/scratch/$USER/${{SLURM_JOB_ID}}}}"; mkdir -p "$TMPDIR"
export CK_UUID="calib" CK_SEED="1" CK_INDEX="0"
export CK_TARGET_PDB="{target_pdb}" CK_RNA="{rna}"
export CK_OUT="$TMPDIR/calib"; mkdir -p "$CK_OUT"
echo "[calibration] model={model['name']} pool={pool['name']} units={units} card={card}G"

# run {units} unit(s) of the model — the probe reads this job's peak VRAM + runtime from sacct
{rc}

echo "[calibration] done; produced $(ls -A "$CK_OUT" 2>/dev/null | wc -l) output item(s)"
'''
    return name, "\n".join(header) + "\n" + body


def write_scripts(plan, campaign, models, inventory, rb_root, outdir, logdir, submit=False):
    """Write one array sbatch per (non-error) workload. Returns [(path, sbatch_cmd, jobid|None)]."""
    os.makedirs(outdir, exist_ok=True)
    os.makedirs(logdir, exist_ok=True)
    results = []
    for w in plan:
        if w.get("error"):
            continue
        model = models.get(w["model"], {"name": w["model"]})
        path = os.path.join(outdir, f"{_slug(w['label'])}.sbatch")
        with open(path, "w") as fh:
            fh.write(render_sbatch(w, campaign, model, inventory, rb_root, logdir))
        jobid = None
        if submit:
            out = subprocess.run(["sbatch", "--parsable", path], capture_output=True, text=True)
            jobid = out.stdout.strip() or f"ERR:{out.stderr.strip()[:80]}"
        results.append((path, f"sbatch {path}", jobid))
    return results
