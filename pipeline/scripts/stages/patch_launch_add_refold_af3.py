#!/usr/bin/env python3
"""
patch_launch_add_refold_af3.py — adds AF3 refold support to launch.py.

Detects config_kind: refold_af3 in YAML configs and orchestrates:
  1. Validate source_run exists and has designs/
  2. Create <source_run>/refold_af3/<label>/  (fails if exists)
  3. Write provenance (resolved_config.yaml, run_meta.json, run.sbatch)
  4. Generate sbatch that calls stages/refold_af3.sh worker
  5. Submit

Backup: launch.py.bak_pre_refold_af3
"""

import sys
import re
from pathlib import Path

LAUNCH = Path("/cluster/work/wicky/lvalsecchi/rna_binder/code/scripts/launch.py")

NEW_BUILDER = '''

def build_refold_af3_sbatch(config, refold_dir, source_run):
    """Build sbatch that calls stages/refold_af3.sh worker.

    refold_dir = <source_run>/refold_af3/<label>/
    """
    slurm = config["slurm"]
    c = config["container"]
    rp = config.get("refolder_params", {})
    rna_seq = config["rna_sequence"]
    model_seeds = rp.get("model_seeds", [1])
    model_seeds_json = json.dumps(model_seeds)

    worker = (
        "/cluster/work/wicky/lvalsecchi/rna_binder/code/scripts/stages/"
        "refold_af3.sh"
    )

    return f"""#!/bin/bash
#SBATCH --job-name={slurm.get('job_name', 'refold_af3')}
#SBATCH --partition={slurm['partition']}
#SBATCH --time={slurm['time']}
#SBATCH --cpus-per-task={slurm['cpus_per_task']}
#SBATCH --mem-per-cpu={slurm['mem_per_cpu']}
#SBATCH --gpus={slurm['gpus']}
#SBATCH --gres={slurm['gres']}
#SBATCH --output={refold_dir}/slurm.log
#SBATCH --error={refold_dir}/slurm.err

set -euo pipefail

export AF3_SIF={c['sif']}
export AF3_PARAMS={c['params_dir']}
export MODEL_SEEDS_JSON='{model_seeds_json}'

bash {worker} \\
    {refold_dir} \\
    {source_run} \\
    {rna_seq}
"""

'''


def main():
    if not LAUNCH.exists():
        sys.exit(f"launch.py not found at {LAUNCH}")

    text = LAUNCH.read_text()

    if "build_refold_af3_sbatch" in text:
        print(f"Already patched: build_refold_af3_sbatch present in {LAUNCH}")
        return

    # ----- Insert the new builder function before BUILDERS dict -----
    marker = "BUILDERS = {"
    if marker not in text:
        sys.exit("Couldn't find 'BUILDERS = {' anchor")
    text = text.replace(marker, NEW_BUILDER + marker, 1)

    # ----- Inject meta-config detection in main() -----
    # New block: detect config_kind == "refold_af3" before existing checks.
    refold_af3_block = '''    # ----- Refold AF3 meta-config -----
    if config.get("config_kind") == "refold_af3":
        # source_run must be provided via --override
        if "source_run" not in config or config["source_run"] == "REQUIRED":
            sys.exit(
                "refold_af3 config requires source_run. Pass it via:\\n"
                "  --override source_run=/cluster/.../path/to/source_run_dir"
            )

        source_run = Path(config["source_run"])
        if not source_run.exists():
            sys.exit(f"source_run does not exist: {source_run}")
        designs_dir = source_run / "designs"
        if not designs_dir.exists():
            sys.exit(f"source_run has no designs/ subdir: {source_run}")

        refold_dir = source_run / "refold_af3" / args.label
        if refold_dir.exists():
            sys.exit(
                f"refold dir already exists: {refold_dir}\\n"
                f"Pick a different --label, or remove the existing dir."
            )

        sbatch_str = build_refold_af3_sbatch(config, refold_dir, source_run)

        if args.dry_run:
            print(f"# Would create: {refold_dir}")
            print(f"# Source designs: {designs_dir}")
            print(f"# Label: {args.label}")
            print("# --- sbatch script: ---")
            print(sbatch_str)
            return

        refold_dir.mkdir(parents=True, exist_ok=False)

        # Write provenance
        (refold_dir / "resolved_config.yaml").write_text(
            yaml.dump(config, default_flow_style=False)
        )
        sbatch_path = refold_dir / "run.sbatch"
        sbatch_path.write_text(sbatch_str)
        sbatch_path.chmod(0o755)

        result = subprocess.run(
            ["sbatch", str(sbatch_path)],
            capture_output=True, text=True, check=True,
        )
        jobid_match = re.search(r"(\\d+)", result.stdout)
        jobid = jobid_match.group(1) if jobid_match else "unknown"

        meta = {
            "jobid": jobid,
            "label": args.label,
            "submitted": datetime.now().isoformat(),
            "host": os.uname().nodename,
            "user": os.environ.get("USER", "unknown"),
            "config_kind": "refold_af3",
            "refolder": "AF3",
            "target": config["target_id"],
            "source_run": str(source_run),
            "config_path": str(args.config_yaml.resolve()),
            "overrides": overrides,
        }
        (refold_dir / "run_meta.json").write_text(json.dumps(meta, indent=2))

        print(f"Submitted AF3 refold: jobid={jobid}")
        print(f"Source:  {source_run}")
        print(f"Refold:  {refold_dir}")
        return

'''

    # Inject right after the multif block. Anchor on the closing of multif's "return"
    # followed by the comment. Use the "Normal (non-meta)" anchor.
    insert_marker = "    # ----- Normal (non-meta) config handling -----"
    if insert_marker not in text:
        sys.exit(f"Couldn't find '{insert_marker}' anchor in main()")
    text = text.replace(insert_marker, refold_af3_block + insert_marker, 1)

    backup = LAUNCH.with_suffix(".py.bak_pre_refold_af3")
    backup.write_text(LAUNCH.read_text())
    LAUNCH.write_text(text)

    print(f"Patched {LAUNCH}")
    print(f"Backup at {backup}")
    print("Verifying syntax...")

    import py_compile
    try:
        py_compile.compile(str(LAUNCH), doraise=True)
        print("Syntax check: OK")
    except py_compile.PyCompileError as e:
        print("SYNTAX ERROR after patch:")
        print(e)
        print(f"Restoring backup from {backup}")
        LAUNCH.write_text(backup.read_text())
        sys.exit(1)


if __name__ == "__main__":
    main()
