#!/usr/bin/env python3
"""
patch_launch_add_refold_b2.py — adds Boltz-2 refold support to launch.py.

Detects config_kind: refold_b2 in YAML configs and orchestrates the same
way as refold_af3: validate source_run, create refold dir, write provenance,
submit sbatch.

Backup: launch.py.bak_pre_refold_b2
"""

import sys
import py_compile
from pathlib import Path

LAUNCH = Path("/cluster/work/wicky/lvalsecchi/rna_binder/code/scripts/launch.py")

NEW_BUILDER = '''

def build_refold_b2_sbatch(config, refold_dir, source_run):
    """Build sbatch that calls stages/refold_b2.sh worker.

    refold_dir = <source_run>/refold_b2/<label>/
    """
    slurm = config["slurm"]
    c = config["container"]
    rp = config.get("refolder_params", {})
    rna_seq = config["rna_sequence"]
    diffusion_samples = rp.get("diffusion_samples", 5)
    output_format = rp.get("output_format", "mmcif")
    write_full_pae = "true" if rp.get("write_full_pae", True) else "false"

    worker = (
        "/cluster/work/wicky/lvalsecchi/rna_binder/code/scripts/stages/"
        "refold_b2.sh"
    )

    return f"""#!/bin/bash
#SBATCH --job-name={slurm.get('job_name', 'refold_b2')}
#SBATCH --partition={slurm['partition']}
#SBATCH --time={slurm['time']}
#SBATCH --nodes={slurm.get('nodes', 1)}
#SBATCH --ntasks-per-node={slurm.get('ntasks_per_node', 4)}
#SBATCH --mem-per-cpu={slurm['mem_per_cpu']}
#SBATCH --gpus={slurm['gpus']}
#SBATCH --gres={slurm['gres']}
#SBATCH --output={refold_dir}/slurm.log
#SBATCH --error={refold_dir}/slurm.err

set -euo pipefail

export B2_SIF={c['sif']}
export BOLTZ_CACHE={c['cache_dir']}
export B2_DIFFUSION_SAMPLES={diffusion_samples}
export B2_OUTPUT_FORMAT={output_format}
export B2_WRITE_FULL_PAE={write_full_pae}

bash {worker} \\
    {refold_dir} \\
    {source_run} \\
    {rna_seq}
"""

'''

NEW_DISPATCH = '''    # ----- Refold B2 meta-config -----
    if config.get("config_kind") == "refold_b2":
        if "source_run" not in config or config["source_run"] == "REQUIRED":
            sys.exit(
                "refold_b2 config requires source_run. Pass it via:\\n"
                "  --override source_run=/cluster/.../path/to/source_run_dir"
            )

        source_run = Path(config["source_run"])
        if not source_run.exists():
            sys.exit(f"source_run does not exist: {source_run}")
        designs_dir = source_run / "designs"
        if not designs_dir.exists():
            sys.exit(f"source_run has no designs/ subdir: {source_run}")

        refold_dir = source_run / "refold_b2" / args.label
        if refold_dir.exists():
            sys.exit(
                f"refold dir already exists: {refold_dir}\\n"
                f"Pick a different --label, or remove the existing dir."
            )

        sbatch_str = build_refold_b2_sbatch(config, refold_dir, source_run)

        if args.dry_run:
            print(f"# Would create: {refold_dir}")
            print(f"# Source designs: {designs_dir}")
            print(f"# Label: {args.label}")
            print("# --- sbatch script: ---")
            print(sbatch_str)
            return

        refold_dir.mkdir(parents=True, exist_ok=False)

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
            "config_kind": "refold_b2",
            "refolder": "B2",
            "target": config["target_id"],
            "source_run": str(source_run),
            "config_path": str(args.config_yaml.resolve()),
            "overrides": overrides,
        }
        (refold_dir / "run_meta.json").write_text(json.dumps(meta, indent=2))

        print(f"Submitted B2 refold: jobid={jobid}")
        print(f"Source:  {source_run}")
        print(f"Refold:  {refold_dir}")
        return

'''


def main():
    if not LAUNCH.exists():
        sys.exit(f"launch.py not found at {LAUNCH}")

    text = LAUNCH.read_text()

    if "build_refold_b2_sbatch" in text:
        print(f"Already patched: build_refold_b2_sbatch present in {LAUNCH}")
        return

    # Insert builder function before BUILDERS dict
    marker = "BUILDERS = {"
    if marker not in text:
        sys.exit("Couldn't find 'BUILDERS = {' anchor")
    text = text.replace(marker, NEW_BUILDER + marker, 1)

    # Insert dispatch block right before "Normal (non-meta)" handling
    insert_marker = "    # ----- Normal (non-meta) config handling -----"
    if insert_marker not in text:
        sys.exit(f"Couldn't find '{insert_marker}' anchor in main()")
    text = text.replace(insert_marker, NEW_DISPATCH + insert_marker, 1)

    backup = LAUNCH.with_suffix(".py.bak_pre_refold_b2")
    backup.write_text(LAUNCH.read_text())
    LAUNCH.write_text(text)

    print(f"Patched {LAUNCH}")
    print(f"Backup at {backup}")
    print("Verifying syntax...")

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
