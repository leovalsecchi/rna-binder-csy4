#!/usr/bin/env python3
"""
patch_launch_add_bd1.py — adds BoltzDesign builder to launch.py.

v2: avoids embedded backslash-newline escape problems.
"""

import sys
import re
from pathlib import Path

LAUNCH = Path("/cluster/work/wicky/lvalsecchi/rna_binder/code/scripts/launch.py")

# Build the new function as a list of lines to avoid escaping hell
BD1_BUILDER_LINES = [
    '',
    'def build_boltzdesign_sbatch(config, run_dir, uuid8):',
    '    """Build a BoltzDesign1 sbatch string."""',
    '    slurm = config["slurm"]',
    '    c = config["container"]',
    '    gp = config["generator_params"]',
    '',
    '    boltz_cache = c["boltz_cache"]',
    '    sif = c["sif"]',
    '',
    '    # Build the --<key> <value> CLI args from generator_params (skip work_dir if present)',
    '    cli_args = []',
    '    for k, v in gp.items():',
    '        if k == "work_dir":',
    '            continue',
    '        cli_args.append(f"--{k} {shlex.quote(str(v))}")',
    '    cli_args_str = " ".join(cli_args)',
    '',
    '    # bind list (static + dynamic)',
    '    bind_sep = " \\\\\\n    "',
    '    static_binds = bind_sep.join(f"--bind {b}" for b in c.get("binds", []))',
    '',
    '    return f"""#!/bin/bash',
    "#SBATCH --job-name={slurm.get('job_name', 'bd1')}_{uuid8}",
    "#SBATCH --time={slurm['time']}",
    "#SBATCH --cpus-per-task={slurm['cpus_per_task']}",
    "#SBATCH --mem-per-cpu={slurm['mem_per_cpu']}",
    "#SBATCH --gpus={slurm['gpus']}",
    "#SBATCH --gres={slurm['gres']}",
    "#SBATCH --tmp={slurm['tmp']}",
    '#SBATCH --output={run_dir}/slurm.log',
    '#SBATCH --error={run_dir}/slurm.err',
    '',
    'set -euo pipefail',
    '',
    'RUN_DIR={run_dir}',
    'SIF={sif}',
    'BOLTZ_CACHE={boltz_cache}',
    '',
    'mkdir -p "${{RUN_DIR}}"',
    '',
    '# Ensure boltz cache exists (idempotent)',
    'mkdir -p "${{BOLTZ_CACHE}}/.boltz"',
    'ln -sfn /cluster/work/wicky/params/boltz/boltz1_conf.ckpt "${{BOLTZ_CACHE}}/.boltz/boltz1_conf.ckpt"',
    'ln -sfn /cluster/work/wicky/params/boltz/ccd.pkl "${{BOLTZ_CACHE}}/.boltz/ccd.pkl"',
    '',
    '# Fake HOME for matplotlib and other writes inside the container',
    'FAKE_HOME="${{TMPDIR}}/fake_home"',
    'mkdir -p "${{FAKE_HOME}}"',
    'export MPLCONFIGDIR="${{FAKE_HOME}}/.config/matplotlib"',
    'mkdir -p "${{MPLCONFIGDIR}}"',
    '',
    'echo "================================================================"',
    'echo " BoltzDesign1 launcher (via launch.py)"',
    'echo " Job ID:        ${{SLURM_JOB_ID}}"',
    'echo " UUID8:         {uuid8}"',
    'echo " Node:          ${{SLURMD_NODENAME}}"',
    'echo " GPU:           ${{CUDA_VISIBLE_DEVICES:-unset}}"',
    'nvidia-smi -L 2>/dev/null | sed \'s/^/ GPU detail:   /\' || true',
    'echo " Run dir:       ${{RUN_DIR}}"',
    'echo " Boltz cache:   ${{BOLTZ_CACHE}}"',
    'echo " Fake HOME:     ${{FAKE_HOME}}"',
    'echo " Started:       $(date)"',
    'echo "================================================================"',
    '',
    'module load eth_proxy',
    '',
    'apptainer exec \\\\',
    '    --nv \\\\',
    '    {static_binds} \\\\',
    '    --bind "${{TMPDIR}}:${{TMPDIR}}" \\\\',
    '    --bind "${{BOLTZ_CACHE}}:/cluster/home/${{USER}}" \\\\',
    '    --env HOME="${{FAKE_HOME}}" \\\\',
    '    --env MPLCONFIGDIR="${{MPLCONFIGDIR}}" \\\\',
    '    "${{SIF}}" \\\\',
    '    bash -c "',
    'set -euo pipefail',
    'cp -r /app/LigandMPNN \'${{RUN_DIR}}/LigandMPNN\'',
    'cd /app',
    'python boltzdesign.py {cli_args_str} --gpu_id 0 --work_dir \'${{RUN_DIR}}\'',
    '"',
    '',
    'EXIT_CODE=$?',
    'echo "================================================================"',
    'echo " Finished:   $(date)"',
    'echo " Exit code:  ${{EXIT_CODE}}"',
    'echo "================================================================"',
    'exit ${{EXIT_CODE}}',
    '"""',
    '',
    '',
]

BD1_BUILDER = "\n".join(BD1_BUILDER_LINES) + "\n"


def main():
    if not LAUNCH.exists():
        sys.exit(f"launch.py not found at {LAUNCH}")

    text = LAUNCH.read_text()

    # Idempotency
    if "build_boltzdesign_sbatch" in text:
        print(f"Already patched: build_boltzdesign_sbatch is present in {LAUNCH}")
        return

    # Insert the new builder function right before the BUILDERS dict
    builders_marker = "BUILDERS = {"
    if builders_marker not in text:
        sys.exit("Couldn't find 'BUILDERS = {' anchor in launch.py")

    text = text.replace(builders_marker, BD1_BUILDER + builders_marker, 1)

    # Add boltzdesign entry to the BUILDERS dict
    bg_line_pattern = r'(\s*"boltzgen":\s+build_boltzgen_sbatch,)'
    bd_line = '\n    "boltzdesign":  build_boltzdesign_sbatch,'

    if re.search(bg_line_pattern, text):
        text = re.sub(bg_line_pattern, r'\1' + bd_line, text, count=1)
    else:
        sys.exit("Couldn't find the boltzgen entry in BUILDERS dict")

    # Back up original
    backup = LAUNCH.with_suffix(".py.bak_pre_bd1")
    backup.write_text(LAUNCH.read_text())

    # Write the patched version
    LAUNCH.write_text(text)

    print(f"Patched {LAUNCH}")
    print(f"Backup at {backup}")
    print("Verifying syntax...")

    import py_compile
    try:
        py_compile.compile(str(LAUNCH), doraise=True)
        print("Syntax check: OK")
    except py_compile.PyCompileError as e:
        print(f"SYNTAX ERROR after patch:")
        print(e)
        print(f"\nRestoring backup from {backup}")
        LAUNCH.write_text(backup.read_text())
        sys.exit(1)


if __name__ == "__main__":
    main()
