#!/usr/bin/env python3
"""patch_launch_add_autorefold.py — add --auto-refold flag to launch.py.

When --auto-refold is passed:
  After the design job is submitted (jobid=N), build refold sbatches for
  the AF3 and B2 refold configs, with `#SBATCH --dependency=afterok:N`,
  and submit them. Refolds will start only when the design job succeeds.

For single-stack generators (rfd3, halludesign, boltzgen, boltzdesign),
this submits 2 dependent jobs (AF3 + B2).

For multif (boltzgen__ligandmpnn + boltzgen__boltzif sharing backbones),
this submits 4 dependent jobs (AF3+B2 for each of the 2 output runs).

Backup: launch.py.bak_pre_autorefold
"""

import py_compile
from pathlib import Path

LAUNCH = Path("/cluster/work/wicky/lvalsecchi/rna_binder/code/scripts/launch.py")


# ============================================================
# PATCH 1: add --auto-refold CLI arg
# ============================================================
OLD_CLI = '''    ap.add_argument("--dry-run", action="store_true",
                    help="Print sbatch but don't submit")
    args = ap.parse_args()'''

NEW_CLI = '''    ap.add_argument("--dry-run", action="store_true",
                    help="Print sbatch but don't submit")
    ap.add_argument("--auto-refold", action="store_true",
                    help="After design submission, auto-queue AF3+B2 refold jobs "
                         "with --dependency=afterok on the design jobid")
    args = ap.parse_args()'''


# ============================================================
# PATCH 2: add the helper function _submit_dependent_refolds()
# Inserted before def main()
# ============================================================
HELPER_FUNC = '''def _submit_dependent_refolds(design_jobid, source_run_dirs, label):
    """Submit AF3 + B2 refold jobs that depend on design_jobid completing.

    Args:
        design_jobid: SLURM jobid of the design job
        source_run_dirs: list of source run dirs to refold (1 for single-stack, 2 for multif)
        label: refold label (reuses the design label by default)

    Returns:
        list of submitted refold jobids
    """
    import re as _re_local, yaml as _yaml_local
    refold_configs = {
        "af3": Path("/cluster/work/wicky/lvalsecchi/rna_binder/code/configs/refold__af3/default.yaml"),
        "b2":  Path("/cluster/work/wicky/lvalsecchi/rna_binder/code/configs/refold__b2/default.yaml"),
    }
    refold_jobids = []
    for source_run in source_run_dirs:
        for refold_name, refold_cfg_path in refold_configs.items():
            refold_cfg = _yaml_local.safe_load(refold_cfg_path.read_text())
            refold_cfg["source_run"] = str(source_run)
            # Build sbatch via the appropriate builder
            kind = refold_cfg["config_kind"]
            if kind == "refold_af3":
                builder = build_refold_af3_sbatch
            elif kind == "refold_b2":
                builder = build_refold_b2_sbatch
            else:
                print(f"  skip: unknown refold kind {kind}")
                continue
            refold_dir = source_run / f"refold_{refold_name}" / label
            sbatch_str = builder(refold_cfg, refold_dir, source_run)
            # Inject --dependency=afterok:<jobid> after the #SBATCH lines
            dep_line = f"#SBATCH --dependency=afterok:{design_jobid}"
            sbatch_str = sbatch_str.replace(
                "set -euo pipefail",
                f"{dep_line}\\nset -euo pipefail",
                1,
            )
            # Create dir + write + submit
            refold_dir.mkdir(parents=True, exist_ok=False)
            (refold_dir / "resolved_config.yaml").write_text(_yaml_local.dump(refold_cfg))
            sbatch_path = refold_dir / "run.sbatch"
            sbatch_path.write_text(sbatch_str)
            sbatch_path.chmod(0o755)
            res = subprocess.run(
                ["sbatch", str(sbatch_path)],
                capture_output=True, text=True, check=True,
            )
            m = _re_local.search(r"(\\d+)", res.stdout)
            jobid = m.group(1) if m else "unknown"
            print(f"  dep refold {refold_name.upper()}: jobid={jobid} → {refold_dir}")
            refold_jobids.append(jobid)
    return refold_jobids


'''


# ============================================================
# PATCH 3: invoke the helper after the design submission
# ============================================================
OLD_POST_SUBMIT = '''    (run_dir / "run_meta.json").write_text(json.dumps(meta, indent=2))
    print(f"Submitted: jobid={jobid}, uuid8={uuid8}")
    print(f"Run dir:   {run_dir}")
if __name__ == "__main__":
    main()'''

NEW_POST_SUBMIT = '''    (run_dir / "run_meta.json").write_text(json.dumps(meta, indent=2))
    print(f"Submitted: jobid={jobid}, uuid8={uuid8}")
    print(f"Run dir:   {run_dir}")
    # Auto-queue refolds with dependency on the design jobid
    if args.auto_refold:
        print(f"--- Queuing refolds dependent on jobid={jobid} ---")
        _submit_dependent_refolds(jobid, [run_dir], args.label)
if __name__ == "__main__":
    main()'''


# ============================================================
# Helper-injection anchor (before `def main():`)
# ============================================================
OLD_BEFORE_MAIN = "def main():\\n    ap = argparse.ArgumentParser"

# Actual literal anchor (no escaping for replacement; use the real \\n)
OLD_BEFORE_MAIN_LITERAL = "def main():"


def main():
    if not LAUNCH.exists():
        raise SystemExit(f"launch.py not found at {LAUNCH}")
    text = LAUNCH.read_text()

    if "--auto-refold" in text:
        print(f"Already patched: --auto-refold flag in launch.py")
        return

    # Apply CLI patch
    if OLD_CLI not in text:
        raise SystemExit(f"Anchor 1 not found:\n{repr(OLD_CLI)}")
    text = text.replace(OLD_CLI, NEW_CLI, 1)

    # Inject helper function before `def main():`
    if OLD_BEFORE_MAIN_LITERAL not in text:
        raise SystemExit("Anchor 'def main():' not found")
    text = text.replace(OLD_BEFORE_MAIN_LITERAL, HELPER_FUNC + OLD_BEFORE_MAIN_LITERAL, 1)

    # Apply post-submit patch
    if OLD_POST_SUBMIT not in text:
        raise SystemExit(f"Anchor 3 not found:\n{repr(OLD_POST_SUBMIT)}")
    text = text.replace(OLD_POST_SUBMIT, NEW_POST_SUBMIT, 1)

    backup = LAUNCH.with_suffix(".py.bak_pre_autorefold")
    backup.write_text(LAUNCH.read_text())
    LAUNCH.write_text(text)

    try:
        py_compile.compile(str(LAUNCH), doraise=True)
        print("Patched launch.py with --auto-refold support")
        print(f"Backup at {backup}")
        print("Syntax check: OK")
    except py_compile.PyCompileError as e:
        print(f"SYNTAX ERROR after patch: {e}")
        print("Restoring backup...")
        LAUNCH.write_text(backup.read_text())
        raise SystemExit(1)


if __name__ == "__main__":
    main()
