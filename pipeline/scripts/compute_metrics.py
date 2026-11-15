#!/usr/bin/env python3
"""compute_metrics.py — append structural + RNA-contact + H-bond metrics.

Step 3 of the analyses/ pipeline (collect -> clean -> compute_metrics).

For each design row it loads the AF3-msaoff refold CIF and computes, via the
project's `rna_binder.analysis.structural.compute_all_structural` (Biotite):

    interface_bsa_approx        buried surface area (Å²)
    ca_rmsd_vs_designed         Cα RMSD of refold vs the designed backbone (Å)
    n_interface_contacts        protein residues within 5 Å of RNA
    radius_of_gyration_ca       Cα radius of gyration (Å)
    helix_fraction / sheet_fraction / loop_fraction      (Biotite annotate_sse)
    sse_helix_fraction / sse_strand_fraction             (pydssp DSSP)
    rna_contact_nt_1..16        per-nucleotide protein↔RNA heavy-atom contact
                                counts for ACUGCCGUAUAGGCAG (chain A binder vs
                                chain B RNA, 4.5 Å)

It also adds the tip-weighted protein↔RNA contact score (rna_binder.analysis.
tip_contact): favorable contacts weighted by proximity to the target hairpin's
loop tip vs. the scaffold-adjacent stem bottom. Columns tip_weighted_score,
tip_bottom_fraction, tip_pass, per-category tip_<hbond|basic_p|base|vdw>_wsum/_n,
etc. By default it scores the all-atom AF3 refold (--tip-source af3) so scores are
comparable across stacks; --tip-source design uses the generator output (note the
LiMPNN stacks emit backbone-only design files, collapsing the sidechain
categories) and --tip-source hybrid falls back to the refold only when the design
lacks sidechains. Disable with --no-tip-contact.

Then it fills `n_h_bonds_interface` by shelling out to the existing
`hbonds_pyrosetta.py` inside pyrosetta.sif (Rosetta's hydrogen-aware, angle-checked
detection — Biotite can't do this on heavy-atom-only AF3/Boltz CIFs). Disable with
--no-hbonds.

CIF paths in the CSV are stable refs (plain path, or `<archive>.tar.zst!<member>`
for packed runs); both are resolved with archive_read.resolve().

Default input is sample3k.csv (fast iteration); pass --full to run on
collected_clean.csv instead. Output defaults to <input_stem>_metrics.csv.

Usage
-----
    python3 compute_metrics.py --dir analyses/<name>            # -> sample3k_metrics.csv
    python3 compute_metrics.py --dir analyses/<name> --full     # -> collected_clean_metrics.csv
    python3 compute_metrics.py --input some.csv -o out.csv --jobs 16
"""
from __future__ import annotations

import argparse
import os
import shutil
import subprocess
import sys
from multiprocessing import Pool
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent          # <root>/code/scripts
RB = SCRIPT_DIR.parents[1]                             # <root>
SRC = RB / "code" / "src"
for _p in (str(SRC), str(SCRIPT_DIR)):
    if _p not in sys.path:
        sys.path.insert(0, _p)

import pandas as pd  # noqa: E402

import archive_read  # noqa: E402
from rna_binder.analysis.structural import compute_all_structural  # noqa: E402
from rna_binder.analysis.tip_contact import (  # noqa: E402
    compute_tip_weighted_contact, TIP_COLUMNS,
    GAMMA as TIP_GAMMA_DEFAULT, BOTTOM_FRAC_CUTOFF as TIP_CUTOFF_DEFAULT)
from rna_binder.analysis.red_contacts import compute_red_contacts, RED_COLUMNS  # noqa: E402
from rna_binder.analysis.seq_complexity import seq_complexity, SEQ_COLUMNS  # noqa: E402

PYROSETTA_SIF = "/cluster/work/wicky/containers/pyrosetta.sif"
HBOND_PATH_COL = "af3_msaoff_cif_path"
HBOND_COL = "n_h_bonds_interface"
# compute_all_structural also returns the refold-derived sequence; we already have
# the design sequence column, so drop the duplicate to avoid confusion.
_DROP_FROM_METRICS = ("binder_sequence",)

# tip-weighted contact score config, set from CLI in main() before the pool forks
# (Linux fork inherits module globals). Keys: enabled, source (af3|design|hybrid),
# gamma, cutoff. "af3" scores the all-atom refold (uniform across stacks); "design"
# the generator output (as-designed, but LiMPNN stacks are backbone-only → the
# sidechain categories collapse); "hybrid" uses design-native when it carries
# sidechains, else falls back to the refold.
_TIP_CFG = {"enabled": True, "source": "af3",
            "gamma": TIP_GAMMA_DEFAULT, "cutoff": TIP_CUTOFF_DEFAULT}
_TIP_KEEP = set(TIP_COLUMNS) | {"tip_contact_error"}

# red (unfavorable) protein-RNA contact pass config, set from CLI before the pool
# forks. red_only = compute ONLY red_* (fast backfill onto an existing metrics CSV).
_RED_CFG = {"enabled": True, "red_only": False}


def _tip_metrics(cif, designed):
    """Tip-weighted contact columns for one design, honoring _TIP_CFG. `cif` and
    `designed` are resolved local paths (or None). Filtered to the scalar schema."""
    src = _TIP_CFG["source"]
    primary = designed if src in ("design", "hybrid") else cif
    if src == "hybrid" and not primary:
        primary = cif
    tm = compute_tip_weighted_contact(
        Path(primary) if primary else None,
        gamma=_TIP_CFG["gamma"], bottom_frac_cutoff=_TIP_CFG["cutoff"])
    # hybrid: design-native lacked sidechains → redo on the all-atom refold
    if (src == "hybrid" and cif and primary != cif
            and tm.get("tip_has_sidechains") is False):
        tm = compute_tip_weighted_contact(
            Path(cif), gamma=_TIP_CFG["gamma"], bottom_frac_cutoff=_TIP_CFG["cutoff"])
    return {k: v for k, v in tm.items() if k in _TIP_KEEP}


def _compute_one(task):
    """(idx, metrics dict) for one design. Never raises — failures become a
    structural_error string so the row is kept with NaN metrics."""
    idx, cif_ref, designed_ref = task
    if not cif_ref or (isinstance(cif_ref, float)):
        return idx, {"structural_error": "no af3_msaoff_cif_path"}

    # red-only backfill: resolve the refold cif and compute ONLY the red_* columns.
    if _RED_CFG["red_only"]:
        try:
            cif = archive_read.resolve(cif_ref)
            return idx, dict(compute_red_contacts(Path(cif)))
        except Exception as e:
            return idx, {"red_contacts_error": f"resolve: {type(e).__name__}: {e}"}

    try:
        cif = archive_read.resolve(cif_ref)
        designed = None
        if designed_ref and not isinstance(designed_ref, float):
            try:
                designed = archive_read.resolve(designed_ref)
            except Exception:
                designed = None
        m = compute_all_structural(Path(cif), Path(designed) if designed else None)
    except Exception as e:
        return idx, {"structural_error": f"{type(e).__name__}: {e}"}
    for k in _DROP_FROM_METRICS:
        m.pop(k, None)
    if _TIP_CFG["enabled"]:
        try:
            m.update(_tip_metrics(cif, designed))
        except Exception as e:  # never let the tip pass sink a row
            m["tip_contact_error"] = f"{type(e).__name__}: {e}"
    if _RED_CFG["enabled"]:
        try:
            m.update(compute_red_contacts(Path(cif)))
        except Exception as e:  # never let the red pass sink a row
            m["red_contacts_error"] = f"{type(e).__name__}: {e}"
    return idx, m


def _run_hbonds(out_csv: Path, sif: str, jobs: int, path_col: str = HBOND_PATH_COL) -> bool:
    """Fill n_h_bonds_interface in-place via hbonds_pyrosetta.py inside pyrosetta.sif.
    Returns True on success; warns and returns False if it can't run."""
    in_container = bool(os.environ.get("APPTAINER_CONTAINER")
                        or os.environ.get("SINGULARITY_CONTAINER"))
    if in_container:
        print("[hbonds] running inside a container; cannot nest apptainer. "
              "Run hbonds_pyrosetta.py separately:", file=sys.stderr)
        print(f"  apptainer exec --bind /cluster/work/wicky {sif} python3 "
              f"{SCRIPT_DIR/'hbonds_pyrosetta.py'} --csv {out_csv} "
              f"--path-col {path_col} --jobs {jobs} --progress", file=sys.stderr)
        return False
    if shutil.which("apptainer") is None:
        print("[hbonds] apptainer not on PATH; skipping H-bond pass.", file=sys.stderr)
        return False
    if not Path(sif).exists():
        print(f"[hbonds] {sif} not found; skipping H-bond pass.", file=sys.stderr)
        return False
    cmd = ["apptainer", "exec", "--bind", "/cluster/work/wicky", "--bind", "/tmp",
           "--env", "MPLCONFIGDIR=/tmp", sif,
           "python3", str(SCRIPT_DIR / "hbonds_pyrosetta.py"),
           "--csv", str(out_csv), "--path-col", path_col,
           "--jobs", str(jobs), "--progress"]
    print(f"[hbonds] {' '.join(cmd)}", file=sys.stderr)
    try:
        subprocess.run(cmd, check=True)
        return True
    except subprocess.CalledProcessError as e:
        print(f"[hbonds] pyrosetta pass failed ({e}); structural CSV still written.",
              file=sys.stderr)
        return False


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--input", "-i", type=Path, default=None,
                    help="input CSV (default: <dir>/sample3k.csv, or "
                         "<dir>/collected_clean.csv with --full)")
    ap.add_argument("--dir", "-d", type=Path, default=Path("."),
                    help="analyses dir used to find the default input/output")
    ap.add_argument("--full", action="store_true",
                    help="default input becomes collected_clean.csv (the whole set)")
    ap.add_argument("--output", "-o", type=Path, default=None,
                    help="output CSV (default: <input_stem>_metrics.csv)")
    ap.add_argument("--jobs", type=int, default=min(16, os.cpu_count() or 1))
    ap.add_argument("--path-col", default=HBOND_PATH_COL,
                    help=f"CIF column to score (default: {HBOND_PATH_COL}); set to "
                         "b2_cif_path / esm_cif_path for per-backend (3-backend) passes")
    ap.add_argument("--no-hbonds", action="store_true",
                    help="skip the PyRosetta H-bond pass")
    ap.add_argument("--pyrosetta-sif", default=PYROSETTA_SIF)
    ap.add_argument("--no-tip-contact", action="store_true",
                    help="skip the tip-weighted protein-RNA contact score")
    ap.add_argument("--tip-source", choices=("af3", "design", "hybrid"), default="af3",
                    help="structure the tip-score runs on: af3 (all-atom refold, "
                         "uniform across stacks; default), design (generator output; "
                         "LiMPNN stacks are backbone-only so basic_p/base/vdw collapse), "
                         "or hybrid (design-native if it has sidechains, else refold)")
    ap.add_argument("--tip-gamma", type=float, default=TIP_GAMMA_DEFAULT,
                    help="tip-preference sharpness: contact weight = t**gamma")
    ap.add_argument("--tip-bottom-frac-cutoff", type=float, default=TIP_CUTOFF_DEFAULT,
                    help="tip_pass = bottom_fraction <= this")
    ap.add_argument("--no-red", action="store_true",
                    help="skip the red (unfavorable) protein-RNA contact pass")
    ap.add_argument("--red-only", action="store_true",
                    help="compute ONLY the red_* columns (skip structural/tip/hbonds) — "
                         "fast backfill onto an existing metrics CSV")
    ap.add_argument("--no-seq", action="store_true",
                    help="skip the sequence low-complexity metrics (seq_max_aa_frac, seq_longest_run)")
    ap.add_argument("--progress", action="store_true")
    args = ap.parse_args()

    _TIP_CFG.update(enabled=not args.no_tip_contact, source=args.tip_source,
                    gamma=args.tip_gamma, cutoff=args.tip_bottom_frac_cutoff)
    _RED_CFG.update(enabled=not args.no_red, red_only=args.red_only)
    if args.red_only:
        args.no_hbonds = True
        _TIP_CFG["enabled"] = False

    inp = args.input or args.dir / ("collected_clean.csv" if args.full else "sample3k.csv")
    if not inp.exists():
        sys.exit(f"input CSV not found: {inp}")
    out = args.output or inp.with_name(inp.stem + "_metrics.csv")

    df = pd.read_csv(inp, low_memory=False)
    if args.path_col not in df.columns:
        sys.exit(f"input CSV missing '{args.path_col}' column: {list(df.columns)}")
    n = len(df)
    print(f"Computing structural metrics for {n} designs from {inp} "
          f"(path-col={args.path_col}, jobs={args.jobs})")

    cifs = df[args.path_col].tolist()
    designed = (df["designed_cif_path"].tolist()
                if "designed_cif_path" in df.columns else [None] * n)
    tasks = list(zip(range(n), cifs, designed))

    results = {}
    n_err = 0
    if args.jobs > 1:
        with Pool(args.jobs) as pool:
            for done, (idx, m) in enumerate(pool.imap_unordered(_compute_one, tasks, chunksize=4)):
                results[idx] = m
                n_err += int("structural_error" in m)
                if args.progress and done % 200 == 0:
                    print(f"  {done}/{n}", file=sys.stderr)
    else:
        for done, t in enumerate(tasks):
            idx, m = _compute_one(t)
            results[idx] = m
            n_err += int("structural_error" in m)
            if args.progress and done % 200 == 0:
                print(f"  {done}/{n}", file=sys.stderr)

    metrics = pd.DataFrame.from_dict(results, orient="index").sort_index()
    # keep input columns first, then the new metric columns (drop any name clashes
    # from the input so the freshly computed values win).
    new_cols = [c for c in metrics.columns if c not in df.columns]
    enriched = pd.concat([df.reset_index(drop=True), metrics[new_cols].reset_index(drop=True)],
                         axis=1)
    # sequence low-complexity metrics — cheap, from the sequence column (no structure)
    if not args.no_seq and "sequence" in enriched.columns:
        sv = [seq_complexity(s) for s in enriched["sequence"]]
        for c in SEQ_COLUMNS:
            enriched[c] = [v[c] for v in sv]
    out.parent.mkdir(parents=True, exist_ok=True)
    enriched.to_csv(out, index=False)
    print(f"Wrote {len(enriched)} rows x {len(enriched.columns)} cols -> {out}")
    print(f"  structural ok: {n - n_err}/{n} ({n_err} failed)")
    if _TIP_CFG["enabled"] and "tip_weighted_score" in enriched.columns:
        tip_ok = int(enriched["tip_weighted_score"].notna().sum())
        n_passed = int((enriched.get("tip_pass") == True).sum())  # noqa: E712
        print(f"  tip-contact ({args.tip_source}) scored: {tip_ok}/{n}"
              f" | tip_pass: {n_passed}/{tip_ok}")
    if _RED_CFG["enabled"] and "red_total_per_ifaceres" in enriched.columns:
        red_ok = int(enriched["red_total_per_ifaceres"].notna().sum())
        print(f"  red-contacts scored: {red_ok}/{n}")

    if not args.no_hbonds:
        if _run_hbonds(out, args.pyrosetta_sif, args.jobs, args.path_col):
            filled = pd.read_csv(out, low_memory=False)[HBOND_COL].notna().sum()
            print(f"  n_h_bonds_interface filled: {int(filled)}/{len(enriched)}")


if __name__ == "__main__":
    main()
