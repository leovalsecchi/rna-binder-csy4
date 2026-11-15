#!/usr/bin/env python3
"""Compute IPSAE interface scores for refold models and emit tiny per-model sidecars.

IPSAE (DunbrackLab, https://github.com/DunbrackLab/IPSAE, MIT) scores a
protein-nucleic-acid interface from the model's *full PAE matrix* -- something
the summary confidence JSONs do not contain. This driver runs the vendored
``ipsae/ipsae.py`` once per refold model, parses its chain-pair output, and
writes ``<stem>_ipsae.json`` next to the canonical ``.cif``. Those sidecars are
tiny (~200 bytes) and get packed into ``designs.tar.zst`` alongside the existing
``*_confidence.json`` / ``*_summary_confidences.json`` sidecars; the large PAE
inputs are consumed here and never persisted.

The aggregator reads the sidecar via ``rna_binder.analysis.refolder_metrics``.

PAE source per backend (``--mode``):
  af3 : ``<stem>_full_data.json`` (+ ``<stem>_summary_confidences.json``) in --designs-dir.
        refold_af3.sh copies AF3's full ``confidences.json`` to canonical as
        ``<stem>_full_data.json`` just before calling us, then deletes it pre-pack.
  b2  : Boltz ``pae_*.npz`` / ``plddt_*.npz`` / ``confidence_*.json`` under --raw-dir
        (the boltz_results_*/predictions/<design>/ tree, still present pre-pack).
  esm : ``pae_<stem>.npz`` / ``plddt_<stem>.npz`` under --raw-dir
        (written per design by esmfold2_infer.py).

Runs on the Euler host python (numpy present); ipsae.py is numpy-only. One
design failing never aborts the batch -- it is logged and skipped, mirroring the
refold workers.

Usage:
  ipsae_sidecar.py --mode af3 --designs-dir DIR [--pae-cutoff 10 --dist-cutoff 15]
  ipsae_sidecar.py --mode b2  --designs-dir DIR --raw-dir OUTPUTS_DIR
  ipsae_sidecar.py --mode esm --designs-dir DIR --raw-dir ESM_PAE_DIR
"""
from __future__ import annotations

import argparse
import json
import re
import shutil
import subprocess
import sys
import tempfile
import traceback
from pathlib import Path

HERE = Path(__file__).resolve().parent
IPSAE_PY = HERE / "ipsae" / "ipsae.py"

# ipsae writes floats for these; ints for these. Column names are the header
# tokens emitted by ipsae.py (whitespace-split).
_FLOAT_COLS = {
    "ipSAE", "ipSAE_d0chn", "ipSAE_d0dom", "ipTM_af", "ipTM_d0chn",
    "pDockQ", "pDockQ2", "LIS", "d0res", "d0chn", "d0dom",
}


def _fmt(pae_cutoff: float, dist_cutoff: float) -> tuple[str, str]:
    """Reproduce ipsae.py's zero-padded cutoff strings used in output filenames."""
    p = str(int(pae_cutoff)); p = "0" + p if pae_cutoff < 10 else p
    d = str(int(dist_cutoff)); d = "0" + d if dist_cutoff < 10 else d
    return p, d


def parse_ipsae_txt(txt_path: Path) -> dict | None:
    """Parse ipsae.py's chain-pair ``.txt`` and return the best-interface metrics.

    For our 2-chain (protein + RNA) complexes there is a single unordered pair
    with two 'asym' rows and one 'max' row; we take the 'max' row (the symmetric
    best-orientation score). For robustness with >2 chains we pick the 'max' row
    with the highest ipSAE, falling back to the best 'asym' row.
    """
    lines = [ln for ln in txt_path.read_text().splitlines() if ln.strip()]
    if not lines:
        return None
    header = lines[0].split()
    idx = {name: i for i, name in enumerate(header)}
    if "ipSAE" not in idx or "Type" not in idx:
        return None
    rows = []
    for ln in lines[1:]:
        parts = ln.split()
        if len(parts) < len(header):
            continue
        row = {}
        for name, i in idx.items():
            val = parts[i]
            if name in _FLOAT_COLS:
                try:
                    val = float(val)
                except ValueError:
                    val = None
            row[name] = val
        rows.append(row)
    if not rows:
        return None

    def _key(r):
        return r["ipSAE"] if isinstance(r.get("ipSAE"), float) else -1.0

    maxrows = [r for r in rows if r.get("Type") == "max"]
    best = max(maxrows, key=_key) if maxrows else max(rows, key=_key)
    return best


def sidecar_from_row(row: dict, pae_cutoff: float, dist_cutoff: float) -> dict:
    """Map an ipsae output row to the compact sidecar schema used by the CSV."""
    return {
        "ipsae": row.get("ipSAE"),            # primary: ipSAE (variable d0, "d0res")
        "ipsae_d0chn": row.get("ipSAE_d0chn"),
        "ipsae_d0dom": row.get("ipSAE_d0dom"),
        "ipsae_iptm": row.get("ipTM_af"),     # the model's own ipTM for the pair
        "ipsae_pdockq": row.get("pDockQ"),
        "ipsae_pdockq2": row.get("pDockQ2"),
        "ipsae_lis": row.get("LIS"),
        "ipsae_chn1": row.get("Chn1"),
        "ipsae_chn2": row.get("Chn2"),
        "ipsae_pae_cutoff": pae_cutoff,
        "ipsae_dist_cutoff": dist_cutoff,
    }


def run_one(pae_arg: Path, cif: Path, out_json: Path,
            pae_cutoff: float, dist_cutoff: float) -> bool:
    """Run ipsae on one model in an isolated temp dir; write the sidecar.

    Files are staged under the constant stem 'm' so ipsae's sibling-file
    derivation (pae->plddt->confidence, or full_data->summary_confidences) and
    output-name construction never collide with characters in the design name.
    """
    work = Path(tempfile.mkdtemp(prefix="ipsae_", dir=out_json.parent))
    try:
        suffix = pae_arg.suffix  # .json (af3) or .npz (boltz/esm)
        if suffix == ".json":
            # AF3 mode: pae_arg is <stem>_full_data.json; ipsae derives the
            # summary as <stem>_summary_confidences.json (co-located).
            staged_pae = work / "m_full_data.json"
            shutil.copy(pae_arg, staged_pae)
            summ = pae_arg.parent / f"{cif.stem}_summary_confidences.json"
            if summ.exists():
                shutil.copy(summ, work / "m_summary_confidences.json")
        else:
            # Boltz mode: pae_arg is pae_*.npz. ipsae derives the sibling files by
            # string-replacing 'pae' in the filename, so staging it as 'pae_m.npz'
            # makes those resolve to 'plddt_m.npz' / 'confidence_m.json'.
            staged_pae = work / "pae_m.npz"
            shutil.copy(pae_arg, staged_pae)
            plddt_src = pae_arg.parent / pae_arg.name.replace("pae", "plddt", 1)
            if plddt_src.exists():
                shutil.copy(plddt_src, work / "plddt_m.npz")
            conf_src = pae_arg.parent / (
                pae_arg.name.replace("pae", "confidence", 1).replace(".npz", ".json"))
            if conf_src.exists():
                # Only useful if it holds a dict-form pair_chains_iptm; harmless otherwise.
                try:
                    cd = json.loads(conf_src.read_text())
                    if isinstance(cd.get("pair_chains_iptm"), dict):
                        shutil.copy(conf_src, work / "confidence_m.json")
                except Exception:
                    pass
        staged_cif = work / "m.cif"
        shutil.copy(cif, staged_cif)

        proc = subprocess.run(
            [sys.executable, str(IPSAE_PY), staged_pae.name, "m.cif",
             str(pae_cutoff), str(dist_cutoff)],
            cwd=work, capture_output=True, text=True, timeout=300,
        )
        p, d = _fmt(pae_cutoff, dist_cutoff)
        txt = work / f"m_{p}_{d}.txt"
        if proc.returncode != 0 or not txt.exists():
            print(f"[FAIL] {cif.name}: ipsae rc={proc.returncode} "
                  f"{proc.stderr.strip()[:200]}", flush=True)
            return False
        row = parse_ipsae_txt(txt)
        if row is None:
            print(f"[FAIL] {cif.name}: could not parse ipsae output", flush=True)
            return False
        out_json.write_text(json.dumps(sidecar_from_row(row, pae_cutoff, dist_cutoff)))
        return True
    finally:
        shutil.rmtree(work, ignore_errors=True)


def _iter_cifs(designs_dir: Path, tag: str):
    """Yield canonical refold .cif paths for a backend tag, skipping done sidecars."""
    for cif in sorted(designs_dir.glob(f"*_{tag}_*.cif")):
        # skip AF3's *_top.cif (no per-seed PAE) and anything already scored
        if cif.stem.endswith("_top"):
            continue
        out_json = cif.parent / f"{cif.stem}_ipsae.json"
        if out_json.exists():
            continue
        yield cif, out_json


def mode_af3(designs_dir: Path, pae_cutoff, dist_cutoff) -> tuple[int, int]:
    ok = fail = 0
    for cif, out_json in _iter_cifs(designs_dir, "AF3"):
        full = designs_dir / f"{cif.stem}_full_data.json"
        if not full.exists():
            print(f"[SKIP] {cif.name}: no full_data.json (full PAE not retained)", flush=True)
            continue
        try:
            ok += 1 if run_one(full, cif, out_json, pae_cutoff, dist_cutoff) else 0
            fail += 0 if out_json.exists() else 1
        except Exception as exc:
            fail += 1
            print(f"[FAIL] {cif.name}: {exc}", flush=True)
            traceback.print_exc()
    return ok, fail


def _boltz_pae_for(raw_dir: Path, design_name: str, sample: str) -> Path | None:
    """Locate Boltz's pae_<design>_model_<sample>.npz under the outputs tree."""
    hits = list(raw_dir.glob(
        f"**/predictions/{design_name}/pae_{design_name}_model_{sample}.npz"))
    if not hits:  # fall back to any pae npz for that design+sample
        hits = list(raw_dir.glob(f"**/{design_name}/pae_*model_{sample}.npz"))
    return hits[0] if hits else None


def mode_b2(designs_dir: Path, raw_dir: Path, pae_cutoff, dist_cutoff) -> tuple[int, int]:
    ok = fail = 0
    for cif, out_json in _iter_cifs(designs_dir, "B2"):
        m = re.match(r"^(?P<name>.+)_B2_sample(?P<s>\d+)$", cif.stem)
        if not m:
            continue
        pae = _boltz_pae_for(raw_dir, m.group("name"), m.group("s"))
        if pae is None:
            print(f"[SKIP] {cif.name}: no pae npz under {raw_dir}", flush=True)
            continue
        try:
            ok += 1 if run_one(pae, cif, out_json, pae_cutoff, dist_cutoff) else 0
            fail += 0 if out_json.exists() else 1
        except Exception as exc:
            fail += 1
            print(f"[FAIL] {cif.name}: {exc}", flush=True)
            traceback.print_exc()
    return ok, fail


def mode_esm(designs_dir: Path, raw_dir: Path, pae_cutoff, dist_cutoff) -> tuple[int, int]:
    ok = fail = 0
    for cif, out_json in _iter_cifs(designs_dir, "ESM"):
        pae = raw_dir / f"pae_{cif.stem}.npz"
        if not pae.exists():
            print(f"[SKIP] {cif.name}: no {pae.name} under {raw_dir}", flush=True)
            continue
        try:
            ok += 1 if run_one(pae, cif, out_json, pae_cutoff, dist_cutoff) else 0
            fail += 0 if out_json.exists() else 1
        except Exception as exc:
            fail += 1
            print(f"[FAIL] {cif.name}: {exc}", flush=True)
            traceback.print_exc()
    return ok, fail


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--mode", required=True, choices=["af3", "b2", "esm"])
    ap.add_argument("--designs-dir", required=True, type=Path,
                    help="canonical designs/ dir holding the refold .cif files")
    ap.add_argument("--raw-dir", type=Path,
                    help="raw PAE source dir (b2: Boltz outputs/, esm: pae npz dir)")
    ap.add_argument("--pae-cutoff", type=float, default=10.0)
    ap.add_argument("--dist-cutoff", type=float, default=15.0)
    args = ap.parse_args()

    if not IPSAE_PY.exists():
        print(f"ERROR: vendored ipsae.py not found at {IPSAE_PY}", file=sys.stderr)
        sys.exit(2)
    if not args.designs_dir.is_dir():
        print(f"ERROR: designs dir not found: {args.designs_dir}", file=sys.stderr)
        sys.exit(2)

    if args.mode == "af3":
        ok, fail = mode_af3(args.designs_dir, args.pae_cutoff, args.dist_cutoff)
    elif args.mode == "b2":
        if not args.raw_dir:
            print("ERROR: --raw-dir required for b2", file=sys.stderr); sys.exit(2)
        ok, fail = mode_b2(args.designs_dir, args.raw_dir, args.pae_cutoff, args.dist_cutoff)
    else:
        if not args.raw_dir:
            print("ERROR: --raw-dir required for esm", file=sys.stderr); sys.exit(2)
        ok, fail = mode_esm(args.designs_dir, args.raw_dir, args.pae_cutoff, args.dist_cutoff)

    print(f"IPSAE sidecars: {ok} written, {fail} failed ({args.mode})", flush=True)
    # Non-fatal: a scoring add-on must never fail the refold job.
    sys.exit(0)


if __name__ == "__main__":
    main()
