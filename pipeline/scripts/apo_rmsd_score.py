#!/usr/bin/env python3
"""apo_rmsd_score.py — CPU pass: Cα RMSD of each apo (protein-only) fold vs its ESM complex
refold, + pLDDT/pTM from the fold's confidence sidecar. Reuses structural.ca_rmsd_vs_designed
(loads both, splits protein, superimposes Cα). biotite host python (no GPU). Clone of the
core/red scoring pattern; per-design tar extraction into an auto-removed tempdir (no inode leak).

Inputs: the apo designs dir (<id>_ESM_sample0.cif + _confidence.json) + the chunk TSV (for each
id's esm_ref) + the packed-archive ref for the apo_cif_path column.
"""
import argparse
import csv
import json
import os
import subprocess
import sys
import tempfile
import warnings
from multiprocessing import Pool

warnings.filterwarnings("ignore")

RB = "/cluster/work/wicky/lvalsecchi/rna_binder"
sys.path.insert(0, os.path.join(RB, "code/src"))
from pathlib import Path
from rna_binder.analysis.structural import ca_rmsd_vs_designed

_DESIGNS = None
_ARCHIVE = None


def _extract(ref, tmpdir):
    if "!" not in ref:
        p = ref if os.path.isabs(ref) else os.path.join(RB, ref)
        return p if os.path.exists(p) else None
    arch, member = ref.split("!", 1)
    arch = arch if os.path.isabs(arch) else os.path.join(RB, arch)
    out = os.path.join(tmpdir, os.path.basename(member))
    with open(out, "wb") as fh:
        r = subprocess.run(["tar", "--zstd", "-xOf", arch, member],
                           stdout=fh, stderr=subprocess.DEVNULL)
    return out if (r.returncode == 0 and os.path.getsize(out) > 0) else None


def _one(task):
    rid, esm_ref, stack = task
    apo = os.path.join(_DESIGNS, f"{rid}_ESM_sample0.cif")
    row = {"id": rid, "stack": stack,
           "apo_cif_path": f"{_ARCHIVE}!{rid}_ESM_sample0.cif"}
    if not os.path.exists(apo):
        row["apo_error"] = "no apo fold"
        return row
    conf = os.path.join(_DESIGNS, f"{rid}_ESM_sample0_confidence.json")
    if os.path.exists(conf):
        try:
            c = json.load(open(conf))
            row["apo_plddt"] = round(float(c.get("plddt")), 4)
            row["apo_ptm"] = round(float(c.get("ptm")), 4)
        except Exception:
            pass
    with tempfile.TemporaryDirectory(prefix="apo_") as td:
        esm = _extract(esm_ref, td)
        if not esm:
            row["apo_error"] = "esm extract failed"
            return row
        r = ca_rmsd_vs_designed(Path(apo), Path(esm))    # apo protein vs ESM-refold protein
    if r is None:
        row["apo_error"] = "rmsd failed (len mismatch?)"
    else:
        row["apo_rmsd_vs_esm"] = round(r, 4)
    return row


def main():
    global _DESIGNS, _ARCHIVE
    ap = argparse.ArgumentParser()
    ap.add_argument("--designs-dir", required=True)
    ap.add_argument("--manifest", required=True)
    ap.add_argument("--archive", required=True, help="RB-relative packed apo archive, e.g. apo_refolds/chunk_00.tar.zst")
    ap.add_argument("-o", "--out", required=True)
    ap.add_argument("-j", "--procs", type=int, default=8)
    a = ap.parse_args()
    _DESIGNS, _ARCHIVE = a.designs_dir, a.archive
    tasks = []
    with open(a.manifest) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            tasks.append((row["id"], row["esm_ref"], row.get("stack", "")))
    sys.stderr.write(f"scoring {len(tasks)} apo folds with {a.procs} procs\n")
    with Pool(a.procs) as p:
        rows = p.map(_one, tasks, chunksize=8)
    cols = ["id", "stack", "apo_rmsd_vs_esm", "apo_plddt", "apo_ptm", "apo_cif_path", "apo_error"]
    with open(a.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, extrasaction="ignore")
        w.writeheader()
        for r in rows:
            w.writerow(r)
    ok = sum(1 for r in rows if "apo_rmsd_vs_esm" in r)
    sys.stderr.write(f"wrote {a.out} ({ok}/{len(rows)} rmsd ok)\n")


if __name__ == "__main__":
    main()
