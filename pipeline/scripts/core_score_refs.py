#!/usr/bin/env python3
"""core_score_refs.py — binder core-packing metrics for a manifest of CIF refs.

Manifest TSV columns: id, ref, stack.  `ref` is an RB-relative "<archive>.tar.zst!<member>"
(or a plain path). Targeted per-member tar extraction into an auto-removed TemporaryDirectory
(no inode leak). Clone of red_score_refs.py — swaps in compute_core_packing.
"""
import argparse
import csv
import os
import subprocess
import sys
import tempfile
import warnings
from multiprocessing import Pool

warnings.filterwarnings("ignore")   # biotite auth_comp_id->label_comp_id fallback spam

RB = "/cluster/work/wicky/lvalsecchi/rna_binder"
sys.path.insert(0, os.path.join(RB, "code/src"))
from rna_binder.analysis.core_packing import compute_core_packing, CORE_COLUMNS


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
    rid, ref, stack = task
    with tempfile.TemporaryDirectory(prefix="core_") as td:
        cif = _extract(ref, td)
        if not cif:
            return {"id": rid, "stack": stack, "core_packing_error": "extract failed"}
        m = compute_core_packing(cif)
    m["id"] = rid
    m["stack"] = stack
    return m


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("manifest")
    ap.add_argument("-o", "--out", required=True)
    ap.add_argument("-j", "--procs", type=int, default=8)
    a = ap.parse_args()
    tasks = []
    with open(a.manifest) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            tasks.append((row["id"], row["ref"], row.get("stack", "")))
    sys.stderr.write(f"scoring {len(tasks)} refs with {a.procs} procs\n")
    with Pool(a.procs) as p:
        rows = p.map(_one, tasks, chunksize=8)
    cols = ["id", "stack"] + list(CORE_COLUMNS) + ["core_packing_error"]
    with open(a.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, extrasaction="ignore")
        w.writeheader()
        for r in rows:
            w.writerow(r)
    ok = sum(1 for r in rows if "core_packing_error" not in r)
    sys.stderr.write(f"wrote {a.out} ({ok}/{len(rows)} ok)\n")


if __name__ == "__main__":
    main()
