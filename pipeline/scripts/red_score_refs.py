#!/usr/bin/env python3
"""red_score_refs.py — red (unfavorable) protein-RNA contacts for a manifest of CIF refs.

Manifest TSV columns: id, ref, stack.  `ref` is an RB-relative
"<archive>.tar.zst!<member>" (or a plain path). Targeted per-member tar extraction
(light for a scattered sample; for the full 100k prefer archive_read whole-archive
caching). Auto chain-detect (design = single binder + one RNA chain).
"""
import argparse
import csv
import os
import subprocess
import sys
import tempfile
from multiprocessing import Pool

RB = "/cluster/work/wicky/lvalsecchi/rna_binder"
sys.path.insert(0, os.path.join(RB, "code/src"))
from rna_binder.analysis.red_contacts import compute_red_contacts, RED_COLUMNS


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
    with tempfile.TemporaryDirectory(prefix="red_") as td:
        cif = _extract(ref, td)
        if not cif:
            return {"id": rid, "stack": stack, "red_contacts_error": "extract failed"}
        m = compute_red_contacts(cif)
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
    cols = ["id", "stack"] + list(RED_COLUMNS) + ["red_contacts_error"]
    with open(a.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, extrasaction="ignore")
        w.writeheader()
        for r in rows:
            w.writerow(r)
    ok = sum(1 for r in rows if "red_contacts_error" not in r)
    sys.stderr.write(f"wrote {a.out} ({ok}/{len(rows)} ok)\n")


if __name__ == "__main__":
    main()
