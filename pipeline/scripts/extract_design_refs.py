#!/usr/bin/env python3
"""extract_design_refs.py — decompress each design archive ONCE and stage the 2000
designed-backbone CIFs the seed scoring needs as loose files.

The scoring array was resolving `<archive>.tar.zst!<member>` per task, so the same
~1522 archives were being decompressed ~29 times over (12.6k decompressions) with
8 tasks contending for Lustre. Staging once turns the scoring into pure DockQ.

Writes <outdir>/<basename>.cif and a design_refs.csv pointing at them.
"""
import argparse
import csv
import os
import subprocess
from collections import defaultdict
from multiprocessing import Pool

SEP = ".tar.zst!"


def extract(item):
    archive, members, outdir = item
    todo = [m for m in members if not os.path.exists(os.path.join(outdir, os.path.basename(m)))]
    if not todo:
        return archive, 0, ""
    r = subprocess.run(["tar", "--zstd", "-xf", archive, "-C", outdir, "--transform",
                        "s|.*/||", *todo], capture_output=True, text=True)
    if r.returncode != 0:
        return archive, 0, r.stderr.strip()[:180]
    return archive, len(todo), ""


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--refs", required=True, help="design_id,designed_cif_path (archive refs)")
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--out-csv", required=True)
    ap.add_argument("--jobs", type=int, default=8)
    a = ap.parse_args()
    os.makedirs(a.outdir, exist_ok=True)

    rows = list(csv.DictReader(open(a.refs)))
    by_arch = defaultdict(set)
    loose = {}
    for r in rows:
        ref = r["designed_cif_path"]
        if SEP in ref:
            arch, member = ref.split(SEP, 1)
            arch = arch + ".tar.zst"
            if not os.path.isabs(arch):
                arch = os.path.join("/cluster/work/wicky/lvalsecchi/rna_binder", arch)
            by_arch[arch].add(member)
            loose[r["design_id"]] = os.path.join(a.outdir, os.path.basename(member))
        else:
            loose[r["design_id"]] = ref
    print(f"{len(rows)} designs, {len(by_arch)} unique archives", flush=True)

    items = [(arch, sorted(m), a.outdir) for arch, m in by_arch.items()]
    n_files, n_err = 0, 0
    with Pool(a.jobs) as pool:
        for i, (arch, k, err) in enumerate(pool.imap_unordered(extract, items, chunksize=4)):
            n_files += k
            if err:
                n_err += 1
                print(f"  ERR {os.path.basename(arch)}: {err}", flush=True)
            if (i + 1) % 200 == 0:
                print(f"  {i + 1}/{len(items)} archives, {n_files} files", flush=True)
    print(f"extracted {n_files} files, {n_err} archive errors", flush=True)

    missing = [d for d, p in loose.items() if not os.path.exists(p)]
    with open(a.out_csv, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["design_id", "designed_cif_path"])
        for d, p in loose.items():
            if os.path.exists(p):
                w.writerow([d, p])
    print(f"wrote {a.out_csv}: {len(loose) - len(missing)} designs "
          f"({len(missing)} missing)", flush=True)


if __name__ == "__main__":
    main()
