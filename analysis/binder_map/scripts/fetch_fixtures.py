#!/usr/bin/env python3
"""Pull a handful of real CIFs from Euler into tests/fixtures/.

Picks structures that exercise the dialect and chain-order differences:
  - RFD3+LiMPNN design  -> protein is chain B (inverted)
  - BG+BIF design       -> different mmCIF dialect, 18-nt RNA
  - BD1+LiMPNN design   -> protein is chain A
  - one af3 refold, one esm refold
Plus a second sequence variant of the RFD3 backbone, to test dedup.

Run from the repo root:  python scripts/fetch_fixtures.py
"""
from __future__ import annotations

import csv
import gzip
import subprocess
import sys
from pathlib import Path

CSV = Path("/Users/leo/Documents/binder_filter/exports/passing.csv")
OUT = Path(__file__).resolve().parent.parent / "tests" / "fixtures"
RB = "/cluster/work/wicky/lvalsecchi/rna_binder"

# design_id -> (column, fixture name)
WANTED = {
    "csy4free_RFD3_LiMPNN_len80_0118c551_design45_seq4": ("design_cif", "rfd3_design"),
    "csy4free_RFD3_LiMPNN_len80_0118c551_design45_seq7": ("design_cif", "rfd3_design_seq7"),
    "csy4free_BG_BIF_len80_c55bf73f_design16_seq8": ("design_cif", "bgbif_design"),
    "csy4free_BD1_LiMPNN_len80_031de1a0_design48_seq1": ("design_cif", "bd1_design"),
    "csy4free_BD1_LiMPNN_len80_031de1a0_design48_seq1_af3": ("af3_cif", "af3_refold"),
    "csy4free_BD1_LiMPNN_len80_031de1a0_design48_seq1_esm": ("esm_cif", "esm_refold"),
}


def rows_by_id(wanted_ids: set[str]) -> dict[str, dict]:
    found: dict[str, dict] = {}
    with CSV.open() as fh:
        fh.readline()
        fh.readline()
        for row in csv.DictReader(fh):
            if row["design_id"] in wanted_ids:
                found[row["design_id"]] = row
            if len(found) == len(wanted_ids):
                break
    return found


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    base_ids = {k.removesuffix("_af3").removesuffix("_esm") for k in WANTED}
    rows = rows_by_id(base_ids)

    missing = base_ids - rows.keys()
    if missing:
        print(f"FATAL: design_ids not in passing.csv: {sorted(missing)}", file=sys.stderr)
        return 1

    for key, (column, name) in WANTED.items():
        design_id = key.removesuffix("_af3").removesuffix("_esm")
        archive, member = rows[design_id][column].split("!", 1)
        cmd = f"cd {RB} && tar --zstd -xOf {archive} {member}"
        proc = subprocess.run(
            ["ssh", "-o", "ConnectTimeout=60", "euler", cmd],
            capture_output=True,
        )
        if proc.returncode != 0 or not proc.stdout.strip():
            print(f"FATAL: could not fetch {name} ({design_id})", file=sys.stderr)
            print(proc.stderr.decode()[:400], file=sys.stderr)
            return 1
        dest = OUT / f"{name}.cif.gz"
        dest.write_bytes(gzip.compress(proc.stdout))
        print(f"  {name:20s} {len(proc.stdout):>9,d} B -> {dest.name}")

    print("fixtures written to", OUT)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
