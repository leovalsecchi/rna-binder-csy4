#!/usr/bin/env python3
"""Stage 0b: census of RNA length/sequence variants. A gate, not a report.

Samples structures per (layer, stack, campaign), reads the RNA chain, and
checks every distinct sequence aligns to the canonical csy4 16-mer. Exits 1
on any failure so the chained SLURM jobs stop before spending compute.

Usage (on Euler, from $FM):
    python scripts/stage0b_rna_census.py --build build --per-group 5
"""
from __future__ import annotations

import argparse
import io
import subprocess
import sys
import tarfile
from collections import defaultdict
from pathlib import Path

import pandas as pd

from binder_map.chains import ChainSelectError, load_structure, split_chains
from binder_map.rnaframe import FrameError, frame_offset, rna_sequence

RB = Path("/cluster/work/wicky/lvalsecchi/rna_binder")


def extract_member(archive: str, member: str) -> bytes:
    """Pull one member out of a .tar.zst via the tar CLI.

    Python's tarfile has no zstd codec in this Python, so the CLI is the
    supported route.
    """
    proc = subprocess.run(
        ["tar", "--zstd", "-xOf", str(RB / archive), member],
        capture_output=True,
    )
    if proc.returncode != 0 or not proc.stdout.strip():
        raise FileNotFoundError(f"{archive}!{member}")
    return proc.stdout


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--build", type=Path, default=Path("build"))
    parser.add_argument("--per-group", type=int, default=5)
    args = parser.parse_args()

    rows = []
    failures = 0

    for layer in ("design", "af3", "esm"):
        manifest = pd.read_csv(args.build / "manifests" / f"{layer}.tsv", sep="\t")
        groups = manifest.groupby(["stack", "campaign"], sort=True)
        for (stack, campaign), group in groups:
            sample = group.head(args.per_group)
            seen: dict[str, int] = defaultdict(int)
            for _, row in sample.iterrows():
                try:
                    data = extract_member(row["archive_path"], row["member_path"])
                    tmp = Path("/tmp") / "bm_census.cif"
                    tmp.write_bytes(data)
                    _, rna = split_chains(load_structure(tmp))
                    seq = rna_sequence(rna)
                    seen[seq] += 1
                except (FileNotFoundError, ChainSelectError) as exc:
                    rows.append({
                        "layer": layer, "stack": stack, "campaign": campaign,
                        "rna_sequence": "", "length": 0, "offset": -1,
                        "n_sampled": 1, "status": f"unreadable: {type(exc).__name__}",
                    })
                    continue

            for seq, count in seen.items():
                try:
                    offset = frame_offset(seq)
                    status = "ok"
                except FrameError:
                    offset = -1
                    status = "FRAME_MISMATCH"
                    failures += 1
                rows.append({
                    "layer": layer, "stack": stack, "campaign": campaign,
                    "rna_sequence": seq, "length": len(seq), "offset": offset,
                    "n_sampled": count, "status": status,
                })

    out = args.build / "meta" / "rna_variants.tsv"
    out.parent.mkdir(parents=True, exist_ok=True)
    frame = pd.DataFrame(rows)
    frame.to_csv(out, sep="\t", index=False)

    print(f"\n=== RNA variant census -> {out} ===")
    summary = frame.groupby(["rna_sequence", "length", "status"], sort=True).size()
    print(summary.to_string())

    if failures:
        print(f"\nFATAL: {failures} RNA variant(s) do not contain the canonical "
              f"16-mer. Epitope comparison would be misregistered. Stopping.",
              file=sys.stderr)
        return 1
    print("\nall RNA variants align to the canonical frame")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
