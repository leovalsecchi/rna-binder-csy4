#!/usr/bin/env python3
"""Split the canonical DB keys into N chunks for the search array.

Query-chunk arrays rather than the module's MPI support: a failed chunk is a
one-line resubmit, whereas a failed MPI rank kills the whole run.
"""
from __future__ import annotations

import argparse
from pathlib import Path

FM = Path("/cluster/work/wicky/lvalsecchi/rna_binder/analyses/foldseek_map")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--layer", required=True)
    parser.add_argument("--n-chunks", type=int, default=32)
    parser.add_argument("--results", type=Path, default=FM / "results")
    args = parser.parse_args()

    out = args.results / args.layer
    keys = (out / "db" / "keep_keys.txt").read_text().split()
    chunk_dir = out / "db" / "chunks"
    chunk_dir.mkdir(parents=True, exist_ok=True)

    for i in range(args.n_chunks):
        part = keys[i::args.n_chunks]
        (chunk_dir / f"q{i:02d}.txt").write_text("\n".join(part) + "\n")
    print(f"{len(keys):,d} queries -> {args.n_chunks} chunks in {chunk_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
