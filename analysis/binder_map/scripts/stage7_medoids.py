#!/usr/bin/env python3
"""Stage 7: bundle every medoid's CA trace into one gzipped JSON.

One file rather than one per medoid: this project has a hard inode cap, and
the explorer wants them all anyway.

Usage:  python scripts/stage7_medoids.py --layer esm
"""
from __future__ import annotations

import argparse
import gzip
import json
import tarfile
from pathlib import Path

import numpy as np
import pandas as pd

FM = Path("/cluster/work/wicky/lvalsecchi/rna_binder/analyses/foldseek_map")


def ca_from_pdb(text: str) -> list[list[float]]:
    out = []
    for line in text.splitlines():
        if line.startswith("ATOM") and line[12:16].strip() == "CA":
            out.append([round(float(line[30:38]), 2),
                        round(float(line[38:46]), 2),
                        round(float(line[46:54]), 2)])
    return out


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--layer", required=True, choices=["design", "af3", "esm"])
    parser.add_argument("--results", type=Path, default=FM / "results")
    args = parser.parse_args()

    out = args.results / args.layer
    wanted: set[str] = set()
    for axis in ("fold", "epitope"):
        path = out / "clusters" / axis / "medoids.parquet"
        if path.exists():
            wanted |= set(pd.read_parquet(path)["struct_name"])
    print(f"{args.layer}: {len(wanted):,d} distinct medoids across axes/resolutions")

    traces: dict[str, list] = {}
    for shard in sorted((out / "shards").glob("shard_*.tar.gz")):
        with tarfile.open(shard, "r:gz") as tar:
            for member in tar:
                name = member.name.removesuffix(".pdb")
                if name not in wanted or name in traces:
                    continue
                fh = tar.extractfile(member)
                if fh is None:
                    continue
                traces[name] = ca_from_pdb(fh.read().decode())
        if len(traces) == len(wanted):
            break

    missing = wanted - traces.keys()
    if missing:
        print(f"  WARNING: {len(missing):,d} medoids had no backbone in the shards")

    dest = out / "medoids"
    dest.mkdir(parents=True, exist_ok=True)
    payload = json.dumps(traces, separators=(",", ":")).encode()
    (dest / "medoids.json.gz").write_bytes(gzip.compress(payload, 9))
    print(f"  wrote {len(traces):,d} CA traces, "
          f"{len(gzip.compress(payload, 9)) / 1e6:.1f} MB gzipped")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
