#!/usr/bin/env python3
"""Re-grid a merged complex blob. Run locally; no cluster round-trip needed.

Extraction writes a 0.1 A grid, which is finer than these backbone models
justify and costs real page weight. Re-gridding trades precision no one can
see for megabytes everyone feels.
"""
from __future__ import annotations

import argparse
import gzip
from pathlib import Path

import numpy as np


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--blob", type=Path, required=True)
    ap.add_argument("--grid", type=float, default=0.2)
    ap.add_argument("--out", type=Path, required=True)
    args = ap.parse_args()

    raw = gzip.decompress(args.blob.read_bytes())
    d = np.frombuffer(raw, dtype="<i2").reshape(-1, 368, 3)
    absolute = np.cumsum(d.astype(np.int32), axis=1) / 10.0   # undo 0.1 A deltas

    q = np.round(absolute / args.grid).astype(np.int32)
    dq = q.copy()
    dq[:, 1:] -= q[:, :-1]
    dq = np.clip(dq, -32767, 32767).astype(np.int16)

    packed = gzip.compress(dq.tobytes(), 9)
    args.out.write_bytes(packed)
    print(f"{args.blob.name}: {len(raw)/1e6:.1f} MB raw -> {len(packed)/1e6:.1f} MB gz "
          f"at {args.grid} A ({len(packed)*4/3/1e6:.1f} MB base64)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
