#!/usr/bin/env python3
"""Stage 8: bake the explorer. Runs locally.

Usage:
    python scripts/stage8_build_html.py --results local_results --build build \
        --out dist/binder_map.html
"""
from __future__ import annotations

import argparse
from pathlib import Path

from binder_map.htmlbuild import bake, build_payload

ROOT = Path(__file__).resolve().parent.parent


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--results", type=Path, default=ROOT / "local_results")
    parser.add_argument("--build", type=Path, default=ROOT / "build")
    parser.add_argument("--out", type=Path, default=ROOT / "dist" / "binder_map.html")
    parser.add_argument("--layers", default="",
                        help="comma-separated subset, e.g. 'design'. Splitting by "
                             "layer keeps each file small enough to share.")
    args = parser.parse_args()

    payload = build_payload(args.results, args.build)
    if args.layers:
        keep = {l.strip() for l in args.layers.split(",") if l.strip()}
        payload["layers"] = {k: v for k, v in payload["layers"].items() if k in keep}
        if not payload["layers"]:
            print(f"no layers matched {sorted(keep)}")
            return 1
    for layer, entry in payload["layers"].items():
        axes = ", ".join(entry["axes"])
        print(f"  {layer}: {entry['n']:,d} nodes, axes: {axes or 'none'}")

    size = bake(payload, ROOT / "web" / "template.html",
                ROOT / "web" / "app.js", ROOT / "web" / "style.css", args.out)
    print(f"wrote {args.out} ({size / 1e6:.1f} MB)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
