#!/usr/bin/env python3
"""Render deck panels in ChimeraX: real cartoons, RNA target, common frame.

ChimeraX has no OpenGL offscreen on this machine, so it runs with its GUI
context and --exit. It still writes the images; it just needs a window.

Every panel is superposed on the RNA (chain B) of a single reference, so the
target sits identically in every image and the panels differ only in how the
binder approaches it.
"""
from __future__ import annotations

import argparse
import json
import subprocess
from pathlib import Path

CX = "/Applications/ChimeraX-1.11.1.app/Contents/bin/ChimeraX"


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--structures", type=Path, default=Path("dist/deck_structures"))
    ap.add_argument("--out", type=Path, default=Path("dist/deck_images"))
    ap.add_argument("--width", type=int, default=1000)
    ap.add_argument("--height", type=int, default=850)
    args = ap.parse_args()

    meta = json.loads((args.structures / "deck_meta.json").read_text())
    panels = meta["panels"]
    args.out.mkdir(parents=True, exist_ok=True)
    ref = args.structures / panels[0]["pdb"]

    lines = [
        "set bgColor white",
        "graphics silhouettes true width 1.2",
        "lighting soft",
        "lighting shadows false",
        f"open {ref.resolve()}",          # #1 is the reference, kept for alignment
        "hide #1 models",
    ]
    for p in panels:
        pdb = (args.structures / p["pdb"]).resolve()
        png = (args.out / f"{Path(p['pdb']).stem}.png").resolve()
        lines += [
            f"open {pdb}",
            # dssp assigns secondary structure from the backbone; these files
            # carry no HELIX/SHEET records of their own.
            "dssp #2",
            # matchmaker, not align: align needs equal atom counts and the RNA
            # is 16 nt in most designs but 18 in BD1+LiMPNN and HD+LiMPNN.
            "matchmaker #2/B to #1/B",
            "hide #2 atoms",
            "show #2 cartoon",
            "cartoon style #2/A arrows true",
            "color #2/A #8d99ae",
            "color #2/A & helix #d1495b",
            "color #2/A & strand #edae49",
            "color #2/B #2a9d8f",
            "nucleotides #2/B tube/slab",
            "view #2",
            "zoom 0.9",
            f"save {png} width {args.width} height {args.height} supersample 3",
            "close #2",
        ]
    lines.append("exit")

    script = args.out / "_render.cxc"
    script.write_text("\n".join(lines) + "\n")
    print(f"rendering {len(panels)} panels via ChimeraX ...")
    r = subprocess.run([CX, "--exit", "--cmd", f"open {script.resolve()}"],
                       capture_output=True, text=True)
    made = sorted(args.out.glob("c*.png"))
    print(f"  {len(made)}/{len(panels)} images written -> {args.out}")
    if len(made) < len(panels):
        tail = [l for l in (r.stdout + r.stderr).splitlines()
                if "rror" in l or "ailed" in l][:6]
        for t in tail:
            print("   ", t)
    return 0 if made else 1


if __name__ == "__main__":
    raise SystemExit(main())
