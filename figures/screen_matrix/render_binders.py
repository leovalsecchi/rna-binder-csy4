#!/usr/bin/env python3
"""Binder cartoons for the screen-matrix header (Fig 5a).

One PNG per candidate column, BINDER ONLY. The header's axis is the binder
series and the target varies down the rows, so drawing the hairpin up here
would put the row variable in the column head.

Orientation is still the Fig 2b shared frame at y120 - the same one Fig 4c's
grid uses - so a reader moving between the two figures sees the same fold from
the same side. Each binder is aligned onto the reference hairpin by RNA C1',
and the hairpin is then hidden: the alignment is what carries the shared
orientation, not the drawn RNA.

The CAMERA, though, is fitted to each binder on its own, so each one fills its
square. The three are therefore NOT at a common scale. That is deliberate -
they are column icons saying *which* binder, and Fig 5 b and Fig 4 c are where
binders are compared.

    python code/screen_matrix/render_binders.py [out_dir]
"""
import argparse
import subprocess
import sys
import time
from pathlib import Path

import pandas as pd

_HERE = Path(__file__).resolve()
REPO = _HERE.parents[2]
sys.path.insert(0, str(REPO / "code"))
sys.path.insert(0, str(REPO / "code" / "Fig2b_map_callouts"))

from common import palette, paths                               # noqa: E402
from render_shared_frame import (CHIMERAX, BASE, canonical_orientation,  # noqa: E402
                                 probe, stack_of)

HITS32 = REPO / "working" / "structures" / "hits32"
OUT = REPO / "working" / "screen_matrix" / "binders"
# the hairpin frame every Fig 2b panel already uses
REF = "csy4free_RFD3_LiMPNN_len80_8b13398e_design50_seq0"
ANGLE = 120          # the Fig 2b / Fig 4c rotation
IDS = (1, 2, 32)     # the header columns, matching COLS in make_figure.py


def _find(design):
    cif = BASE / "cifs" / f"{design}.cif"
    return cif if cif.exists() else HITS32 / f"{design}.pdb"


def candidates():
    """Ordered-binder ID -> design, from the flow table the panel b plot uses.

    The ID is the plate position, P<plate>_C<col> with six columns a plate,
    exactly as Fig5a_wetlab_on_off reads it - so "Candidate 32" here and
    "Ordered binder 32" there are the same well.
    """
    d = pd.read_csv(paths.INPUTS / "lv2_34_1_flow_on_off.csv")
    pc = d.binder.str.extract(r"P(\d+)_C(\d+)").astype(int)
    d["ID"] = ((pc[0] - 1) * 6 + pc[1]).values
    m = d.set_index("ID").design_id
    return [(i, m[i]) for i in IDS]


def build(out_dir, width):
    picks = candidates()
    ids = [REF] + [d for _, d in picks]
    files = [_find(d) for d in ids]
    meta = [probe(p) for p in files]                 # (pro, rna, res0, res1)

    L = [f"open {p}" for p in files]
    L += ["set bgColor white", "hide atoms", "cartoon"]

    rpro, rrna, r0, r1 = meta[0]
    for i in range(2, len(files) + 1):
        _, rna, a0, a1 = meta[i - 1]
        L.append(f"align #{i}/{rna}:{a0}-{a1}@C1' toAtoms "
                 f"#1/{rrna}:{r0}-{r1}@C1' move structures")

    # Every RNA hidden, reference included, and the reference binder with it.
    # No surface is made at all - this is the one place in the thesis where the
    # shared hairpin is aligned to and then not drawn.
    L += [f"hide #1/{rrna} cartoons", f"hide #1/{rrna} atoms",
          f"hide #1/{rpro} cartoons"]
    for i in range(2, len(files) + 1):
        rna = meta[i - 1][1]
        L += [f"hide #{i}/{rna} cartoons", f"hide #{i}/{rna} atoms"]

    for i in range(2, len(files) + 1):
        L.append(f"color #{i}/{meta[i-1][0]} "
                 f"{palette.stack_color(stack_of(ids[i-1]))}")

    centre, v, ang = canonical_orientation(files[0])
    allm = f"#1-{len(files)}"
    for ax, d in zip("xyz", -centre):
        L.append(f"move {ax} {d:.4f} models {allm}")
    L.append(f"turn {v[0]:.6f},{v[1]:.6f},{v[2]:.6f} {ang:.4f} center 0,0,0 "
             f"models {allm}")
    L.append(f"turn y {ANGLE}")

    L += ["lighting soft",
          f"graphics silhouettes true width 1.5 color {palette.NEUTRAL['700']}"]

    binders = [f"#{i}/{meta[i-1][0]}" for i in range(2, len(files) + 1)]
    for k, (cid, _) in enumerate(picks):
        # `view` re-aims and re-zooms the camera WITHOUT touching the
        # orientation set above, which is what lets each binder fill its own
        # frame while all three keep the shared rotation.
        L.append("hide " + " ".join(binders) + " cartoons")
        L.append(f"show {binders[k]} cartoons")
        L.append(f"view {binders[k]}")
        L.append("zoom 0.92")
        L.append(f"save {out_dir}/candidate{cid:02d}.png width {width} "
                 "supersample 3 transparentBackground true")
    L.append("exit")
    return "\n".join(L), picks


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("out_dir", nargs="?", default=str(OUT))
    ap.add_argument("--width", type=int, default=600,
                    help="px. The icon prints about 9.5 mm wide, so 600 is "
                         "already well past 600 dpi; it is downsampled again "
                         "when make_figure.py embeds it.")
    ap.add_argument("--timeout", type=int, default=300)
    a = ap.parse_args()
    out = Path(a.out_dir).resolve()
    out.mkdir(parents=True, exist_ok=True)
    timeout = a.timeout
    script_text, picks = build(out, a.width)
    script = out / "render_binders.cxc"
    script.write_text(script_text + "\n")
    for cid, design in picks:
        print(f"  candidate {cid:2d}  {design}")
    print("script ->", script)

    # macOS has no OSMesa, so ChimeraX cannot render offscreen: a window must
    # exist or `save` raises. Running the binary straight from subprocess gets
    # no usable GL context and hangs before the first save. `open` gives it a
    # real app launch; `-g` keeps it behind whatever the user is looking at and
    # `-n` starts an instance of its own, so a ChimeraX already on screen is
    # neither disturbed nor handed this script.
    targets = [out / f"candidate{cid:02d}.png" for cid, _ in picks]
    for t in targets:                      # never wait on a stale render
        t.unlink(missing_ok=True)
    app = str(Path(CHIMERAX).parents[2])
    subprocess.run(["open", "-n", "-g", "-a", app, "--args",
                    "--exit", "--cmd", f"open {script}"], check=True)

    # `open` returns at once, so wait on the FILES - and on their size settling,
    # since a PNG appears before it is finished being written.
    deadline, seen = time.time() + timeout, {}
    while time.time() < deadline:
        time.sleep(2)
        sizes = {t: t.stat().st_size if t.exists() else -1 for t in targets}
        if all(v > 0 for v in sizes.values()) and sizes == seen:
            break
        seen = sizes
    pngs = sorted(t for t in targets if t.exists())
    print(f"wrote {len(pngs)} png (expect {len(picks)}): "
          + ", ".join(p.name for p in pngs))
    if len(pngs) != len(picks):
        sys.exit(f"timed out after {timeout}s - is a ChimeraX window blocked "
                 f"on a dialog?")


if __name__ == "__main__":
    main()
