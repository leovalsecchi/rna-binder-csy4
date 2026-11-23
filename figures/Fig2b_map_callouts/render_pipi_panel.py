#!/usr/bin/env python3
"""The pi-stacking panel: an overview in the same house style and the same
shared RNA frame as the groove panels, plus a clean zoom from the SAME camera
angle showing the interaction.

Overview follows code/common/CHIMERAX_RENDERING.md exactly (design binder in
stack colour, AF3 refold in neutral 400 at 35%, RNA as a light surface, soft
lighting, neutral-700 silhouettes). The only departure is that the two
pi-stacking partners, TRP40 and U10, are drawn as sticks in saturated colour so
the viewer can see which pair the zoom is about.

Zoom is RNAcon-inspired (preset AutoRNAproteincontactscore): no molecular
surface, `nucleotides fill` for filled base rings, backbone and sidechains as
sticks, contacts drawn as pseudobonds. Everything that is not the stacking pair
is faded, so the pair reads without losing context.

The camera is set ONCE. The zoom uses `view <spec>`, which recentres and
rescales without rotating, so both images are from the identical angle.

The hairpin is placed in the shared frame by aligning onto the SAME reference
RNA the groove panels use (panel a's design), then applying the canonical
upright orientation: helix axis vertical, open end down, viewed side on.
"""
import argparse
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from common import palette                                      # noqa: E402
from render_shared_frame import (BASE, probe, canonical_orientation,  # noqa: E402
                                 stack_of)

CHIMERAX = "/Applications/ChimeraX-1.11.1.app/Contents/MacOS/ChimeraX"
REF = "csy4free_RFD3_LiMPNN_len80_8b13398e_design50_seq0"       # frame reference
DID = "csy4free_RFD3_LiMPNN_len80_e2a4d6b9_design57_seq0"
TRP, NT = 40, 10                       # the stacking pair, measured not guessed
ANGLE = 130                            # groove-panel view, 5 deg the other way


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("out_dir")
    ap.add_argument("--width", type=int, default=3000)
    a = ap.parse_args()
    out = Path(a.out_dir).resolve()
    out.mkdir(parents=True, exist_ok=True)

    ref = BASE / "cifs" / f"{REF}.cif"
    dsn = BASE / "cifs" / f"{DID}.cif"
    rfd = BASE / "af3" / f"{DID}_AF3_seed1_sample0.cif"
    rpro, rrna, r0, r1 = probe(ref)
    dpro, drna, d0, d1 = probe(dsn)
    fpro, frna, f0, f1 = probe(rfd)
    col = palette.stack_color(stack_of(DID))
    ink = palette.NEUTRAL["700"]
    centre, v, ang = canonical_orientation(ref)

    L = [f"open {ref}", f"open {dsn}", f"open {rfd}",
         "set bgColor white", "delete H", "hide atoms", "cartoon",
         f"align #2/{drna}:{d0}-{d1}@C1' toAtoms #1/{rrna}:{r0}-{r1}@C1' move structures",
         f"align #3/{frna}:{f0}-{f1}@C1' toAtoms #1/{rrna}:{r0}-{r1}@C1' move structures",
         f"hide #1 models",                       # reference is frame only
         ]
    for ax, dd in zip("xyz", -centre):
        L.append(f"move {ax} {dd:.4f} models #1-3")
    L.append(f"turn {v[0]:.6f},{v[1]:.6f},{v[2]:.6f} {ang:.4f} center 0,0,0 models #1-3")

    # ---- overview, house convention + highlighted stacking pair
    L += [
        f"color #2/{dpro} {col}",
        f"color #3/{fpro} {palette.NEUTRAL['400']}",
        f"transparency #3/{fpro} 35 target c",
        f"surface #2/{drna}:{d0}-{d1}",
        f"hide #3/{frna} cartoons",
        f"hide #3/{frna} atoms",
        f"color #2/{drna} #d7dcdf target s",
        f"transparency #2/{drna} 45 target s",
        f"hide #2/{drna} cartoons",
        f"show #2/{drna}:{NT} atoms",
        f"style #2/{drna}:{NT} stick",
        f"size #2/{drna} stickRadius 0.30",
        f"color #2/{drna}:{NT} {palette.THEME['600']}",
        f"show #3/{fpro}:{TRP} atoms",
        f"style #3/{fpro}:{TRP} stick",
        f"size #3/{fpro} stickRadius 0.30",
        f"color #3/{fpro}:{TRP} {col}",
        f"transparency #3/{fpro}:{TRP} 0 target a",
        "lighting soft",
        f"graphics silhouettes true width 1.5 color {ink}",
        "view", "zoom 1.35",
        f"turn y {ANGLE}",
        f"save {out}/pipi_overview_y{ANGLE:03d}.png width {a.width} "
        "supersample 3 transparentBackground true",
    ]

    # ---- zoom, RNAcon style, SAME camera angle (view recentres, never rotates)
    L += [
        f"~surface #2/{drna}",
        f"hide #2 models",
        f"show #3/{frna} atoms",
        f"style #3/{frna} stick",
        f"size #3/{frna} stickRadius 0.14",
        f"nucleotides #3/{frna} fill",
        f"color #3/{frna} {palette.THEME['300']}",
        f"show #3/{fpro} & protein & sidechain atoms",
        f"style #3/{fpro} stick",
        f"size #3/{fpro} stickRadius 0.14",
        f"cartoon #3/{fpro}",
        f"color #3/{fpro} {col}",
        f"transparency #3/{fpro} 80 target cas",
        f"transparency #3/{frna} 70 target as",
        f"transparency #3/{frna}:{NT} 0 target as",
        f"transparency #3/{fpro}:{TRP} 0 target as",
        f"color #3/{frna}:{NT} {palette.THEME['700']}",
        f"color #3/{fpro}:{TRP} {col}",
        # Annotation colours come from palette.json (COLORS.md: never hand-pick
        # a hex in a figure script). Both are the palette ink: black reads as
        # annotation and does not compete with the petrol RNA or the binder.
        f"contacts #3/{fpro}:{TRP} restrict #3/{frna}:{NT} overlapCutoff -0.6 "
        f"hbondAllowance 0.0 reveal true color {palette.ANNOTATION['interaction']} "
        "radius 0.075 dashes 0",
        f"hbonds #3/{fpro} restrict #3/{frna} reveal true "
        f"color {palette.ANNOTATION['hbond']} radius 0.075 dashes 0",
        f"view #3/{fpro}:{TRP} #3/{frna}:{NT}",
        "zoom 1.15",
        f"save {out}/pipi_zoom_y{ANGLE:03d}.png width {a.width} "
        "supersample 3 transparentBackground true",
        # Second pass from the IDENTICAL camera with everything except the
        # stacking pair hidden. Compositing this into a disc on the full zoom
        # clears the structures in front of and behind the interaction without
        # touching the surrounding context.
        f"hide #3/{frna} atoms",
        f"hide #3/{fpro} atoms",
        f"hide #3/{fpro} cartoons",
        f"show #3/{frna}:{NT} atoms",
        f"show #3/{fpro}:{TRP} atoms",
        f"save {out}/pipi_paironly_y{ANGLE:03d}.png width {a.width} "
        "supersample 3 transparentBackground true",
        f"show #3/{frna} atoms",
        f"show #3/{fpro} & protein & sidechain atoms",
        f"cartoon #3/{fpro}",
        # Square source for the circular cutout. `view` on the pair alone
        # centres it, and a square canvas means a centred circular mask lands
        # on the interaction without having to hunt for it in pixels.
        f"view #3/{fpro}:{TRP} #3/{frna}:{NT}",
        "zoom 1.9",
        f"save {out}/pipi_cutout_src_y{ANGLE:03d}.png width 2400 height 2400 "
        "supersample 3 transparentBackground true",
        "exit",
    ]

    script = out / "render_pipi_panel.cxc"
    script.write_text("\n".join(L) + "\n")
    print("script ->", script)
    r = subprocess.run([CHIMERAX, "--exit", "--cmd", f"open {script}"],
                       capture_output=True, text=True)
    for ln in (r.stdout + r.stderr).strip().splitlines()[-15:]:
        if "TSMSend" not in ln:
            print("  cx|", ln)
    print([p.name for p in sorted(out.glob("pipi_*.png"))
           if not p.name.startswith("._")])


if __name__ == "__main__":
    main()
