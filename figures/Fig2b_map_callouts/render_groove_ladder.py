#!/usr/bin/env python3
"""Re-render the two groove callouts (hd_b, bd1_b) as NUCLEOTIDE LADDERS.

Leo, 2026-11-17: "in 2b, we want to have the minor and major groove
representations, which are now surface representations, to be normal array
ladder representations like the one in the binding to base not ribose field."

So hd_b ("major groove") and bd1_b ("minor groove") drop the molecular surface
and take the same `nucleotides tube/slab` display as rfd3_b ("binding to base
not ribose"): ribbon backbone, bases as flat slabs, binder a 90%-transparent
ghost with the six nearest sidechains as sticks.

HOW THE GROOVE SURVIVES THE SURFACE GOING AWAY. On a surface the groove was a
tint on the skin, coloured from the groove-edge ATOMS. There is no skin now, so
the claim moves onto the SLABS: every contacted nucleotide is painted by the
face the binder reads, in the same pair of colours the surface patch used, so
the two conventions agree and the old panels stay comparable.

    thesis_petrol_800   major-groove edge   N7 O6 N6 N4 O4 C5 C6
    thesis_grey_800     minor-groove edge   N3 N2 C2 O2 O2'
    thesis_grey_300     neither, or a tie

petrol_800 and grey_800 are a matched pair by design (L* 26 against 22, dE 20.6,
23.2 simulated deuteranope); petrol_600 at L* 47 read as a highlight instead of
as one of two faces. Classification is MAJOR/MINOR from analyse_examples, per
base type, never a pooled atom set: C6 is a major-groove atom in U and C only
and C2 a minor-groove atom in the purines only, so pooling would paint purine C6
"major" and make a panel disagree with its own caption number.

Both faces are drawn in BOTH panels, which is what makes the caption checkable
rather than asserted: hd_b comes out 11 petrol slabs to 1 grey, bd1_b 10 grey to
3 petrol.

WHY THIS IS A SEPARATE SCRIPT. The eight callout panels are built by
output/Fig2b_map_composite/symmetric/zoom/make_zooms.py, which is where this
belongs - as a `kind in ("major", "minor")` branch on its nucleotide-readout
path, replacing the surface branch. That file is outside the directory this
change was allowed to touch, so the two panels are re-rendered from here
instead, reusing make_zooms' own geometry (site, near_residues, target_mask,
view_groove) and writing the same files, at the same size, into the same place:

    output/Fig2b_map_composite/symmetric/zoom/zoom_{hd_b,bd1_b}.png
    (width 2000, supersample 3, transparent background)

make_paired.py picks them up unchanged. Note that re-running make_zooms.py for
either panel puts the surface version back.

Usage:  python render_groove_ladder.py [--dry]      (--dry writes only the .cxc)
"""
import argparse
import subprocess
import sys
import time
import csv
import glob
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
ZOOM = ROOT / "output/Fig2b_map_composite/symmetric/zoom"
# analyse_examples lives in THIS directory; make_zooms.py still points its
# sys.path at the old name (code/fig02b_fold_stacks), so prepend the real one
# before importing it or the import fails.
sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(ZOOM))
import make_zooms as M                                          # noqa: E402
from analyse_examples import MAJOR, MINOR, CUT                  # noqa: E402

PANELS = ("hd_b", "bd1_b")


def groove_slabs(rna, pro, pos):
    """Contacted nucleotides split into the two groove faces, per base type.

    Returns (major, minor) residue-id lists for colouring slabs. One nucleotide
    is one slab and a slab takes its WHOLE residue's colour - an atom-subset
    selection colours a slab not at all - so a base cannot be painted half
    major and half minor: each contacted nucleotide goes to whichever face the
    binder touches more atoms of. Ties stay grey; a tie is not evidence for
    either face, and awarding it to the panel's own caption would be assuming
    the answer.
    """
    nm = np.array([str(x).strip() for x in rna.atom_name])
    rn = np.array([str(x).strip() for x in rna.res_name])
    d = ((rna.coord[:, None, :].astype(float)
          - pro.coord[None, :, :].astype(float)) ** 2).sum(-1).min(1) ** 0.5
    score = {}
    for i in np.flatnonzero(d <= CUT):
        rid = int(rna.res_id[i])
        if not 1 <= pos.get(rid, 0) <= 16:        # outside the canonical 16-mer
            continue
        mj, mn = score.get(rid, (0, 0))
        if nm[i] in MAJOR[rn[i]]:
            mj += 1
        elif nm[i] in MINOR[rn[i]]:
            mn += 1
        score[rid] = (mj, mn)
    major = sorted(r for r, (a, b) in score.items() if a > b)
    minor = sorted(r for r, (a, b) in score.items() if b > a)
    assert not set(major) & set(minor), "a nucleotide cannot read both faces"
    return major, minor


def script_for(panel, row):
    """The .cxc for one groove panel, and a one-line summary of what it shows."""
    cif = glob.glob(f"{M.AF3}/{row['design']}_AF3*.cif")
    if not cif:
        raise SystemExit(f"no AF3 structure for {panel}")
    kind, label = M.FEATURE[panel]
    assert kind in ("major", "minor"), f"{panel} is not a groove panel ({kind})"
    pres, rres, st, rch, pch, pos, rna, pro, rnames, rids = M.site(Path(cif[0]), kind)
    assert st is None
    col = M.STACKCOL[row["stack"]]
    # the six binder residues nearest THIS panel's groove face, and the segment
    # that lies in it: both are make_zooms' own definitions, so the camera and
    # the sticks match the other readout panels.
    nearsel = ",".join(str(x) for x in
                       M.near_residues(rna, pro, M.target_mask(kind, rna, pos)))
    lo, hi = pres[0] - 2, pres[-1] + 2
    maj, mnr = groove_slabs(rna, pro, pos)

    cxc = [f"cd {ZOOM}", f"open {ROOT/'config'/'thesis_colors.cxc'}",
           "set bgColor white",
           "graphics silhouettes true width 1.5 color thesis_grey_700",
           "lighting soft",
           f"open {cif[0]}",
           "hide #1 atoms", "show #1 cartoon",
           f"color #1/{pch} {col}", f"color #1/{rch} thesis_grey_300",
           "style #1 stick", "size stickRadius 0.24",
           # --- the readout display, identical to rfd3_b ---
           # The backbone tube is drawn by the cartoon machinery AND the slabs
           # need the atoms shown: with either hidden, tube/slab silently falls
           # back to a thicket of backbone sticks. `nucleotides` then hides the
           # atoms it replaces.
           "hide #1 atoms", "hide #1 cartoon",
           f"show #1/{rch} cartoon", f"show #1/{rch} atoms",
           f"style #1/{rch} stick", "size stickRadius 0.16",
           f"nucleotides #1/{rch} tube/slab",
           f"color #1/{rch} thesis_grey_300",
           f"show #1/{pch} cartoon", f"color #1/{pch} {col}",
           f"transparency #1/{pch} 90 cartoons",
           f"show #1/{pch}:{nearsel} atoms",
           f"style #1/{pch}:{nearsel} stick",
           # byhetero RESTRICTED to the binder: a bare `color byhetero` would
           # recolour the RNA heteroatoms this panel depends on being grey.
           f"color #1/{pch}:{nearsel} {col}",
           f"color #1/{pch}:{nearsel} byhetero"]
    # --- the groove, on the slabs ---
    for ids, colour in ((maj, "thesis_petrol_800"), (mnr, "thesis_grey_800")):
        if ids:                     # an empty spec `#1/B:` aborts the whole .cxc
            cxc.append(f"color #1/{rch}:{','.join(str(x) for x in ids)} {colour}")
    cxc += [
        # Ribbon back to pale grey AFTER the slabs: `color` paints ribbon and
        # slab together and `target c` reverts the ribbon only. The backbone
        # carries no groove claim, so it must not keep the groove colour.
        f"color #1/{rch} thesis_grey_300 target c",
        # Same quiet H-bond convention as the other readout panels: the slab
        # colour carries the groove, the dashes only show where the polar
        # contacts land.
        f"hbonds (#1/{pch}) restrict (#1/{rch}) {M.HB} "
        f"color thesis_grey_900 radius 0.07 dashes 6 name hb",
        # helix axis up the page, occupied groove towards the camera
        M.view_groove(rna, pos, pro, lo, hi),
        f"view #1/{rch} #1/{pch}:{nearsel}",
        "zoom 0.95",
        f"save {ZOOM}/zoom_{panel}.png width 2000 supersample 3 "
        "transparentBackground true",
        "exit",
    ]
    nt = lambda ids: ",".join(str(pos[i]) for i in ids)
    return cxc, (f"{panel:6s} {label:13s} segment {pres[0]}-{pres[-1]}  "
                 f"major slabs {len(maj)} (nt {nt(maj)})  "
                 f"minor slabs {len(mnr)} (nt {nt(mnr)})  sticks {nearsel}")


def render(panel, script):
    """One ChimeraX run, exactly as make_zooms does it.

    ChimeraX cannot render offscreen on macOS (no OSMesa), so a window must
    exist; `open -g` keeps it off the front, and `-n` forces a new instance -
    without it macOS just activates a running ChimeraX and DROPS the --cmd, so
    nothing renders and the wait below times out.
    """
    out = ZOOM / f"zoom_{panel}.png"
    if out.exists():
        out.unlink()
    for _ in range(900):                       # never two instances at once
        if subprocess.run(["pgrep", "-x", "ChimeraX"],
                          capture_output=True).returncode != 0:
            break
        time.sleep(0.1)
    subprocess.run(["open", "-g", "-n", "-a", M.CHIMERAX_APP, "--args",
                    "--exit", "--cmd", f"open {script}"], check=False)
    for _ in range(600):                       # up to 60 s
        time.sleep(0.1)
        if subprocess.run(["pgrep", "-x", "ChimeraX"],
                          capture_output=True).returncode != 0 and out.exists():
            break
    ok = out.exists()
    if not ok:                                 # do not let it cascade
        subprocess.run(["pkill", "-x", "ChimeraX"], check=False)
    return ok


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dry", action="store_true",
                    help="write the .cxc files, run nothing")
    a = ap.parse_args()
    rows = {r["panel"]: r for r in
            csv.DictReader(open(ROOT / "inputs/map8_callouts.csv"))}
    for panel in PANELS:
        cxc, summary = script_for(panel, rows[panel])
        script = ZOOM / f"render_{panel}.cxc"
        script.write_text("\n".join(cxc) + "\n")
        print(" ", summary)
        if a.dry:
            print(f"  wrote {script} (dry run)")
            continue
        print(f"  {panel}: {'OK' if render(panel, script) else 'FAILED'} "
              f"-> {ZOOM/f'zoom_{panel}.png'}")


if __name__ == "__main__":
    main()
