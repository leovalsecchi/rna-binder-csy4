#!/usr/bin/env python3
"""Flat, icon-styled render of the real Csy4 crRNA hairpin for the Fig2a slot.

The pipeline panel is drawn as flat two-tone glyphs, so the structure has to be
rendered the same way rather than as a lit cartoon: `lighting flat` (ambient
only, no specular, no shadows), one dark petrol for the backbone and one light
petrol for the bases, transparent background. That gives a real structure with
poster colours instead of a photographic one.

Geometry is the campaign target itself, the 16-nt hairpin held fixed in every
design, taken from the Fig2b shared-frame reference so this render and the
interaction panels show the same molecule in the same canonical pose (loop up,
open end down, flipped-out U10 toward the camera).

    python code/Fig2a_pipeline/render_hairpin.py            # variant sheet
    python code/Fig2a_pipeline/render_hairpin.py --pick slab_y000
"""
import argparse
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from common import palette, paths                               # noqa: E402
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "Fig2b_map_callouts"))
from render_shared_frame import (CHIMERAX, BASE, canonical_orientation,  # noqa: E402
                                 probe, _atoms, RNA_RES, CANON, PAIRS, RING)

REF = BASE / "cifs" / "csy4free_RFD3_LiMPNN_len80_8b13398e_design50_seq0.cif"
BACKBONE = "thesis_petrol_700"   # the colour the glyph it replaces already used
BASES = "thesis_petrol_300"
EDGE = "thesis_petrol_900"
STYLES = {"slab": "tube/slab", "ladder": "ladder",
          "fill": "fill/slab", "stubs": "stubs"}


def _rungs():
    """The six Watson-Crick C1'-C1' vectors in the canonical frame.

    Same frame `canonical_orientation` establishes, recomputed here rather than
    returned by it, so that function keeps doing one thing. These vectors are
    the rungs: their angle on screen is what "horizontal bases" means.
    """
    import numpy as np
    import biotite.structure as struc
    a = _atoms(REF)
    rc = [str(c) for c in np.unique(a.chain_id)
          if set(map(str, np.unique(a[a.chain_id == c].res_name))) <= RNA_RES][0]
    rna = a[a.chain_id == rc]
    ids, names = struc.get_residues(rna)
    off = "".join(map(str, names)).find(CANON)
    cen, c1 = {}, {}
    for k in range(1, 17):
        r = rna[rna.res_id == ids[off + k - 1]]
        cen[k] = r[np.isin(r.atom_name, RING[CANON[k - 1]])].coord.mean(0)
        c1[k] = r[r.atom_name == "C1'"].coord[0]
    mids = np.array([(cen[i] + cen[j]) / 2 for i, j in PAIRS])
    centre = np.array([cen[k] for k in range(1, 17)]).mean(0)
    axis = np.linalg.svd(mids - mids.mean(0))[2][0]
    if np.dot(cen[9] - centre, axis) < 0:
        axis = -axis
    f = cen[10] - centre
    f -= np.dot(f, axis) * axis
    f /= np.linalg.norm(f)
    R = np.vstack([np.cross(axis, f), axis, f])
    V = np.array([R @ (c1[j] - c1[i]) for i, j in PAIRS])
    return V / np.linalg.norm(V, axis=1)[:, None]


def level_roll(tilt, az, _V=[]):
    """The in-plane roll that levels the rungs at this tilt and azimuth.

    The hairpin's own A-form twist and bend leave the six pairs spread over
    about 35 degrees, so no camera makes them all horizontal; this levels their
    MEAN, and the stem leans by the same amount. Positive roll is clockwise.
    """
    import numpy as np
    if not _V:
        _V.append(_rungs())
    rx, ry = np.radians(tilt), np.radians(az)
    Rx = np.array([[1, 0, 0], [0, np.cos(rx), -np.sin(rx)], [0, np.sin(rx), np.cos(rx)]])
    Ry = np.array([[np.cos(ry), 0, np.sin(ry)], [0, 1, 0], [-np.sin(ry), 0, np.cos(ry)]])
    W = (Ry @ Rx @ _V[0].T).T
    ang = np.degrees(np.arctan2(W[:, 1], W[:, 0]))
    ang = np.where(ang > 90, ang - 180, np.where(ang < -90, ang + 180, ang))
    return -float(ang.mean())


def build(out, width, angles, styles, tilts, roll, spins, probe_nt=None):
    pro, rna, r0, r1 = probe(REF)
    spec = f"#1/{rna}:{r0}-{r1}"
    # A second render of the SAME frame showing only these target nucleotides,
    # so a caller can measure where they land on the published PNG instead of
    # eyeballing it. Fig 2a uses it to put a bracket over the three flipped
    # loop bases; its alpha bbox and the full render's are in one pixel frame.
    pspec = (f"#1/{rna}:{r0 + probe_nt[0] - 1}-{r0 + probe_nt[-1] - 1}"
             if probe_nt else None)
    L = [f"open {REF}", f"open {paths.ROOT}/config/thesis_colors.cxc",
         "set bgColor white", "hide atoms", "hide cartoons",
         f"cartoon {spec}", f"show {spec} atoms"]

    centre, v, ang = canonical_orientation(REF)
    for ax, d in zip("xyz", -centre):
        L.append(f"move {ax} {d:.4f} models #1")
    L.append(f"turn {v[0]:.6f},{v[1]:.6f},{v[2]:.6f} {ang:.4f} center 0,0,0 models #1")

    # Orthographic, so the rungs stay parallel instead of converging: a glyph
    # has no vanishing point, and perspective is what made the first pass read
    # as a photograph of a molecule.
    L += ["lighting flat", f"graphics silhouettes true width 3 color {EDGE}",
          "camera ortho", "view", "zoom 0.95"]
    for st in styles:
        L.append(f"nucleotides {spec} {STYLES[st]}")
        L.append(f"color {spec} {BACKBONE}")
        L.append(f"color {spec} {BASES} target a")
        # tilt about x tips the base planes toward the camera, which is what
        # turns edge-on rungs into visible faces; azimuth then spins about the
        # hairpin axis. Both are undone before the next tilt so every frame
        # starts from the canonical pose.
        for t in tilts:
            if t:
                L.append(f"turn x {t}")
            # absolute azimuth, so a single-angle run is that angle and not
            # zero: turns are relative, so each frame moves by the step from
            # the previous one and the first frame moves by the angle itself.
            for i, a in enumerate(angles):
                step = a - (angles[i - 1] if i else 0)
                if step:
                    L.append(f"turn y {step}")
                # roll is per frame: the angle that levels the rungs depends on
                # the azimuth, so "auto" is a different number in every frame
                # and has to be undone before the next turn y.
                r = level_roll(t, a) if roll == "auto" else float(roll)
                tag = f"_r{r:+.1f}" if r else ""
                if r:
                    L.append(f"turn z {r:.2f}")
                # spin about SCREEN x, which after the roll IS the rung
                # direction: rotating about a horizontal screen axis leaves
                # horizontal vectors horizontal, so the rungs stay level while
                # the hairpin turns through the base-pair plane and the faces
                # come round toward the camera.
                for sp in spins:
                    if sp:
                        L.append(f"turn x {sp}")
                    stag = f"_s{sp:03d}" if sp else ""
                    L.append(f"save {out}/hairpin_{st}_x{t:+04d}_y{a:03d}"
                             f"{tag}{stag}.png width {width} supersample 3 "
                             "transparentBackground true")
                    if pspec:
                        L += [f"hide {spec} atoms", f"hide {spec} cartoons",
                              f"show {pspec} atoms",
                              f"save {out}/hairpin_probe.png width {width} "
                              "supersample 3 transparentBackground true",
                              f"show {spec} atoms", f"cartoon {spec}"]
                    if sp:
                        L.append(f"turn x {-sp}")
                if r:
                    L.append(f"turn z {-r:.2f}")
            if angles[-1] % 360:
                L.append(f"turn y {-angles[-1] % 360}")
            if t:
                L.append(f"turn x {-t}")
    L.append("exit")
    return "\n".join(L)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=None)
    ap.add_argument("--width", type=int, default=1400)
    ap.add_argument("--angles", default="0,30,60,90")
    ap.add_argument("--styles", default="slab,ladder")
    ap.add_argument("--tilts", default="0")
    ap.add_argument("--spins", default="0",
                    help="rotation about the rung axis, in degrees, applied "
                         "after the roll; the rungs stay horizontal throughout")
    ap.add_argument("--probe", default=None,
                    help="target nucleotides (e.g. 8-10) to render ALONE into "
                         "hairpin_probe.png, same camera, for measuring where "
                         "they sit on the published render")
    ap.add_argument("--roll", default="0",
                    help="in-plane roll in degrees, or 'auto' to level the "
                         "rungs at each azimuth")
    a = ap.parse_args()
    out = Path(a.out) if a.out else paths.fig_output("Fig2a_pipeline") / "hairpin"
    out.mkdir(parents=True, exist_ok=True)
    angles = [int(x) for x in a.angles.split(",")]
    script = out / "render_hairpin.cxc"
    tilts = [int(x) for x in a.tilts.split(",")]
    pn = [int(x) for x in a.probe.split("-")] if a.probe else None
    script.write_text(build(out, a.width, angles, a.styles.split(","), tilts, a.roll,
                            [int(x) for x in a.spins.split(",")], pn) + "\n")
    r = subprocess.run([CHIMERAX, "--exit", "--cmd", f"open {script}"],
                       capture_output=True, text=True)
    for ln in (r.stdout + r.stderr).strip().splitlines()[-15:]:
        print("  cx|", ln)
    pngs = sorted(p for p in out.glob("hairpin_*.png") if not p.name.startswith("._"))
    print(f"wrote {len(pngs)} png (expect "
          f"{len(angles)*len(tilts)*len(a.spins.split(','))*len(a.styles.split(','))})")


if __name__ == "__main__":
    main()
