"""Figure 5 panel c: four binders on the Csy4 hairpin, as SURFACES.

Both the RNA and the binder are molecular surfaces, all four in one shared
hairpin frame, in one of two cameras:

  --frame hits (default)  EXACTLY the four cameras of code/Fig5c_best_hits, so
      each surface can be set directly beside, or under, the cartoon close-up of
      the same design and the two show the same thing. The per-panel `turn y`
      angles, the residues and nucleotides each view is framed on, and the
      reference model are imported from that script, not copied.
  --frame fig4            the Fig 2b shared RNA frame plus `turn y 120`, which
      is what Fig 3b and the 32-hit grid of Fig 4 are rendered in.
  --frame axis            looking straight DOWN the hairpin's own helical axis
      from the loop end, so the head is face-on. Only useful with --cut.

and in one of two cuts:

  (default)   whole surfaces, nothing clipped.
  --cut       a SLAB A section perpendicular to the line of sight, so the
      reader looks into a cut through the complex and sees the interface inside
      it. Protein residues whose cross-section falls in the slice AND which
      touch the RNA are painted in one highlight colour.

Four decisions behind the panel, settled with Leo 2026-11-18:

1. WHICH FOUR. Flow fold-change ranks 1 to 4 from LV_32Candidates/index.tsv,
   the same four as code/Fig5c_best_hits:

       twist 23  FC 1.56  BG+LiMPNN     e8cb617b_design26_seq2
       twist 15  FC 1.54  BG+BIF        25a69495_design18_seq0
       twist 16  FC 1.47  RFD3+LiMPNN   73d23e78_design49_seq4
       twist 07  FC 1.41  BD1+LiMPNN    c5072321_design5_seq1

   CAVEAT the caption must carry: twist 15 fails the self-consistency fidelity
   gate (best-refold DockQ 0.020, the worst of the 32). Its modelled interface
   is not reproduced by either refolder, so what is drawn for it is a structure
   no predictor stands behind. Leo kept it in for the FC ranking.

2. WHICH MODEL. Each candidate PDB holds three superposed models,
   1 = design, 2 = AF3 refold, 3 = ESMFold2 refold. The one with the lowest
   self-consistency RMSD is used, never AF3 by default - AF3 is the worse model
   for three of the four. Same choice as Fig5c_best_hits, so the two panels
   agree: twist 23 -> AF3, twist 15/16/07 -> ESM.

3. WHERE THE CUT SITS (--cut only). ONE shared plane for all four. It is placed
   on the MEASURED interface: the mean depth of every protein heavy atom within
   CONTACT A of RNA, pooled over the four. Centring it on the hairpin head
   instead puts it ~8 A too far out, where three of the four panels come back
   nearly empty.

4. WHAT IS HIGHLIGHTED. Nothing, by default. `--hilite` paints the residues
   that touch the RNA (and, with --cut, only those the slice also passes
   through) in HILITE; it was tried on the panel and dropped - a second colour
   over the stack colour reads as a chemistry claim the panel cannot support,
   and the shape is the point.

Each panel shows its OWN RNA copy, never the reference one: pairing binder i
against a hairpin it was not modelled with would invent clashes and gaps at the
interface, which is the one thing this panel is read for. Only the FRAME is
shared.

House rendering recipe: code/common/CHIMERAX_RENDERING.md. Colours come from
config/thesis_colors.cxc by name; no hex is written here.

Run (as used for the figure):

    ~/.venvs/binder_map/bin/python render.py ../../output/Fig5c_cutaway_surface
    ~/.venvs/binder_map/bin/python render.py ../../output/Fig5c_cutaway_surface \
        --cut --frame axis

The first writes w23/w15/w16/w07.png, the second c23/.../c07.png, into the same
folder. --probe renders twist 23 alone at four slab depths, to re-pick SHIFT.
"""
import argparse
import html
import re
import subprocess
import sys
from pathlib import Path

import numpy as np
import biotite.structure as struc
import biotite.structure.io.pdb as pdbio

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "code"))
sys.path.insert(0, str(ROOT / "code" / "Fig2b_map_callouts"))
from render_shared_frame import (CANON, PAIRS, RING, RNA_RES,  # noqa: E402
                                 _atoms, stack_of)


def _load(name, path):
    """Import a module by FILE PATH. `import render` would pick up this file -
    the panel script next door is also called render.py - so the sibling has to
    be loaded explicitly."""
    import importlib.util
    s = importlib.util.spec_from_file_location(name, path)
    m = importlib.util.module_from_spec(s)
    s.loader.exec_module(m)
    return m


# The existing panel c's own per-panel cameras, imported rather than copied, so
# the two panels can never drift apart: (twist, model, residues, nucleotides,
# turn y, FC, dashes).
HITS_PANELS = _load("fig5c_best_hits_render",
                    ROOT / "code/Fig5c_best_hits/render.py").PANELS
HITS_ZOOM = 1.05          # ONE zoom for the whole row (see panel() below): the
                          # global `view` frames the reference hairpin, and this
                          # pulls back far enough that the widest binder, pushed
                          # off-centre by its own panel's move, still fits.
                          # Measured, not guessed: at 0.52 the four used 46 % of
                          # the frame's width and 39 % of its height, so this is
                          # the tightest the worst panel allows, less 8 %.

CHIMERAX = "/Applications/ChimeraX-1.11.1.app/Contents/MacOS/ChimeraX"
CAND = ROOT.parent / "LV_BinderDesign" / "LV_32Candidates"

# THE reference hairpin of the whole book: Fig 2b panel a's design. Every
# structure panel - Fig 2b, Fig 3b, the 32-hit grid in Fig 4 - aligns onto this
# RNA by C1' and takes its camera from it, so using anything else here would put
# these four in a frame of their own (Leo, 2026-11-18: "align them exactly like
# in the other figure 4").
REF = (ROOT / "working/structures/fig02b_examples/cifs"
       / "csy4free_RFD3_LiMPNN_len80_8b13398e_design50_seq0.cif")
FIG4_TURN = 120.0   # the `turn y` Fig 3b and the Fig 4 grid are rendered at

DESIGN, AF3, ESM = 1, 2, 3
# twist, model, flow FC
PANELS = [("23", AF3, 1.56), ("15", ESM, 1.54), ("16", ESM, 1.47),
          ("07", ESM, 1.41)]

CXC = {"BG+BIF": "thesis_bg_bif", "BG+LiMPNN": "thesis_bg_limpnn",
       "RFD3+LiMPNN": "thesis_rfd3_limpnn",
       "BD1+LiMPNN": "thesis_boltzdesign_limpnn",
       "HD+LiMPNN": "thesis_hallu_limpnn"}
RNA_S = "thesis_grey_400"       # RNA surface, neutral (COLORS.md rule 3). 300 is
                                # so pale the hairpin all but vanishes into the page.
HILITE = "thesis_petrol_800"    # contacting residues in the slice; a dark theme
                                # step, so it separates from all four stack
                                # colours by lightness as well as by hue
SIL = "thesis_grey_700"

# --- as-run geometry -------------------------------------------------------
# HEAD is the three bases that actually leave the stack, U8/A9/U10, not the
# nominal five-nucleotide loop: only those three are flipped out (ring-centroid
# stacking breaks, measured 2026-11-14). It is reported, not used to place the
# cut - see decision 3 in the docstring.
HEAD = (8, 9, 10)
SLAB = 10.0        # A, slab thickness along the view axis
SHIFT = 0.0        # A, slab centre offset from the interface centroid
CONTACT = 4.5      # A, heavy-atom cutoff for "touches the RNA"
# Surfaces are SOLID. Transparency was tried at 40/55/70/82 % on 2026-11-18 and
# rejected: past ~55 % the stack colour stops being readable and the grey RNA
# vanishes into the white page, which costs the panel both things it encodes.
# The knob stays, on --alpha.
SURF_T = 0         # % transparency of both surfaces
HILITE_T = 0       # % transparency of the contacting residues
ZOOM = {"hits": HITS_ZOOM, "fig4": 0.80, "axis": 0.62}  # `view` frames ATOMS; a
                                      # surface stands ~3 A proud and the binder
                                      # reaches past it
PROBE_SHIFTS = (-6.0, -2.0, 2.0, 6.0)

# --- the X-ray layer (--xray) ---------------------------------------------
# Panel c's OWN cartoon recipe, rendered in THIS panel's cameras, so the two
# renders register pixel for pixel and the cartoon can be dropped into a window
# cut in the surface (Leo, 2026-11-18: "you use the surface one as the
# background ... and then you just crop the figure C style one, have it on top
# ... then it looks kind of like an X-ray into the figure").
#
# It started as panel c's recipe copied setting for setting (copied, not
# imported: they are literals inside that script's build()). Two of Leo's
# calls, 2026-11-18, have moved it off c since:
#
#   "remove the dark aminoacids. They should be the same color."  c paints the
#   contacted BASES a darker grey and leaves them opaque while the rest of the
#   RNA is 60 % transparent, which inside a 13 mm window reads as dark blobs,
#   not as marked bases. There is now ONE grey across the whole hairpin. (They
#   are nucleotides, not amino acids - no side chain has been drawn in this
#   panel since the sticks came off - but the call is the same either way:
#   nothing in the window gets a second colour.)
#
#   "make the hairpin and protein in the image a bit more color intensive" -
#   every transparency cut, roughly in half. c is read at 45 mm across; this
#   layer is read at 13 mm through a window, where c's values wash out.
#
# The depth cue stays c's (on; this panel turns it off for the surfaces).
XR_RNA = "thesis_grey_400"   # = palette.NEUTRAL["400"], the hex panel c writes
XR_THICK = 0.4               # protein cartoon thickness
XR_T_PRO = 40                # % transparency, protein cartoon   (c: 65)
XR_T_RNA = 10                # % transparency, RNA cartoon       (c: 25)
XR_T_NUC = 35                # % transparency, the nucleotides   (c: 60)


def atoms(path, model):
    a = pdbio.PDBFile.read(str(path)).get_structure(model=model)
    return a[a.element != "H"]


def probe_chains(a, name):
    """Chain roles BY COMPOSITION, never by chain ID: RFD3+LiMPNN design files
    have the RNA as chain A and the protein as chain B, the other way round
    from every refold, and that is 53 % of the campaign."""
    rna_c = [str(c) for c in np.unique(a.chain_id)
             if set(map(str, np.unique(a[a.chain_id == c].res_name))) <= RNA_RES]
    pro_c = [str(c) for c in np.unique(a.chain_id) if str(c) not in rna_c]
    if len(rna_c) != 1 or len(pro_c) != 1:
        raise ValueError(f"{name}: rna={rna_c} prot={pro_c}")
    rna = a[a.chain_id == rna_c[0]]
    ids, names = struc.get_residues(rna)
    off = "".join(map(str, names)).find(CANON)
    if off < 0:
        raise ValueError(f"{name}: canonical 16-nt frame not found")
    return pro_c[0], rna_c[0], int(ids[off]), int(ids[off + 15])


def ring_centres(rna, a0):
    """Base ring centroid per canonical position 1..16."""
    return {k: rna[(rna.res_id == a0 + k - 1)
                   & np.isin(rna.atom_name, RING[CANON[k - 1]])].coord.mean(0)
            for k in range(1, 17)}


def rot_y(deg):
    c, s = np.cos(np.radians(deg)), np.sin(np.radians(deg))
    return np.array([[c, 0, s], [0, 1.0, 0], [-s, 0, c]])


def frame(rna, a0, mode, fig4_turn=None):
    """Scene -> screen transform, as the matching ChimeraX turns would leave it.

    `R` is the house upright frame, identical to
    Fig2b_map_callouts/render_shared_frame.canonical_orientation: the helix axis
    (first principal direction of the six Watson-Crick pair midpoints, signed so
    the loop is positive) onto +Y, the flipped-out U10 onto +Z. Then

      hits   nothing: every panel turns to its own angle, from Fig5c_best_hits.
      fig4   turn y 120, the rotation Fig 3b and the Fig 4 32-hit grid use.
      axis   turn x 90 tips the hairpin back so its axis points at the camera
             and the head is face-on; turn z 180 spins in-plane so U10 is up.

    Verify the axis by its RISE (~2.8-3.4 A per bp), never by which spread is
    longest: this stem is 6 bp, ~20 A long and ~18 A across, so "the longest
    axis" is not a safe test.
    """
    cen = ring_centres(rna, a0)
    mids = np.array([(cen[i] + cen[j]) / 2 for i, j in PAIRS])
    centre = np.array([cen[k] for k in range(1, 17)]).mean(0)
    axis = np.linalg.svd(mids - mids.mean(0))[2][0]
    if np.dot(cen[9] - centre, axis) < 0:
        axis = -axis
    f = cen[10] - centre
    f -= np.dot(f, axis) * axis
    f /= np.linalg.norm(f)
    R = np.vstack([np.cross(axis, f), axis, f])
    rise = float(np.ptp(mids @ axis)) / (len(PAIRS) - 1)
    if mode == "hits":
        extra, turns = np.eye(3), []       # each panel sets its own turn y
    elif mode == "fig4":
        ty = FIG4_TURN if fig4_turn is None else fig4_turn
        extra, turns = rot_y(ty), [f"turn y {ty:g}"]
    else:
        rx90 = np.array([[1.0, 0, 0], [0, 0, -1.0], [0, 1.0, 0]])
        rz180 = np.array([[-1.0, 0, 0], [0, -1.0, 0], [0, 0, 1.0]])
        extra, turns = rz180 @ rx90, ["turn x 90", "turn z 180"]
    return centre, R, extra @ R, turns, rise


def axis_angle(R):
    """Rotation matrix -> (unit axis, degrees). Quaternion route, so it stays
    stable at 180 degrees where the antisymmetric-part shortcut degenerates."""
    t = np.trace(R)
    if t > -0.9:
        w = np.sqrt(max(t + 1.0, 0.0)) / 2.0
        v = np.array([R[2, 1] - R[1, 2], R[0, 2] - R[2, 0],
                      R[1, 0] - R[0, 1]]) / (4.0 * w)
    else:                                   # near 180: take the largest column
        i = int(np.argmax(np.diag(R)))
        j, k = (i + 1) % 3, (i + 2) % 3
        s = np.sqrt(max(R[i, i] - R[j, j] - R[k, k] + 1.0, 0.0))
        v = np.zeros(3)
        v[i], v[j], v[k] = s / 2, (R[j, i] + R[i, j]) / (2 * s), \
            (R[k, i] + R[i, k]) / (2 * s)
        w = (R[k, j] - R[j, k]) / (2 * s)
    n = np.linalg.norm(v)
    if n < 1e-9:
        return np.array([0.0, 1.0, 0.0]), 0.0
    return v / n, float(np.degrees(2 * np.arctan2(n, w)))


def kabsch(P, Q):
    """Rigid transform taking P onto Q, as (R, t) with x -> R@x + t. Hand-rolled
    so this does not depend on which biotite release is installed; ChimeraX's
    own `align` does the same least-squares fit on the same ordered C1' pairs,
    so the render and these numbers agree."""
    pc, qc = P.mean(0), Q.mean(0)
    U, _, Vt = np.linalg.svd((P - pc).T @ (Q - qc))
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    R = Vt.T @ np.diag([1.0, 1.0, d]) @ U.T
    return R, qc - R @ pc


def c1(a, chain, a0):
    return np.array([a[(a.chain_id == chain) & (a.res_id == a0 + k)
                       & (a.atom_name == "C1'")].coord[0] for k in range(16)])


def build(out_dir, width, mode="hits", cut=False, probe=False, shift=SHIFT,
          zoom=None, slab=SLAB, alpha=SURF_T, hilite=False, turn=None,
          recentre=False, height=None, xray=False):
    files = [next(CAND.glob(f"Q1_{t}_*.pdb")) for t, *_ in PANELS]
    mdl = [m for _, m, _ in PANELS]
    arr = [atoms(p, m) for p, m in zip(files, mdl)]
    geo = [probe_chains(a, p.name) for a, p in zip(arr, files)]
    did = ["csy4free_" + p.name.split("csy4free_")[1][:-4] for p in files]
    zoom = ZOOM[mode] if zoom is None else zoom

    # --- the frame ----------------------------------------------------------
    # hits:      twist 23's own chosen model is the alignment target and its
    #            model 1 supplies the camera, exactly as Fig5c_best_hits does.
    # fig4/axis: the book's reference hairpin, Fig 2b panel a's design.
    if mode == "hits":
        ref_a, ref_model, cam_a = arr[0], f"#2.{mdl[0]}", atoms(files[0], DESIGN)
    else:
        ref_a, ref_model, cam_a = _atoms(REF), "#1", _atoms(REF)
    rp, rr, r0, r1 = probe_chains(ref_a, "reference")
    ref_rna = ref_a[ref_a.chain_id == rr]
    _, cc, c0, _ = probe_chains(cam_a, "camera")
    centre, R, M, turns, rise = frame(cam_a[cam_a.chain_id == cc], c0, mode, turn)
    print(f"frame {mode}: stem rise {rise:.2f} A/bp (expect 2.8-3.4)")
    ref_c1 = c1(ref_rna, rr, r0)

    # --- every complex into that frame, contacts measured there ------------
    pz, contact, xyz_all = [], [], []
    for i, (a, (pro, rna, a0, a1)) in enumerate(zip(arr, geo)):
        mob = c1(a, rna, a0)
        rot, tr = kabsch(mob, ref_c1)
        rms = np.sqrt((((rot @ mob.T).T + tr - ref_c1) ** 2).sum(1).mean())
        print(f"  twist {PANELS[i][0]}: C1' fit to the reference hairpin "
              f"{rms:.3f} A")
        xyz = (M @ (((rot @ a.coord.T).T + tr) - centre).T).T
        prot, rnam = a.chain_id == pro, a.chain_id == rna
        d = np.linalg.norm(xyz[prot][:, None, :] - xyz[rnam][None, :, :], axis=2)
        pz.append(xyz[prot])
        contact.append(d.min(1) <= CONTACT)
        xyz_all.append(xyz)

    # --- where each panel's camera is aimed, in the `hits` frame -----------
    # One scale for the whole row is set by the GLOBAL camera below; all a panel
    # gets of its own is this offset, which slides its interface to the middle
    # of the frame. See panel().
    #
    # The selection is the one the cartoon panel frames on - that design's
    # highlighted residues plus the nucleotides they read - and the offset is
    # its bounding-box centre, which is what `view` would have centred.
    ctr = []
    if mode == "hits":
        for i, (a, (pro, rna, a0, a1)) in enumerate(zip(arr, geo)):
            sel, nts = HITS_PANELS[i][2], HITS_PANELS[i][3]
            m = (((a.chain_id == pro) & np.isin(a.res_id, sel))
                 | ((a.chain_id == rna)
                    & np.isin(a.res_id, [n + a0 - 1 for n in nts])))
            # In the frame the camera will see it in: the panel loop turns the
            # models by this panel's own `turn y` before it is saved.
            q = (rot_y(HITS_PANELS[i][4]) @ xyz_all[i][m].T).T
            ctr.append(((q[:, 0].min() + q[:, 0].max()) / 2,
                        (q[:, 1].min() + q[:, 1].max()) / 2))
        print("  interface offsets " + ", ".join(
            f"twist {PANELS[i][0]} ({x:+.1f}, {y:+.1f}) A"
            for i, (x, y) in enumerate(ctr)))

    # --- where the cut sits -------------------------------------------------
    cen = ring_centres(ref_rna, r0)
    z_head = float(np.mean([(M @ (cen[k] - centre))[2] for k in HEAD]))
    z_mid = float(np.concatenate([z[c, 2] for z, c in zip(pz, contact)]).mean())
    z_lo, z_hi = z_mid + shift - slab / 2, z_mid + shift + slab / 2
    # Which residues get the highlight. With a cut it is the ones the slice
    # passes through AND that touch RNA; without one there is no slice, so it is
    # every contacting residue - the binder's whole footprint on the hairpin.
    hits = []
    for i, (a, (pro, _, _, _)) in enumerate(zip(arr, geo)):
        inslab = (pz[i][:, 2] >= z_lo) & (pz[i][:, 2] <= z_hi)
        rid = a.res_id[a.chain_id == pro]
        keep = contact[i] & inslab if cut else contact[i]
        hits.append(sorted({int(r) for r in rid[keep]}))
    if cut or probe:
        print(f"head centroid z {z_head:+.2f} A, interface centroid z "
              f"{z_mid:+.2f} A  ->  slab [{z_lo:+.2f}, {z_hi:+.2f}] "
              f"(shift {shift:+.1f})")
        for i, h in enumerate(hits):
            print(f"  twist {PANELS[i][0]}: {len(h)} residues in the slice and "
                  f"within {CONTACT} A of RNA: {h}")

    # --- the ChimeraX script ------------------------------------------------
    # #1 is the reference hairpin. It is never drawn; it is open only so the
    # four are aligned and framed exactly as every other structure panel is.
    L = [f"open {ROOT / 'config' / 'thesis_colors.cxc'}", f"open {REF}"]
    L += [f"open {p}" for p in files]
    n = len(files) + 1
    allm = f"#1-{n}"
    L += ["set bgColor white", "hide atoms", "hide cartoons", "hide surfaces",
          "~label"]
    ref = f"{ref_model}/{rr}:{r0}-{r1}"
    for i, (_, rna, a0, a1) in enumerate(geo, start=2):
        if f"#{i}.{mdl[i-2]}" == ref_model:
            continue                       # the target does not align to itself
        L.append(f"align #{i}.{mdl[i-2]}/{rna}:{a0}-{a1}@C1' toAtoms {ref}@C1' "
                 "move structures")
    for ax, d in zip("xyz", -centre):
        L.append(f"move {ax} {d:.4f} models {allm}")
    v, ang = axis_angle(R)
    L.append(f"turn {v[0]:.6f},{v[1]:.6f},{v[2]:.6f} {ang:.4f} center 0,0,0 "
             f"models {allm}")
    L += [f"{t} center 0,0,0 models {allm}" for t in turns]
    L += ["lighting soft",
          # Depth cueing fades whatever lies deeper toward the background, which
          # on a white page turns the far half of every SURFACE into a pale
          # smear and reads as a rendering fault rather than as depth. Panel c
          # leaves it on, and --xray matches c, since that layer has to look
          # like it was rendered there.
          f"lighting depthCue {'true' if xray else 'false'}",
          f"graphics silhouettes true width 1.5 color {SIL}",
          "surface cap true",
          # One camera for all four, framed on the reference hairpin, so the
          # hairpin lands on the same pixels in every panel.
          #
          # SHOW THE ATOMS BEFORE FRAMING, then hide them again. `view <spec>`
          # fails with "No displayed objects specified." when the spec has no
          # displayed ATOMS - a molecular surface over them does not count - and
          # that failure ABORTS the whole command file silently: exit status 0,
          # no message on stdout (it goes to the GUI log), every PNG unwritten.
          # It is why the --surface mode of code/Fig5c_best_hits/render.py
          # reported "wrote 0 png".
          f"show {ref} atoms", f"view {ref}", f"hide {ref} atoms",
          f"zoom {zoom:.3f}"]
    if recentre:
        # `view` centres on the reference hairpin. At a close zoom that leaves
        # the interface off to one side, so slide every model (together, so no
        # relative geometry changes) until the pooled contact centroid is in the
        # middle of the frame. x and y only: z is the cut's own axis.
        ix, iy = (np.concatenate([z[c] for z, c in zip(pz, contact)])[:, :2]
                  .mean(0))
        L += [f"move x {-ix:.4f} models {allm}", f"move y {-iy:.4f} models {allm}"]
        print(f"recentred on the interface, screen offset "
              f"({ix:+.2f}, {iy:+.2f}) A")

    # Clip to the slab. The near plane hides what lies between it and the
    # camera, the far plane what lies beyond; both are measured from the centre
    # of rotation, which `move` put at the origin, along the line of sight.
    # ChimeraX measures them AWAY from the camera, hence the sign flip from
    # screen z.
    def clip_cmd(lo, hi):
        return f"clip near {-hi:.3f} far {-lo:.3f} position 0,0,0"
    clip = clip_cmd(z_lo, z_hi) if cut else "clip off"

    def panel(i, tag, cl):
        _, mi, _ = PANELS[i]
        pro, rna, a0, a1 = geo[i]
        me = f"#{i + 2}.{mi}"
        rspec, pspec = f"{me}/{rna}:{a0}-{a1}", f"{me}/{pro}"
        cmds = [f"hide {allm} surfaces"]
        if xray:
            cmds += [f"hide {allm} cartoons", f"hide {allm} atoms", "~label"]
        cmds += [f"surface {rspec}", f"surface {pspec}",
                 f"color {rspec} {RNA_S} target s",
                 f"color {pspec} {CXC[stack_of(did[i])]} target s"]
        if alpha:
            cmds += [f"transparency {rspec} {alpha} target s",
                     f"transparency {pspec} {alpha} target s"]
        if hits[i] and hilite:
            hspec = f"{pspec}:{','.join(map(str, hits[i]))}"
            cmds += [f"color {hspec} {HILITE} target s",
                     # after the blanket transparency above, so the contacts
                     # come back solid inside the glass
                     f"transparency {hspec} {HILITE_T} target s"]
        # A zoomed panel is a full-bleed CROP, so its aspect is the saved
        # frame's, not the molecule's: --height sets how short the row is.
        hgt = f"height {height} " if height else ""
        back = []
        if mode == "hits":
            # MOVE, never re-frame. The camera is the global one, framed on the
            # reference hairpin, so the hairpin - the same molecule in all four
            # - comes out the same size in all four. The panel only slides its
            # own interface to the middle of the frame, where make_figure.py
            # cuts the X-ray window, and slides it back afterwards so the next
            # panel's `turn y center 0,0,0` still turns about the origin.
            #
            # This replaced a per-panel `view` on that design's framing
            # selection, corrected by the selection's radius. It could not hold
            # one scale: `view` fits a bounding BOX and the correction divided
            # out a centroid RADIUS, and the two disagree by up to 20 % - the
            # hairpin came out a fifth smaller in twist 7's panel than in twist
            # 23's (Leo, 2026-11-18: "check that they are visually all about the
            # same size").
            ix, iy = ctr[i]
            cmds += [f"move x {-ix:.4f} models {allm}",
                     f"move y {-iy:.4f} models {allm}"]
            back = [f"move x {ix:.4f} models {allm}",
                    f"move y {iy:.4f} models {allm}"]
        if xray:
            # Panel c's representation, swapped in AFTER the framing above, so
            # the camera is the surface run's camera to the last digit and the
            # two layers land on the same pixels. The surfaces are built,
            # coloured and framed on, then hidden: that costs one surface
            # calculation per panel and buys exact registration, which
            # re-framing on the cartoon would NOT give - `view` fits whatever
            # geometry is displayed and a surface stands ~3 A proud of it.
            #
            # No labels, and no second grey on the contacted bases: panel c
            # labels and darkens them, but c is read at 45 mm and this is read
            # at 13 mm through a window, where both turn into blobs. c, two
            # rows up, still carries the labels.
            cmds += [f"hide {allm} surfaces",
                     f"cartoon {rspec}", f"cartoon {pspec}",
                     f"cartoon style {pspec} thickness {XR_THICK}",
                     f"color {rspec} {XR_RNA}",
                     f"color {pspec} {CXC[stack_of(did[i])]}",
                     f"transparency {pspec} {XR_T_PRO} target c",
                     f"transparency {rspec} {XR_T_RNA} target c",
                     f"show {rspec} atoms", f"nucleotides {rspec} fill",
                     f"color {rspec} {XR_RNA} target ac",
                     f"transparency {rspec} {XR_T_NUC} target a"]
        return cmds + [cl, f"save {out_dir}/{tag}.png width {width} {hgt}"
                           "supersample 3 transparentBackground true"] + back

    if probe:
        for s in PROBE_SHIFTS:
            lo, hi = z_mid + s - slab / 2, z_mid + s + slab / 2
            nm = "probe_shift" + f"{s:+.0f}".replace("+", "p").replace("-", "m")
            L += panel(0, nm, clip_cmd(lo, hi))
    else:
        pre = "x" if xray else "c" if cut else "z" if recentre else "w"
        prev = 0.0
        for i, (tw, _, _) in enumerate(PANELS):
            if mode == "hits":           # each panel its own turn y, as there
                th = HITS_PANELS[i][4]
                if th != prev:
                    L.append(f"turn y {th - prev:g} center 0,0,0 models {allm}")
                    prev = th
            L += panel(i, f"{pre}{tw}", clip)
        if prev:
            L.append(f"turn y {-prev:g} center 0,0,0 models {allm}")
    L.append("exit")
    return "\n".join(L)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("out_dir")
    ap.add_argument("--width", type=int, default=2000)
    ap.add_argument("--frame", choices=("hits", "fig4", "axis"),
                    default="hits",
                    help="hits = the four cameras of Fig5c_best_hits (default); "
                         "fig4 = the Fig 3b / Fig 4-grid camera; "
                         "axis = down the hairpin's helical axis")
    ap.add_argument("--cut", action="store_true",
                    help="clip to a slab and highlight the contacts in it")
    ap.add_argument("--probe", action="store_true",
                    help="twist 23 only, at four slab depths, to re-pick SHIFT")
    ap.add_argument("--shift", type=float, default=SHIFT)
    ap.add_argument("--zoom", type=float, default=None)
    ap.add_argument("--slab", type=float, default=SLAB)
    ap.add_argument("--alpha", type=float, default=SURF_T,
                    help="%% transparency of both surfaces")
    ap.add_argument("--height", type=int, default=None,
                    help="saved pixel height; sets a cropped panel's aspect")
    ap.add_argument("--recentre", action="store_true",
                    help="centre the frame on the pooled interface, not on\n                          the hairpin; for close zooms")
    ap.add_argument("--turn", type=float, default=None,
                    help="override the fig4 `turn y` angle")
    ap.add_argument("--xray", action="store_true",
                    help="render panel c's see-through CARTOON in this panel's "
                         "cameras, as x*.png, for make_figure.py to drop into "
                         "a window cut in the surface")
    ap.add_argument("--hilite", action="store_true",
                    help="paint the RNA-contacting residues in HILITE. OFF by "
                         "default: Leo, 2026-11-18, 'remove the blue coloring "
                         "of the aminoacids' - the panel is about the shape "
                         "the binder makes on the hairpin, and a second colour "
                         "on top of the stack colour was reading as chemistry "
                         "the panel does not support.")
    a = ap.parse_args()
    out = Path(a.out_dir).resolve()
    out.mkdir(parents=True, exist_ok=True)
    tag = ("probe" if a.probe else "xray" if a.xray else
           "cutaway" if a.cut else "zoom" if a.recentre else "whole")
    script = out / f"render_{tag}.cxc"
    script.write_text(build(out, a.width, a.frame, a.cut, a.probe, a.shift,
                            a.zoom, a.slab, a.alpha, a.hilite,
                            a.turn, a.recentre, a.height, a.xray) + "\n")
    print("script ->", script)
    # Clear only THIS variant's images, so a cut-away run does not wipe the
    # uncut set sitting beside it in the same folder.
    pre = ("probe_" if a.probe else "x" if a.xray else
           "c" if a.cut else "z" if a.recentre else "w")
    for p in out.glob(f"{pre}*.png"):
        if not p.name.startswith("._"):
            p.unlink()
    # A command-file error goes to the GUI log, NOT to stdout, and it aborts the
    # file while the process still exits 0. The only way to see it from a script
    # is to dump the log afterwards as a second --cmd: that one runs even when
    # the first aborted.
    log = out / "cx_log.html"
    r = subprocess.run([CHIMERAX, "--exit", "--cmd", f"open {script}",
                        "--cmd", f"log save {log}"],
                       capture_output=True, text=True)
    for ln in (r.stdout + r.stderr).splitlines():
        if "rror" in ln:
            print("  cx|", ln)
    if log.exists():
        txt = html.unescape(re.sub(r"<[^>]+>", "\n", log.read_text()))
        bad = [l.strip() for l in txt.splitlines()
               if re.search(r"No displayed|rror|[Uu]nknown|[Ii]nvalid|"
                            r"[Mm]issing|expected", l)]
        for l in dict.fromkeys(bad):
            print("  cx|", l)
    print("wrote", len([p for p in out.glob(f"{pre}*.png")
                        if not p.name.startswith("._")]), "png")


if __name__ == "__main__":
    main()
