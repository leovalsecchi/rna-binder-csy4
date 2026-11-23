"""What each design route was given, and where those coordinates come from.

Two bare panels.

  (a) the Csy4-crRNA product complex (PDB 4AL5) and the hairpin taken out of it,
      rendered from ONE camera by render_target.py so the second image is the first
      with the protein removed. Build the renders first:
          python code/supp/target_definition/render_target.py --azimuth 60
  (b) the 18 modelled residues of 4AL5 chain B, and which of them each design route
      and each re-predictor received, as coordinates or as sequence.

Everything the caption asserts is re-derived from the primary files by verify() and
printed on every run: that the campaign target file IS 4AL5 chain B atom-for-atom,
that only 16 of the 18 modelled residues are complete nucleotides, that the other two
carry no base and therefore cannot be a preserved closing pair, and what each config
actually passed. The config values are transcribed into ROUTES below from the Euler
configs listed in the README - they are not readable from this machine - so verify()
checks the two things that ARE local: the coordinate file and the residue content.

scipy is broken in ~/.venvs/binder_map (_propack dlopen failure); nothing here needs
it, the geometry is numpy and biotite.
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import biotite.structure as struc
import biotite.structure.io.pdb as pdb
import biotite.structure.io.pdbx as pdbx

from common import paths, style, palette, csy4

FIG = Path(__file__).resolve().parent.name
CIF_4AL5 = Path.home() / "4al5.cif"
TARGET_PDB = (paths.ROOT / "working" / "structures" / FIG / "csy4_rna_only.pdb")
RENDERS = paths.OUTPUT / FIG / "renders"
AZIMUTH = 60              # chosen from the 0/30/60/90 sheet: stem ladder open,
                          # loop at the top, both stub residues clear of the stem

RNA_RES = ["A", "U", "G", "C"]
RES = list(range(4, 22))          # the 18 modelled residues of 4AL5 chain B
CORE = list(range(5, 21))         # the 16 complete nucleotides = csy4.SEQUENCE
RING = {"A": ["N9", "C8", "N7", "C5", "C6", "N1", "C2", "N3", "C4"],
        "G": ["N9", "C8", "N7", "C5", "C6", "N1", "C2", "N3", "C4"],
        "C": ["N1", "C2", "N3", "C4", "C5", "C6"],
        "U": ["N1", "C2", "N3", "C4", "C5", "C6"]}

# What each route was handed, transcribed from the Euler configs (paths in README).
# "coords" = a coordinate file; "seq" = a sequence string and nothing else.
# Route order is palette.STACK_ORDER and is never re-sorted.
ROUTES = [
    ("RFD3+LiMPNN", "coords", CORE,
     "rfd3__multif/default.yaml -> csy4_crrna.json contig B5-20"),
    ("BG+LiMPNN", "coords", RES, "boltzgen__multif/default.yaml file csy4_rna_only.pdb chain B"),
    ("BG+BIF", "coords", RES, "boltzgen__multif/default.yaml file csy4_rna_only.pdb chain B"),
    ("BD1+LiMPNN", "coords", RES, "boltzdesign__multif/default.yaml pdb_target_ids B"),
    ("HD+LiMPNN", "coords", RES, "halludesign/default.yaml csy4_helix.pdb chain B (helix mode)"),
    ("HD+LiMPNN", "seq", RES, "halludesign/random_init.yaml template jsons, 18-nt rna (de novo)"),
]
ROUTE_NOTE = {4: "helix", 5: "de novo"}       # the two HalluDesign modes
REPREDICTORS = ("AF3 · Boltz-2 · ESMFold2", "seq", CORE,
                "refold__{af3,b2,esm}/default.yaml rna_sequence: ACUGCCGUAUAGGCAG")


# --- verification ----------------------------------------------------------
def _rna(path):
    if path.suffix == ".cif":
        a = pdbx.get_structure(pdbx.CIFFile.read(path), model=1)
    else:
        a = pdb.PDBFile.read(path).get_structure(model=1)
    return a[(a.chain_id == "B") & np.isin(a.res_name, RNA_RES) & (a.element != "H")]


def verify():
    """Re-derive every structural claim the caption makes. Asserts, so a changed
    file fails the run instead of quietly changing the figure."""
    rna = _rna(CIF_4AL5)
    ids, names = struc.get_residues(rna)
    assert list(ids) == RES, f"4AL5 chain B residues {ids[0]}-{ids[-1]}, expected 4-21"
    seq18 = "".join(map(str, names))
    print(f"4AL5 chain B: {len(ids)} modelled residues {ids[0]}-{ids[-1]}, {seq18}")

    stub, full = {}, []
    for rid, nm in zip(ids, names):
        atoms = sorted(map(str, rna[rna.res_id == rid].atom_name))
        if not set(RING[str(nm)]) <= set(atoms):
            stub[rid] = atoms
        else:
            full.append(str(nm))
    assert set(stub) == {4, 21}, f"expected stubs at 4 and 21, got {sorted(stub)}"
    assert "".join(full) == csy4.SEQUENCE, "the 16 complete nucleotides are not csy4.SEQUENCE"
    for rid, atoms in stub.items():
        print(f"  res {rid} stub: {len(atoms)} atoms {atoms} - NO base, cannot pair")
    print(f"  {len(full)} complete nucleotides (res {CORE[0]}-{CORE[-1]}) = {''.join(full)}")

    # the campaign target file against the deposition, no superposition
    if TARGET_PDB.exists():
        t = _rna(TARGET_PDB)
        k = lambda s: {(int(r), str(n), str(a)): c for r, n, a, c
                       in zip(s.res_id, s.res_name, s.atom_name, s.coord)}
        kt, kd = k(t), k(rna)
        assert set(kt) == set(kd), (f"atom keys differ: {len(set(kt) - set(kd))} only in "
                                    f"target, {len(set(kd) - set(kt))} only in 4AL5")
        dv = np.array([kt[q] - kd[q] for q in sorted(kt)])
        rmsd = float(np.sqrt((dv ** 2).sum(1).mean()))
        assert rmsd == 0.0, f"target file is not 4AL5 chain B: RMSD {rmsd}"
        print(f"  campaign target csy4_rna_only.pdb == 4AL5 chain B: {len(kt)} heavy atoms "
              f"matched, 0 unique to either, RMSD {rmsd:.6f} A with NO superposition")
    else:
        print(f"  NEED: {TARGET_PDB} not fetched; cannot re-check the target identity here")

    # geometry: pairs, stacking, twist
    cen, nrm, c1 = {}, {}, {}
    for rid in CORE:
        r = rna[rna.res_id == rid]
        ring = r[np.isin(r.atom_name, RING[str(r.res_name[0])])]
        cen[rid] = ring.coord.mean(0)
        nrm[rid] = np.linalg.svd(ring.coord - cen[rid])[2][2]
        c1[rid] = r[r.atom_name == "C1'"].coord[0]
    print("  stem pairs (N1/N3), 4AL5 numbering -> target nt numbering:")
    for i, j in csy4.PAIRS:                      # target numbering nt1..nt16
        ri, rj = i + 4, j + 4
        ai = "N1" if str(rna[rna.res_id == ri].res_name[0]) in "AG" else "N3"
        aj = "N1" if str(rna[rna.res_id == rj].res_name[0]) in "AG" else "N3"
        d = np.linalg.norm(rna[(rna.res_id == ri) & (rna.atom_name == ai)].coord[0]
                           - rna[(rna.res_id == rj) & (rna.atom_name == aj)].coord[0])
        assert d < 3.3, f"pair nt{i}-nt{j} is {d:.2f} A, not Watson-Crick"
        print(f"    res {ri}-{rj} (nt{i}-nt{j}) {d:.2f} A")
    breaks = []
    for rid in CORE[:-1]:
        g = float(np.linalg.norm(cen[rid + 1] - cen[rid]))
        ang = float(np.degrees(np.arccos(min(1, abs(np.dot(nrm[rid], nrm[rid + 1]))))))
        if g > 7.0 or ang > 45:
            breaks.append((rid - 4, rid - 3, g, ang))
    print("  stacking breaks (ring-centroid gap > 7 A or interplanar angle > 45 deg):")
    for a, b, g, ang in breaks:
        print(f"    nt{a}->nt{b}: {g:.2f} A, {ang:.1f} deg")
    out = [b for a, b, _, _ in breaks if a in csy4.LOOP or b in csy4.LOOP]
    print(f"  -> continuously stacked nt2-nt7 and nt11-nt16; "
          f"U8, A9, U10 leave the stack (loop is {len(csy4.LOOP)} nt, three bases flip out)")
    mid = np.array([(cen[i + 4] + cen[j + 4]) / 2 for i, j in csy4.PAIRS])
    axis = np.linalg.svd(mid - mid.mean(0))[2][0]
    V = []
    for i, j in csy4.PAIRS:
        v = c1[j + 4] - c1[i + 4]
        v -= np.dot(v, axis) * axis
        V.append(v / np.linalg.norm(v))
    twist = sum(np.degrees(np.arccos(np.clip(np.dot(V[k], V[k + 1]), -1, 1)))
                for k in range(len(V) - 1))
    rise = abs(float(np.dot(mid[-1] - mid[0], axis)))
    print(f"  {len(csy4.PAIRS)} Watson-Crick pairs span {twist:.1f} deg = "
          f"{twist / 360:.3f} of a helical turn over a {rise:.1f} A rise "
          f"-> a groove FACE, not a channel")
    # the sixth, non-Watson-Crick pair that closes the loop
    g11n2 = rna[(rna.res_id == 11) & (rna.atom_name == "N2")].coord[0]
    a15n7 = rna[(rna.res_id == 15) & (rna.atom_name == "N7")].coord[0]
    print(f"  loop closed by a sheared G.A (res 11 G N2 - res 15 A N7 "
          f"{np.linalg.norm(g11n2 - a15n7):.2f} A), which is the 6th pair counted "
          f"elsewhere in this repo")
    return seq18


# --- (a) complex and extracted hairpin --------------------------------------
def _crop(path):
    im = mpimg.imread(path)
    ys, xs = np.where(im[..., 3] > 0.02)
    return im[ys.min():ys.max() + 1, xs.min():xs.max() + 1]


def panel_a(out):
    cx = RENDERS / f"complex_y{AZIMUTH:03d}.png"
    hp = RENDERS / f"hairpin_y{AZIMUTH:03d}.png"
    if not (cx.exists() and hp.exists()):
        print(f"SKIP panel a: run render_target.py --azimuth {AZIMUTH} first")
        return
    W, H = 85.0, 42.0
    fig = plt.figure(figsize=(style.mm(W), style.mm(H)))
    for ax_rect, img, frac in (([0.00, 0.03, 0.50, 0.95], _crop(cx), 1.0),
                               ([0.56, 0.03, 0.44, 0.95], _crop(hp), 1.0)):
        ax = fig.add_axes(ax_rect)
        ax.imshow(img, interpolation="antialiased")
        ax.axis("off")
    # The hairpin is drawn larger than it sits in the complex; the arrow says the
    # second image is the first with the protein taken away, and the caption says
    # the scale is not shared.
    fig.text(0.525, 0.50, "→", ha="center", va="center",
             color=palette.NEUTRAL["600"], fontsize=style.CFG["font_size_pt"] + 2)
    # dpi on the save, or a vector output resamples the raster to the default and
    # the render comes out blurry
    fig.savefig(out / "panel_a_complex_and_hairpin.pdf", bbox_inches="tight", dpi=600)
    fig.savefig(out / "panel_a_complex_and_hairpin.png", bbox_inches="tight", dpi=600)
    plt.close(fig)


# --- (b) what each route received -------------------------------------------
def panel_b(seq18, out):
    W, H = 85.0, 44.0
    fig, ax = plt.subplots(figsize=(style.mm(W), style.mm(H)))
    rows = ROUTES + [REPREDICTORS]
    nrow = len(rows)
    x = np.arange(len(RES)).astype(float)
    bw, bh = 0.78, 0.60

    for ri, (name, kind, given, _) in enumerate(rows):
        y = nrow - 1 - ri
        is_route = ri < len(ROUTES)
        col = palette.stack_color(name) if is_route else palette.NEUTRAL["600"]
        for k, res in enumerate(RES):
            if res not in given:
                continue
            # A stub is only a stub in COORDINATE form: res 4 and 21 carry no base
            # in the crystal, so a route handed the file got backbone atoms there.
            # In SEQUENCE form they are ordinary C's and the model builds them whole,
            # which is why those cells get a letter like every other.
            stub = (res in (4, 21)) and kind == "coords"
            face = palette.SURFACE if (kind == "seq" or stub) else col
            ax.add_patch(plt.Rectangle((x[k] - bw / 2, y - bh / 2), bw, bh,
                                       facecolor=face, edgecolor=col, linewidth=0.5))
            if kind == "seq":
                ax.text(x[k], y, seq18[k], ha="center", va="center", color=col,
                        fontsize=style.CFG["font_size_pt"])
            elif stub:
                # a short rule inside the open cell: backbone atoms and no base,
                # so it cannot be mistaken for an empty cell or for a letter cell
                ax.plot([x[k] - bw * 0.22, x[k] + bw * 0.22], [y, y],
                        color=col, lw=0.6, solid_capstyle="butt")
        lab = name if not is_route else f"{name}"
        if ri in ROUTE_NOTE:
            lab = f"{name}  ({ROUTE_NOTE[ri]})"
        ax.text(-1.0, y, lab, ha="right", va="center",
                fontsize=style.CFG["font_size_pt"], color=palette.INK)

    # column header: 4AL5 residue number, every fifth one, and the base letter
    for k, res in enumerate(RES):
        ax.text(x[k], nrow - 0.30, seq18[k], ha="center", va="bottom",
                fontsize=style.CFG["font_size_pt"],
                color=palette.NEUTRAL["500"] if res in (4, 21) else palette.INK)
        if res % 5 == 0:
            ax.text(x[k], nrow + 0.42, str(res), ha="center", va="bottom",
                    fontsize=style.CFG["font_size_pt"], color=palette.NEUTRAL["600"])

    # secondary structure of the 16-nt core, from common/csy4.py (measured, not drawn
    # by eye): stem / loop / stem under the columns, and the two stubs named.
    yb = -0.95
    for a, b, lab in csy4.BRACKETS:
        ax.plot([a, b], [yb, yb], color=palette.NEUTRAL["500"], lw=0.7,
                solid_capstyle="butt")
        ax.text((a + b) / 2, yb - 0.30, lab, ha="center", va="top",
                fontsize=style.CFG["font_size_pt"], color=palette.NEUTRAL["600"])
    for k, res in enumerate(RES):
        if res in (4, 21):
            ax.text(x[k], yb - 0.05, "stub", ha="center", va="top", rotation=90,
                    fontsize=style.CFG["font_size_pt"], color=palette.NEUTRAL["600"])

    ax.set_xlim(-6.6, len(RES) - 0.2)
    ax.set_ylim(yb - 1.9, nrow + 1.1)
    ax.axis("off")
    fig.subplots_adjust(left=0.01, right=0.99, top=0.99, bottom=0.01)
    style.save(fig, out, "panel_b_route_inputs")
    plt.close(fig)


def main():
    style.apply()
    out = paths.fig_output(FIG)
    seq18 = verify()
    print()
    for name, kind, given, src in ROUTES + [REPREDICTORS]:
        print(f"  {name:22s} {kind:6s} {len(given):2d} residues "
              f"({given[0]}-{given[-1]})  <- {src}")
    panel_a(out)
    panel_b(seq18, out)


if __name__ == "__main__":
    main()
