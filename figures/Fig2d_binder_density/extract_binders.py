#!/usr/bin/env python
"""
extract_binders.py
==================

Stage 1 of Laetitia Lehmann's binder-density pipeline, rebuilt for the csy4 RNA
campaign.  Her original (same name, run on Euler) is not on the NAS and her Euler
home is not readable, so this writes the same .npz contract her downstream
scripts read:

    binder_density_map.py      -> ball clouds (.pdb with dummy atoms)
    binder_density_volume.py   -> .mrc density maps, contour level == percentage

What it does
------------
For every design: superpose the design's RNA TARGET onto the reference target,
apply that same transform to the design's BINDER, and keep the binder's heavy
atoms.  All binders of a group therefore end up in one common frame - the
reference target's own frame - which is the whole point: the density is only
meaningful if every binder is measured in the same coordinate system.

The RNA core is located BY SEQUENCE (str.find of the 16-mer), not by residue
number, because the stacks disagree on numbering: RFD3+LiMPNN 1-16,
BG+* 2-17, BD1/HD+LiMPNN 1-18 (an 18-mer that contains the 16-mer at offset 1).
Chains are picked BY COMPOSITION via rna_binder's split_chains, never by chain
ID: RFD3+LiMPNN design CIFs have chain A = RNA, inverted from every other stack.

Output naming
-------------
Her campaign parser wants <pid>_<accession>_<domain>_v<n>, and only accepts
D<n>, "t" or "f" as the domain token.  One target with several stacks therefore
becomes one pid with one D<n> per stack:

    1_csy4_D1_v1.npz   RFD3+LiMPNN        (whatever order the groups come in)
    1_csy4_D2_v1.npz   BG+LiMPNN
    ...

The D<n> -> group-name mapping is written to legend.csv beside the .npz files,
because nothing downstream knows the real names.

Usage
-----
    python extract_binders.py \
        --designs passing_v6_final.csv --path-col design_cif --group-col stack \
        --ref /path/to/rna-binder-csy4/targets/csy4_crrna_free/csy4_rna_only.pdb \
        --out ./npz --max-per-group 3000

    python extract_binders.py --self-check      # no data needed
"""

import argparse
import csv
import importlib.util
import os
import random
import sys

import numpy as np

# The csy4 target as every refold and every design carries it (16 nt).
CORE_SEQ = "ACUGCCGUAUAGGCAG"

# Canonical stack order and colours, copied from LV_ThesisFigures/config/palette.json
# v0.5 ("stacks").  That file is the source of truth; if it changes, change this.
# Groups are numbered D1, D2, ... in THIS order, so the figure colour and the
# domain label always agree.  Short aliases because CSVs spell stacks either way.
STACK_ORDER = [
    ("RFD3 + LigandMPNN",
     ("RFD3", "RFD3+LiMPNN", "RFD3_LiMPNN", "rfd3"), "#a83902"),
    ("BoltzGen + LigandMPNN",
     ("BG+LiMPNN", "BG_LiMPNN", "BoltzGen", "bg", "BG"), "#a376e9"),
    ("BoltzGen + BoltzIF",
     ("BG+BIF", "BG_BIF", "bif"), "#6938a7"),
    ("BoltzDesign + LigandMPNN",
     ("BD1", "BD1+LiMPNN", "BD1_LiMPNN", "BoltzDesign", "bd1"), "#bd8703"),
    ("HalluDesign + LigandMPNN",
     ("HD", "HD+LiMPNN", "HD_LiMPNN", "HalluDesign", "hd"), "#d35182"),
]
UNKNOWN_COLOUR = "#919699"          # palette.json neutral.500

# Design files are NOT comparable atom for atom: RFD3+LiMPNN, BG+LiMPNN and BD1
# store the backbone alone (N, CA, C, O - 4 atoms/residue), while BG+BIF (7.8)
# and HD (7.2) carry side chains.  A full-atom cloud is simply fatter than a
# backbone one, so mixing them would compare file formats rather than binder
# placement.  Backbone atoms are the common denominator: that is the default.
BACKBONE = ("N", "CA", "C", "O")


def rank_and_colour(name):
    """(canonical position, colour) for a group name; unknown names sort last."""
    for i, (canon, aliases, hexcol) in enumerate(STACK_ORDER):
        if name == canon or name in aliases:
            return i, hexcol
    return len(STACK_ORDER), UNKNOWN_COLOUR

DEFAULT_SRC = os.path.expanduser(
    "~/rna-binder-csy4/pipeline/src/rna_binder/analysis/structural.py")


def load_helpers(src):
    """load_structure + split_chains from rna_binder, by file path.

    By path rather than by package import so that rna_binder.analysis.__init__
    (pandas, the whole aggregate stack) is not dragged in.
    """
    spec = importlib.util.spec_from_file_location("_rb_structural", src)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod.load_structure, mod.split_chains


# --------------------------------------------------------------------------
# geometry
# --------------------------------------------------------------------------

def kabsch(mobile, fixed):
    """Rotation+translation taking `mobile` onto `fixed`: x -> R @ x + t."""
    cm, cf = mobile.mean(0), fixed.mean(0)
    u, _, vt = np.linalg.svd((mobile - cm).T @ (fixed - cf))
    d = np.sign(np.linalg.det(vt.T @ u.T))
    r = vt.T @ np.diag([1.0, 1.0, d]) @ u.T
    return r, cf - r @ cm


def apply_rt(r, t, coords):
    return coords @ r.T + t


def rmsd(a, b):
    return float(np.sqrt(np.mean(np.sum((a - b) ** 2, axis=1))))


# --------------------------------------------------------------------------
# RNA core location
# --------------------------------------------------------------------------

def residue_order(chain):
    """(sequence string, residue ids) for one RNA chain, in residue order."""
    ids = sorted(set(int(r) for r in chain.res_id))
    seq = []
    for r in ids:
        names = chain.res_name[chain.res_id == r]
        seq.append(str(names[0]))
    return "".join(seq), ids


def core_residues(chain, core_seq=CORE_SEQ):
    """Residue ids of the core 16-mer inside an RNA chain, or None."""
    seq, ids = residue_order(chain)
    off = seq.find(core_seq)
    if off < 0:
        return None
    return ids[off:off + len(core_seq)]


def core_atoms(chain, ids):
    """{(core position, atom name): coord} for the core residues of a chain."""
    out = {}
    for pos, rid in enumerate(ids):
        m = chain.res_id == rid
        for name, xyz in zip(chain.atom_name[m], chain.coord[m]):
            out.setdefault((pos, str(name)), np.asarray(xyz, dtype=np.float64))
    return out


def matched(a, b):
    """Two coordinate arrays over the atoms present in both, same order."""
    keys = sorted(set(a) & set(b))
    if not keys:
        return None, None, keys
    return (np.array([a[k] for k in keys]),
            np.array([b[k] for k in keys]), keys)


def pick_atoms(atoms, mode):
    """Backbone atoms (the default, and comparable across stacks) or all heavy."""
    if mode == "backbone":
        return atoms[np.isin(atoms.atom_name, BACKBONE)]
    return atoms[np.array([str(e).upper() != "H" for e in atoms.element])]


# --------------------------------------------------------------------------
# one design
# --------------------------------------------------------------------------

def binder_in_reference_frame(path, ref_core, load_structure, split_chains,
                              mode="backbone"):
    """Binder atom coords moved into the reference target's frame.

    Returns (coords, fit_rmsd, n_matched) or raises ValueError with a reason.
    """
    atoms = load_structure(path)
    protein, rna = split_chains(atoms)
    if protein is None or rna is None:
        raise ValueError("no protein/RNA pair in file")
    ids = core_residues(rna)
    if ids is None:
        seq, _ = residue_order(rna)
        raise ValueError("core 16-mer not found in RNA chain (%s)" % seq)
    mob, fix, keys = matched(core_atoms(rna, ids), ref_core)
    if mob is None or len(keys) < 3:
        raise ValueError("fewer than 3 target atoms shared with the reference")
    r, t = kabsch(mob, fix)
    fit = rmsd(apply_rt(r, t, mob), fix)
    coords = apply_rt(r, t, pick_atoms(protein, mode).coord.astype(np.float64))
    if len(coords) == 0:
        raise ValueError("binder has no %s atoms" % mode)
    return coords, fit, len(keys)


# --------------------------------------------------------------------------
# a group -> one .npz
# --------------------------------------------------------------------------

def write_npz(path, coords_list, fits, ref_core, ref_core_ids):
    """The .npz contract binder_density_map.py / _volume.py read.

    domain_* describe the target, written as the REFERENCE's own atoms, so her
    frame check (worst atom <= 0.01 A against the raw reference file) is exact.
    """
    offsets = np.zeros(len(coords_list) + 1, dtype=np.int64)
    offsets[1:] = np.cumsum([len(c) for c in coords_list])
    keys = sorted(ref_core)
    fits = np.asarray(fits, dtype=np.float64)
    np.savez(
        path,
        coords=np.concatenate(coords_list).astype(np.float64),
        offsets=offsets,
        rmsd=fits,
        rmsd_whole_domain=fits,
        domain_resseq=np.asarray(ref_core_ids, dtype=np.int64),
        domain_heavy=np.array([ref_core[k] for k in keys], dtype=np.float64),
        domain_heavy_resseq=np.array([ref_core_ids[k[0]] for k in keys],
                                     dtype=np.int64),
        domain_heavy_name=np.array([k[1] for k in keys], dtype="S4"),
        fit_mode=np.array("global"),
        fit_radius=np.array(0.0),
    )


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--designs", help="CSV of designs")
    ap.add_argument("--path-col", default="design_cif",
                    help="column holding the design structure path")
    ap.add_argument("--group-col", default="stack",
                    help="column to group by; one .npz per value")
    ap.add_argument("--ref", help="reference target structure (.pdb/.cif)")
    ap.add_argument("--out", help="output directory for the .npz files")
    ap.add_argument("--pid", default="1", help="target number in the campaign name")
    ap.add_argument("--accession", default="csy4")
    ap.add_argument("--max-per-group", type=int, default=3000,
                    help="random subsample; 0 = all. A density needs thousands, "
                         "not millions, and percentage_field loops per binder.")
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--atoms", choices=("backbone", "heavy"), default="backbone",
                    help="backbone (N,CA,C,O) is the default and the only "
                         "cross-stack-comparable choice; see BACKBONE above")
    ap.add_argument("--src", default=DEFAULT_SRC,
                    help="path to rna_binder/analysis/structural.py")
    ap.add_argument("--self-check", action="store_true")
    args = ap.parse_args(argv)

    if args.self_check:
        return self_check()
    for req in ("designs", "ref", "out"):
        if not getattr(args, req):
            ap.error("--%s is required (or use --self-check)" % req)

    load_structure, split_chains = load_helpers(args.src)

    ref = load_structure(args.ref)
    _, ref_rna = split_chains(ref)
    if ref_rna is None:
        ap.error("no RNA chain in the reference %s" % args.ref)
    ref_ids = core_residues(ref_rna)
    if ref_ids is None:
        seq, _ = residue_order(ref_rna)
        ap.error("core 16-mer %s not in the reference (%s)" % (CORE_SEQ, seq))
    ref_core = core_atoms(ref_rna, ref_ids)
    print("reference %s: core residues %d-%d, %d atoms"
          % (os.path.basename(args.ref), ref_ids[0], ref_ids[-1], len(ref_core)))

    groups = {}
    with open(args.designs, newline="") as fh:
        for row in csv.DictReader(fh):
            if args.path_col not in row:
                ap.error("no column %r in %s" % (args.path_col, args.designs))
            groups.setdefault(row.get(args.group_col) or "all", []).append(
                row[args.path_col])

    os.makedirs(args.out, exist_ok=True)
    rng = random.Random(args.seed)
    legend = []
    ordered = sorted(groups, key=lambda g: (rank_and_colour(g)[0], g))
    for n, name in enumerate(ordered, start=1):
        paths = groups[name]
        if args.max_per_group and len(paths) > args.max_per_group:
            paths = rng.sample(paths, args.max_per_group)
        coords_list, fits, bad = [], [], {}
        for p in paths:
            try:
                c, fit, _ = binder_in_reference_frame(
                    p, ref_core, load_structure, split_chains, args.atoms)
            except Exception as exc:           # one bad file must not stop a run
                bad[str(exc)[:60]] = bad.get(str(exc)[:60], 0) + 1
                continue
            coords_list.append(c)
            fits.append(fit)
        label = "D%d" % n
        if not coords_list:
            print("%-28s %s: NOTHING extracted - %s" % (name, label, bad))
            continue
        stem = "%s_%s_%s_v1" % (args.pid, args.accession, label)
        write_npz(os.path.join(args.out, stem + ".npz"),
                  coords_list, fits, ref_core, ref_ids)
        f = np.asarray(fits)
        print("%-28s %s  %6d binders  %9d atoms  fit RMSD mean %.3f max %.3f%s"
              % (name, label, len(coords_list), sum(len(c) for c in coords_list),
                 f.mean(), f.max(), "  skipped %d" % sum(bad.values()) if bad else ""))
        if bad:
            for why, k in sorted(bad.items(), key=lambda kv: -kv[1]):
                print("      %5d x %s" % (k, why))
        legend.append((label, name, rank_and_colour(name)[1], len(coords_list),
                       # binder_density_volume.py names maps <pid>_<accession>_<label>_density.mrc
                       "%s_%s_%s_density.mrc" % (args.pid, args.accession, label),
                       stem + ".npz"))

    with open(os.path.join(args.out, "legend.csv"), "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["domain_label", "group", "colour", "n_binders",
                    "expected_map", "npz"])
        w.writerows(legend)
    print("wrote %d npz + legend.csv to %s" % (len(legend), args.out))
    return 0


# --------------------------------------------------------------------------
# self-check: synthetic target + binder, no data files needed
# --------------------------------------------------------------------------

def self_check():
    rng = np.random.default_rng(0)

    # kabsch recovers a known rigid motion exactly
    truth_r, _ = np.linalg.qr(rng.normal(size=(3, 3)))
    if np.linalg.det(truth_r) < 0:
        truth_r[:, 0] *= -1
    truth_t = np.array([13.0, -7.0, 2.5])
    x = rng.normal(size=(40, 3)) * 8
    y = apply_rt(truth_r, truth_t, x)
    r, t = kabsch(x, y)
    assert rmsd(apply_rt(r, t, x), y) < 1e-8, "kabsch failed on an exact motion"
    assert np.allclose(r, truth_r, atol=1e-8) and np.allclose(t, truth_t, atol=1e-8)

    # a reflection must NOT be fitted (det +1 enforced)
    z = x * np.array([1.0, 1.0, -1.0])
    r, t = kabsch(x, z)
    assert np.linalg.det(r) > 0, "kabsch returned a reflection"

    # core located by sequence, not by numbering: an 18-mer offset by 1,
    # numbered from 1, must still match a reference numbered from 5
    # `pad` = how many residues precede the core, so that the core's GEOMETRY is
    # identical in both chains and only the numbering differs.  (If the two cores
    # differ by a translation as well, the binder legitimately lands translated
    # too - the transform preserves the binder's position relative to its own
    # target, which is all it promises.  The distance-matrix check below is the
    # invariant that holds either way.)
    class FakeChain(object):
        def __init__(self, seq, first, pad=0):
            ids, names, coord, anames = [], [], [], []
            for i, ch in enumerate(seq):
                j = i - pad
                for an in ("P", "C1'", "N1"):
                    ids.append(first + i)
                    names.append(ch)
                    anames.append(an)
                    coord.append([j * 5.0 + len(an), j * 2.0, float(len(an))])
            self.res_id = np.array(ids)
            self.res_name = np.array(names)
            self.atom_name = np.array(anames)
            self.coord = np.array(coord, dtype=np.float64)

    ref_chain = FakeChain(CORE_SEQ, 5)
    ref_ids = core_residues(ref_chain)
    assert ref_ids == list(range(5, 21)), ref_ids
    design_chain = FakeChain("C" + CORE_SEQ + "C", 1, pad=1)   # the BD1/HD 18-mer
    d_ids = core_residues(design_chain)
    assert d_ids == list(range(2, 18)), d_ids

    # end to end: move the design target+binder by a known motion, check the
    # binder comes back to where it started
    ref_core = core_atoms(ref_chain, ref_ids)
    binder0 = rng.normal(size=(25, 3)) * 10
    moved_core = {k: apply_rt(truth_r, truth_t, v.reshape(1, 3))[0]
                  for k, v in core_atoms(design_chain, d_ids).items()}
    mob, fix, keys = matched(moved_core, ref_core)
    assert len(keys) == 16 * 3, len(keys)
    r, t = kabsch(mob, fix)
    back = apply_rt(r, t, apply_rt(truth_r, truth_t, binder0))
    assert rmsd(back, binder0) < 1e-6, "binder did not return to the ref frame"
    assert rmsd(apply_rt(r, t, mob), fix) < 1e-6, "target fit is not exact"

    # the invariant that matters: the binder keeps its position relative to the
    # target it was designed against, whatever the two frames were
    def dmat(a, b):
        return np.linalg.norm(a[:, None, :] - b[None, :, :], axis=2)
    before = dmat(binder0, np.array([core_atoms(design_chain, d_ids)[k]
                                     for k in sorted(keys)]))
    after = dmat(back, fix)
    assert np.abs(before - after).max() < 1e-6, "binder-to-target geometry moved"

    # the atom filter keeps exactly the backbone, whatever else is in the file
    class FakeProt(object):
        def __init__(self, names, elems):
            self.atom_name = np.array(names)
            self.element = np.array(elems)
            self.coord = np.zeros((len(names), 3))

        def __getitem__(self, m):
            return FakeProt(self.atom_name[m], self.element[m])

        def __len__(self):
            return len(self.atom_name)

    prot = FakeProt(["N", "CA", "C", "O", "CB", "CG", "OXT", "HA"],
                    ["N", "C", "C", "O", "C", "C", "O", "H"])
    assert sorted(pick_atoms(prot, "backbone").atom_name) == ["C", "CA", "N", "O"]
    assert len(pick_atoms(prot, "heavy")) == 7, "heavy must drop H and keep the rest"

    # the .npz round-trips under allow_pickle=False, the way her scripts read it
    import tempfile
    with tempfile.TemporaryDirectory() as td:
        p = os.path.join(td, "1_csy4_D1_v1.npz")
        write_npz(p, [binder0, binder0 + 1.0], [0.1, 0.2], ref_core, ref_ids)
        d = np.load(p, allow_pickle=False)
        assert list(d["offsets"]) == [0, 25, 50], list(d["offsets"])
        assert len(d["coords"]) == 50 and d["coords"].shape[1] == 3
        assert str(d["fit_mode"]) == "global"
        names = [n.decode() for n in d["domain_heavy_name"]]
        assert set(names) == {"P", "C1'", "N1"}, set(names)
        # her frame check: every domain atom sits on the reference, exactly
        idx = {(int(rs), nm): xyz for rs, nm, xyz
               in zip(d["domain_heavy_resseq"], names, d["domain_heavy"])}
        worst = 0.0
        for (pos, nm), xyz in ref_core.items():
            got = idx[(ref_ids[pos], nm)]
            worst = max(worst, float(np.linalg.norm(got - xyz)))
        assert worst == 0.0, "frame deviation against the reference is %g" % worst
        assert sorted(d["domain_resseq"]) == ref_ids

    print("self-check passed")
    return 0


if __name__ == "__main__":
    sys.exit(main())
