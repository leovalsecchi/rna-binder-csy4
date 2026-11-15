#!/usr/bin/env python3
"""scaffold_fit.py — does a designed binder fit the csy4 hairpin *within the full 121-nt
scaffold*? Pure geometry (biotite Kabsch + KDTree), no folding.

For each design: take its binder + 16-nt epitope from the AF3 refold, superpose that epitope
onto a scaffold hairpin (the scaffold has TWO), move the binder with it, then measure:
  (1) interface ON the hairpin  -> scaffold_interface_frac (contacts on hairpin / all RNA contacts)
  (2) NO clashes with the rest  -> scaffold_clashes (binder heavy atoms <CLASH_D of non-hairpin RNA)
Keeps the better-fitting hairpin. Emits a sidecar CSV.

  python3 scaffold_fit.py --input <csv w/ design_id,af3_msaoff_cif_path> --out <sidecar.csv> [--jobs 16] [--limit N]
"""
from __future__ import annotations
import argparse, csv, os, shlex, subprocess, sys, tempfile, warnings
warnings.filterwarnings("ignore")  # biotite auth_atom_id fallback warnings, one per file
from functools import partial
from multiprocessing import Pool
from pathlib import Path

import numpy as np

RB = "/cluster/work/wicky/lvalsecchi/rna_binder"
sys.path.insert(0, f"{RB}/code/src")
from rna_binder.analysis.structural import load_structure, split_chains  # noqa: E402
import biotite.structure as struc  # noqa: E402
from scipy.spatial import cKDTree  # noqa: E402

SCAFFOLD_CIF = f"{RB}/analyses/scaffold_fit/scaffold_csy4_dualrepeat.BEST_Boltz2.cif"
CORE = "CUGCCGUAUAGGCAG"    # the conserved hairpin motif shared by epitope + scaffold (res 29 & 77)
CLASH_D = 2.5              # heavy-atom overlap distance = a steric clash
CONTACT_D = 4.5           # heavy-atom interface contact distance
# scaffold_fits = BOTH of your criteria: interface on the hairpin AND no clashes with the rest
FIT_MAX_CLASHES = 3        # tolerate a couple of borderline atom-pairs, fail real overlaps
MIN_HAIRPIN_CONTACTS = 10  # the binder must actually engage the hairpin (else it "fits" by touching nothing)
INTERFACE_FRAC_MIN = 0.80  # its RNA interface must be predominantly on the hairpin, not spilling onto the rest


def _extract(ref, scratch):
    """Plain path or '<archive.tar.zst>!member' -> local file. None on miss."""
    if not isinstance(ref, str) or not ref.strip():
        return None
    if "!" not in ref:
        p = ref if ref.startswith("/") else f"{RB}/{ref}"
        return p if os.path.exists(p) else None
    arch, member = ref.split("!", 1)
    if not arch.startswith("/"):
        arch = f"{RB}/{arch}"
    out = os.path.join(scratch, os.path.basename(member))
    subprocess.run(f"tar --zstd -xOf {shlex.quote(arch)} {shlex.quote(member.lstrip('/'))} "
                   f"> {shlex.quote(out)}", shell=True, stderr=subprocess.DEVNULL)
    return out if os.path.exists(out) and os.path.getsize(out) > 0 else None


def _res_order(rna):
    """(ordered list of res_id, sequence string) for an RNA chain."""
    ids, seq, seen = [], [], set()
    for i in range(rna.array_length()):
        r = int(rna.res_id[i])
        if r not in seen:
            seen.add(r); ids.append(r); seq.append(str(rna.res_name[i]).strip())
    return ids, "".join(seq)


def _c1_array(rna, res_ids):
    """C1' atoms of the given residues, ordered by res_id (matches motif order)."""
    an = np.char.replace(rna.atom_name.astype(str), "*", "'")
    mask = np.isin(rna.res_id, res_ids) & (an == "C1'")
    sub = rna[mask]
    return sub[np.argsort(sub.res_id)]


def _heavy(atoms):
    return atoms[atoms.element != "H"]


# ---- scaffold: load once, find the two hairpins, cache non-hairpin RNA heavy coords ----
def scaffold_context():
    scaf = load_structure(Path(SCAFFOLD_CIF))
    _, rna = split_chains(scaf)
    if rna is None:
        rna = scaf  # the scaffold is RNA-only
    ids, seq = _res_order(rna)
    starts = [i for i in range(len(seq)) if seq[i:i + len(CORE)] == CORE]
    hairpins = [[ids[k] for k in range(s, s + len(CORE))] for s in starts]  # res_ids per hairpin
    return {"rna": rna, "ids": ids, "seq": seq, "hairpins": hairpins}


def _fit_one_hairpin(binder_h, epi_c1_by_core, scaf, hp_res):
    """Superpose the epitope's core C1' onto this hairpin's core C1', move the binder, score."""
    rna = scaf["rna"]
    scaf_c1 = _c1_array(rna, hp_res)                       # 15 C1' (fixed)
    if scaf_c1.array_length() != epi_c1_by_core.array_length() or scaf_c1.array_length() < 6:
        return None
    fitted, transform = struc.superimpose(scaf_c1, epi_c1_by_core)   # (fixed, mobile)
    rmsd = float(struc.rmsd(scaf_c1, fitted))
    binder_moved = transform.apply(binder_h)
    bcoord = binder_moved.coord
    # scaffold RNA heavy atoms split into docked-hairpin vs the rest
    hp_set = set(int(r) for r in hp_res)
    rna_h = _heavy(rna)
    is_hp = np.isin(rna_h.res_id, list(hp_set))
    hp_coord, rest_coord = rna_h.coord[is_hp], rna_h.coord[~is_hp]
    btree = cKDTree(bcoord)
    n_hp = sum(len(x) for x in btree.query_ball_point(hp_coord, CONTACT_D)) if len(hp_coord) else 0
    n_rest = sum(len(x) for x in btree.query_ball_point(rest_coord, CONTACT_D)) if len(rest_coord) else 0
    clashes = sum(len(x) for x in btree.query_ball_point(rest_coord, CLASH_D)) if len(rest_coord) else 0
    ifrac = (n_hp / (n_hp + n_rest)) if (n_hp + n_rest) else 0.0
    return {"clashes": int(clashes), "interface_frac": round(ifrac, 4),
            "n_hairpin_contacts": int(n_hp), "align_rmsd": round(rmsd, 3)}


def score_design(args, scaf):
    design_id, ref = args
    scratch = tempfile.mkdtemp(prefix="scaffold_fit_")
    try:
        p = _extract(ref, scratch)
        if not p:
            return design_id, None
        atoms = load_structure(Path(p))
        prot, epi = split_chains(atoms)
        if prot is None or epi is None:
            return design_id, None
        eids, eseq = _res_order(epi)
        es = eseq.find(CORE)
        if es < 0:
            return design_id, None
        epi_core_res = [eids[k] for k in range(es, es + len(CORE))]
        epi_c1 = _c1_array(epi, epi_core_res)
        binder_h = _heavy(prot)
        scored = []
        for hi, hp_res in enumerate(scaf["hairpins"], 1):
            r = _fit_one_hairpin(binder_h, epi_c1, scaf, hp_res)
            if r:
                r["hairpin"] = hi; scored.append(r)
        if not scored:
            return design_id, None
        # best = a hairpin the binder actually engages, then fewest clashes (binding either is valid)
        best = min(scored, key=lambda r: (r["n_hairpin_contacts"] < MIN_HAIRPIN_CONTACTS, r["clashes"]))
        return design_id, best
    except Exception as e:  # noqa: BLE001
        return design_id, {"error": f"{type(e).__name__}: {e}"}
    finally:
        import shutil; shutil.rmtree(scratch, ignore_errors=True)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--jobs", type=int, default=16)
    ap.add_argument("--limit", type=int, default=0)
    ap.add_argument("--progress", action="store_true")
    a = ap.parse_args()

    scaf = scaffold_context()
    print(f"scaffold: {len(scaf['ids'])} nt, {len(scaf['hairpins'])} hairpins at core-starts "
          f"{[h[0] for h in scaf['hairpins']]}", file=sys.stderr)

    with open(a.input, newline="") as f:
        rows = list(csv.DictReader(f))
    pairs = [(r["design_id"], r.get("af3_msaoff_cif_path", "")) for r in rows if r.get("design_id")]
    pairs = list(dict(pairs).items())          # unique design_id
    if a.limit:
        pairs = pairs[:a.limit]
    print(f"{len(pairs)} designs, jobs={a.jobs}", file=sys.stderr)

    results = {}
    worker = partial(score_design, scaf=scaf)
    with Pool(a.jobs) as pool:
        for i, (did, res) in enumerate(pool.imap_unordered(worker, pairs, chunksize=16)):
            results[did] = res
            if a.progress and i % 500 == 0:
                print(f"  {i}/{len(pairs)}", file=sys.stderr)

    cols = ["design_id", "scaffold_clashes", "scaffold_fits", "scaffold_interface_frac",
            "scaffold_best_hairpin", "scaffold_align_rmsd", "scaffold_n_hairpin_contacts"]
    n_ok = n_fit = 0
    with open(a.out, "w", newline="") as f:
        w = csv.writer(f); w.writerow(cols)
        for did, r in results.items():
            if not r or "error" in (r or {}):
                w.writerow([did] + [""] * (len(cols) - 1)); continue
            fits = int(r["clashes"] <= FIT_MAX_CLASHES                 # (2) no clashes with the rest
                       and r["interface_frac"] >= INTERFACE_FRAC_MIN   # (1) interface on the hairpin...
                       and r["n_hairpin_contacts"] >= MIN_HAIRPIN_CONTACTS)  # ...and actually engaging it
            n_ok += 1; n_fit += fits
            w.writerow([did, r["clashes"], fits, r["interface_frac"], r["hairpin"],
                        r["align_rmsd"], r["n_hairpin_contacts"]])
    print(f"wrote {a.out}: {n_ok} scored, {n_fit} fit cleanly (<= {FIT_MAX_CLASHES} clashes)",
          file=sys.stderr)


if __name__ == "__main__":
    main()
