#!/usr/bin/env python3
"""Annotate a structural-metrics CSV with DockQ (refold vs DESIGNED = self-consistency).

Runs on the Euler HOST python3 (3.12) with DockQ+gemmi vendored in code/pylibs_dockq
(boltzgen.sif is py3.11 and DockQ ships no cp311 wheel; host py3.12 does). For each row
it computes DockQ between the refold structure (--path-col, default abs_path) and its
designed reference (--ref-col, default designed_source), both archive-aware (.tar.zst!member,
extracted via the bundled zstd). Adds columns: dockq_self, dockq_LRMSD, dockq_iRMSD, dockq_fnat.
Mirrors hbonds_pyrosetta.py (one compute per unique structure, atomic CSV rewrite). Stdlib +
DockQ/gemmi only.

Usage (Euler host, NOT a container):
    PYTHONPATH=code/pylibs_dockq:code/scratch_dockq PATH=code/bin:$PATH \\
        python3 code/scripts/dockq_self.py --csv <structural.csv> [--limit N] [--jobs 8] --progress
"""
from __future__ import annotations
import argparse, csv, os, shlex, shutil, subprocess, sys, tempfile
from multiprocessing import Pool

RB = "/cluster/work/wicky/lvalsecchi/rna_binder"
sys.path.insert(0, f"{RB}/code/pylibs_dockq")
sys.path.insert(0, f"{RB}/code/scratch_dockq")   # ligand_rmsd._kind
COLS = ["dockq_self", "dockq_LRMSD", "dockq_iRMSD", "dockq_fnat"]

# {ref('archive!member' or loose path) -> local file}, built once by _batch_extract
# before the Pool forks so every worker inherits the resolved paths.
_LOCAL = {}


def _batch_extract(items, scratch):
    """Decompress each source .tar.zst ONCE and extract all needed members in a single
    tar pass, vs the old per-design full-archive decompress (a leaf's ~640 designs each
    re-decompressed the same ~100MB archive). Returns {'archive!member' -> local path}."""
    by_arch = {}
    for path, ref in items:
        for x in (path, ref):
            if x and "!" in x:
                arch, member = x.split("!", 1)
                by_arch.setdefault(arch, set()).add(member)
    local = {}
    for ai, (arch, members) in enumerate(by_arch.items()):
        d = os.path.join(scratch, f"a{ai:06d}")
        os.makedirs(d, exist_ok=True)
        mlist = " ".join(shlex.quote(m) for m in members)
        subprocess.run(f"zstd -dc {shlex.quote(arch)} | tar -x -C {shlex.quote(d)} {mlist}",
                       shell=True, stderr=subprocess.DEVNULL)
        for m in members:
            p = os.path.join(d, m)
            if os.path.exists(p) and os.path.getsize(p) > 0:
                local[f"{arch}!{m}"] = p
    return local


def _resolve(ref):
    """Map a loose path or 'archive!member' to a local file via the batch cache."""
    if not ref:
        return None
    if "!" not in ref:
        return ref if os.path.exists(ref) else None
    return _LOCAL.get(ref)


def _chain_types(path):
    import gemmi
    from ligand_rmsd import _kind
    st = gemmi.read_structure(str(path)); out = {}
    for c in st[0]:
        kinds = {_kind(r.name) for r in c} - {"other"}
        if "protein" in kinds: out[c.name] = "protein"
        elif "na" in kinds: out[c.name] = "na"
    return out


def _chain_map(model_path, native_path):
    mt = _chain_types(model_path); n_by_type = {t: ch for ch, t in _chain_types(native_path).items()}
    return {mc: n_by_type[t] for mc, t in mt.items() if t in n_by_type}



def _load_pdb_safe(path):
    """load_PDB, retrying on a gemmi-normalized copy for minimal mmCIFs (e.g. ESMFold2
    refolds) that omit _atom_site.occupancy and crash DockQ's Biopython mmCIF parser.
    AF3/B2 cifs parse first try, so they never hit the gemmi path."""
    from DockQ.DockQ import load_PDB
    try:
        return load_PDB(path)
    except Exception:
        import gemmi
        st = gemmi.read_structure(path); st.setup_entities()
        norm = path + ".norm.cif"
        st.make_mmcif_document().write_file(norm)
        return load_PDB(norm)


def _dockq_one(args):
    path, ref = args
    from DockQ.DockQ import run_on_all_native_interfaces
    m = _resolve(path); n = _resolve(ref)   # pre-extracted by _batch_extract
    if not m or not n:
        return path, None
    try:
        res, _ = run_on_all_native_interfaces(_load_pdb_safe(m), _load_pdb_safe(n), chain_map=_chain_map(m, n))
        if not res:
            return path, None
        v = next(iter(res.values()))
        return path, (round(v["DockQ"], 4), round(v["LRMSD"], 3), round(v["iRMSD"], 3), round(v["fnat"], 4))
    except Exception:
        return path, None


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--csv", required=True)
    ap.add_argument("--output", default=None, help="write here instead of overwriting --csv")
    ap.add_argument("--path-col", default="abs_path")
    ap.add_argument("--ref-col", default="designed_source")
    ap.add_argument("--limit", type=int, default=0, help="process only first N unique structures")
    ap.add_argument("--jobs", type=int, default=min(8, os.cpu_count() or 1))
    ap.add_argument("--progress", action="store_true")
    args = ap.parse_args()

    with open(args.csv, newline="") as f:
        reader = csv.DictReader(f); rows = list(reader); fieldnames = list(reader.fieldnames or [])
    for c in (args.path_col, args.ref_col):
        if c not in fieldnames:
            sys.exit(f"column {c!r} not in CSV")
    for c in COLS:
        if c not in fieldnames:
            fieldnames.append(c)
    # one compute per unique (path, ref)
    pairs = {}
    for r in rows:
        p = r.get(args.path_col, "")
        if p and p not in pairs:
            pairs[p] = r.get(args.ref_col, "")
    items = list(pairs.items())
    if args.limit:
        items = items[:args.limit]
    print(f"{len(rows)} rows, {len(items)} unique structures, jobs={args.jobs}", file=sys.stderr)

    # batch-extract every source archive ONCE (before Pool forks so workers inherit _LOCAL)
    scratch = tempfile.mkdtemp(prefix="dockq_self_batch_")
    global _LOCAL
    n_arch = len({x.split("!", 1)[0] for pr in items for x in pr if x and "!" in x})
    if args.progress:
        print(f"  batch-extracting {n_arch} archives (one decompress each)...", file=sys.stderr)
    _LOCAL = _batch_extract(items, scratch)
    if args.progress:
        print(f"  extracted {len(_LOCAL)} members from {n_arch} archives", file=sys.stderr)

    results = {}; nerr = 0
    try:
        with Pool(args.jobs) as pool:
            for i, (path, val) in enumerate(pool.imap_unordered(_dockq_one, items, chunksize=8)):
                results[path] = val; nerr += int(val is None)
                if args.progress and i % 200 == 0:
                    print(f"  {i}/{len(items)}", file=sys.stderr)
    finally:
        shutil.rmtree(scratch, ignore_errors=True)

    nfill = 0
    for r in rows:
        v = results.get(r.get(args.path_col, ""))
        if v:
            r["dockq_self"], r["dockq_LRMSD"], r["dockq_iRMSD"], r["dockq_fnat"] = v; nfill += 1
        else:
            for c in COLS:
                r.setdefault(c, "")
    out = args.output or args.csv; tmp = out + ".tmp"
    with open(tmp, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames); w.writeheader(); w.writerows(rows)
    os.replace(tmp, out)
    print(f"DockQ filled for {nfill}/{len(rows)} rows ({nerr} structures failed) -> {out}", file=sys.stderr)


if __name__ == "__main__":
    main()
