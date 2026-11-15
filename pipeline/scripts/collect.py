#!/usr/bin/env python3
"""collect.py — one row per design with its AF3-msaoff / B2 / ESM refold metrics.

Step 1 of the analyses/ pipeline (collect -> clean -> compute_metrics).

Walks the campaign run dirs given as CLI args and reuses the project aggregator
(`rna_binder.analysis.aggregate_metrics`) to discover every refold sample, parse
its canonical filename, and read its confidence JSON sidecar. That aggregator is
archive-aware (.tar.zst), HD-aware (recycle/sample/pack axes) and --num-runs-aware
(per-uuid design slices), so we get the discovery + parsing + JSON extraction for
free instead of re-walking the tree by hand.

The aggregator emits ONE ROW PER REFOLD SAMPLE (long form). This script PIVOTS that
to ONE ROW PER DESIGN (wide form): for each refolder we keep that design's
best-scoring sample (by the refolder's own iPTM) and lay AF3/B2/ESM side by side.

Output columns (in order):
    design_id   stack   uuid8   sequence
    af3_msaoff_iptm  af3_msaoff_ptm  b2_iptm  esm_iptm
    has_af3_msaoff   has_b2   has_esm
    af3_msaoff_cif_path  b2_cif_path  esm_cif_path
    + gen, invfold, length          (handy for per-stack grouping)
    + designed_cif_path, source_run (so compute_metrics.py can do Cα-RMSD-vs-designed)

Notes
-----
* AF3 is restricted to the **msaoff** MSA variant (msaon was retired campaign-wide);
  B2 and ESM have no MSA variants so all their samples are eligible.
* `design_id` is the canonical design stem (e.g.
  `csy4free_RFD3_LiMPNN_len80_b4ed9efe_design0`) — globally unique, and the exact
  basename of the source design file, so it joins straight back to designs/.
* CIF paths are the aggregator's *stable refs*: a plain path for live runs, or
  `<archive>.tar.zst!<member>` for packed runs. compute_metrics.py resolves both
  via archive_read.resolve().

Usage
-----
    python3 collect.py RUN_DIR [RUN_DIR ...] -o collected_raw.csv [--progress]

Example:
    python3 collect.py \\
      campaigns/csy4_crrna_free/rfd3__ligandmpnn/midrun_v2_* \\
      campaigns/csy4_crrna_free/boltzgen__boltzif/midrun_v2_* \\
      -o analyses/midrun_v2_pipeline/collected_raw.csv --progress
"""
from __future__ import annotations

import argparse
import sys
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent          # <root>/code/scripts
RB = SCRIPT_DIR.parents[1]                             # <root>
SRC = RB / "code" / "src"
for _p in (str(SRC), str(SCRIPT_DIR)):
    if _p not in sys.path:
        sys.path.insert(0, _p)

import pandas as pd  # noqa: E402

from rna_binder.analysis import aggregate_metrics  # noqa: E402
from rna_binder.analysis.aggregate import _find_source_design, _stable_ref  # noqa: E402
from rna_binder.analysis.parsing import parse_filename  # noqa: E402

# Per-refolder: (refolder tag, output prefix, the iPTM column to rank samples by).
REFOLDERS = [
    ("AF3", "af3_msaoff", "af3_iptm"),
    ("B2", "b2", "b2_iptm"),
    ("ESM", "esm", "esm_iptm"),
]

OUT_COLS = [
    "design_id", "stack", "uuid8", "sequence",
    "af3_msaoff_iptm", "af3_msaoff_ptm", "b2_iptm", "esm_iptm",
    "has_af3_msaoff", "has_b2", "has_esm",
    "af3_msaoff_cif_path", "b2_cif_path", "esm_cif_path",
    "gen", "invfold", "length",
    "designed_cif_path", "source_run",
]


def design_stem(r) -> str:
    """Canonical design stem (no refolder suffix, no extension) for a parsed row.

    Mirrors the source-design filename so the id joins back to designs/.
    """
    parts = [r["target"], r["gen"], r["invfold"], f"len{int(r['length'])}", r["uuid8"]]
    if pd.notna(r.get("cycle")):                       # HalluDesign axes
        parts += [f"recycle{int(r['cycle'])}",
                  f"sample{int(r['hd_sample'])}",
                  f"pack{int(r['pack'])}"]
    else:                                              # design<D>[_seq<S>]
        parts.append(f"design{int(r['design_id'])}")
        if pd.notna(r.get("seq_id")):
            parts.append(f"seq{int(r['seq_id'])}")
    return "_".join(parts)


def pick_best(sub: pd.DataFrame, iptm_col: str):
    """Best sample row from `sub` by `iptm_col` (max). Falls back to the first row
    when every iPTM is NaN, so a refold with a CIF but no score still counts as
    'available'. Returns None when `sub` is empty."""
    if sub.empty:
        return None
    s = sub[iptm_col] if iptm_col in sub.columns else None
    if s is not None and s.notna().any():
        return sub.loc[s.idxmax()]
    return sub.iloc[0]


def _designed_ref(af3_row) -> str | None:
    """Stable ref to the source design file backing this AF3 sample (for RMSD)."""
    parsed = parse_filename(af3_row["file"])
    if parsed is None:
        return None
    refold_run = Path(af3_row["refold_run"])
    designed = _find_source_design(refold_run, parsed)
    return _stable_ref(refold_run.parent.parent / "designs", designed)


def _leaf_labels(runs, refolders):
    """All distinct refold-leaf names (e.g. 'msaoff_<uuid8>') across the run dirs —
    the unit we shard parallel aggregation over. A uuid's AF3/B2/ESM leaves share a
    name, so chunking by name keeps a design's three refolds in the same worker."""
    labels = set()
    for run in runs:
        for rf in refolders:
            root = Path(run) / f"refold_{rf.lower()}"
            if root.is_dir():
                labels.update(d.name for d in root.iterdir() if d.is_dir())
    return sorted(labels)


def _agg_chunk(task):
    """Worker (own process): aggregate one chunk of leaf labels -> long df."""
    runs, refolders, labels = task
    from rna_binder.analysis import aggregate_metrics as _agg
    return _agg(source_runs=runs, refolders=tuple(refolders), labels=labels,
                structural=False, sequences=True, progress=False)


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("runs", nargs="+", type=Path, help="campaign run directories")
    ap.add_argument("-o", "--output", type=Path, default=Path("collected_raw.csv"))
    ap.add_argument("--refolders", nargs="+", default=["AF3", "B2", "ESM"],
                    choices=["AF3", "B2", "ESM"])
    ap.add_argument("--progress", action="store_true",
                    help="print one line per refold file the aggregator reads")
    ap.add_argument("--jobs", type=int, default=1,
                    help="parallel worker processes (default 1). Shards refold leaves "
                         "across processes; the written CSV is identical to --jobs 1.")
    args = ap.parse_args()

    runs = [p for p in args.runs if p.exists()]
    if not runs:
        sys.exit(f"none of the {len(args.runs)} run dirs exist")

    if args.jobs > 1:
        labels = _leaf_labels(runs, args.refolders)
        nchunks = min(len(labels), max(args.jobs * 4, args.jobs)) if labels else 0
        chunks = [c for c in ([labels[i::nchunks] for i in range(nchunks)] if nchunks else []) if c]
        print(f"[collect] {len(labels)} leaves -> {len(chunks)} chunks x {args.jobs} workers",
              flush=True)
        tasks = [(runs, args.refolders, c) for c in chunks]
        parts, done = [], 0
        with ProcessPoolExecutor(max_workers=args.jobs) as ex:
            for df_part in ex.map(_agg_chunk, tasks):
                parts.append(df_part); done += 1
                print(f"  chunk {done}/{len(chunks)} done ({len(df_part)} rows)", flush=True)
        long = pd.concat(parts, ignore_index=True) if parts else pd.DataFrame()
    else:
        long = aggregate_metrics(
            source_runs=runs,
            refolders=tuple(args.refolders),
            structural=False,
            sequences=True,
            progress=args.progress,
        )
    if long.empty:
        sys.exit(f"no refold samples found across {len(runs)} run dir(s)")

    # AF3 -> msaoff only (msaon retired). B2/ESM have no MSA variant.
    if "msa_variant" in long.columns:
        af3_bad = long["refolder"].eq("AF3") & long["msa_variant"].ne("msaoff")
        long = long[~af3_bad].copy()

    long["design_id"] = long.apply(design_stem, axis=1)

    records = []
    for stem, g in long.groupby("design_id", sort=False):
        base = g.iloc[0]
        seqs = g["binder_sequence"].dropna() if "binder_sequence" in g else pd.Series([], dtype=object)
        rec = {
            "design_id": stem,
            "stack": base["stack"],
            "uuid8": base["uuid8"],
            "sequence": seqs.iloc[0] if not seqs.empty else None,
            "gen": base["gen"],
            "invfold": base["invfold"],
            "length": int(base["length"]),
            "designed_cif_path": None,
            "source_run": base.get("source_run"),
        }
        for tag, pfx, iptm_col in REFOLDERS:
            best = pick_best(g[g["refolder"].eq(tag)], iptm_col)
            rec[f"has_{pfx}"] = best is not None
            rec[f"{pfx}_cif_path"] = None if best is None else best.get("abs_path")
            rec[f"{pfx}_iptm"] = None if best is None else best.get(iptm_col)
            if tag == "AF3":
                rec["af3_msaoff_ptm"] = None if best is None else best.get("af3_ptm")
                if best is not None:
                    rec["designed_cif_path"] = _designed_ref(best)
        records.append(rec)

    df = pd.DataFrame.from_records(records)
    df = df.reindex(columns=OUT_COLS)
    df = df.sort_values(["stack", "uuid8", "design_id"]).reset_index(drop=True)

    args.output.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(args.output, index=False)

    n_all3 = int((df["has_af3_msaoff"] & df["has_b2"] & df["has_esm"]).sum())
    print(f"Wrote {len(df)} designs x {len(df.columns)} cols -> {args.output}")
    print(f"  per stack: {df.groupby('stack').size().to_dict()}")
    print(f"  with AF3-msaoff: {int(df['has_af3_msaoff'].sum())} | "
          f"B2: {int(df['has_b2'].sum())} | ESM: {int(df['has_esm'].sum())} | "
          f"all three: {n_all3}")


if __name__ == "__main__":
    main()
