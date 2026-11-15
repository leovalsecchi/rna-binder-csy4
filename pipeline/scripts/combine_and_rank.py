#!/usr/bin/env python3
"""combine_and_rank.py — build the cross-campaign table and top-N export.

Concatenates the v3 (HD + RFD3) and v2 (BG + BD) analyses produced by the SAME
collect->clean->compute pipeline, so all five stacks share one uniform schema.
Ranks by an AF3 & B2 *consensus* iPTM (a design only ranks high if BOTH refolders
like it) and exports the winners.

Consensus score:
    consensus_iptm = min(af3_msaoff_iptm, b2_iptm)      # both-must-be-high
    mean_iptm      = mean(af3_msaoff_iptm, b2_iptm)      # reported alongside
Only designs with BOTH an AF3-msaoff and a B2 refold are eligible.

Outputs under analyses/<out>/:
    combined_clean.csv              — v3+v2 full clean set (ranking population)
    combined_sample10k_metrics.csv   — v3+v2 600/stack sample (plot input)
    top10_across_stacks.csv         — global top 10 by consensus
    topN_per_stack.csv              — top N per stack by consensus
    top_designs/                    — extracted designed + AF3 refold CIFs (top 10 global)
    README.md

Usage:
    combine_and_rank.py --v3 analyses/midrun_v3 --v2 analyses/cross_v2_bg_bd \
        --out analyses/cross_campaign [--per-stack 10] [--top 10]
"""
from __future__ import annotations

import argparse
import shutil
import sys
from pathlib import Path

import pandas as pd

RB = Path("/cluster/work/wicky/lvalsecchi/rna_binder")
sys.path.insert(0, str(RB / "code" / "scripts"))
import archive_read  # noqa: E402  (resolve archive refs -> local paths)

AF3, B2, ESM = "af3_msaoff_iptm", "b2_iptm", "esm_iptm"


def load_many(dirlist, fname: str) -> pd.DataFrame:
    """Concat the same file from every analysis dir, tagging campaign."""
    frames = []
    for camp, d in dirlist:
        f = d / fname
        if not f.exists():
            sys.exit(f"missing {f}")
        df = pd.read_csv(f, low_memory=False)
        df["campaign"] = camp
        frames.append(df)
        print(f"  {camp}: {len(df):>8,} rows  <- {f.name}")
    out = pd.concat(frames, ignore_index=True, sort=False)
    print(f"  combined: {len(out):>8,} rows")
    return out


def add_consensus(df: pd.DataFrame) -> pd.DataFrame:
    for c in (AF3, B2):
        if c not in df.columns:
            sys.exit(f"column {c} absent — cannot rank on AF3&B2 consensus")
    both = df[AF3].notna() & df[B2].notna()
    df = df.copy()
    df["has_consensus"] = both
    df["consensus_iptm"] = df[[AF3, B2]].min(axis=1).where(both)
    df["mean_iptm"] = df[[AF3, B2]].mean(axis=1).where(both)
    return df


MANIFEST_COLS = ["rank", "stack", "campaign", "design_id", "uuid8",
                 "consensus_iptm", "mean_iptm", AF3, B2, ESM, "length",
                 "designed_cif_path", "af3_msaoff_cif_path", "b2_cif_path"]


def manifest(df: pd.DataFrame) -> pd.DataFrame:
    cols = [c for c in MANIFEST_COLS if c in df.columns or c == "rank"]
    m = df.copy()
    m.insert(0, "rank", range(1, len(m) + 1))
    return m[[c for c in cols if c in m.columns]]


def export_cifs(top: pd.DataFrame, dst: Path) -> None:
    """Extract designed + AF3-refold CIFs for the top designs (archive-aware)."""
    if dst.exists():
        shutil.rmtree(dst)      # clear stale top-N CIFs from prior runs before re-extracting
    dst.mkdir(parents=True, exist_ok=True)
    ok = 0
    for _, r in top.iterrows():
        tag = f"rank{int(r['rank']):02d}_{r['stack'].replace('+','_')}_{r['design_id']}"
        for kind, col in (("designed", "designed_cif_path"),
                          ("af3", "af3_msaoff_cif_path")):
            ref = r.get(col)
            if not isinstance(ref, str) or not ref:
                continue
            try:
                local = archive_read.resolve(ref)
                out = dst / f"{tag}.{kind}.cif"
                shutil.copy(local, out)
                ok += 1
            except Exception as e:  # noqa: BLE001  (best-effort export)
                print(f"    (skip {kind} for {r['design_id']}: {e})")
    print(f"  extracted {ok} CIFs -> {dst}")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--dirs", nargs="+",
                    help="campaign inputs as LABEL:PATH (e.g. v2:analyses/cross_v2_bg_bd "
                         "v3:analyses/midrun_v3 v4:analyses/midrun_v4 v5:analyses/midrun_v5)")
    ap.add_argument("--v3", help="(legacy) equivalent to --dirs v3:<path>")
    ap.add_argument("--v2", help="(legacy) equivalent to --dirs v2:<path>")
    ap.add_argument("--out", required=True)
    ap.add_argument("--top", type=int, default=10, help="global top-N across stacks")
    ap.add_argument("--per-stack", type=int, default=10, help="top-N per stack")
    args = ap.parse_args()

    if args.dirs:
        dirlist = []
        for spec in args.dirs:
            lab, p = spec.split(":", 1) if ":" in spec else (Path(spec).name, spec)
            dirlist.append((lab, Path(p)))
    elif args.v3 and args.v2:
        dirlist = [("v3", Path(args.v3)), ("v2", Path(args.v2))]
    else:
        ap.error("provide --dirs LABEL:PATH ...  (or legacy --v3 --v2)")

    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    print("campaigns: " + ", ".join(f"{lab}={d}" for lab, d in dirlist))

    print("=== combine clean (ranking population) ===")
    clean = add_consensus(load_many(dirlist, "collected_clean.csv"))
    clean.to_csv(out / "combined_clean.csv", index=False)

    print("=== combine sample10k_metrics (plot input) ===")
    sample = add_consensus(load_many(dirlist, "sample10k_metrics.csv"))
    sample.to_csv(out / "combined_sample10k_metrics.csv", index=False)

    # combined WHOLE-SET metrics (only if EVERY campaign ran compute_metrics --full)
    if all((d / "collected_clean_metrics.csv").exists() for _, d in dirlist):
        print("=== combine collected_clean_metrics (FULL, whole set) ===")
        full = add_consensus(load_many(dirlist, "collected_clean_metrics.csv"))
        full.to_csv(out / "combined_clean_metrics.csv", index=False)
    else:
        miss = [lab for lab, d in dirlist if not (d / "collected_clean_metrics.csv").exists()]
        print(f"=== (skip combined_clean_metrics FULL: no --full metrics for {miss}) ===")

    elig = clean[clean["has_consensus"]].copy()
    print(f"\n{len(elig):,}/{len(clean):,} designs have BOTH AF3+B2 (rankable)")
    per_stack_counts = elig["stack"].value_counts()
    print("rankable per stack:\n" + per_stack_counts.to_string())

    # global top-N
    g = elig.sort_values("consensus_iptm", ascending=False).head(args.top)
    gm = manifest(g)
    gm.to_csv(out / "top10_across_stacks.csv", index=False)
    print(f"\n=== global top {args.top} (consensus) ===")
    print(gm[["rank", "stack", "design_id", "consensus_iptm", AF3, B2]].to_string(index=False))

    # top-N per stack
    parts = []
    for stk, grp in elig.groupby("stack"):
        parts.append(grp.sort_values("consensus_iptm", ascending=False).head(args.per_stack))
    ps = pd.concat(parts, ignore_index=True).sort_values(
        ["stack", "consensus_iptm"], ascending=[True, False])
    psm = manifest(ps)
    psm.to_csv(out / "topN_per_stack.csv", index=False)
    print(f"\n=== top {args.per_stack} per stack written ({len(psm)} rows) ===")

    print("\n=== extract top-global CIFs ===")
    export_cifs(gm, out / "top_designs")

    (out / "README.md").write_text(
        f"# Cross-campaign sweep\n\n"
        f"Campaigns combined: **{', '.join(lab for lab, _ in dirlist)}**, "
        f"all re-run through the same collect->clean->compute pipeline.\n\n"
        f"Ranking = **AF3 & B2 consensus** = min(af3_msaoff_iptm, b2_iptm); "
        f"only designs with both refolds are eligible.\n\n"
        f"- `combined_clean.csv` — {len(clean):,} clean designs (ranking population)\n"
        f"- `combined_sample10k_metrics.csv` — {len(sample):,} sampled designs (plot input)\n"
        f"- `top10_across_stacks.csv` — global top {args.top}\n"
        f"- `topN_per_stack.csv` — top {args.per_stack} per stack\n"
        f"- `top_designs/` — extracted designed + AF3-refold CIFs of the global top {args.top}\n"
        f"- `plots/` — interactive HTML + contact heatmaps (plot_cross.py)\n")
    print(f"\nDONE -> {out}")


if __name__ == "__main__":
    main()
