#!/usr/bin/env python3
"""Validate the full run. Every check either passes or explains itself.

Usage:  python scripts/validate.py --results results --build build
"""
from __future__ import annotations

import argparse
import json
import subprocess
from pathlib import Path

import numpy as np
import pandas as pd

FM = Path("/cluster/work/wicky/lvalsecchi/rna_binder/analyses/foldseek_map")
RB = Path("/cluster/work/wicky/lvalsecchi/rna_binder")
LAYERS = ("design", "af3", "esm")

# Campaign tables that carry rna_contact_nt_1..16 (checked 2026-08-21).
CONTACT_TABLES = {
    "v3": RB / "analyses/midrun_v3/collected_clean_metrics.csv",
    "v4": RB / "analyses/midrun_v4/collected_clean_metrics.csv",
    "v5": RB / "analyses/midrun_v5/collected_clean_metrics.csv",
    "v9": RB / "analyses/midrun_v9/collected_clean_metrics.csv",
}


def check(name: str, ok: bool, detail: str = "") -> bool:
    print(f"  [{'PASS' if ok else 'FAIL'}] {name}" + (f" — {detail}" if detail else ""))
    return ok


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--results", type=Path, default=FM / "results")
    parser.add_argument("--build", type=Path, default=FM / "build")
    parser.add_argument("--skip-contact-crosscheck", action="store_true")
    args = parser.parse_args()

    ok = True
    manifest: dict = {"layers": {}, "checks": {}}

    manifest_counts = json.loads((args.build / "meta" / "counts.json").read_text())
    print("=== per-layer reconciliation ===")
    for layer in LAYERS:
        out = args.results / layer
        if not (out / "dedup" / "nodes.tsv").exists():
            ok = check(f"{layer}: results present", False, "missing dedup/nodes.tsv") and ok
            continue
        nodes = pd.read_csv(out / "dedup" / "nodes.tsv", sep="\t")
        hashes = pd.concat([pd.read_csv(f, sep="\t")
                            for f in sorted((out / "hashes").glob("shard_*.tsv"))],
                           ignore_index=True)
        missing = pd.concat([pd.read_csv(f, sep="\t")
                             for f in sorted((out / "missing").glob("shard_*.tsv"))],
                            ignore_index=True)
        expected = manifest_counts["layer_counts"][layer]
        seen = len(hashes) + len(missing)
        ok = check(f"{layer}: every design accounted for",
                   seen == expected, f"{seen:,d} of {expected:,d}") and ok

        # Chain tally: RFD3 design structures MUST be inverted.
        if layer == "design":
            man = pd.read_csv(args.build / "manifests" / "design.tsv", sep="\t")
            joined = hashes.merge(man[["design_id", "stack"]], on="design_id")
            tally = joined.groupby(["stack", "protein_chain"]).size().unstack(fill_value=0)
            print(tally.to_string())
            rfd3 = tally.loc["RFD3+LiMPNN"] if "RFD3+LiMPNN" in tally.index else None
            ok = check("design: RFD3+LiMPNN protein is chain B",
                       rfd3 is not None and rfd3.get("B", 0) > 0
                       and rfd3.get("A", 0) == 0,
                       "if this flips, chain selection is broken") and ok
        else:
            ok = check(f"{layer}: protein always chain A",
                       set(hashes["protein_chain"]) == {"A"},
                       str(hashes['protein_chain'].value_counts().to_dict())) and ok

        ok = check(f"{layer}: protein always 80 residues",
                   set(hashes["n_res"]) == {80}) and ok

        manifest["layers"][layer] = {
            "structures": int(len(hashes)),
            "unique_backbones": int(len(nodes)),
            "missing": int(len(missing)),
            "missing_reasons": missing["reason"].value_counts().to_dict()
                               if len(missing) else {},
        }

    # Self-hits must score 1.0
    print("\n=== fold sanity ===")
    for layer in LAYERS:
        files = sorted((args.results / layer / "hits").glob("hits_*.tsv"))
        if not files:
            continue
        head = pd.read_csv(files[0], sep="\t", header=None, nrows=200000,
                           names=["query", "target", "qtm", "ttm", "alntm",
                                  "lddt", "fident", "alnlen", "evalue", "bits"])
        self_hits = head[head["query"] == head["target"]]
        if len(self_hits):
            ok = check(f"{layer}: self-hits score 1.0",
                       bool((self_hits["qtm"] > 0.99).all()),
                       f"min {self_hits['qtm'].min():.3f}") and ok

    # Epitope cross-check against the stored rna_contact_nt_* columns
    if not args.skip_contact_crosscheck:
        print("\n=== epitope cross-check vs stored rna_contact_nt_* ===")
        per_nt = [f"per_nt_{i}" for i in range(1, 17)]
        stored_cols = [f"rna_contact_nt_{i}" for i in range(1, 17)]
        for campaign, path in CONTACT_TABLES.items():
            if not path.exists():
                print(f"  [skip] {campaign}: {path} not found")
                continue
            stored = pd.read_csv(path, usecols=["design_id"] + stored_cols,
                                 low_memory=False).set_index("design_id")
            hit_any = False
            for layer in LAYERS:
                fp_path = args.results / layer / "fp" / "canonical.parquet"
                if not fp_path.exists():
                    continue
                fp = pd.read_parquet(fp_path).set_index("design_id")
                common = fp.index.intersection(stored.index)
                if len(common) < 20:
                    print(f"  [skip] {campaign}/{layer}: only {len(common)} shared "
                          f"design_ids (dedup elects one id per backbone, so the "
                          f"design layer overlaps less than the refold layers)")
                    continue
                hit_any = True
                a = fp.loc[common, per_nt].to_numpy(dtype=float)
                b = stored.loc[common, stored_cols].to_numpy(dtype=float)
                good = np.isfinite(b).all(axis=1)
                if good.sum() < 20:
                    continue
                corr = np.corrcoef(a[good].ravel(), b[good].ravel())[0, 1]
                print(f"  {campaign}/{layer}: n={int(good.sum()):,d} "
                      f"r={corr:.3f}")
                ok = check(f"{campaign}/{layer}: recomputed contacts track stored",
                           corr > 0.9,
                           "a low correlation means the cutoff, frame, or chain "
                           "selection diverged from rna_binder") and ok
            if not hit_any:
                print(f"  [skip] {campaign}: no overlapping design_ids")

    manifest["checks"]["all_passed"] = ok
    try:
        version = subprocess.run(
            ["bash", "-lc", "module load stack/2024-06 gcc/12.2.0 openmpi/4.1.6 "
             "foldseek/8-ef4e960 >/dev/null 2>&1; foldseek version"],
            capture_output=True, text=True).stdout.strip()
    except Exception:
        version = "unknown"
    manifest["tools"] = {"foldseek": version or "8-ef4e960"}
    manifest["input"] = manifest_counts

    (args.results / "manifest.json").write_text(json.dumps(manifest, indent=2))
    print(f"\nmanifest -> {args.results / 'manifest.json'}")
    print("ALL CHECKS PASSED" if ok else "SOME CHECKS FAILED")
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
