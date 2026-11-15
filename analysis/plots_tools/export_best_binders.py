#!/usr/bin/env python3
"""export_best_binders.py — extract the best designs' CIFs UNCOMPRESSED into
analyses/best_binders/, named by ranking (rest of the design name preserved).

Two sets (same consensus ipTM = min(AF3,B2) ranking as plots 01_1 / 01_2):
    best_binders/top<N>_across_all/   — global top-N              (rank01..)
    best_binders/best<N>_per_stack/   — best-N of each stack      (rank01.. per stack)

For each design it writes every available structure — the DESIGNED model plus each refold
(AF3 / Boltz-2 / ESMFold2) — as  rank<NN>_<design_id>.<source>.cif  so all four can be
opened/superposed for viewing. A manifest CSV is written per set. Run standalone or via
run_all.py (which calls export() after the plots).

  python3 export_best_binders.py [--master CSV] [--stacks s1 s2 ...]
"""
from __future__ import annotations
import argparse
import subprocess
from pathlib import Path

import config as cfg
import data_loader as D

# structure sources to extract per design: (filename tag, wide-CSV path column)
SOURCES = [
    ("designed", "designed_cif_path"),
    ("af3", "af3_msaoff_cif_path"),
    ("b2", "b2_cif_path"),
    ("esm", "esm_cif_path"),
]


def _extract(spec, out: Path) -> bool:
    """Resolve a plain path or '<archive.tar.zst>!member' and write it UNCOMPRESSED to out."""
    if not isinstance(spec, str) or not spec.strip():
        return False
    src = spec if spec.startswith("/") else str(cfg.RB / spec)
    if "!" in src:
        arch, member = src.split("!", 1)
        r = subprocess.run(["tar", "--zstd", "-xOf", arch, member.lstrip("/")],
                           capture_output=True)
        if r.returncode != 0 or not r.stdout:
            return False
        out.write_bytes(r.stdout)
        return True
    p = Path(src)
    if p.exists():
        out.write_bytes(p.read_bytes())
        return True
    return False


def _export_set(ranked, rank_field, subdir):
    d = cfg.BEST_BINDERS_DIR / subdir
    d.mkdir(parents=True, exist_ok=True)
    # clear stale rank*.cif so a re-run doesn't leave orphans from a previous ranking
    for old in d.glob("rank*.cif"):
        old.unlink()
    n_designs = 0
    n_files = 0
    for _, x in ranked.iterrows():
        rk = int(x[rank_field])
        did = str(x["design_id"])
        got = 0
        for tag, col in SOURCES:
            out = d / f"rank{rk:02d}_{did}.{tag}.cif"
            if _extract(x.get(col), out):
                got += 1
        n_designs += 1
        n_files += got
    return d, n_designs, n_files


def export(master_csv=None, stacks=None):
    r = D.ranked_wide(master_csv=master_csv, stacks=stacks)
    if r.empty:
        print("[best_binders] no eligible designs (need both AF3 + B2 refolds)")
        return
    cfg.BEST_BINDERS_DIR.mkdir(parents=True, exist_ok=True)
    man_cols = [c for c in ["rank", "rank_in_stack", "stack", "design_id", "consensus_iptm",
                            "mean_iptm", "af3_msaoff_iptm", "b2_iptm", "esm_iptm", "length",
                            "dockq_af3_msa_off", "dockq_b2_msa_off", "dockq_esm_msa_off"]
                if c in r.columns]

    top = r.head(cfg.TOP_N_ACROSS)
    d1, nd1, nf1 = _export_set(top, "rank", "top%d_across_all" % cfg.TOP_N_ACROSS)
    top[man_cols].to_csv(cfg.BEST_BINDERS_DIR / "top_across_manifest.csv", index=False)

    per = r[r["rank_in_stack"] <= cfg.TOP_N_PER_STACK].sort_values(["stack", "rank_in_stack"])
    d2, nd2, nf2 = _export_set(per, "rank_in_stack", "best%d_per_stack" % cfg.TOP_N_PER_STACK)
    per[man_cols].to_csv(cfg.BEST_BINDERS_DIR / "best_per_stack_manifest.csv", index=False)

    print(f"[best_binders] {d1.name}: {nd1} designs / {nf1} cifs   |   "
          f"{d2.name}: {nd2} designs / {nf2} cifs   -> {cfg.BEST_BINDERS_DIR}")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--master", default=None)
    ap.add_argument("--stacks", nargs="*", default=None)
    a = ap.parse_args()
    export(master_csv=a.master, stacks=a.stacks)
