### Task 12: RNA-variant census (the frame gate)

Runs on Euler before any compute. A misregistered frame would invalidate every epitope result while looking entirely normal, so this fails loudly rather than proceeding.

**Files:**
- Create: `scripts/stage0b_rna_census.py`
- Create: `slurm/stage_single.sbatch`

**Interfaces:**
- Consumes: `binder_map.chains.load_structure`, `split_chains`; `binder_map.rnaframe.rna_sequence`, `frame_offset`, `FrameError`.
- Produces: `$FM/build/meta/rna_variants.tsv` — columns `layer, stack, campaign, rna_sequence, length, offset, n_sampled, status`. Exit code 1 if any variant fails to align.

- [ ] **Step 1: Write the parameterised single-node sbatch**

Create `slurm/stage_single.sbatch`. Stages 2, 4, 5, 6, 6b, 7 and this census all share it.

```bash
#!/bin/bash
# Parameterised single-node stage runner.
#   sbatch --job-name=bm_s4_esm slurm/stage_single.sbatch <script> [args...]
#SBATCH --partition=normal.4h
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu=4G
#SBATCH --output=/cluster/work/wicky/lvalsecchi/rna_binder/analyses/foldseek_map/logs/%x_%j.log
set -uo pipefail            # NOT -e: report failures rather than vanishing

FM=/cluster/work/wicky/lvalsecchi/rna_binder/analyses/foldseek_map
PY=/cluster/work/wicky/lvalsecchi/binder_map_venv/bin/python
export PYTHONPATH="$FM/src:${PYTHONPATH:-}"

mkdir -p "$FM/logs"
echo "=== $SLURM_JOB_NAME on $(hostname) at $(date) ==="
echo "cmd: $PY $*"

cd "$FM" || exit 1
"$PY" "$@"
rc=$?
echo "=== exit $rc at $(date) ==="
exit $rc
```

- [ ] **Step 2: Implement `scripts/stage0b_rna_census.py`**

```python
#!/usr/bin/env python3
"""Stage 0b: census of RNA length/sequence variants. A gate, not a report.

Samples structures per (layer, stack, campaign), reads the RNA chain, and
checks every distinct sequence aligns to the canonical csy4 16-mer. Exits 1
on any failure so the chained SLURM jobs stop before spending compute.

Usage (on Euler, from $FM):
    python scripts/stage0b_rna_census.py --build build --per-group 5
"""
from __future__ import annotations

import argparse
import io
import subprocess
import sys
import tarfile
from collections import defaultdict
from pathlib import Path

import pandas as pd

from binder_map.chains import ChainSelectError, load_structure, split_chains
from binder_map.rnaframe import FrameError, frame_offset, rna_sequence

RB = Path("/cluster/work/wicky/lvalsecchi/rna_binder")


def extract_member(archive: str, member: str) -> bytes:
    """Pull one member out of a .tar.zst via the tar CLI.

    Python's tarfile has no zstd codec in this Python, so the CLI is the
    supported route.
    """
    proc = subprocess.run(
        ["tar", "--zstd", "-xOf", str(RB / archive), member],
        capture_output=True,
    )
    if proc.returncode != 0 or not proc.stdout.strip():
        raise FileNotFoundError(f"{archive}!{member}")
    return proc.stdout


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--build", type=Path, default=Path("build"))
    parser.add_argument("--per-group", type=int, default=5)
    args = parser.parse_args()

    rows = []
    failures = 0

    for layer in ("design", "af3", "esm"):
        manifest = pd.read_csv(args.build / "manifests" / f"{layer}.tsv", sep="\t")
        groups = manifest.groupby(["stack", "campaign"], sort=True)
        for (stack, campaign), group in groups:
            sample = group.head(args.per_group)
            seen: dict[str, int] = defaultdict(int)
            for _, row in sample.iterrows():
                try:
                    data = extract_member(row["archive_path"], row["member_path"])
                    tmp = Path("/tmp") / "bm_census.cif"
                    tmp.write_bytes(data)
                    _, rna = split_chains(load_structure(tmp))
                    seq = rna_sequence(rna)
                    seen[seq] += 1
                except (FileNotFoundError, ChainSelectError) as exc:
                    rows.append({
                        "layer": layer, "stack": stack, "campaign": campaign,
                        "rna_sequence": "", "length": 0, "offset": -1,
                        "n_sampled": 1, "status": f"unreadable: {type(exc).__name__}",
                    })
                    continue

            for seq, count in seen.items():
                try:
                    offset = frame_offset(seq)
                    status = "ok"
                except FrameError:
                    offset = -1
                    status = "FRAME_MISMATCH"
                    failures += 1
                rows.append({
                    "layer": layer, "stack": stack, "campaign": campaign,
                    "rna_sequence": seq, "length": len(seq), "offset": offset,
                    "n_sampled": count, "status": status,
                })

    out = args.build / "meta" / "rna_variants.tsv"
    out.parent.mkdir(parents=True, exist_ok=True)
    frame = pd.DataFrame(rows)
    frame.to_csv(out, sep="\t", index=False)

    print(f"\n=== RNA variant census -> {out} ===")
    summary = frame.groupby(["rna_sequence", "length", "status"], sort=True).size()
    print(summary.to_string())

    if failures:
        print(f"\nFATAL: {failures} RNA variant(s) do not contain the canonical "
              f"16-mer. Epitope comparison would be misregistered. Stopping.",
              file=sys.stderr)
        return 1
    print("\nall RNA variants align to the canonical frame")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
```

- [ ] **Step 3: Sync the repo to Euler**

```bash
cd /Users/leo/Documents/binder_map
ssh euler 'mkdir -p /cluster/work/wicky/lvalsecchi/rna_binder/analyses/foldseek_map/logs'
rsync -av --exclude '.git' --exclude 'tests/fixtures' --exclude 'build' \
  ./ euler:/cluster/work/wicky/lvalsecchi/rna_binder/analyses/foldseek_map/
rsync -av build/ euler:/cluster/work/wicky/lvalsecchi/rna_binder/analyses/foldseek_map/build/
```

- [ ] **Step 4: Run the census**

```bash
ssh euler 'cd /cluster/work/wicky/lvalsecchi/rna_binder/analyses/foldseek_map && \
  PYTHONPATH=src /cluster/work/wicky/lvalsecchi/binder_map_venv/bin/python \
  scripts/stage0b_rna_census.py --build build --per-group 3'
```

Expected: a table dominated by `ACUGCCGUAUAGGCAG` (16 nt, offset 0) plus the `CACUGCCGUAUAGGCAGC` variant (18 nt, offset 1), final line `all RNA variants align to the canonical frame`, exit 0.

If any row reports `FRAME_MISMATCH`, stop. Read the offending sequence and decide whether it is a genuine second target or a parsing artefact before touching anything else.

- [ ] **Step 5: Commit**

```bash
cd /Users/leo/Documents/binder_map
git add scripts/stage0b_rna_census.py slurm/stage_single.sbatch
git commit -m "Add RNA-variant census as a hard gate before compute"
```

---

