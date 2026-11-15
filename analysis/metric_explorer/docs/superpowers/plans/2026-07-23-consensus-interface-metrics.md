# Consensus Interface Metrics Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Give the AF3-only Rosetta+Biotite interface metrics a worst-case 3-backend (AF3/B2/ESM) consensus, following the existing house convention, as pipeline infrastructure plus a one-time backfill of the current master.

**Architecture:** Mirror the proven `red`/`hbond` consensus pattern. Recompute the interface suite (`contacts_sidecar.py` + `interface_rosetta.py`, both `--path-col`-aware) on the B2 and ESM refold CIFs into per-backend sidecars via a sharded SLURM array; a new `interface_consensus_integrate.py` joins them onto the master and overwrites each bare column with `min`/`max(af3,b2,esm)`, preserving the AF3 value as `<col>_af3`. Wire into `run_analysis.sh` for go-forward; backfill the current master; fix the stale `refresh_data.sh` so the explorer inherits it.

**Tech Stack:** Python 3 (pandas, stdlib csv, multiprocessing), Bash + SLURM arrays, Apptainer (`pyrosetta.sif`, `boltzgen.sif`), Euler cluster.

## Global Constraints

- **House convention (verbatim):** bare column `M` = the consensus; `M_af3` = preserved original AF3; `M_b2`/`M_esm` = per-backend. **No `_consensus` suffix.** Matches `hbond_consensus_integrate.py` / `red_consensus_integrate.py`.
- **Direction map — MIN (higher=better):** `dsasa_total`, `dsasa_polar`, `dsasa_hphobic`, `sc_value`, `packstat`, `n_heavy_contacts_4A`, `n_heavy_contacts_4p5A`, `n_contact_res_prot`, `n_contact_res_rna`, `n_contacts_base`, `n_contacts_backbone`, `n_contacts_phosphate`, `n_contacts_sugar`, `frac_base_contacts`.
- **Direction map — MAX (lower=better):** `delta_unsat_hbonds`, `n_interface_clashes`, `clash_overlap_sum`, `worst_clash`.
- **SKIP (untouched, no per-backend cols):** `dsasa_hphobic_frac`, `aa_base_*` (20), `aa_bb_*` (20), and everything from the structural pass (`n_interface_contacts`, `interface_bsa_approx`, `radius_of_gyration_ca`, `helix/sheet/loop_fraction`, `sse_*`).
- **Master (source of truth):** `/cluster/work/wicky/lvalsecchi/rna_binder/analyses/cross_campaign/app_metrics_100k.csv` (`$RB` = `/cluster/work/wicky/lvalsecchi/rna_binder`).
- **Join key:** `design_id` (unique).
- **Idempotent + backed-up:** every mutating step drops prior `_af3/_b2/_esm` for managed cols and writes a `.bak_pre_iface_consensus` first.
- **Producers accept `--path-col`** (default `af3_msaoff_cif_path`) and `--output`; do NOT modify them.
- **Containers:** Biotite pass → `boltzgen.sif` with `PYTHONPATH=code/src:code/scripts`; Rosetta pass → `pyrosetta.sif`. Bind `--bind /cluster/work/wicky --bind /tmp --env MPLCONFIGDIR=/tmp`, `export TMPDIR=/tmp` (fresh node-local `/tmp` per array task — the `/local`-scratch-leak fix).
- **No commit trailer:** plain commit messages, no `Co-Authored-By`.

## File Structure

New (author locally in `~/scratch_scoring/scripts/`, unit-test locally, deploy to `$RB/code/scripts/` on Euler):
- `interface_consensus_integrate.py` — join sidecars + build consensus (core logic; pure pandas).
- `build_iface_be_manifests.py` — slim `design_id,<cif_col>` manifests for the B2/ESM passes.
- `interface_be_parallel.sh` — sharded per-backend recompute driver → `iface_be_<be>.csv` sidecars.
- `validate_consensus.py` — post-backfill invariant checker.
- `test_interface_consensus.py` — pytest for the two Python scripts (local only).

Modify:
- `$RB/code/scripts/run_analysis.sh` — add B2/ESM passes + integrate after the AF3 interface pass.
- `~/Documents/metric_explorer/refresh_data.sh` — fix stale `SRC`.
- `~/Documents/metric_explorer/app.py` — update METRIC descriptions (bare = consensus).

---

### Task 1: `interface_consensus_integrate.py` (core consensus logic)

**Files:**
- Create: `~/scratch_scoring/scripts/interface_consensus_integrate.py`
- Test: `~/scratch_scoring/scripts/test_interface_consensus.py`

**Interfaces:**
- Consumes: a master CSV with `design_id` + bare managed cols (AF3 values); a `be_dir` containing `iface_be_b2.csv` and `iface_be_esm.csv`, each with `design_id` + recomputed managed cols.
- Produces: CLI `python3 interface_consensus_integrate.py <master.csv> <be_dir> [--dry-run]`. Overwrites master in place (backup `.bak_pre_iface_consensus`); adds `M_af3/M_b2/M_esm` for managed cols, sets bare `M`=min/max, adds `consensus_n_structures` (0–3). Module constants `MIN_COLS`, `MAX_COLS`, `MANAGED` importable by tests.

- [ ] **Step 1: Write the failing test**

```python
# ~/scratch_scoring/scripts/test_interface_consensus.py
import os, subprocess, sys, textwrap
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
SCRIPT = os.path.join(HERE, "interface_consensus_integrate.py")

def _write(tmp):
    # master: 3 designs, AF3 values in bare cols + one SKIP col + one unrelated col
    master = pd.DataFrame({
        "design_id": ["d1", "d2", "d3"],
        "sc_value":      [0.70, 0.60, 0.50],   # MIN
        "n_contacts_base":[100,   80,   60],   # MIN
        "worst_clash":   [0.5,  0.5,  0.5],    # MAX
        "aa_base_ARG":   [5,     4,    3],     # SKIP (must stay untouched)
        "iptm":          [0.9,  0.8,  0.7],    # unrelated (untouched)
    })
    mpath = os.path.join(tmp, "master.csv"); master.to_csv(mpath, index=False)
    be = os.path.join(tmp, "be"); os.makedirs(be, exist_ok=True)
    pd.DataFrame({"design_id": ["d1","d2","d3"],
                  "sc_value":[0.65,0.55,0.80], "n_contacts_base":[90,120,50],
                  "worst_clash":[0.9,0.4,0.6], "aa_base_ARG":[9,9,9]}
                 ).to_csv(os.path.join(be,"iface_be_b2.csv"), index=False)
    pd.DataFrame({"design_id": ["d1","d2","d3"],
                  "sc_value":[0.60,0.90,0.55], "n_contacts_base":[110,70,55],
                  "worst_clash":[0.7,0.3,0.55], "aa_base_ARG":[1,1,1]}
                 ).to_csv(os.path.join(be,"iface_be_esm.csv"), index=False)
    return mpath, be

def test_consensus_min_max_and_preservation(tmp_path):
    tmp = str(tmp_path); mpath, be = _write(tmp)
    subprocess.run([sys.executable, SCRIPT, mpath, be], check=True)
    out = pd.read_csv(mpath)
    # MIN: bare = min(af3,b2,esm); _af3 preserved
    assert list(out["sc_value_af3"]) == [0.70, 0.60, 0.50]
    assert list(out["sc_value"]) == [0.60, 0.55, 0.50]
    assert list(out["n_contacts_base"]) == [90, 70, 50]
    # MAX: bare = max
    assert list(out["worst_clash"]) == [0.9, 0.5, 0.6]
    # SKIP untouched, no per-backend cols added
    assert list(out["aa_base_ARG"]) == [5, 4, 3]
    assert "aa_base_ARG_b2" not in out.columns
    # n_structures = 3 everywhere
    assert list(out["consensus_n_structures"]) == [3, 3, 3]

def test_idempotent(tmp_path):
    tmp = str(tmp_path); mpath, be = _write(tmp)
    subprocess.run([sys.executable, SCRIPT, mpath, be], check=True)
    first = pd.read_csv(mpath)
    subprocess.run([sys.executable, SCRIPT, mpath, be], check=True)  # run again
    second = pd.read_csv(mpath)
    assert list(first["sc_value"]) == list(second["sc_value"])
    assert sum(c.endswith("_af3") for c in second.columns) == sum(c.endswith("_af3") for c in first.columns)
```

- [ ] **Step 2: Run test to verify it fails**

Run: `cd ~/scratch_scoring/scripts && python3 -m pytest test_interface_consensus.py -q`
Expected: FAIL (`interface_consensus_integrate.py` does not exist / cannot run).

- [ ] **Step 3: Write minimal implementation**

```python
# ~/scratch_scoring/scripts/interface_consensus_integrate.py
#!/usr/bin/env python3
"""Join B2 + ESM interface-suite sidecars onto the master and build the worst-case
3-backend consensus, following hbond/red_consensus_integrate: bare column BECOMES the
consensus (min for higher-is-better, max for lower-is-better), original AF3 preserved as
<col>_af3, per-backend as <col>_b2/_esm. SKIP metrics untouched. Backs up; idempotent.

    python3 interface_consensus_integrate.py <master.csv> <be_dir> [--dry-run]

<be_dir> holds iface_be_b2.csv and iface_be_esm.csv, each with 'design_id' + the
recomputed interface columns for that backend.
"""
import argparse, os, shutil, sys
import pandas as pd

MIN_COLS = ["dsasa_total", "dsasa_polar", "dsasa_hphobic", "sc_value", "packstat",
            "n_heavy_contacts_4A", "n_heavy_contacts_4p5A", "n_contact_res_prot",
            "n_contact_res_rna", "n_contacts_base", "n_contacts_backbone",
            "n_contacts_phosphate", "n_contacts_sugar", "frac_base_contacts"]
MAX_COLS = ["delta_unsat_hbonds", "n_interface_clashes", "clash_overlap_sum", "worst_clash"]
MANAGED = MIN_COLS + MAX_COLS
SENTINEL = "sc_value"   # presence => that backend produced a Rosetta score for the design


def load_sidecar(be_dir, be):
    path = os.path.join(be_dir, f"iface_be_{be}.csv")
    if not os.path.exists(path):
        sys.exit(f"missing sidecar {path}")
    d = pd.read_csv(path, low_memory=False).drop_duplicates("design_id")
    have = [c for c in MANAGED if c in d.columns]
    d = d[["design_id"] + have].rename(columns={c: f"{c}_{be}" for c in have})
    print(f"{be}: {len(d)} rows, {len(have)} managed cols")
    return d


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("master")
    ap.add_argument("be_dir")
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args()

    m = pd.read_csv(args.master, low_memory=False)
    drop = [f"{c}_{s}" for c in MANAGED for s in ("af3", "b2", "esm")
            if f"{c}_{s}" in m.columns]
    if drop:
        m = m.drop(columns=drop)
    m = (m.merge(load_sidecar(args.be_dir, "b2"), on="design_id", how="left")
          .merge(load_sidecar(args.be_dir, "esm"), on="design_id", how="left"))

    for col in MANAGED:
        if col not in m.columns:
            print(f"WARN: master missing managed col {col}; skipping")
            continue
        m[f"{col}_af3"] = m[col]
        trio = [f"{col}_af3", f"{col}_b2", f"{col}_esm"]
        for t in trio:
            if t not in m.columns:
                m[t] = pd.NA
        vals = m[trio].apply(pd.to_numeric, errors="coerce")
        m[col] = vals.min(axis=1) if col in MIN_COLS else vals.max(axis=1)

    sent = [f"{SENTINEL}_{s}" for s in ("af3", "b2", "esm") if f"{SENTINEL}_{s}" in m.columns]
    m["consensus_n_structures"] = (m[sent].apply(pd.to_numeric, errors="coerce")
                                   .notna().sum(axis=1)) if sent else 0

    if args.dry_run:
        for col in ("sc_value", "n_contacts_base", "worst_clash"):
            if f"{col}_af3" in m.columns:
                print(f"  {col}: af3_med={m[f'{col}_af3'].median():.4g} "
                      f"consensus_med={m[col].median():.4g}")
        print("consensus_n_structures:\n"
              f"{m['consensus_n_structures'].value_counts().sort_index()}")
        return

    shutil.copy(args.master, args.master + ".bak_pre_iface_consensus")
    m.to_csv(args.master, index=False)
    print(f"wrote {args.master}: {len(m)} rows x {len(m.columns)} cols")


if __name__ == "__main__":
    main()
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `cd ~/scratch_scoring/scripts && python3 -m pytest test_interface_consensus.py -q`
Expected: PASS (2 passed).

- [ ] **Step 5: Commit**

```bash
cd ~/scratch_scoring && git add scripts/interface_consensus_integrate.py scripts/test_interface_consensus.py && git commit -m "interface consensus integrate: 3-backend worst-case, house convention"
```
(If `~/scratch_scoring` is not a git repo, skip the commit; the deploy in Task 4 is the checkpoint.)

---

### Task 2: `build_iface_be_manifests.py` (slim per-backend manifests)

**Files:**
- Create: `~/scratch_scoring/scripts/build_iface_be_manifests.py`
- Test: append to `~/scratch_scoring/scripts/test_interface_consensus.py`

**Interfaces:**
- Consumes: master CSV with `design_id`, `b2_cif_path`, `esm_cif_path`.
- Produces: CLI `python3 build_iface_be_manifests.py <master.csv> <outdir>` → `<outdir>/iface_manifest_b2.csv` (cols `design_id,b2_cif_path`) and `iface_manifest_esm.csv` (cols `design_id,esm_cif_path`). Column name matches the `--path-col` the driver passes.

- [ ] **Step 1: Write the failing test** (append)

```python
def test_manifests(tmp_path):
    import pandas as pd, subprocess, sys, os
    tmp = str(tmp_path)
    m = pd.DataFrame({"design_id":["d1","d2"], "b2_cif_path":["a!x.cif","b!y.cif"],
                      "esm_cif_path":["c!x.cif", ""], "stack":["s","s"]})
    mp = os.path.join(tmp,"m.csv"); m.to_csv(mp, index=False)
    BUILD = os.path.join(HERE, "build_iface_be_manifests.py")
    subprocess.run([sys.executable, BUILD, mp, tmp], check=True)
    b2 = pd.read_csv(os.path.join(tmp,"iface_manifest_b2.csv"))
    assert list(b2.columns) == ["design_id","b2_cif_path"] and len(b2) == 2
    esm = pd.read_csv(os.path.join(tmp,"iface_manifest_esm.csv"))
    assert len(esm) == 1  # empty esm ref dropped
```

- [ ] **Step 2: Run to verify it fails**

Run: `cd ~/scratch_scoring/scripts && python3 -m pytest test_interface_consensus.py::test_manifests -q`
Expected: FAIL (script missing).

- [ ] **Step 3: Write minimal implementation**

```python
# ~/scratch_scoring/scripts/build_iface_be_manifests.py
#!/usr/bin/env python3
"""Build slim per-backend interface manifests (design_id, <backend>_cif_path) from the
master for the B2/ESM interface recompute. The producers match rows by --path-col, so each
manifest keeps design_id + that backend's CIF ref column.

    python3 build_iface_be_manifests.py <master.csv> <outdir>
"""
import csv, os, sys
csv.field_size_limit(10 ** 7)

SRC, OUT = sys.argv[1], sys.argv[2]
os.makedirs(OUT, exist_ok=True)
rows = list(csv.DictReader(open(SRC)))
for be, col in [("b2", "b2_cif_path"), ("esm", "esm_cif_path")]:
    if not rows or col not in rows[0]:
        sys.exit(f"missing column {col}")
    recs = [(r["design_id"], r[col]) for r in rows if r.get(col)]
    out = os.path.join(OUT, f"iface_manifest_{be}.csv")
    with open(out, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["design_id", col])
        w.writerows(recs)
    print(f"{be}: {len(recs)} refs -> {out}")
```

- [ ] **Step 4: Run to verify it passes**

Run: `cd ~/scratch_scoring/scripts && python3 -m pytest test_interface_consensus.py -q`
Expected: PASS (3 passed).

- [ ] **Step 5: Commit**

```bash
cd ~/scratch_scoring && git add scripts/build_iface_be_manifests.py scripts/test_interface_consensus.py && git commit -m "build slim B2/ESM interface manifests"
```

---

### Task 3: `interface_be_parallel.sh` (sharded per-backend recompute driver)

**Files:**
- Create: `~/scratch_scoring/scripts/interface_be_parallel.sh`

**Interfaces:**
- Consumes: a slim manifest `iface_manifest_<be>.csv` (`design_id,<col>`), the path-col name, a backend tag, N chunks.
- Produces: shards the manifest, array-runs `contacts_sidecar.py` + `interface_rosetta.py --resume --path-col <col>` per chunk, merges → `<CK>/iface_be_<be>.csv` (design_id + col + all managed cols). Adapted from `interface_metrics_parallel.sh` (the AF3 driver) with `--path-col` threaded and non-in-place output.

- [ ] **Step 1: Write the script**

```bash
#!/bin/bash
# interface_be_parallel.sh — recompute the interface suite (Biotite contacts +
# Rosetta dSASA/sc/packstat) on ONE backend's CIFs into a per-backend sidecar,
# sharded across a SLURM array with fresh node-local /tmp per task. Does NOT touch
# the master. Mirrors interface_metrics_parallel.sh but threads --path-col and writes
# a backend sidecar.
#   sbatch interface_be_parallel.sh <iface_manifest_<be>.csv> <path_col> <be_tag> [N=40]
#SBATCH --job-name=ifbe_prep
#SBATCH --partition=normal.4h
#SBATCH --time=00:40:00
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=6G
#SBATCH --output=/cluster/work/wicky/lvalsecchi/rna_binder/analyses/_ifbe_prep_%j.log
set -euo pipefail
RB=/cluster/work/wicky/lvalsecchi/rna_binder
cd "$RB"
MAN="${1:?usage: interface_be_parallel.sh <manifest.csv> <path_col> <be_tag> [N]}"
COL="${2:?need path-col name}"
BE="${3:?need backend tag}"
N="${4:-40}"
CK="${MAN%.csv}_${BE}_chunks"
[ -f "$MAN" ] || { echo "FATAL: no $MAN"; exit 1; }

rm -rf "$CK"; mkdir -p "$CK/logs"
python3 - "$MAN" "$CK" "$N" <<'PY'
import sys, math, pandas as pd
inp, ck, N = sys.argv[1], sys.argv[2], int(sys.argv[3])
df = pd.read_csv(inp, low_memory=False)
per = math.ceil(len(df) / N)
for i in range(N):
    df.iloc[i*per:(i+1)*per].to_csv(f"{ck}/chunk_{i:03d}.csv", index=False)
print(f"{len(df):,} rows -> {N} chunks (~{per} each)")
PY

cat > "$CK/_compute.sbatch" <<EOF
#!/bin/bash
#SBATCH --job-name=ifbe_${BE}
#SBATCH --partition=normal.4h
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=3G
#SBATCH --requeue
#SBATCH --output=$CK/logs/%A_%a.log
set -uo pipefail
cd $RB
IDX=\$(printf "%03d" "\$SLURM_ARRAY_TASK_ID")
IN="$CK/chunk_\${IDX}.csv"
[ -s "\$IN" ] || { echo "empty \$IN"; exit 0; }
export TMPDIR=/tmp
J="\${SLURM_CPUS_PER_TASK:-8}"
echo "chunk \$IDX contacts ($BE) \$(date)"
apptainer exec --bind /cluster/work/wicky --bind /tmp --env MPLCONFIGDIR=/tmp \
  /cluster/work/wicky/containers/boltzgen.sif \
  env PYTHONPATH=code/src:code/scripts python3 \
  code/scripts/contacts_sidecar.py --csv "\$IN" --path-col "$COL" --jobs "\$J" --progress
echo "chunk \$IDX rosetta ($BE) \$(date)"
apptainer exec --bind /cluster/work/wicky --bind /tmp --env MPLCONFIGDIR=/tmp \
  /cluster/work/wicky/containers/pyrosetta.sif \
  python3 code/scripts/interface_rosetta.py --csv "\$IN" --path-col "$COL" --resume --jobs "\$J" --progress
echo "chunk \$IDX DONE \$(date)"
EOF

cat > "$CK/_merge.sbatch" <<EOF
#!/bin/bash
#SBATCH --job-name=ifbe_${BE}_merge
#SBATCH --partition=normal.4h
#SBATCH --time=00:20:00
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=6G
#SBATCH --output=$CK/_merge_%j.log
set -euo pipefail
cd $RB
OUT="$CK/iface_be_${BE}.csv"
FIRST=\$(ls $CK/chunk_*.csv | sort | head -1)
head -1 "\$FIRST" > "\$OUT"
for f in \$(ls $CK/chunk_*.csv | sort); do tail -n +2 "\$f" >> "\$OUT"; done
echo "merged -> \$OUT (\$(( \$(wc -l < "\$OUT") - 1 )) rows)"
python3 - "\$OUT" <<'PYX'
import csv, sys
csv.field_size_limit(10**7)
rows=list(csv.DictReader(open(sys.argv[1])))
for c in ["dsasa_total","sc_value","n_heavy_contacts_4A","frac_base_contacts","worst_clash"]:
    n=sum(1 for r in rows if str(r.get(c,"")).strip() not in ("","nan"))
    print(f"  {c:20s} {n}/{len(rows)}")
PYX
EOF

ARR=\$(sbatch --parsable --array=0-$((N-1)) "$CK/_compute.sbatch")
echo "compute array ($BE): $ARR"
MRG=\$(sbatch --parsable --dependency=afterok:"$ARR" "$CK/_merge.sbatch")
echo "merge ($BE): $MRG"
```

- [ ] **Step 2: Static-check the script**

Run: `bash -n ~/scratch_scoring/scripts/interface_be_parallel.sh` (syntax) and `shellcheck ~/scratch_scoring/scripts/interface_be_parallel.sh || true`
Expected: `bash -n` produces no output (valid syntax); review any shellcheck notes.

- [ ] **Step 3: Commit**

```bash
cd ~/scratch_scoring && git add scripts/interface_be_parallel.sh && git commit -m "sharded per-backend interface recompute driver"
```

---

### Task 4: Deploy to Euler + smoke-test on 5 designs

**Files:** none new; deploys Task 1–3 artifacts to `$RB/code/scripts/`.

- [ ] **Step 1: Deploy the four scripts to Euler**

```bash
cd ~/scratch_scoring/scripts
scp interface_consensus_integrate.py build_iface_be_manifests.py interface_be_parallel.sh \
    euler:/cluster/work/wicky/lvalsecchi/rna_binder/code/scripts/
```
Expected: 3 files copied (validator deploys in Task 6).

- [ ] **Step 2: Build a 5-row smoke manifest from the master**

```bash
ssh euler 'cd /cluster/work/wicky/lvalsecchi/rna_binder && \
  head -6 analyses/cross_campaign/app_metrics_100k.csv > /tmp/mini_master.csv && \
  python3 code/scripts/build_iface_be_manifests.py /tmp/mini_master.csv /tmp/ifbe_smoke && \
  head -3 /tmp/ifbe_smoke/iface_manifest_b2.csv'
```
Expected: prints `b2: 5 refs -> …` and `esm: 5 refs -> …`, then a 2-column header `design_id,b2_cif_path`.

- [ ] **Step 3: Run the B2 pass on the 5-row manifest (N=1) and confirm the sidecar fills**

```bash
ssh euler 'cd /cluster/work/wicky/lvalsecchi/rna_binder && \
  sbatch code/scripts/interface_be_parallel.sh /tmp/ifbe_smoke/iface_manifest_b2.csv b2_cif_path b2 1'
```
Then poll: `ssh euler 'squeue -u lvalsecchi -n ifbe_b2,ifbe_b2_merge'` until clear, and:
```bash
ssh euler 'F=/tmp/ifbe_smoke/iface_manifest_b2_b2_chunks/iface_be_b2.csv; wc -l $F; \
  head -1 $F | tr "," "\n" | grep -E "^(dsasa_total|sc_value|frac_base_contacts|worst_clash)$"'
```
Expected: 6 lines (5 rows + header); the four managed columns present. If clashes/sc are blank on all 5, stop and inspect a `logs/*.log` (a container/path problem, not a plan problem).

---

### Task 5: Backfill the current master (B2 + ESM passes + integrate + validate)

**Files:**
- Create: `~/scratch_scoring/scripts/validate_consensus.py` (+ deploy)

- [ ] **Step 1: Write `validate_consensus.py`**

```python
# ~/scratch_scoring/scripts/validate_consensus.py
#!/usr/bin/env python3
"""Assert the consensus invariants on a backfilled master.
    python3 validate_consensus.py <master.csv>
Exits non-zero on any violation."""
import sys
import pandas as pd
from interface_consensus_integrate import MIN_COLS, MAX_COLS

m = pd.read_csv(sys.argv[1], low_memory=False)
bad = 0
for col in MIN_COLS + MAX_COLS:
    trio = [f"{col}_af3", f"{col}_b2", f"{col}_esm"]
    if col not in m.columns or not all(t in m.columns for t in trio):
        print(f"MISSING cols for {col}"); bad += 1; continue
    vals = m[trio].apply(pd.to_numeric, errors="coerce")
    want = vals.min(axis=1) if col in MIN_COLS else vals.max(axis=1)
    got = pd.to_numeric(m[col], errors="coerce")
    mism = int((~((got - want).abs() < 1e-6) & want.notna()).sum())
    if mism:
        print(f"VIOLATION {col}: {mism} rows bare != {'min' if col in MIN_COLS else 'max'}(trio)")
        bad += 1
print(f"consensus_n_structures value counts:\n{m['consensus_n_structures'].value_counts().sort_index()}")
sys.exit(1 if bad else 0)
```

- [ ] **Step 2: Deploy validator + build full manifests**

```bash
scp ~/scratch_scoring/scripts/validate_consensus.py euler:/cluster/work/wicky/lvalsecchi/rna_binder/code/scripts/
ssh euler 'cd /cluster/work/wicky/lvalsecchi/rna_binder && \
  python3 code/scripts/build_iface_be_manifests.py \
    analyses/cross_campaign/app_metrics_100k.csv analyses/cross_campaign/ifbe'
```
Expected: `b2: 100000 refs -> …`, `esm: 100000 refs -> …`.

- [ ] **Step 3: Launch both backend passes (N=40 each)**

```bash
ssh euler 'cd /cluster/work/wicky/lvalsecchi/rna_binder && \
  sbatch code/scripts/interface_be_parallel.sh analyses/cross_campaign/ifbe/iface_manifest_b2.csv  b2_cif_path  b2  40 && \
  sbatch code/scripts/interface_be_parallel.sh analyses/cross_campaign/ifbe/iface_manifest_esm.csv esm_cif_path esm 40'
```
Then wait for both merges: `ssh euler 'squeue -u lvalsecchi -n ifbe_b2,ifbe_esm,ifbe_b2_merge,ifbe_esm_merge'` until empty. Verify fill fractions in the two merge logs are ~100% (allow a small pathological-CIF tail, matching the natural `rosetta_interface_ok` behavior).

- [ ] **Step 4: Stage sidecars into a single be_dir and integrate**

```bash
ssh euler 'cd /cluster/work/wicky/lvalsecchi/rna_binder/analyses/cross_campaign && \
  mkdir -p ifbe/be && \
  cp ifbe/iface_manifest_b2_b2_chunks/iface_be_b2.csv   ifbe/be/ && \
  cp ifbe/iface_manifest_esm_esm_chunks/iface_be_esm.csv ifbe/be/ && \
  cd /cluster/work/wicky/lvalsecchi/rna_binder && \
  python3 code/scripts/interface_consensus_integrate.py \
    analyses/cross_campaign/app_metrics_100k.csv \
    analyses/cross_campaign/ifbe/be --dry-run'
```
Expected: dry-run prints `sc_value af3_med=… consensus_med=…` (consensus ≤ AF3), `worst_clash` consensus ≥ AF3, and `consensus_n_structures` mostly 3. Then re-run **without** `--dry-run` to write.

- [ ] **Step 5: Validate invariants on the written master**

```bash
ssh euler 'cd /cluster/work/wicky/lvalsecchi/rna_binder/code/scripts && \
  python3 validate_consensus.py ../../analyses/cross_campaign/app_metrics_100k.csv'
```
Expected: no `VIOLATION`/`MISSING` lines; exit 0; `consensus_n_structures` mostly `3`.

- [ ] **Step 6: Commit deploy note**

```bash
cd ~/scratch_scoring && git add scripts/validate_consensus.py && git commit -m "consensus invariant validator + backfill run"
```

---

### Task 6: Wire into `run_analysis.sh` (go-forward infra)

**Files:**
- Modify: `$RB/code/scripts/run_analysis.sh` (after the interface pass; edit on Euler, keep a local copy in `~/scratch_scoring/scripts/run_analysis.sh` for reference)

- [ ] **Step 1: Read the current tail of `run_analysis.sh` to find the interface-pass insertion point**

Run: `ssh euler 'sed -n "1,120p" /cluster/work/wicky/lvalsecchi/rna_binder/code/scripts/run_analysis.sh'`
Expected: locate where AF3 interface metrics are added (or confirm they are added by a separate `interface_metrics_parallel.sh` call downstream — if `run_analysis.sh` does not itself call the interface pass, add the block after `compute_metrics.py`).

- [ ] **Step 2: Add the B2/ESM + integrate block** (exact lines depend on Step 1; insert after the AF3 interface pass, before final aggregation)

```bash
# --- 3-backend interface consensus (B2/ESM recompute + integrate) ---
python3 "$S/build_iface_be_manifests.py" "$OUT/${MASTER}" "$OUT/ifbe"
IB2=$(sbatch --parsable "$S/interface_be_parallel.sh" "$OUT/ifbe/iface_manifest_b2.csv"  b2_cif_path  b2  40)
IES=$(sbatch --parsable "$S/interface_be_parallel.sh" "$OUT/ifbe/iface_manifest_esm.csv" esm_cif_path esm 40)
# integrate runs after both backend merges complete (chain via afterok on their merge jobs)
```
(Adapt variable names `$OUT`/`$MASTER`/`$S` to those already in `run_analysis.sh`; the integrate call mirrors Task 5 Step 4 and should be gated on both merges via `--dependency=afterok`.)

- [ ] **Step 3: Dry-run validate the edit**

Run: `ssh euler 'bash -n /cluster/work/wicky/lvalsecchi/rna_binder/code/scripts/run_analysis.sh'`
Expected: no syntax errors.

- [ ] **Step 4: Commit** (local reference copy)

```bash
cd ~/scratch_scoring && git add scripts/run_analysis.sh && git commit -m "run_analysis: emit 3-backend interface consensus for future campaigns"
```

---

### Task 7: Fix `refresh_data.sh` + update `app.py` descriptions + refresh

**Files:**
- Modify: `~/Documents/metric_explorer/refresh_data.sh`
- Modify: `~/Documents/metric_explorer/app.py`

- [ ] **Step 1: Fix the stale `SRC` in `refresh_data.sh`**

```bash
# in ~/Documents/metric_explorer/refresh_data.sh, change:
#   SRC=…/analyses/cross_campaign/combined_sample10k_metrics.csv
# to:
SRC=/cluster/work/wicky/lvalsecchi/rna_binder/analyses/cross_campaign/app_metrics_100k.csv
```

- [ ] **Step 2: Pull the refreshed master**

Run: `cd ~/Documents/metric_explorer && bash refresh_data.sh`
Expected: `refreshed: 100000 rows`; the "live metric columns" grep lists `dsasa_total`, `sc_value`, `frac_base_contacts` (now consensus) and the header now also contains `sc_value_af3` etc.

- [ ] **Step 3: Update the affected METRIC descriptions in `app.py`**

For each consensus'd metric already in the `METRICS` list (e.g. `dsasa_total`, `sc_value`, `frac_base_contacts`, clash metrics), append to its description: `" — worst-of-3 consensus across AF3/B2/ESM (single-structure AF3 value in <col>_af3)."` Match the wording style of the existing `n_h_bonds_interface` entry.

- [ ] **Step 4: Verify the app loads and shows consensus columns**

Run: `cd ~/Documents/metric_explorer && python3 -c "import pandas as pd; d=pd.read_csv('data/metrics.csv', nrows=5); print([c for c in d.columns if c in ('sc_value','sc_value_af3','consensus_n_structures')])"`
Expected: `['sc_value', 'sc_value_af3', 'consensus_n_structures']`.

- [ ] **Step 5: Commit**

```bash
cd ~/Documents/metric_explorer && git add refresh_data.sh app.py && git commit -m "explorer: pull app_metrics_100k master; describe interface metrics as 3-backend consensus"
```

---

### Task 8: Regression — natural-vs-design comparison on `_af3`

**Files:**
- Create: `~/scratch_scoring/scripts/compare_nat_vs_design_af3.py` (adapt the existing comparison to read `<col>_af3` for the design side).

- [ ] **Step 1: Adapt the comparison to use single-structure AF3 columns**

Copy the earlier `compare_nat_vs_design.py` logic; for any managed metric, read `<col>_af3` from the design master (fall back to bare only if `_af3` absent). Compare against `metrics_natural_59.csv` as before.

- [ ] **Step 2: Run and record the deltas**

Run: `python3 ~/scratch_scoring/scripts/compare_nat_vs_design_af3.py`
Expected: the base-reading (`frac_base_contacts`) and aromatic-depletion signals persist; the H-bond and clash gaps are now single-structure-fair (no longer distorted by consensus). Save the output table alongside the plan.

- [ ] **Step 3: Commit**

```bash
cd ~/scratch_scoring && git add scripts/compare_nat_vs_design_af3.py && git commit -m "natural-vs-design comparison on single-structure AF3 columns"
```

---

## Self-Review

**Spec coverage:** house convention → Task 1; direction map → Global Constraints + Task 1 constants; per-backend recompute → Tasks 2–5; go-forward infra → Task 6; backfill → Task 5; `refresh_data.sh` fix → Task 7; app descriptions → Task 7; validation (invariants, fixture, spot-check via smoke, regression) → Tasks 1/4/5/8; caveat (compare on `_af3`) → Task 8. All spec sections map to a task.

**Placeholder scan:** No TBD/TODO. Task 6 Steps 1–2 depend on reading the live `run_analysis.sh` (its exact variable names aren't knowable offline) — the step gives the concrete block and names the adaptation explicitly rather than hand-waving; this is a genuine read-first dependency, not a placeholder.

**Type/name consistency:** `MIN_COLS`/`MAX_COLS`/`MANAGED` defined in Task 1 and imported by the Task 5 validator; sidecar filenames `iface_be_<be>.csv` consistent across Tasks 3/5; manifest column names (`design_id`,`b2_cif_path`/`esm_cif_path`) consistent Tasks 2/3/5; `design_id` join key throughout.

---

## As-built (executed 2026-07-24)

Executed inline via the executing-plans skill. Deviations from the plan above, and why:

1. **Tests are stdlib-runnable, not pytest.** `pytest` isn't installed locally; rewrote the
   tests as a plain `python3 test_interface_consensus.py` harness (tempfile + asserts).
   Same coverage (consensus min/max + `_af3` preservation, idempotency, manifest builder) —
   3/3 pass.

2. **Driver hardened for the `/local` leak (Task 3).** Before submitting, checked the
   metrics-`/local`-scratch-leak memory + the hardened `chunk.sh`. Added to
   `interface_be_parallel.sh`'s compute task: `#SBATCH --tmp=680000`,
   `export ARCH_READ_CACHE_MAX=3`, and **removed** the `TMPDIR=/tmp` override (so extraction
   lands in SLURM's reserved node-local `$TMPDIR`, e.g. `/scratch/tmp.<job>`). Each B2/ESM
   ref extracts a ~12 GB `.tar.zst` via `archive_read`, so this is mandatory. Verified live:
   `MinTmpDiskNode=680000`, `tmpdir=/scratch/tmp.<job>`.

3. **Fixed a bug in the plan's driver code.** The final `sbatch` submission lines used
   escaped `\$(sbatch …)` (copied from heredoc context) which would fail in the main script
   body. The deployed script uses the correct `interface_metrics_parallel.sh` form. **The
   deployed `$RB/code/scripts/*.sh|*.py` are the source of truth**, mirrored in
   `~/scratch_scoring/scripts/`.

4. **N=48 shards** (matching the proven AF3 `interface_metrics_parallel.sh`), not 40.

5. **Task 6 = orchestrator, not `run_analysis.sh`.** `run_analysis.sh` does **not** run the
   interface pass (interface metrics are applied to the cross-campaign master separately), so
   wiring there was wrong. Built `interface_consensus_all.sh <master> [N]` instead —
   one command: manifests → B2+ESM arrays → integrate+validate as an `afterok`-dependent job.
   Use it on any future master to keep consensus current.

6. **Producers:** the actual design-side producers are `interface_rosetta.py` (Rosetta) and
   `contacts_sidecar.py` (Biotite `interface_geom`), both `--path-col`-aware. `contacts_sidecar`
   emits `CONTACT_COLUMNS` (contacts/aa/clashes/frac_base) but **not** `n_interface_contacts`/
   `interface_bsa_approx` (structural pass) — those stayed out of the consensus set.

**Outcome:** master `app_metrics_100k.csv` 244→299 cols, all 100k `consensus_n_structures=3`,
`validate_consensus.py` = OK, backup `.bak_pre_iface_consensus`. Explorer refreshed +
relabeled. Fair natural-vs-design uses `_af3`: n_h_bonds nat 14 vs design 2 (d −9.4).
