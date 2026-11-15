# Seed Self-Consistency Test — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Measure how reproducible each refold backend is by folding 2,000 designs (400 random per stack) 5 times each on AF3, Boltz-2 and ESMFold2, then quantifying pose agreement (pairwise DockQ) and confidence spread (iPTM/pTM/pLDDT).

**Architecture:** A self-contained directory on Euler. A manifest is drawn from the 100k table; fold inputs are built directly from the binder `sequence` column (no archive extraction, since the 2,000 designs are scattered over ~2,000 archives). Three **local copies** of the production refold drivers run each backend in its native multi-sample mode (5 structures per design per pass). A chunked SLURM array publishes one packed archive + one CSV per chunk. A scoring script emits one row per design × backend.

**Tech Stack:** Python 3 + pandas/numpy (Euler `module load stack/2024-06 python/3.12.8`), apptainer containers (`af3_protenix_esm.sif`, `boltz2.sif`, `esmc.sif`), SLURM arrays, DockQ via the repo's existing `dockq_self.py`.

## Global Constraints

- `$RB` = `/cluster/work/wicky/lvalsecchi/rna_binder`. Everything for this experiment lives under `$RB/seedtest_2026-07-30/`.
- **NEVER modify production drivers.** `$RB/code/scripts/stages/refold_af3.sh`, `refold_b2.sh`, `refold_esm.sh`, `esmfold2_infer.py` are READ-ONLY for this work. All edits go into copies under `seedtest_2026-07-30/drivers/`. Verify with `git`-less checksum comparison before launch (Task 9).
- Subsample: **400 designs per stack, 5 stacks, 2,000 total**, `random_state=20260730`.
- Stack order everywhere: `RFD3+LiMPNN, HD+LiMPNN, BG+LiMPNN, BG+BIF, BD1+LiMPNN`.
- **5 draws per design per backend.** AF3 `MODEL_SEEDS_JSON="[1,2,3,4,5]"` + `AF3_NUM_DIFFUSION_SAMPLES=1`; B2 `B2_DIFFUSION_SAMPLES=5` + `B2_SEED=1`; ESM `ESM_NUM_DIFFUSION_SAMPLES=5` + `ESM_SEED=1`.
- MSA **off** (`msaoff`) for all backends, matching the existing 100k refolds.
- All heavy IO node-local (`$TMPDIR`); each chunk publishes exactly ONE packed archive + ONE csv to Lustre. Never leave loose CIFs on Lustre (inode quota).
- In pandas, the stack column MUST be accessed as `df["stack"]` — `df.stack` is the DataFrame method and raises `KeyError: False`.
- Local dev/test of pure-python pieces uses `/Users/leo/miniconda3/bin/python3` on the Mac; cluster execution uses Euler's module python.

---

### Task 1: Subsample manifest builder

**Files:**
- Create: `$RB/seedtest_2026-07-30/scripts/build_seedtest_manifest.py`
- Create: `$RB/seedtest_2026-07-30/scripts/test_build_manifest.py`
- Reads: `/Users/leo/Documents/metric_explorer/data/metrics_100k.csv` (local dev) or `$RB/analyses/cross_campaign/app_metrics_100k.csv` (cluster)

**Interfaces:**
- Consumes: nothing.
- Produces: `build_manifest(df, target_rna, n_per_stack=400, seed=20260730) -> pandas.DataFrame` with exactly the columns `design_id, stack, sequence, rna_target, uuid8` and 400 rows per stack. Written as TSV to `seedtest_manifest.tsv`.

- [ ] **Step 1: Write the failing test**

```python
# test_build_manifest.py
import pandas as pd
from build_seedtest_manifest import build_manifest

STACKS = ["RFD3+LiMPNN", "HD+LiMPNN", "BG+LiMPNN", "BG+BIF", "BD1+LiMPNN"]

def _fake_100k():
    rows = []
    for s in STACKS:
        for i in range(1000):
            rows.append({"design_id": f"{s}_{i}", "stack": s,
                         "sequence": "A" * 80, "uuid8": f"{i:08x}"})
    return pd.DataFrame(rows)

def test_400_per_stack_and_columns():
    m = build_manifest(_fake_100k(), "GUUCACUGCCGUAUAGGCAG", n_per_stack=400, seed=20260730)
    assert len(m) == 2000
    assert list(m.columns) == ["design_id", "stack", "sequence", "rna_target", "uuid8"]
    assert (m.groupby("stack").size() == 400).all()
    assert m["design_id"].is_unique
    assert (m["rna_target"] == "GUUCACUGCCGUAUAGGCAG").all()

def test_deterministic_for_same_seed():
    a = build_manifest(_fake_100k(), "GUUC", seed=20260730)
    b = build_manifest(_fake_100k(), "GUUC", seed=20260730)
    assert a["design_id"].tolist() == b["design_id"].tolist()

def test_different_seed_gives_different_sample():
    a = build_manifest(_fake_100k(), "GUUC", seed=1)
    b = build_manifest(_fake_100k(), "GUUC", seed=2)
    assert a["design_id"].tolist() != b["design_id"].tolist()

def test_rejects_null_or_short_sequences():
    df = _fake_100k()
    df.loc[df.index[:900], "sequence"] = None      # too few usable in stack 0
    try:
        build_manifest(df, "GUUC", n_per_stack=400, seed=1)
    except ValueError as e:
        assert "usable" in str(e).lower()
    else:
        raise AssertionError("should raise when a stack lacks 400 usable designs")
```

- [ ] **Step 2: Run test to verify it fails**

Run: `cd $RB/seedtest_2026-07-30/scripts && python3 -m pytest test_build_manifest.py -v`
Expected: FAIL — `ModuleNotFoundError: No module named 'build_seedtest_manifest'`

- [ ] **Step 3: Write minimal implementation**

```python
#!/usr/bin/env python3
"""Draw the seed-self-consistency subsample: 400 random designs per stack."""
import argparse
import pandas as pd

STACKS = ["RFD3+LiMPNN", "HD+LiMPNN", "BG+LiMPNN", "BG+BIF", "BD1+LiMPNN"]
COLUMNS = ["design_id", "stack", "sequence", "rna_target", "uuid8"]


def build_manifest(df, target_rna, n_per_stack=400, seed=20260730):
    out = []
    for s in STACKS:
        sub = df[df["stack"] == s]                      # df.stack is a METHOD - never use it
        sub = sub[sub["sequence"].notna() & (sub["sequence"].str.len() > 0)]
        if len(sub) < n_per_stack:
            raise ValueError(f"stack {s}: only {len(sub)} usable designs, need {n_per_stack}")
        out.append(sub.sample(n=n_per_stack, random_state=seed))
    m = pd.concat(out, ignore_index=True)
    m["rna_target"] = target_rna
    if "uuid8" not in m.columns:
        m["uuid8"] = ""
    return m[COLUMNS].reset_index(drop=True)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--metrics-csv", required=True)
    ap.add_argument("--rna-target", required=True)
    ap.add_argument("--n-per-stack", type=int, default=400)
    ap.add_argument("--seed", type=int, default=20260730)
    ap.add_argument("-o", "--out", required=True)
    a = ap.parse_args()
    df = pd.read_csv(a.metrics_csv, usecols=lambda c: c in
                     {"design_id", "stack", "sequence", "uuid8"}, low_memory=False)
    m = build_manifest(df, a.rna_target, a.n_per_stack, a.seed)
    m.to_csv(a.out, sep="\t", index=False)
    print(f"wrote {a.out}: {len(m)} designs, {m.groupby('stack').size().to_dict()}")


if __name__ == "__main__":
    main()
```

- [ ] **Step 4: Run test to verify it passes**

Run: `python3 -m pytest test_build_manifest.py -v`
Expected: PASS (4 passed)

- [ ] **Step 5: Resolve the real RNA target, then generate the real manifest**

The target must be the SAME 16-nt Csy4 target the original refolds used. Find it:

```bash
ssh euler 'cd /cluster/work/wicky/lvalsecchi/rna_binder && ls targets/ && head -40 targets/*.yaml 2>/dev/null | head -60'
```

Then generate and verify:

```bash
ssh euler 'cd /cluster/work/wicky/lvalsecchi/rna_binder/seedtest_2026-07-30 && \
  module load stack/2024-06 python/3.12.8 && \
  python3 scripts/build_seedtest_manifest.py \
    --metrics-csv ../analyses/cross_campaign/app_metrics_100k.csv \
    --rna-target "<RESOLVED_RNA_SEQ>" \
    -o seedtest_manifest.tsv && \
  awk -F"\t" "NR>1{c[\$2]++} END{for(k in c) print k, c[k]}" seedtest_manifest.tsv'
```

Expected: 2000 rows total, each stack printing exactly `400`.

- [ ] **Step 6: Commit**

```bash
cd $RB/seedtest_2026-07-30 && tar -cf - scripts | (cd . && echo staged)
# repo code/ is not git-tracked; snapshot instead:
cp scripts/build_seedtest_manifest.py scripts/build_seedtest_manifest.py.v1
```

---

### Task 2: Fold-input builder

**Files:**
- Create: `$RB/seedtest_2026-07-30/scripts/build_seedtest_inputs.py`
- Create: `$RB/seedtest_2026-07-30/scripts/test_build_inputs.py`
- Reference (READ ONLY, to copy the exact format): `$RB/code/scripts/stages/refold_af3.sh` lines ~120-215 (AF3 JSON), and any existing `inputs/*.yaml` under a `refold_b2` dir.

**Interfaces:**
- Consumes: `seedtest_manifest.tsv` from Task 1.
- Produces: `write_b2_yaml(design_id, protein_seq, rna_seq, out_dir) -> Path` and `write_af3_json(design_id, protein_seq, rna_seq, seeds, out_dir) -> Path`. ESM consumes the same YAML as B2 (that is what `esmfold2_infer.py` parses).

- [ ] **Step 1: Write the failing test**

```python
# test_build_inputs.py
import json, yaml, tempfile
from pathlib import Path
from build_seedtest_inputs import write_b2_yaml, write_af3_json

P, R = "MKTAYIAKQ", "GUUCACUGCC"

def test_b2_yaml_has_protein_A_and_rna_B():
    with tempfile.TemporaryDirectory() as td:
        p = write_b2_yaml("d1", P, R, Path(td))
        d = yaml.safe_load(p.read_text())
        seqs = d["sequences"]
        prot = [s["protein"] for s in seqs if "protein" in s][0]
        rna  = [s["rna"] for s in seqs if "rna" in s][0]
        assert prot["id"] == "A" and prot["sequence"] == P
        assert rna["id"] == "B" and rna["sequence"] == R

def test_af3_json_carries_all_five_seeds():
    with tempfile.TemporaryDirectory() as td:
        p = write_af3_json("d1", P, R, [1,2,3,4,5], Path(td))
        d = json.loads(p.read_text())
        assert d["modelSeeds"] == [1,2,3,4,5]
        chains = {c.get("protein", c.get("rna", {})).get("id") for c in d["sequences"]}
        assert chains == {"A", "B"}

def test_filenames_are_design_id_based():
    with tempfile.TemporaryDirectory() as td:
        assert write_b2_yaml("abc", P, R, Path(td)).name == "abc.yaml"
        assert write_af3_json("abc", P, R, [1], Path(td)).name == "abc.json"
```

- [ ] **Step 2: Run test to verify it fails**

Run: `python3 -m pytest test_build_inputs.py -v`
Expected: FAIL — `ModuleNotFoundError: No module named 'build_seedtest_inputs'`

- [ ] **Step 3: Write minimal implementation**

```python
#!/usr/bin/env python3
"""Build per-design fold inputs (protein chain A + RNA chain B) for all 3 backends."""
import argparse, json
from pathlib import Path
import pandas as pd
import yaml


def write_b2_yaml(design_id, protein_seq, rna_seq, out_dir):
    """Boltz-style YAML - consumed by BOTH refold_b2.sh and esmfold2_infer.py."""
    out_dir = Path(out_dir); out_dir.mkdir(parents=True, exist_ok=True)
    doc = {"version": 1,
           "sequences": [{"protein": {"id": "A", "sequence": protein_seq}},
                         {"rna": {"id": "B", "sequence": rna_seq}}]}
    p = out_dir / f"{design_id}.yaml"
    p.write_text(yaml.safe_dump(doc, sort_keys=False))
    return p


def write_af3_json(design_id, protein_seq, rna_seq, seeds, out_dir):
    out_dir = Path(out_dir); out_dir.mkdir(parents=True, exist_ok=True)
    doc = {"name": design_id,
           "modelSeeds": list(seeds),
           "sequences": [{"protein": {"id": "A", "sequence": protein_seq}},
                         {"rna": {"id": "B", "sequence": rna_seq}}],
           "dialect": "alphafold3", "version": 1}
    p = out_dir / f"{design_id}.json"
    p.write_text(json.dumps(doc, indent=2))
    return p


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--manifest", required=True)
    ap.add_argument("--out-root", required=True)
    ap.add_argument("--seeds", default="1,2,3,4,5")
    a = ap.parse_args()
    seeds = [int(x) for x in a.seeds.split(",")]
    m = pd.read_csv(a.manifest, sep="\t")
    root = Path(a.out_root)
    n = 0
    for _, r in m.iterrows():
        write_b2_yaml(r.design_id, r.sequence, r.rna_target, root / "b2")
        write_b2_yaml(r.design_id, r.sequence, r.rna_target, root / "esm")
        write_af3_json(r.design_id, r.sequence, r.rna_target, seeds, root / "af3")
        n += 1
    print(f"wrote inputs for {n} designs under {root}/(af3|b2|esm)")


if __name__ == "__main__":
    main()
```

- [ ] **Step 4: Run test to verify it passes**

Run: `python3 -m pytest test_build_inputs.py -v`
Expected: PASS (3 passed)

- [ ] **Step 5: Validate the AF3 JSON against the real production format**

The AF3 JSON must match what `refold_af3.sh` writes, or AF3 will reject it. Diff one generated file against a production-built one:

```bash
ssh euler 'cd /cluster/work/wicky/lvalsecchi/rna_binder && sed -n "120,215p" code/scripts/stages/refold_af3.sh'
```

Compare key-for-key with the generated JSON. If production sets extra keys (e.g. `userCCD`, different `version`), add them to `write_af3_json` and re-run the tests.

- [ ] **Step 6: Generate the real inputs and verify counts**

```bash
ssh euler 'cd /cluster/work/wicky/lvalsecchi/rna_binder/seedtest_2026-07-30 && \
  module load stack/2024-06 python/3.12.8 && \
  python3 scripts/build_seedtest_inputs.py --manifest seedtest_manifest.tsv --out-root inputs && \
  for b in af3 b2 esm; do echo "$b: $(ls inputs/$b | wc -l)"; done'
```

Expected: `af3: 2000`, `b2: 2000`, `esm: 2000`.

---

### Task 3: Local driver copies (production untouched)

**Files:**
- Create (copies): `$RB/seedtest_2026-07-30/drivers/refold_af3_seedtest.sh`, `refold_b2_seedtest.sh`, `refold_esm_seedtest.sh`, `esmfold2_infer_seedtest.py`
- Create: `$RB/seedtest_2026-07-30/drivers/CHECKSUMS.production.txt`
- READ ONLY: `$RB/code/scripts/stages/refold_{af3,b2,esm}.sh`, `esmfold2_infer.py`

**Interfaces:**
- Consumes: `inputs/{af3,b2,esm}/` from Task 2.
- Produces: three driver scripts invoked by the array in Task 5, each writing canonical CIFs into a `designs/` dir: AF3 `<id>_AF3_seed<S>_sample0.cif` (S=1..5), B2 `<id>_B2_sample<N>.cif` (N=0..4), ESM `<id>_ESM_sample<N>.cif` (N=0..4).

- [ ] **Step 1: Record production checksums so we can prove we did not touch them**

```bash
ssh euler 'cd /cluster/work/wicky/lvalsecchi/rna_binder && mkdir -p seedtest_2026-07-30/drivers && \
  md5sum code/scripts/stages/refold_af3.sh code/scripts/stages/refold_b2.sh \
         code/scripts/stages/refold_esm.sh code/scripts/stages/esmfold2_infer.py \
    > seedtest_2026-07-30/drivers/CHECKSUMS.production.txt && \
  cat seedtest_2026-07-30/drivers/CHECKSUMS.production.txt'
```

- [ ] **Step 2: Copy the four production files into `drivers/`**

```bash
ssh euler 'cd /cluster/work/wicky/lvalsecchi/rna_binder/seedtest_2026-07-30/drivers && \
  cp ../../code/scripts/stages/refold_af3.sh   refold_af3_seedtest.sh && \
  cp ../../code/scripts/stages/refold_b2.sh    refold_b2_seedtest.sh && \
  cp ../../code/scripts/stages/refold_esm.sh   refold_esm_seedtest.sh && \
  cp ../../code/scripts/stages/esmfold2_infer.py esmfold2_infer_seedtest.py && ls -la'
```

- [ ] **Step 3: Patch the B2 copy to pass `--seed`**

In `refold_b2_seedtest.sh`, find the `boltz predict` invocation (around line 220) and add the seed line:

```bash
    boltz predict "$INPUTS_DIR" \
        --out_dir "$OUTPUTS_DIR" \
        --cache "$BOLTZ_CACHE" \
        --output_format "$B2_OUTPUT_FORMAT" \
        --diffusion_samples "$B2_DIFFUSION_SAMPLES" \
        --seed "${B2_SEED:-1}" \
        $WRITE_FULL_PAE_FLAG
```

Verify only that one line changed:

Run: `diff ../../code/scripts/stages/refold_b2.sh refold_b2_seedtest.sh`
Expected: exactly one added line containing `--seed`.

- [ ] **Step 4: Patch the ESM copy to keep ALL diffusion samples**

In `esmfold2_infer_seedtest.py`, the current code writes only `sample0` and skips on `sample0` existing. Replace the per-design block so it iterates over every returned structure:

```python
        spi = parse_yaml_to_spi(y)
        n_s = args.num_diffusion_samples
        # resume-skip MUST key on the LAST expected sample, else a requeued job
        # treats a partially-written design as complete
        last = cdir / f"{name}_{TAG}_sample{n_s - 1}.cif"
        if last.exists():
            n_skip += 1
            print(f"[SKIP] {name} (all {n_s} samples present)", flush=True)
            continue
        result = builder.fold(model, spi,
                              num_loops=args.num_loops,
                              num_sampling_steps=args.num_sampling_steps,
                              num_diffusion_samples=n_s,
                              seed=args.seed)
        structures = getattr(result, "structures", None) or [result]
        for si, res_i in enumerate(structures):
            out_cif = cdir / f"{name}_{TAG}_sample{si}.cif"
            tmp = out_cif.with_suffix(".cif.part")
            tmp.write_text(res_i.to_cif() if hasattr(res_i, "to_cif") else res_i.cif)
            tmp.replace(out_cif)
            (cdir / f"{name}_{TAG}_sample{si}_confidence.json").write_text(json.dumps({
                "plddt": float(res_i.plddt.mean()),
                "ptm": float(res_i.ptm),
                "iptm": float(res_i.iptm),
                "seed": args.seed, "sample": si}))
```

**The exact attribute holding the 5 structures is unverified** — `result.structures` is the assumption. Confirm it in Step 5 and correct the code if the API differs.

- [ ] **Step 5: Probe the ESM multi-sample API on ONE design (this is the riskiest unknown)**

```bash
ssh euler 'cd /cluster/work/wicky/lvalsecchi/rna_binder/seedtest_2026-07-30 && \
  mkdir -p _probe/in _probe/out && cp inputs/esm/$(ls inputs/esm | head -1) _probe/in/ && \
  sbatch --wrap="apptainer exec --nv /cluster/work/wicky/containers/esmc.sif python3 \
    drivers/esmfold2_infer_seedtest.py --inputs _probe/in --canonical-dir _probe/out \
    --num-diffusion-samples 5 --seed 1" --gpus=1 --time=00:30:00 -o _probe/esm_probe.log'
```

Then: `ls _probe/out/*.cif` — Expected: **5** files, `*_ESM_sample0.cif` … `*_ESM_sample4.cif`.

If only 1 file appears, the API does not return multiple structures from one call. **Fallback:** revert the loop and instead run the driver 5× with `--seed 1..5`, writing `sample0` each time into 5 separate output dirs, then rename to `sample0..4`. ESM is ~4 s/design so this costs little. Record which path was taken.

- [ ] **Step 6: Confirm production files are byte-identical**

```bash
ssh euler 'cd /cluster/work/wicky/lvalsecchi/rna_binder && md5sum -c seedtest_2026-07-30/drivers/CHECKSUMS.production.txt'
```

Expected: four lines all ending `: OK`. If any says FAILED, restore that file immediately from a sibling copy and re-verify — production must not drift.

---

### Task 4: Calibration probe (measure before scaling)

**Files:**
- Create: `$RB/seedtest_2026-07-30/probe/probe_manifest.tsv` (20 designs)
- Create: `$RB/seedtest_2026-07-30/probe/PROBE_RESULTS.md`

**Interfaces:**
- Consumes: drivers from Task 3, inputs from Task 2.
- Produces: measured seconds-per-design and peak GPU memory per backend, which set the chunk size and walltime tier in Task 5.

- [ ] **Step 1: Build a 20-design probe manifest (4 per stack)**

```bash
ssh euler 'cd /cluster/work/wicky/lvalsecchi/rna_binder/seedtest_2026-07-30 && mkdir -p probe && \
  head -1 seedtest_manifest.tsv > probe/probe_manifest.tsv && \
  awk -F"\t" "NR>1{if(++c[\$2]<=4) print}" seedtest_manifest.tsv >> probe/probe_manifest.tsv && \
  wc -l probe/probe_manifest.tsv'
```

Expected: 21 lines (header + 20).

- [ ] **Step 2: Run the probe for each backend, one small GPU job each**

Copy the sbatch directives from a known-good sibling rather than hand-writing them:

```bash
ssh euler 'cd /cluster/work/wicky/lvalsecchi/rna_binder && sed -n "1,20p" apo_fold_array_v2v7.sh'
```

Submit one probe job per backend using those directives with `--time=02:00:00`, `--gpus=1`, no `gpumem` (>11g forces Features=gpuhe), `--nodelist=eu-g3-[062-077,083-084]`.

- [ ] **Step 3: Verify each probe produced 5 structures per design**

```bash
ssh euler 'cd /cluster/work/wicky/lvalsecchi/rna_binder/seedtest_2026-07-30/probe && \
  for b in af3 b2 esm; do echo "$b: $(ls out/$b/designs/*.cif 2>/dev/null | wc -l) cifs (expect 100)"; done'
```

Expected: 100 CIFs each (20 designs × 5). **If any backend yields 20, its multi-sample flag did not take effect — stop and fix before scaling.**

- [ ] **Step 4: Record runtime and VRAM**

```bash
ssh euler 'sacct -j <JOBID> --format=JobID,Elapsed,State,MaxRSS -P; \
  sacct -j <JOBID>.batch --format=JobID,TRESUsageInMax -P'
```

Write `probe/PROBE_RESULTS.md` with, per backend: total elapsed, seconds/design at 5 samples, peak `gres/gpumem`, and the derived chunk size for a ~3 h target chunk.

- [ ] **Step 5: Commit the probe results**

```bash
ssh euler 'cd /cluster/work/wicky/lvalsecchi/rna_binder/seedtest_2026-07-30 && cat probe/PROBE_RESULTS.md'
```

---

### Task 5: Chunked SLURM array

**Files:**
- Create: `$RB/seedtest_2026-07-30/seedtest_array.sh`
- Create: `$RB/seedtest_2026-07-30/scripts/make_chunks.py`
- Template (READ ONLY): `$RB/apo_fold_array_v2v7.sh`

**Interfaces:**
- Consumes: `seedtest_manifest.tsv`, `inputs/`, `drivers/`, chunk size from Task 4.
- Produces: per chunk `outputs/<backend>/chunk_<NNN>.tar.zst` + `outputs/<backend>/chunk_<NNN>_conf.csv` (design_id, backend, sample, seed, iptm, ptm, plddt).

- [ ] **Step 1: Write the chunk splitter**

```python
#!/usr/bin/env python3
"""Split the manifest into N-design chunks, one TSV per chunk."""
import argparse
from pathlib import Path
import pandas as pd

ap = argparse.ArgumentParser()
ap.add_argument("--manifest", required=True)
ap.add_argument("--chunk-size", type=int, required=True)
ap.add_argument("--out-dir", required=True)
a = ap.parse_args()
m = pd.read_csv(a.manifest, sep="\t")
out = Path(a.out_dir); out.mkdir(parents=True, exist_ok=True)
n = 0
for i in range(0, len(m), a.chunk_size):
    m.iloc[i:i + a.chunk_size].to_csv(out / f"chunk_{n:03d}.tsv", sep="\t", index=False)
    n += 1
print(f"{n} chunks of <= {a.chunk_size} designs -> {out}")
```

Run it with the size chosen in Task 4, then `ls chunks/ | wc -l` to confirm.

- [ ] **Step 2: Write the array script, copying directives from the apo array**

`seedtest_array.sh` must:
1. Take `BACKEND` (af3|b2|esm) from an env var and the chunk id from `$SLURM_ARRAY_TASK_ID`.
2. Be **idempotent**: `[ -f "$OUT_ARCHIVE" ] && [ -f "$OUT_CSV" ] && exit 0` at the top.
3. Stage that chunk's inputs into `$TMPDIR` (node-local), run the matching driver from `drivers/` with the backend's 5-sample env, and keep ALL output node-local.
4. Publish exactly one `.tar.zst` of the `designs/` dir and one `_conf.csv` to Lustre.
5. `#SBATCH --requeue` and `--open-mode=append` (preemption safety on bulk).

- [ ] **Step 3: Dry-run ONE array task before submitting the array**

```bash
ssh euler 'cd /cluster/work/wicky/lvalsecchi/rna_binder/seedtest_2026-07-30 && \
  BACKEND=esm sbatch --array=0-0 --time=02:00:00 --gpus=1 seedtest_array.sh'
```

Verify: `outputs/esm/chunk_000.tar.zst` and `chunk_000_conf.csv` both exist, and the archive contains 5 × chunk-size CIFs (`tar --zstd -tf ... | grep -c '\.cif$'`).

- [ ] **Step 4: Verify idempotency by re-running the same task**

Re-submit the identical command. Expected: the job exits almost immediately with the "already done" message, and the archive's mtime is unchanged.

- [ ] **Step 5: Submit the full arrays, throttled**

One array per backend, `%N` throttle to respect the 250-job scheduler window and stay neighbourly with the running v2v7 apo campaign. Use the walltime tier that matches the probe measurement (remember: walltime ROUTES to tier and overrides `--partition`; `gpu.bulk` submits get remapped, so reach bulk via `--hold` → `scontrol update partition/qos` → `scontrol release`).

- [ ] **Step 6: Monitor to completion**

```bash
ssh euler 'squeue -u lvalsecchi -o "%.18i %.30j %.8T %.10M" | grep -i seedtest | head -20; \
  cd /cluster/work/wicky/lvalsecchi/rna_binder/seedtest_2026-07-30 && \
  for b in af3 b2 esm; do echo "$b: $(ls outputs/$b/*.tar.zst 2>/dev/null | wc -l) chunks done"; done'
```

---

### Task 6: Scoring — pairwise DockQ + confidence spread

**Files:**
- Create: `$RB/seedtest_2026-07-30/scripts/score_seedtest.py`
- Create: `$RB/seedtest_2026-07-30/scripts/test_score_seedtest.py`
- Reference (READ ONLY): the repo's existing `dockq_self.py`

**Interfaces:**
- Consumes: `outputs/<backend>/chunk_*.tar.zst` and `chunk_*_conf.csv`.
- Produces: `seedtest_results.csv`, one row per design × backend, with columns
  `design_id, stack, backend, n_structures_found, n_pairs, dockq_pair_mean, dockq_pair_median, dockq_pair_min, dockq_pair_max, dockq_pair_std, iptm_mean, iptm_std, iptm_min, iptm_max, iptm_range, ptm_mean, ptm_std, ptm_min, ptm_max, ptm_range, plddt_mean, plddt_std, plddt_min, plddt_max, plddt_range`.

- [ ] **Step 1: Write the failing test for the pure-math parts**

```python
# test_score_seedtest.py
import math
from score_seedtest import pair_stats, spread_stats

def test_pair_stats_counts_ten_pairs_for_five_structures():
    vals = [0.5] * 10
    s = pair_stats(vals)
    assert s["n_pairs"] == 10
    assert math.isclose(s["dockq_pair_mean"], 0.5)
    assert math.isclose(s["dockq_pair_std"], 0.0)

def test_pair_stats_min_max():
    s = pair_stats([0.1, 0.9, 0.5])
    assert s["dockq_pair_min"] == 0.1
    assert s["dockq_pair_max"] == 0.9
    assert math.isclose(s["dockq_pair_median"], 0.5)

def test_pair_stats_empty_is_all_none():
    s = pair_stats([])
    assert s["n_pairs"] == 0 and s["dockq_pair_mean"] is None

def test_spread_stats_prefixes_and_range():
    s = spread_stats([0.2, 0.4, 0.6], "iptm")
    assert math.isclose(s["iptm_mean"], 0.4)
    assert math.isclose(s["iptm_min"], 0.2)
    assert math.isclose(s["iptm_max"], 0.6)
    assert math.isclose(s["iptm_range"], 0.4)

def test_spread_stats_single_value_has_zero_std():
    s = spread_stats([0.5], "ptm")
    assert s["ptm_std"] == 0.0 and s["ptm_range"] == 0.0
```

- [ ] **Step 2: Run test to verify it fails**

Run: `python3 -m pytest test_score_seedtest.py -v`
Expected: FAIL — `ModuleNotFoundError: No module named 'score_seedtest'`

- [ ] **Step 3: Write minimal implementation of the stats helpers**

```python
#!/usr/bin/env python3
"""Score seed self-consistency: pairwise DockQ + confidence spread per design x backend."""
import itertools
import statistics as st


def pair_stats(vals):
    vals = [v for v in vals if v is not None]
    if not vals:
        return {"n_pairs": 0, "dockq_pair_mean": None, "dockq_pair_median": None,
                "dockq_pair_min": None, "dockq_pair_max": None, "dockq_pair_std": None}
    return {"n_pairs": len(vals),
            "dockq_pair_mean": st.fmean(vals),
            "dockq_pair_median": st.median(vals),
            "dockq_pair_min": min(vals),
            "dockq_pair_max": max(vals),
            "dockq_pair_std": st.pstdev(vals) if len(vals) > 1 else 0.0}


def spread_stats(vals, prefix):
    vals = [v for v in vals if v is not None]
    if not vals:
        return {f"{prefix}_{k}": None for k in ("mean", "std", "min", "max", "range")}
    return {f"{prefix}_mean": st.fmean(vals),
            f"{prefix}_std": st.pstdev(vals) if len(vals) > 1 else 0.0,
            f"{prefix}_min": min(vals),
            f"{prefix}_max": max(vals),
            f"{prefix}_range": max(vals) - min(vals)}


def all_pairs(structures):
    """The 10 unordered pairs of 5 structures."""
    return list(itertools.combinations(structures, 2))
```

- [ ] **Step 4: Run test to verify it passes**

Run: `python3 -m pytest test_score_seedtest.py -v`
Expected: PASS (5 passed)

- [ ] **Step 5: Add the driver that walks the archives and calls DockQ**

Extend `score_seedtest.py` with a `main()` that, per backend and per chunk archive:
extract to `$TMPDIR` (node-local, bounded), group CIFs by `design_id`, compute DockQ for each of `all_pairs(...)` via the existing `dockq_self.py` code path, read the per-sample confidence JSONs, and emit one row per design × backend through `pair_stats` + `spread_stats` for `iptm`/`ptm`/`plddt`.

Run it as a CPU job (`normal.4h`, `--mem-per-cpu`), NOT on the login node — the login node kills long single-threaded jobs.

- [ ] **Step 6: Verify the output shape**

```bash
ssh euler 'cd /cluster/work/wicky/lvalsecchi/rna_binder/seedtest_2026-07-30 && \
  module load stack/2024-06 python/3.12.8 && python3 -c "
import pandas as pd
d=pd.read_csv(\"seedtest_results.csv\")
print(\"rows:\",len(d),\"(expect 6000 = 2000 x 3)\")
print(d.groupby([\"backend\"]).size())
print(\"n_structures_found value counts:\"); print(d.n_structures_found.value_counts())
print(\"designs with all 5:\", (d.n_structures_found==5).mean()*100, \"%\")"'
```

Expected: 6,000 rows, 2,000 per backend, the overwhelming majority with `n_structures_found == 5`.

---

### Task 7: Analysis + report

**Files:**
- Create: `$RB/seedtest_2026-07-30/REPORT.md`
- Create: `$RB/seedtest_2026-07-30/scripts/analyse_seedtest.py`

**Interfaces:**
- Consumes: `seedtest_results.csv`, plus the 100k table joined on `design_id`.
- Produces: the four answers listed in the spec.

- [ ] **Step 1: Compute the headline numbers**

Per backend and per stack: median `dockq_pair_mean`, median `iptm_std`, and the fraction of designs whose 5 seeds agree well (`dockq_pair_min > 0.8`) versus disagree badly (`dockq_pair_min < 0.23`).

- [ ] **Step 2: Answer "does seed-stability track quality?"**

Join to the 100k table on `design_id` and report Spearman correlation between `dockq_pair_mean` (stability) and the existing `consensus_iptm` / `dockq_consensus` (quality), overall and per stack. State N.

- [ ] **Step 3: Quantify how much consensus-metric spread is noise**

Compare the within-design seed spread (`iptm_std`) against the between-design spread of `consensus_iptm` in the 100k. If within-design ≈ between-design, the metric is mostly noise at the individual-design level — an important caveat for any ranking built on it.

- [ ] **Step 4: Write `REPORT.md`**

Lead with the single number per backend (median pairwise DockQ across seeds), then the per-stack table, the quality-correlation result, the noise-floor result, and a Caveats section (predicted structures; 5 draws is a small sample per design; MSA off).

- [ ] **Step 5: Verify every number in the report reproduces**

Re-run the analysis script and diff its printed values against the numbers written in `REPORT.md`. Any mismatch is a bug in the report, not the script.

---

## Self-Review

**Spec coverage:** subsample → Task 1; input building → Task 2; per-backend 5-draw config and the no-production-edits constraint → Task 3; probe → Task 4; chunked inode-safe execution → Task 5; pairwise DockQ + confidence spread → Task 6; the four analysis questions → Task 7. No spec section is unaddressed.

**Placeholder scan:** the two genuinely unverified items are called out explicitly rather than hidden — the ESM multi-sample attribute (Task 3 Step 5, with a concrete fallback) and the AF3 JSON key set (Task 2 Step 5, with the exact command to check). Chunk size and walltime are deliberately deferred to the probe measurement rather than guessed.

**Type consistency:** `build_manifest` emits exactly the five columns Task 2 reads; `write_b2_yaml` is reused for both B2 and ESM (matching how `esmfold2_infer.py` parses input); `pair_stats`/`spread_stats` key names match the `seedtest_results.csv` column list in Task 6's Interfaces block, which Task 7 consumes.
