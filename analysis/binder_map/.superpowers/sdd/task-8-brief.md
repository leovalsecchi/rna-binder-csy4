### Task 8: Scaling probe, then both runs

The ~2,100 CPU-h stage-3 figure is an N² extrapolation from the v6 run
(31,763 nodes → 8.5 CPU-h), not a measurement. Measure it before committing a
128-task array.

**Files:**
- Create: `probe_scaling.sh`
- Create: `submit_maps.sh`

**Interfaces:**
- Consumes: everything above
- Produces: nothing importable; these are drivers.

- [ ] **Step 1: Write the probe**

Create `probe_scaling.sh`:

```bash
#!/bin/bash
# Measure how Foldseek search actually scales on this node set before
# committing the full run. Builds two sub-DBs from the already-built stage-2
# database and times one query chunk against each.
set -uo pipefail
FM=/cluster/work/wicky/lvalsecchi/rna_binder/analyses/foldseek_map
RESULTS=${RESULTS:-results_full}
cd "$FM" || exit 1

for N in 50000 150000; do
  sbatch --job-name="probe_$N" --partition=normal.4h --time=04:00:00 \
    --cpus-per-task=16 --mem-per-cpu=4G \
    --output="$FM/logs/probe_${N}_%j.log" \
    --wrap "bash -lc '
      set -e
      module load stack/2024-06 gcc/12.2.0 openmpi/4.1.6 foldseek/8-ef4e960
      OUT=$FM/$RESULTS/design
      WORK=\${TMPDIR:-/tmp}/probe_$N
      mkdir -p \$WORK; trap \"rm -rf \$WORK\" EXIT
      head -n $N \$OUT/db/uniq.lookup | cut -f1 > \$WORK/keys.txt
      head -n 2000 \$WORK/keys.txt > \$WORK/qkeys.txt
      for suffix in \"\" _ss _ca; do
        foldseek createsubdb \$WORK/keys.txt  \$OUT/db/uniq\$suffix \$WORK/sub\$suffix
        foldseek createsubdb \$WORK/qkeys.txt \$OUT/db/uniq\$suffix \$WORK/q\$suffix
      done
      echo \"=== N=$N : 2000 queries against a $N-target DB ===\"
      /usr/bin/time -v foldseek search \$WORK/q \$WORK/sub \$WORK/aln \$WORK/tmp \
        -c 0.8 --cov-mode 0 -e 10 --max-seqs 300 --alignment-type 1 -a --threads 16
    '"
done
echo "Two probes submitted. When both finish:"
echo "  grep -H 'Elapsed (wall clock)' $FM/logs/probe_*.log"
echo "p = log(t150/t50)/log(3). Full map-1 search ~ t50 * (500000/50000)^p * (500000/2000)."
```

- [ ] **Step 2: Write the submit driver**

Create `submit_maps.sh`:

```bash
#!/bin/bash
# Three-map build: one shared 683 GB extraction, split per map, then stages 2-5
# per map, then the comparison. Map 2 (results_sw_k30_f50_mut) already exists
# and is never rebuilt.
# Parameters are frozen to match results_sw_k30_f50_mut exactly.
set -uo pipefail
FM=/cluster/work/wicky/lvalsecchi/rna_binder/analyses/foldseek_map
PY=/cluster/work/wicky/lvalsecchi/binder_map_venv/bin/python
cd "$FM" || exit 1
LAYER=design

# --- stage 1: one task per outer archive, over the UNION manifest -------------
J1=$(sbatch --parsable --job-name=cwfm_s1 --array=0-34%35 \
     --export=ALL,N_SHARDS=35,BUILD=build_shared,RESULTS=results_shared \
     slurm/stage1_extract_nested.sbatch) || exit 1
echo "stage1 extract (shared)  $J1"

# --- stage 1b: split into the two per-map trees ------------------------------
J1B=$(sbatch --parsable --job-name=cwfm_s1b --dependency=afterok:$J1 \
      --partition=normal.4h --time=02:00:00 --cpus-per-task=4 --mem-per-cpu=4G \
      --output="$FM/logs/cwfm_s1b_%j.log" \
      --wrap "bash -lc 'cd $FM && PYTHONPATH=src $PY scripts/stage1b_split.py \
        --shared $FM/results_shared \
        --map full:build_full:results_full \
        --map lib12k:build_lib12k:results_lib12k'") || exit 1
echo "stage1b split            $J1B"

# --- stages 2-5, once per map ------------------------------------------------
submit_map () {
  local TAG=$1 RESULTS=$2 CHUNKS=$3 THROTTLE=$4 DEP=$5
  local J2 J3 J4 J5
  J2=$(sbatch --parsable --job-name=${TAG}_s2 --dependency=afterok:$DEP \
       --wrap "bash -lc 'cd $FM && RESULTS=$RESULTS bash scripts/stage2_builddb.sh $LAYER && \
         PYTHONPATH=src $PY scripts/stage3_split_queries.py --layer $LAYER \
         --n-chunks $CHUNKS --results $FM/$RESULTS'" \
       --partition=normal.4h --time=04:00:00 --cpus-per-task=8 --mem-per-cpu=8G \
       --output="$FM/logs/${TAG}_s2_%j.log") || exit 1
  J3=$(sbatch --parsable --job-name=${TAG}_s3 --dependency=afterok:$J2 \
       --array=0-$((CHUNKS-1))%${THROTTLE} \
       --export=ALL,ALN_TYPE=1,RESULTS=$RESULTS \
       slurm/stage3_search.sbatch $LAYER) || exit 1
  J4=$(sbatch --parsable --job-name=${TAG}_s4 --dependency=afterok:$J3 \
       --partition=normal.4h --time=04:00:00 --cpus-per-task=16 --mem-per-cpu=8G \
       --output="$FM/logs/${TAG}_s4_%j.log" \
       --wrap "bash -lc 'cd $FM && PYTHONPATH=src $PY scripts/stage4_graph.py \
         --layer $LAYER --tm-floor 0.5 --epi-floor 0.6 --k 30 --mutual \
         --results $FM/$RESULTS'") || exit 1
  J5=$(sbatch --parsable --job-name=${TAG}_s5 --dependency=afterok:$J4 \
       slurm/stage_single.sbatch scripts/stage5_cluster.py \
         --layer $LAYER --axis fold --results $FM/$RESULTS) || exit 1
  echo "$TAG: s2=$J2 s3=$J3 s4=$J4 s5=$J5" >&2
  echo "$J5"
}

J5_LIB=$(submit_map lib12k results_lib12k 4   4  $J1B)
J5_FULL=$(submit_map full  results_full  128 64 $J1B)

# --- stage 9: compare all three ---------------------------------------------
J9=$(sbatch --parsable --job-name=cwfm_s9 \
     --dependency=afterok:$J5_FULL:$J5_LIB \
     --partition=normal.4h --time=01:00:00 --cpus-per-task=8 --mem-per-cpu=4G \
     --output="$FM/logs/cwfm_s9_%j.log" \
     --wrap "bash -lc 'cd $FM && PYTHONPATH=src $PY scripts/stage9_compare.py'") || exit 1
echo "stage9 compare           $J9"
echo "$J9" > $FM/logs/last_job_cwfm.txt
```

- [ ] **Step 3: Run stage 0, then launch**

```bash
cd $FM && PYTHONPATH=src:scripts $PY scripts/stage0_sample.py --n 500000 --seed 42
cd $FM && bash submit_maps.sh
```

Expected: seven job ids printed, `map3_n` in `build_shared/meta/sample.json`
exactly 9,759.

- [ ] **Step 4: Gate on hash identity before paying for the search**

When `cwfm_s1b` completes — **before** `full_s3` starts — verify the repack
rewrite did not change what gets extracted:

```bash
cd $FM && PYTHONPATH=src $PY - <<'EOF'
import glob
import pandas as pd
FM = "/cluster/work/wicky/lvalsecchi/rna_binder/analyses/foldseek_map"

def hashes(results):
    files = sorted(glob.glob(f"{FM}/{results}/design/hashes/shard_*.tsv"))
    return pd.concat([pd.read_csv(f, sep="\t") for f in files]).set_index("design_id")["coord_hash"]

new, old = hashes("results_shared"), hashes("results_v6")
shared = new.index.intersection(old.index)
print(f"design_ids in both runs: {len(shared):,d}")
assert len(shared) > 5000, "too few shared ids to be a real check"
mismatch = int((new.loc[shared] != old.loc[shared]).sum())
print(f"coord_hash mismatches: {mismatch}")
assert mismatch == 0, "extraction changed -- STOP, do not run the search"
print("PASS: extraction is byte-identical to the v6 run")
EOF
```

Expected: `coord_hash mismatches: 0` and `PASS`. If it fails, `scancel` the
`full_s3`/`s4`/`s5` ids printed by `submit_maps.sh` before ~2,000 CPU-h is spent.

- [ ] **Step 4b: Check the shard member names match the manifest ids**

Stage 0 writes the manifests; stage 1 writes the PDB member names; stage 1b
matches one against the other by stripping `.pdb`. If those two ever disagree on
id format by a single character, **every shard produces a well-formed, EMPTY
output tar and nothing raises** — the map would come out valid and empty. Check
it before the search, not after:

```bash
cd $FM && PYTHONPATH=src $PY - <<'EOF'
import glob, tarfile
import pandas as pd
FM = "/cluster/work/wicky/lvalsecchi/rna_binder/analyses/foldseek_map"

man = set(pd.read_csv(f"{FM}/build_shared/manifests/design.tsv",
                      sep="\t", usecols=["design_id"])["design_id"])
shard = sorted(glob.glob(f"{FM}/results_shared/design/shards/shard_*.tar.gz"))[0]
with tarfile.open(shard, "r:gz") as tar:
    stems = [m.name[:-4] for m in tar if m.isfile() and m.name.endswith(".pdb")]
print(f"shard {shard.split('/')[-1]}: {len(stems):,d} members")
assert stems, "shard has no .pdb members at all"
hit = sum(s in man for s in stems)
print(f"members whose stem is a manifest design_id: {hit:,d} / {len(stems):,d}")
print("examples:", stems[:2])
assert hit == len(stems), "ID FORMAT MISMATCH -- stage 1b would emit empty tars"

# and the split itself must not be silently empty
for tree in ("results_full", "results_lib12k"):
    tars = sorted(glob.glob(f"{FM}/{tree}/design/shards/shard_*.tar.gz"))
    assert tars, f"{tree}: no shards written"
    with tarfile.open(tars[0], "r:gz") as tar:
        n = sum(1 for m in tar if m.isfile())
    print(f"{tree}: first shard holds {n:,d} members")
    assert n > 0, f"{tree}: split produced an EMPTY tar"
print("PASS")
EOF
```

Expected: every stem is a manifest `design_id`, and both split trees have
non-empty shards. A mismatch here is the one failure that produces a complete,
plausible, empty result — fix it before spending the search.

- [ ] **Step 5: Check the library map is a complete census**

```bash
cd $FM && $PY -c "
import glob, pandas as pd
FM='/cluster/work/wicky/lvalsecchi/rna_binder/analyses/foldseek_map'
n=pd.read_csv(f'{FM}/results_lib12k/design/dedup/nodes.tsv',sep='\t')
miss=pd.concat([pd.read_csv(f,sep='\t') for f in
                sorted(glob.glob(f'{FM}/results_lib12k/design/missing/shard_*.tsv'))])
print('library nodes:', len(n), 'missing:', len(miss))
if len(miss): print(miss.reason.value_counts())
assert len(n) + 0 >= 9700, 'library census is short -- investigate before comparing'
print('PASS')"
```

Expected: `library nodes: 9759` (or a shortfall fully explained by `missing/`)
and `PASS`. Map 3 is a census, not a sample — a silent shortfall here biases
every comparison that follows.

- [ ] **Step 6: Spot-check TM identity against the v6 hits**

```bash
cd $FM && PYTHONPATH=src $PY - <<'EOF'
import pandas as pd
FM = "/cluster/work/wicky/lvalsecchi/rna_binder/analyses/foldseek_map"
COLS = ["query","target","qtmscore","ttmscore","alntmscore",
        "lddt","fident","alnlen","evalue","bits"]
read = lambda f: pd.read_csv(f, sep="\t", names=COLS, header=None,
                             usecols=["query","target","qtmscore","ttmscore"])
new = pd.concat(read(f"{FM}/results_lib12k/design/hits/hits_0{i}.tsv") for i in range(4))
old = pd.concat(read(f"{FM}/results_v6/design/hits/hits_{i:02d}.tsv") for i in range(8))
m = new.merge(old, on=["query","target"], suffixes=("_new","_old"))
print(f"shared query/target pairs: {len(m):,d}")
assert len(m) > 100, "too few shared pairs -- widen the chunks compared"
d = (m["qtmscore_new"] - m["qtmscore_old"]).abs()
print(f"max |dTM| = {d.max():.6f}")
assert d.max() < 1e-3, "Foldseek scores drifted -- STOP"
print("PASS: TM scores reproduce the v6 run")
EOF
```

Expected: `max |dTM|` at or near 0 and `PASS`. The library map is used here
because it finishes first and its node set overlaps map 2 almost entirely.

- [ ] **Step 7: Read the comparison**

When `cwfm_s9` completes:

```bash
cd $FM && $PY -c "
import pandas as pd
s=pd.read_parquet('compare/map_summary.parquet'); print(s.to_string(index=False))
c=pd.read_parquet('compare/crosswalk.parquet');  print(c.to_string(index=False))
assert len(s)==3, 'expected all three maps'
print('PASS')"
```

Expected: three rows — `sample500k`, `passing60k`, `library12k` — with node
counts 500,000 / 31,763 / 9,759, and three pairwise ARI/AMI rows.

- [ ] **Step 8: Commit**

```bash
cd ~/Documents/binder_map
for f in probe_scaling.sh submit_maps.sh; do
  done
git add probe_scaling.sh submit_maps.sh
git commit -m "Add scaling probe and the three-map submit driver

One shared extraction feeds both new maps; stages 2-5 then run once per map with
frozen parameters. The chain gates on a coord-hash identity check against
results_v6 before the 500k search is paid for."
```

---

## Notes for whoever executes this

- **Stages 2, 3 and 5 are not modified.** If a change to them looks necessary,
  stop and re-read the spec — the clustering parameters are the one thing this
  work must not perturb, because comparability with map 2 is the whole point.
- **Map 2 is never rebuilt.** `results_sw_k30_f50_mut` is the reference.
- **`map3_n` must be exactly 9,760**, which stage 2 then dedups to 9,759 graph
  nodes. These are different stages; comparing them is an error the spec made once.
  The library-selection export has its header on **line 1** (no comment preamble),
  unlike the binder_filter exports which have two — count leading `#` lines rather
  than hardcoding `skiprows`.
- **Watch `missing/` after stage 1.** A large `member_missing` count means the
  manifest's nested paths disagree with the archives; check `nested_ref` against
  a real `.list` sidecar before blaming the data.
- **Do not "fix" the dedup key.** It over-splits BG+BIF ~9.8× (rigid-body offset,
  0.94 Å superposed RMSD) and HD ~10× (0.022 Å packer jitter). Map 2 has the same
  defect; matching it is what makes the three maps comparable.
- **Cluster ids are not comparable across maps.** Three independent Leiden runs on
  three node sets. Compare size distributions and composition; the crosswalk's
  ARI/AMI is a result to report, not a correspondence to rely on.
