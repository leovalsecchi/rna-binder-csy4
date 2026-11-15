#!/usr/bin/env python3
"""Reapply HD local patches onto the 1b08312 worktree (3 upstream-rewritten files).
Strict: asserts every precondition; writes NOTHING unless ALL edits succeed.
Drops debug probes (PRE-PREDICT PROBE, _metrics_probe) by omission; keeps traceback.print_exc."""
import sys, os

UP = "/cluster/work/wicky/lvalsecchi/halludesign_upgrade"

def read(rel):
    with open(os.path.join(UP, rel)) as f:
        return f.readlines()

def patch_models(lines):
    # --- patch2 site 293 (idx 292): af3_op, 8-space indent ---
    i = 292
    assert 'copied_file = metrics[0]["packed_path"]' in lines[i], f"293 content: {lines[i]!r}"
    assert lines[i][8] == 'c', f"293 indent (want 8 sp): {lines[i]!r}"
    lines[i] = lines[i].replace('metrics[0]["packed_path"]', 'source_file')
    # --- patch2 site 849 (idx 848): cross_model Protenix branch, 16-space, near ProtenixInferrer ---
    i = 848
    assert 'copied_file = metrics[0]["packed_path"]' in lines[i], f"849 content: {lines[i]!r}"
    assert lines[i][16] == 'c', f"849 indent (want 16 sp): {lines[i]!r}"
    assert any('ProtenixInferrer' in lines[j] for j in range(848, 855)), "849 not near ProtenixInferrer"
    lines[i] = lines[i].replace('metrics[0]["packed_path"]', 'source_file')
    # --- patch3 cif_path 1054 (idx 1053): first tag.replace -> metrics[0]["file_name"] ---
    i = 1053
    assert 'tag.replace(".pdb", ""), tag.replace(".pdb", "")' in lines[i], f"1054: {lines[i]!r}"
    lines[i] = lines[i].replace('tag.replace(".pdb", "")', 'metrics[0]["file_name"]', 1)
    # --- regex fix (unique): _recycle_.*_model -> _recycle_.* ---
    idx = [k for k, l in enumerate(lines) if 're.sub(r"_recycle_.*_model"' in l]
    assert len(idx) == 1, f"regex sites: {idx}"
    lines[idx[0]] = lines[idx[0]].replace('r"_recycle_.*_model"', 'r"_recycle_.*"')
    # --- traceback (unique, keep): prepend traceback.print_exc to the file-error print ---
    idx = [k for k, l in enumerate(lines) if 'print(f"file error: {str(e)}")' in l]
    assert len(idx) == 1, f"traceback sites: {idx}"
    lines[idx[0]] = lines[idx[0]].replace('print(f"file error',
                                          'import traceback; traceback.print_exc(); print(f"file error')
    return "".join(lines)

def patch_eval(lines):
    # --- patch1 line 917 (idx 916): op_ ptm .item() -> guarded ---
    i = 916
    assert "metrics[f'op_{label}_ptm']" in lines[i] and ".item()" in lines[i], f"917: {lines[i]!r}"
    lines[i] = lines[i].replace("dict_result['chain_ptm'][i].item()",
                                "dict_result['chain_ptm'][i] if i < len(dict_result['chain_ptm']) else None")
    # --- patch1 line 1065 (idx 1064): op_ ptm [i] (no .item) -> guarded ---
    i = 1064
    assert "metrics[f'op_{label}_ptm']" in lines[i] and ".item()" not in lines[i], f"1065: {lines[i]!r}"
    lines[i] = lines[i].replace("dict_result['chain_ptm'][i]",
                                "dict_result['chain_ptm'][i] if i < len(dict_result['chain_ptm']) else None")
    return "".join(lines)

def patch_run(lines):
    # --- patch6a: insert --mpnn_num_packs after the temperature help line ---
    idx = [k for k, l in enumerate(lines) if "help='mpnn temperature to use')" in l]
    assert len(idx) == 1, f"temp-help sites: {idx}"
    ins = ["    parser.add_argument('--mpnn_num_packs', type=int, default=1,\n",
           "                       help='number_of_packs_per_design for the MPNN side-chain packer (sidechain pack count per design)')\n"]
    lines[idx[0] + 1:idx[0] + 1] = ins
    # --- patch6b: thread args.mpnn_num_packs into all 5 MPNNModel blocks ---
    content = "".join(lines)
    n = content.count("number_of_packs_per_design=1,")
    assert n == 5, f"packs sites (want 5): {n}"
    content = content.replace("number_of_packs_per_design=1,", "number_of_packs_per_design=args.mpnn_num_packs,")
    return content

jobs = [("models_utility.py", patch_models),
        ("eval/evaluation.py", patch_eval),
        ("HalluDesign_run.py", patch_run)]

results = []
for rel, fn in jobs:
    try:
        results.append((rel, fn(read(rel))))
    except Exception as e:
        print(f"FAIL {rel}: {e}")
        sys.exit(1)

for rel, content in results:
    with open(os.path.join(UP, rel), "w") as f:
        f.write(content)
    print(f"WROTE {rel}")
print("ALL PATCHES APPLIED OK")
