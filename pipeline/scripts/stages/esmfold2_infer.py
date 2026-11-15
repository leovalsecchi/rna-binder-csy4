#!/usr/bin/env python3
"""ESMFold2 batch refold worker.

Folds each Boltz-style YAML input (protein + RNA) with the local biohub/ESMFold2
model and writes a canonical `<design>_ESM_sample0.cif` plus a
`<design>_ESM_sample0_confidence.json` sidecar (pLDDT / pTM / ipTM) into the
canonical designs dir.

The model is loaded once and reused across the whole batch (load is ~20s, so a
per-design relaunch would dominate runtime -- this mirrors the Boltz-2 worker,
not the AF3 per-design loop). Resume-safe: a design whose canonical .cif already
exists is skipped, so a requeued/preempted job continues where it stopped.

Invoked inside esmc.sif (see refold_esm.sh). Run env supplies the HF cache.
"""
import argparse
import json
import sys
import traceback
from pathlib import Path

import yaml

try:
    import numpy as np
except Exception:  # numpy ships in esmc.sif; degrade gracefully if ever absent
    np = None

from esm.models.esmfold2 import (
    ESMFold2InputBuilder,
    ProteinInput,
    RNAInput,
    StructurePredictionInput,
)
from transformers.models.esmfold2.modeling_esmfold2 import ESMFold2Model

# Optional entity types (present in this esm build, but guard for portability).
try:
    from esm.models.esmfold2 import DNAInput
except Exception:  # pragma: no cover
    DNAInput = None

TAG = "ESM"


def _to_numpy(t):
    """Torch tensor (any device) or array-like -> float32 numpy, batch dim squeezed."""
    if hasattr(t, "detach"):
        t = t.detach().cpu()
    arr = np.asarray(t, dtype="float32")
    if arr.ndim >= 3 and arr.shape[0] == 1:  # drop a leading singleton sample dim
        arr = arr[0]
    return arr


def parse_yaml_to_spi(path):
    """Boltz-style YAML -> StructurePredictionInput (protein + rna [+ dna])."""
    data = yaml.safe_load(path.read_text())
    if not data or "sequences" not in data:
        raise ValueError("no 'sequences' key")
    seqs = []
    for record in data["sequences"]:
        (etype, ent), = record.items()
        etype = etype.lower()
        ids = ent.get("id")
        ids = ids if isinstance(ids, (list, tuple)) else [ids]
        for chain_id in ids:
            if etype == "protein":
                seqs.append(ProteinInput(id=chain_id, sequence=ent["sequence"]))
            elif etype == "rna":
                seqs.append(RNAInput(id=chain_id, sequence=ent["sequence"]))
            elif etype == "dna" and DNAInput is not None:
                seqs.append(DNAInput(id=chain_id, sequence=ent["sequence"]))
            else:
                raise ValueError(f"unsupported entity type {etype!r}")
    if not seqs:
        raise ValueError("parsed zero sequences")
    return StructurePredictionInput(sequences=seqs)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("inputs_dir", help="dir of *.yaml inputs")
    ap.add_argument("canonical_dir", help="dir to write <name>_ESM_sample0.cif into")
    ap.add_argument("--num-loops", type=int, default=20)
    ap.add_argument("--num-sampling-steps", type=int, default=100)
    ap.add_argument("--num-diffusion-samples", type=int, default=1)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--model", default="biohub/ESMFold2")
    ap.add_argument("--pae-dir", default=None,
                    help="if set, dump per-model pae_<stem>.npz / plddt_<stem>.npz "
                         "here (Boltz layout) for downstream IPSAE scoring")
    args = ap.parse_args()

    inputs = sorted(p for p in Path(args.inputs_dir).iterdir()
                    if p.suffix.lower() in {".yaml", ".yml"})
    cdir = Path(args.canonical_dir)
    cdir.mkdir(parents=True, exist_ok=True)
    pae_dir = Path(args.pae_dir) if args.pae_dir else None
    if pae_dir is not None:
        pae_dir.mkdir(parents=True, exist_ok=True)
    if not inputs:
        print(f"ERROR: no YAML inputs in {args.inputs_dir}", file=sys.stderr)
        sys.exit(1)

    print(f"Loading {args.model} ...", flush=True)
    model = ESMFold2Model.from_pretrained(args.model).cuda().eval()
    builder = ESMFold2InputBuilder()

    n_ok = n_fail = n_skip = 0
    for y in inputs:
        name = y.stem
        out_cif = cdir / f"{name}_{TAG}_sample0.cif"
        if out_cif.exists():
            n_skip += 1
            print(f"[SKIP] {name} (already folded)", flush=True)
            continue
        try:
            spi = parse_yaml_to_spi(y)
            result = builder.fold(
                model, spi,
                num_loops=args.num_loops,
                num_sampling_steps=args.num_sampling_steps,
                num_diffusion_samples=args.num_diffusion_samples,
                seed=args.seed,
            )
            plddt = float(result.plddt.mean())
            ptm = float(result.ptm)
            iptm = float(result.iptm)
            # Chain-pair iPTM matrix (chain A=0 protein, B=1 RNA). Mirrors B2's
            # pair_chains_iptm so the aggregator can derive esm_protein_rna_iptm.
            pci = getattr(result, "pair_chains_iptm", None)
            if pci is not None and hasattr(pci, "tolist"):
                pci = pci.tolist()
            # Write to a temp name then rename so a preemption mid-write never
            # leaves a half-written .cif that the resume-skip would treat as done.
            tmp_cif = cdir / f".{name}_{TAG}_sample0.cif.tmp"
            tmp_cif.write_text(result.complex.to_mmcif())
            tmp_cif.rename(out_cif)
            conf = {"iptm": iptm, "ptm": ptm, "plddt": plddt, "complex_plddt": plddt}
            if pci is not None:
                conf["pair_chains_iptm"] = pci
            (cdir / f"{name}_{TAG}_sample0_confidence.json").write_text(json.dumps(conf))
            # Dump the full PAE (+ per-residue pLDDT) in Boltz-npz layout so the
            # finalize step can score this model with IPSAE. Non-fatal on failure.
            if pae_dir is not None and np is not None:
                pae = getattr(result, "pae", None)
                if pae is not None:
                    stem = f"{name}_{TAG}_sample0"
                    try:
                        np.savez(pae_dir / f"pae_{stem}.npz", pae=_to_numpy(pae))
                        np.savez(pae_dir / f"plddt_{stem}.npz", plddt=_to_numpy(result.plddt))
                    except Exception as exc:
                        print(f"[warn] {name}: PAE dump failed: {exc}", flush=True)
            n_ok += 1
            print(f"[OK]   {name}: pLDDT {plddt:.3f} pTM {ptm:.3f} ipTM {iptm:.3f}",
                  flush=True)
        except Exception as exc:  # one bad design shouldn't kill the batch
            n_fail += 1
            print(f"[FAIL] {name}: {exc}", flush=True)
            traceback.print_exc()

    print(f"Done: {n_ok} ok, {n_fail} failed, {n_skip} skipped of {len(inputs)}",
          flush=True)


if __name__ == "__main__":
    main()
