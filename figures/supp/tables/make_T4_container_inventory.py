#!/usr/bin/env python3
"""T4 - container and version inventory. Writes T4_container_inventory.{csv,md}.

Run from the repo root:
    ~/.venvs/binder_map/bin/python3 code/supp/tables/make_T4_container_inventory.py

Digests are NEVER retyped: they are read from sources/container_digests.txt, so a
transcription slip is impossible. The tool/version column is the hand-curated
part; every value in VERSIONS below was read out of the image itself with the
command recorded in the `probe` field, and all nine sha256 digests were
independently re-computed on Euler (see README.md, T4).
"""
import csv
import re
from pathlib import Path

SRC = Path(__file__).resolve().parent / "sources"
OUT = Path(__file__).resolve().parent

UNRECORDED = "not recorded - the image carries no version string"

# Keyed by image basename. `probe` is the command whose output `version` came
# from, run as: ssh euler 'apptainer exec <image> bash -lc "<probe>"'
VERSIONS = {
    "foundry_0.1.9.sif": dict(
        tool="rc-foundry (carries RFdiffusion3; Institute for Protein Design)",
        version="0.1.9",
        probe="cat /app/foundry/.venv/lib/python*/site-packages/"
              "rc_foundry-*.dist-info/METADATA",
        extra="torch 2.7.1. The upstream git commit cannot be recovered: the "
              "image carries no git metadata and no `git` binary, so this "
              "digest is its only identifier. Configs reference this image via "
              "the symlink /cluster/work/wicky/containers/foundry.sif.",
        used_by="RFdiffusion3+LigandMPNN (generation)"),
    "boltzgen.sif": dict(
        tool="BoltzGen",
        version="0.3.2",
        probe="pip list | grep -i boltzgen",
        extra="git commit 31d9d9b9c72245b4ed6fe8742d6fbf4e1a3552a0 "
              "(`cd /app && git log -1 --format=%H`); torch 2.12.0.",
        used_by="BoltzGen+LigandMPNN and BoltzGen+BoltzIF (shared generation)"),
    "boltz_design1.sif": dict(
        tool="BoltzDesign1",
        version=UNRECORDED + " (no version string in the image)",
        probe="cd /app && git log -1 --format=%H",
        extra="git commit 31aeedd7f187ba88dd9e92514070aee13381cb93 is recorded "
              "and is the only identifier this image carries.",
        used_by="BoltzDesign1+LigandMPNN (generation)"),
    "LigandMPNN.sif": dict(
        tool="LigandMPNN",
        version=UNRECORDED + " (no version string in the image)",
        probe="cd /app/LigandMPNN && git log -1 --format=%H",
        extra="git commit 091ab1ff5fb4d13854cf6a7c41ec531e1d9d3e67; torch 2.2.1. "
              "Checkpoint ligandmpnn_v_32_010_25.pt, loaded from outside the "
              "image (/cluster/work/wicky/params/foundry/), so it is not covered "
              "by this digest.",
        used_by="inverse folding for RFdiffusion3, BoltzGen and BoltzDesign1 "
                "routes (NOT HalluDesign - see note)"),
    "alphafold3.sif": dict(
        tool="AlphaFold3",
        version="3.0.0",
        probe="pip list | grep -i alphafold3",
        extra="jax 0.4.34. The weights blob "
              "/cluster/work/wicky/params/af3_params/af3.bin sits outside the "
              "image and carries no version string, so it cannot be "
              "versioned beyond this path.",
        used_by="refold backend 1 (all five routes)"),
    "boltz2.sif": dict(
        tool="Boltz-2",
        version="2.2.1",
        probe="pip list | grep -iE '^boltz '",
        extra="torch 2.8.0.",
        used_by="refold backend 2 (all five routes)"),
    "esmc.sif": dict(
        tool="ESMFold2 (ESMC-6B backbone)",
        version="esm 3.3.0",
        probe="pip list | grep -iE '^esm |^transformers '",
        extra="transformers 4.57.6, torch 2.6.0+cu124. Hugging Face weights "
              "cache staged outside the image "
              "(rna_binder/params/esmfold2_hf_cache); its checkpoint revision "
              "is not recorded. Image built from a docker-archive tar that "
              "no longer "
              "exists, so there is no build recipe.",
        used_by="refold backend 3 (all five routes)"),
    "pyrosetta.sif": dict(
        tool="PyRosetta",
        version="2025.37+release.df75a9c48e",
        probe="pip list | grep -i pyrosetta",
        extra="Not referenced anywhere in code/configs/ - it is an analysis "
              "image (interface H-bond and packing metrics), not a design or "
              "refold image.",
        used_by="post-hoc interface metrics (not in the design pipeline)"),
    "af3_protenix_esm.sif": dict(
        tool="HalluDesign predictor image (AlphaFold3 + Protenix + ESM)",
        version="alphafold3 3.0.0, fair-esm 2.0.0",
        probe="pip list | grep -iE 'protenix|alphafold|esm|^torch '",
        extra="torch 2.6.0+cu126. The Protenix version and checkpoint cannot "
              "be recovered: there is no pip package and no version string; "
              "the weights ship inside the "
              "image, so this digest is the only identifier for them. "
              "HalluDesign's own source is NOT in this image - see note.",
        used_by="HalluDesign+LigandMPNN (generation)"),
}


def main():
    text = (SRC / "container_digests.txt").read_text()
    entries = re.findall(r"^([0-9a-f]{64})\s+(\S+)$", text, re.M)
    # VERIFY the reported count of nine sha256 digests.
    assert len(entries) == 9, f"expected 9 digests, found {len(entries)}"
    assert len({d for d, _ in entries}) == 9, "duplicate digest"
    assert set(Path(p).name for _, p in entries) == set(VERSIONS), \
        "image set in the digest file does not match VERSIONS"

    rows = []
    for digest, path in entries:
        v = VERSIONS[Path(path).name]
        rows.append(dict(
            image=path,
            sha256_digest=digest,
            tool=v["tool"],
            version=v["version"],
            used_by=v["used_by"],
            version_probe=v["probe"],
            notes=v["extra"],
        ))

    cols = list(rows[0])
    with (OUT / "T4_container_inventory.csv").open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols)
        w.writeheader()
        w.writerows(rows)

    hdr = ["Image", "sha256 digest", "Tool", "Version", "Used by"]
    md = ["# T4 - container and version inventory", "",
          "Generated by `code/supp/tables/make_T4_container_inventory.py` from "
          "`sources/container_digests.txt`; digests are copied programmatically, "
          "never retyped. All nine were independently re-computed on Euler and "
          "reproduce exactly.", "",
          "These are `sha256sum` digests **of the `.sif` files on disk**, not OCI "
          "image digests.", "",
          "| " + " | ".join(hdr) + " |",
          "|" + "|".join("---" for _ in hdr) + "|"]
    for r in rows:
        md.append("| `{}` | `{}` | {} | {} | {} |".format(
            r["image"], r["sha256_digest"], r["tool"], r["version"], r["used_by"]))
    md += ["", "## Per-image notes", ""]
    md += [f"- **`{Path(r['image']).name}`** - {r['notes']}" for r in rows]
    md += ["", "See `README.md` for the two coverage gaps this inventory does "
               "not close (HalluDesign's own source tree, and the externally "
               "mounted model weights)."]
    (OUT / "T4_container_inventory.md").write_text("\n".join(md) + "\n")
    print(f"wrote T4_container_inventory.csv / .md  "
          f"({len(rows)} images, VERIFY: count is nine)")
    unrec = [Path(r["image"]).name for r in rows if UNRECORDED in r["version"]]
    print(f"  versions not recorded in the image: {unrec or 'none'}")


if __name__ == "__main__":
    main()
