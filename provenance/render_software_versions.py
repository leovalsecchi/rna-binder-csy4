#!/usr/bin/env python3
"""versions.tsv -> the software and container table, written beside it as
software_versions.md.

The README and the Methods cite that table from several places, so it is generated
rather than retyped. It was called S1.md until 2026-11-23, after an earlier
supplement numbering; the supplement's own container inventory is a narrower cut.
Rows whose status is not_recoverable are kept and marked; an omitted row would
read as an absent dependency rather than an unrecoverable one.
"""
import csv
import pathlib

HERE = pathlib.Path(__file__).parent
COLS = [
    ("role", "Stage"), ("tool", "Tool"), ("version", "Version"),
    ("commit", "Commit"), ("container", "Container image"),
    ("image_date", "Image date"), ("sha256", "Image sha256"),
    ("weights_params", "Weights / parameters"), ("config", "Configuration as run"),
]


def rows():
    with (HERE / "versions.tsv").open() as fh:
        for row in csv.DictReader((l for l in fh if not l.startswith("#")), delimiter="\t"):
            yield row


def cell(row, key):
    v = (row.get(key) or "").strip()
    if v in ("", "n/a"):
        return "--"
    if v == "not_recoverable":
        return "*not recoverable*"
    return v.replace("|", r"\|")


def main():
    rs = list(rows())
    assert rs, "versions.tsv parsed to zero rows"
    out = ["# Software, containers and parameter sets", "",
           "One row per pipeline stage: tool, version, commit, container image and its",
           "sha256 digest, the weights it loaded and the configuration it ran under. The",
           "supplement's container inventory (Supplementary Table S4) is a narrower,",
           "per-image cut of the same material; this file is the fuller per-stage version",
           "and is the one the archive cites. It was called `S1.md` until 2026-11-23,",
           "after an earlier supplement numbering.", ""]
    out.append("| " + " | ".join(h for _, h in COLS) + " |")
    out.append("|" + "|".join(["---"] * len(COLS)) + "|")
    for r in rs:
        out.append("| " + " | ".join(cell(r, k) for k, _ in COLS) + " |")

    missing = [r for r in rs if r["status"] == "not_recoverable"]
    if missing:
        out += ["", "## Components with no recoverable version", ""]
        for r in missing:
            out.append(f"- **{r['tool']}** — {r['evidence']}")

    out += ["", "Image dates are filesystem modification times: they record when an image was",
            "placed in the container store, not when the upstream release was cut.", ""]
    (HERE / "software_versions.md").write_text("\n".join(out))
    print(f"software_versions.md: {len(rs)} rows, {len(missing)} not recoverable")


if __name__ == "__main__":
    main()
