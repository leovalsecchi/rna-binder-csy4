#!/usr/bin/env python3
"""versions.tsv -> Supplementary Table S1 (markdown), written next to it as S1.md.

Methods cites S1 from three places, so it is generated rather than retyped.
Rows whose status is not_recoverable are kept and marked; an omitted row would
read as an absent dependency rather than an unrecoverable one.
"""
import csv
import pathlib

HERE = pathlib.Path(__file__).parent
COLS = [
    ("role", "Stage"), ("tool", "Tool"), ("version", "Version"),
    ("commit", "Commit"), ("container", "Container image"),
    ("image_date", "Image date"), ("weights_params", "Weights / parameters"),
    ("config", "Configuration as run"),
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
    out = ["# Supplementary Table S1 — software, containers and parameter sets", ""]
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
    (HERE / "S1.md").write_text("\n".join(out))
    print(f"S1.md: {len(rs)} rows, {len(missing)} not recoverable")


if __name__ == "__main__":
    main()
