#!/usr/bin/env python3
"""Render T2 from its CSV.

T3-T6 generate their .md beside their .csv; T2 is a transcription with no
generator, so this renders it the same way. Provenance stays out of the table
file: it lives in T2_findings.md, which is not part of this archive.
"""
import csv
from pathlib import Path

OUT = Path(__file__).resolve().parent

TABLES = {
    "T2_primers": dict(
        title="T2 - primers",
        intro=["Every primer used in this work, as ordered, 5'->3'. Oligo-pool "
               "amplification, read-out RT-PCR, Sanger verification.",
               "Where a primer anneals on the minus strand the sequence given is the "
               "oligonucleotide, not the plus-strand site it reads; the site is in the "
               "position column. Recovery of each sequence is recorded in "
               "`T2_findings.md`."],
        heads=["Primer", "Sequence", "Position", "Use", "Source"]),
}


def render(slug, spec):
    rows = list(csv.reader((OUT / f"{slug}.csv").open()))
    cols, data = rows[0], rows[1:]
    assert len(spec["heads"]) == len(cols), (slug, len(spec["heads"]), len(cols))
    md = [f"# {spec['title']}", ""]
    for p in spec["intro"]:
        md += [p, ""]
    md += ["| " + " | ".join(spec["heads"]) + " |",
           "|" + "---|" * len(cols)]
    for r in data:
        md.append("| " + " | ".join(c.replace("|", r"\|").replace("\n", " ") or "-"
                                   for c in r) + " |")
    md += ["", f"{len(data)} rows."]
    (OUT / f"{slug}.md").write_text("\n".join(md) + "\n")
    return len(data)


if __name__ == "__main__":
    for slug, spec in TABLES.items():
        print(f"wrote {slug}.md ({render(slug, spec)} rows)")
