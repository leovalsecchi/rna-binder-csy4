#!/usr/bin/env python3
"""Render T1 and T2 from their CSVs.

T3-T6 generate their .md beside their .csv; T1 and T2 are transcriptions with no
generator, so this renders them the same way. Provenance stays out: it lives in
T1_findings.md and T2_findings.md.
"""
import csv
from pathlib import Path

OUT = Path(__file__).resolve().parent

TABLES = {
    "T1_construct_inventory": dict(
        title="T1 - construct inventory",
        intro=["Every plasmid built or used in this work: the entry and scaffold "
               "plasmids, the screening plasmids, the 32 ordered candidates and the "
               "helper plasmids used but not built here. Ordered synthesis fragments "
               "are **not** listed; the plasmids assembled from them are.",
               "`size_bp` is the full plasmid length throughout. Every value is "
               "measured from the NAS DNA record, never from a journal label. Scope, "
               "per-construct provenance and the contradictions found while building "
               "it are in `T1_findings.md`."],
        heads=["pART", "pLV", "Insert", "PaqCI", "BsmBI", "BsaI", "BbsI", "Size (bp)"]),
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
