"""Read the lab's own sequence files: SnapGene `.dna` and GenBank `.gb`.

Three supplementary figures (`scaffold_design`, `scaffold_selection`,
`negative_controls`) re-derive their numbers from the construct files on the
NAS rather than from any written record, so the reader lives here once.

Biopython is NOT in `~/.venvs/binder_map`, and a SnapGene file is a short
tag-length-value stream, so the parser is six lines rather than a dependency:
each packet is one type byte, a four-byte big-endian length, then the payload;
type 0 is the DNA sequence, whose first payload byte is a flags byte.
"""
import re
import struct
from pathlib import Path

# The Csy4 extended repeat as it appears in the construct files (DNA alphabet).
CSY4 = "GTTCACTGCCGTATAGGCAGCTAAGAAA"
# What the negative controls put in its place at every hairpin position.
NEG_LOOP = "ATTTCGATCAATGGATTAT"


def _snapgene(path):
    b = Path(path).read_bytes()
    i, seq = 0, None
    while i + 5 <= len(b):
        t = b[i]
        n = struct.unpack(">I", b[i + 1:i + 5])[0]
        if t == 0:
            seq = b[i + 6:i + 5 + n].decode("ascii", "replace")
        i += 5 + n
    if seq is None:
        raise ValueError(f"no DNA packet in {path}")
    return seq


def _genbank(path):
    out, on = [], False
    for ln in Path(path).read_text(errors="replace").splitlines():
        if ln.startswith("ORIGIN"):
            on = True
        elif ln.startswith("//"):
            on = False
        elif on:
            out.append(re.sub(r"[^A-Za-z]", "", ln))
    return "".join(out)


def read(path):
    """Sequence of a `.dna`, `.gb`/`.gbk` or `.fa`/`.fasta` file, upper case."""
    p = Path(path)
    if p.suffix == ".dna":
        s = _snapgene(p)
    elif p.suffix in (".gb", ".gbk"):
        s = _genbank(p)
    else:
        s = "".join(l.strip() for l in p.read_text(errors="replace").splitlines()
                    if not l.startswith(">"))
    return s.upper()


def find_all(seq, motif):
    """Every start index of `motif`, overlaps included."""
    out, i = [], seq.find(motif)
    while i >= 0:
        out.append(i)
        i = seq.find(motif, i + 1)
    return out


def gc(seq):
    """GC percent."""
    return 100.0 * sum(c in "GC" for c in seq) / len(seq)


def _selfcheck():
    """The two motifs are what the construct files actually contain.

    Checked against the four pLV00x files, which are the authoritative sizes
    (journal cross-check A4: the pART bp figures in the journal are working
    notes against intermediates, the NAS DNA files are authoritative).
    """
    assert len(CSY4) == 28 and len(NEG_LOOP) == 19
    assert round(gc(NEG_LOOP), 1) == 26.3
    assert round(gc(CSY4), 1) == 46.4


_selfcheck()
