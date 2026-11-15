"""seq_complexity.py — cheap sequence low-complexity metrics for the binder sequence.

No structure needed — computed from the protein sequence string alone, so the whole
set scores in seconds. Flags repetitive / poly-X sequences (e.g. a LigandMPNN
sequence collapse, which shows up as one residue dominating the chain).

    seq_max_aa_frac  : largest single-amino-acid share of the sequence (0-1).
                       High = poly-X (one residue dominates). Lower = better.
    seq_longest_run  : longest run of the same amino acid in a row.
                       High = a homopolymer stretch. Lower = better.
    seq_max_kmer     : length of the longest substring (k-mer) that occurs at least
                       twice in the sequence (classic "longest repeated substring",
                       occurrences may overlap). Catches PERIODIC / dispersed
                       multi-residue repeats (e.g. ELELEL, KTKTKT) that the single-AA
                       metrics above miss. A random 80-mer over 20 AAs sits ~3-4;
                       a repetitive design runs much longer. Lower = better. Filter
                       out designs whose largest repeated k-mer exceeds some k.

Public entry point: ``seq_complexity(seq) -> dict`` (never raises; empty/None input
yields None values). ``SEQ_COLUMNS`` is the column schema.
"""
from __future__ import annotations

from collections import Counter

SEQ_COLUMNS = ("seq_max_aa_frac", "seq_longest_run", "seq_max_kmer")


def _longest_repeated_kmer(s: str) -> int:
    """Length of the longest substring occurring >= 2 times in ``s`` (overlaps
    allowed). Binary search on length: a repeat of length L implies a repeat of
    every length < L (its prefix), so 'a duplicate of length L exists' is monotone
    in L. O(n log n) with a per-length substring set — trivial for protein lengths."""
    n = len(s)
    if n < 2:
        return 0

    def has_dup(length: int) -> bool:
        seen = set()
        for i in range(n - length + 1):
            sub = s[i:i + length]
            if sub in seen:
                return True
            seen.add(sub)
        return False

    lo, hi, best = 1, n - 1, 0
    while lo <= hi:
        mid = (lo + hi) // 2
        if has_dup(mid):
            best, lo = mid, mid + 1
        else:
            hi = mid - 1
    return best


def seq_complexity(seq) -> dict:
    """Low-complexity metrics for one protein sequence string."""
    none = {"seq_max_aa_frac": None, "seq_longest_run": None, "seq_max_kmer": None}
    if not isinstance(seq, str):
        return dict(none)
    s = "".join(c for c in seq.strip().upper() if c.isalpha())
    if not s:
        return dict(none)
    max_aa = max(Counter(s).values()) / len(s)
    longest = run = 1
    for i in range(1, len(s)):
        run = run + 1 if s[i] == s[i - 1] else 1
        if run > longest:
            longest = run
    return {
        "seq_max_aa_frac": round(max_aa, 4),
        "seq_longest_run": longest,
        "seq_max_kmer": _longest_repeated_kmer(s),
    }
