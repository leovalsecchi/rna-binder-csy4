"""Shared constants for the binder map pipeline.

Values are pinned to match rna_binder.analysis.structural so that recomputed
per-nucleotide contacts can be cross-checked against the stored
rna_contact_nt_* columns.
"""

# The csy4 target RNA, 16 nt. This is the canonical epitope frame: every
# fingerprint position refers to an index into this sequence.
RNA_CANONICAL = "ACUGCCGUAUAGGCAG"
RNA_LEN = 16

# Angstroms, heavy-atom. Must equal rna_binder's RNA_CONTACT_CUTOFF or the
# cross-check in Task 19 is meaningless.
RNA_CONTACT_CUTOFF = 4.5

# Every stochastic step (Leiden, UMAP) uses this.
SEED = 42

# Moieties the per-nucleotide contact counts are split across.
MOIETIES = ("phosphate", "sugar", "base")

# 16 nucleotides x 3 moieties
FINGERPRINT_DIM = RNA_LEN * len(MOIETIES)

assert len(RNA_CANONICAL) == RNA_LEN
