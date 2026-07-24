"""Region model.

The construct's parts, a 5' adaptor, the variable (coding) region, and a 3'
adaptor, are kept as **explicit columns** on the library, so region boundaries are
never inferred from character case. The single place case carries meaning is the
uSort-M wire format (`usortm_sequence`), applied only at that export boundary
because uSort-M strips flanking by case.
"""
from __future__ import annotations

# Full IUPAC complement (upper- and lowercase). A superset of plain ACGT/N, so
# shipped inputs (enzyme sites, screened primers, WT-derived overhangs) are
# unaffected, while degenerate custom motifs/primers complement correctly instead
# of silently passing through unchanged.
_COMPLEMENT = str.maketrans(
    "ACGTRYSWKMBDHVNacgtryswkmbdhvn",
    "TGCAYRSWMKVHDBNtgcayrswmkvhdbn",
)


def assemble(adaptor_5: str, variable: str, adaptor_3: str) -> str:
    """Full construct as plain uppercase DNA (no case-encoded regions)."""
    return f"{adaptor_5}{variable}{adaptor_3}".upper()


def usortm_sequence(adaptor_5: str, variable: str, adaptor_3: str) -> str:
    """uSort-M wire format: flanking lowercase, variable UPPERCASE. Used only when
    writing ``variants.csv``, derived from the explicit region columns, not parsed
    back, since uSort-M distinguishes the variable region from flanking by case."""
    return f"{adaptor_5.lower()}{variable.upper()}{adaptor_3.lower()}"


def reverse_complement(seq: str) -> str:
    """Reverse complement, case-preserving, over the full IUPAC alphabet. Characters
    outside that alphabet pass through unchanged."""
    return seq.translate(_COMPLEMENT)[::-1]
