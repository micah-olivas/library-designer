"""The translation round-trip check, run on the variable region alone because the adaptors
are not coding.
"""
from __future__ import annotations

from dnachisel import translate


def translates_to(variable_dna: str, expected_protein: str) -> bool:
    """True if the variable-region DNA translates to the intended protein
    (an internal stop codon translates to ``'*'``, matching amber variants)."""
    return translate(variable_dna) == expected_protein
