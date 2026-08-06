"""Reusable DNA Chisel constraint builders, one place to define sequence rules,
so the optimizer does not keep two near-identical constraint lists.
"""
from __future__ import annotations

from dnachisel import AvoidPattern, EnforceSequence, EnforceTranslation

from ..spec import LibrarySpec


def sequence_rules(spec: LibrarySpec) -> list:
    """Location-independent constraints every optimized CDS must satisfy: enforce translation,
    avoid forbidden motifs, avoid restriction sites, and keep single-base runs within
    ``spec.max_homopolymer`` when one is set.

    ``AvoidPattern("<enzyme>_site")`` is enzyme-aware and checks both strands.
    """
    rules: list = [EnforceTranslation()]
    rules += [AvoidPattern(pattern) for pattern in spec.avoid_patterns]
    rules += [AvoidPattern(f"{enzyme}_site") for enzyme in spec.avoid_enzymes]
    # One per base, so a run longer than the limit cannot survive optimization.
    rules += [AvoidPattern(pattern) for pattern in spec.homopolymer_patterns]
    return rules


def protect_codon(dna: str, aa_index: int) -> EnforceSequence:
    """Lock the codon at a given amino-acid index (e.g. a hand-placed stop)."""
    start = aa_index * 3
    return EnforceSequence(sequence=dna[start:start + 3], location=(start, start + 3))
