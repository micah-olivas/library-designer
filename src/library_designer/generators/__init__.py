"""Generators, the objects that turn a spec into a variant table.

Each one takes a ``LibrarySpec`` and hands back a ``Library`` from ``generate()``. The
generator you pick is what sets how the members relate to each other, and so how they are
codon-optimized. ``SubstitutionScan`` stamps single mutants onto one frozen reference,
while ``SequenceSet`` optimizes every member on its own.
"""
from .base import Generator
from .scanning import SubstitutionScan
from .sequence_set import SequenceSet

__all__ = ["Generator", "SubstitutionScan", "SequenceSet"]
