"""Codon optimization. ``codon_optimize`` reverse-translates a protein under the spec's
sequence rules, and ``backbone`` freezes one WT reference that way and stamps each
variant's codon onto it.
"""
from . import constraints
from .codon import codon_optimize

__all__ = ["codon_optimize", "constraints"]
