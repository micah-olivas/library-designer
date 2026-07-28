"""QC. ``check_library`` runs everything the design supports and returns a ``CheckReport``.

The per-topic checks live in the sibling modules (``motifs``, ``translation``, ``tiled``,
``vector``, ``assembly``) and are called from ``report`` rather than by hand.
"""
from .report import CheckReport, check_library

__all__ = ["CheckReport", "check_library"]
