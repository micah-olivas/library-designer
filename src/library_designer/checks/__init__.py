"""QC. ``check_library`` runs every check that applies to a library and returns a ``CheckReport``.

The per-topic checks live in the sibling modules (``motifs``, ``translation``, ``tiled``,
``vector``, ``overhangs``, ``mispriming``, ``assembly``) and are called from ``report``
rather than by hand.
"""
from .report import CheckReport, check_library

__all__ = ["CheckReport", "check_library"]
