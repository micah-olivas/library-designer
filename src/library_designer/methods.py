"""Synthesis-platform helpers.

library_designer emits two order-form shapes: a *pooled* oligo order (one shared pool,
Twist oligo-pool style) and an *arrayed* order (one row per named construct, IDT eBlock
/ gene-fragment style). ``spec.platform`` selects between them, either the bare type
``"pooled"`` / ``"arrayed"`` or one of the provider slugs below that map to a shape.

The package ships no synthesis-limit registry: set the length gate yourself with
``spec.max_oligo_length`` (or, for tiled assembly, ``tiled.oligo_budget``).
"""
from __future__ import annotations

__all__ = ["platform_type", "PLATFORM_TYPES"]

POOLED = "pooled"
ARRAYED = "arrayed"

# Provider slugs mapped to their order-form shape. The two bare types map to
# themselves. Extend this with your own vendor slugs as needed; an unknown platform
# raises so a typo cannot silently pick the wrong order form.
PLATFORM_TYPES: dict[str, str] = {
    POOLED: POOLED,
    ARRAYED: ARRAYED,
    "twist_oligo_pools": POOLED,
    "twist_gene_fragments": ARRAYED,
    "idt_oligo_pools": POOLED,
    "idt_eblocks": ARRAYED,
    "agilent_oligo_pools": POOLED,
}


def platform_type(platform: str) -> str:
    """Resolve a platform to its order-form shape, ``"pooled"`` or ``"arrayed"``."""
    try:
        return PLATFORM_TYPES[platform]
    except KeyError:
        raise ValueError(
            f"Unknown platform {platform!r}. Use 'pooled' or 'arrayed', or one of "
            f"{sorted(PLATFORM_TYPES)}; or pass method= to to_vendor() directly."
        ) from None
