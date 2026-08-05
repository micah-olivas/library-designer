"""QC specific to tiled-assembly libraries, run in addition to the standard checks.

Each assembled oligo must carry exactly the two intended Golden Gate sites and no
extras: adding the flanking primers, sites, and spacers can
create an unintended recognition site spanning a junction, which would misdirect
digestion. We also confirm every oligo fits the budget, and hand the fused overhangs to
``checks/overhangs.py``, which scores how much homology every pair of them shares (so the
fragment inserts directionally and the cut vector cannot re-close on itself).
"""
from __future__ import annotations

from .motifs import ENZYME_SITES  # noqa: F401  (kept for discoverability)


def check_tiled(library) -> dict:
    """Return tiled-specific findings as lists of offending variant/tile labels."""
    from ..layout.tiled import extra_sites

    spec = library.spec
    # The params the library was laid out with, not spec.tiled: an explicit
    # tile(params) must be checked against the layout it actually produced.
    params = library.tiled_params
    df = library.df
    tiles = library.tiles
    ref = library.reference

    extra: list[str] = []
    over_budget: list[str] = []
    for name, oligo, ti in zip(df["name"], df["oligo"], df["tile"]):
        if not isinstance(oligo, str):
            continue
        t = tiles[int(ti)]
        # lead/trail include the pad, so the two intended sites are where they belong even
        # when the pool was evened out to one length.
        if extra_sites(oligo, len(t.lead), len(t.trail), params.enzyme):
            extra.append(str(name))
        if len(oligo) > params.oligo_budget:
            over_budget.append(str(name))

    # Overhangs come straight from the resolved params, whose terminal contexts
    # tile_library already derived from the backbone when a starting vector is set.
    # checks/overhangs.py reads them with the same helper the oligos were built with, so QC
    # cannot drift off the layout, and it scores every pair rather than only exact matches.
    from .overhangs import overhang_findings

    overhang_issues, overhang_advisories = overhang_findings(library)

    # Unplaced non-WT variants would signal a codon split across tile boundaries
    # (codon-aligned tiling prevents this). The global WT row is unplaced by design, it has
    # no position; the per-tile WT_Tile_<i> controls are placed and checked like any oligo.
    unplaced = [n for n in getattr(library, "_unplaced", []) if n != "WT"]

    # Stray enzyme sites on the assembled destination vector (computed at tile time).
    # In use_vector_cds mode these are flagged here rather than raised, so the user who
    # chose to keep their own CDS still sees them.
    vector_extra_sites = list(getattr(library, "_vector_extra_sites", []))

    # Advisories: things kept verbatim from a chosen reference (SD sites, avoid-motifs,
    # a native BsaI). Informational, they do not fail the report, because the user opted
    # into this sequence rather than letting us recode it.
    from .report import verbatim_advisories

    vec = spec.resolve_vector(params)
    verbatim = spec.cds is not None or (vec is not None and vec.use_vector_cds)
    advisories = verbatim_advisories(spec, ref) if verbatim and ref else []

    return {
        "oligo_extra_sites": extra,
        "oligo_over_budget": over_budget,
        "overhang_issues": overhang_issues,
        "overhang_advisories": overhang_advisories,
        "unplaced": unplaced,
        "vector_extra_sites": vector_extra_sites,
        "reference_advisories": advisories,
    }
