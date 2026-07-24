"""QC specific to tiled-assembly libraries, run in addition to the standard checks.

The load-bearing check is that each assembled oligo carries *exactly* the two intended
Golden Gate sites and no extras: adding the flanking primers, sites, and spacers can
create an unintended recognition site spanning a junction, which would misdirect
digestion. We also confirm every oligo fits the budget and that each tile's two fused
overhangs are distinct and non-palindromic (so the fragment inserts directionally and
cannot self-ligate).
"""
from __future__ import annotations

import re

from ..regions import reverse_complement
from .motifs import ENZYME_SITES, count_enzyme_sites  # noqa: F401  (ENZYME_SITES kept for discoverability)


def check_tiled(library) -> dict:
    """Return tiled-specific findings as lists of offending variant/tile labels."""
    from ..layout.tiled import extra_sites, tile_contexts

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
        if extra_sites(oligo, len(t.fwd), len(t.rev), params.enzyme):
            extra.append(str(name))
        if len(oligo) > params.oligo_budget:
            over_budget.append(str(name))

    # Overhangs come straight from the resolved params, whose terminal contexts
    # tile_library already derived from the backbone when a starting vector is set. Using
    # the same helper the oligos were built with keeps QC from drifting off the layout.
    overhang_issues: list[str] = []
    for t in tiles:
        c5, c3 = tile_contexts(ref, t.start, t.end, params)
        if (
            c5 == c3                          # 5' and 3' overhangs identical, so not directional
            or c5 == reverse_complement(c5)    # palindromic, so it self-ligates
            or c3 == reverse_complement(c3)
            or c5 == reverse_complement(c3)    # complementary, so the fragment can flip
        ):
            overhang_issues.append(f"tile{t.index}")

    # Unplaced non-WT variants would signal a codon split across tile boundaries
    # (codon-aligned tiling prevents this); the WT control is unplaced by design.
    unplaced = [n for n in getattr(library, "_unplaced", []) if n != "WT"]

    # Stray enzyme sites on the assembled destination vector (computed at tile time).
    # In use_vector_cds mode these are flagged here rather than raised, so the user who
    # chose to keep their own CDS still sees them.
    vector_extra_sites = list(getattr(library, "_vector_extra_sites", []))

    # Advisories: things kept verbatim from a chosen reference (SD sites, avoid-motifs,
    # a native BsaI). Informational, they do not fail the report, because the user opted
    # into this sequence rather than letting us recode it.
    advisories: list[str] = []
    verbatim = spec.cds is not None or (params is not None and params.use_vector_cds)
    if verbatim and ref:
        for e in spec.avoid_enzymes:
            n = count_enzyme_sites(ref, e)
            if n:
                advisories.append(
                    f"reference carries {n} internal {e} site(s), kept verbatim "
                    "(a Golden Gate hazard; see the destination-vector check)"
                )
        for p in spec.avoid_patterns:
            n = len(re.findall(p, ref))
            if n:
                advisories.append(f"reference carries {n} /{p}/ motif(s), kept verbatim (not recoded)")

    return {
        "oligo_extra_sites": extra,
        "oligo_over_budget": over_budget,
        "overhang_issues": overhang_issues,
        "unplaced": unplaced,
        "vector_extra_sites": vector_extra_sites,
        "reference_advisories": advisories,
    }
