"""Restriction-enzyme recognition sites, Type IIS cut geometry, and the sequence lookups
that read them.

Two tables carry the data. ``ENZYME_SITES`` maps an enzyme name to its recognition
sequence, and ``ENZYME_CUTS`` maps it to where it cuts, which is what the tiling layout
needs to place a fused overhang. Every lookup goes through ``recognition_site`` or
``cut_geometry`` so an enzyme the package does not know fails with a readable error
instead of a bare ``KeyError`` partway through a run.
"""
from __future__ import annotations

from ..regions import reverse_complement

# Recognition sequences for the restriction enzymes we support checking.
ENZYME_SITES: dict[str, str] = {
    "BsaI": "GGTCTC",
    "BsmBI": "CGTCTC",
    "Esp3I": "CGTCTC",
    "BbsI": "GAAGAC",
    "SapI": "GCTCTTC",
}



# Where each Type IIS enzyme cuts, as (spacer, overhang_len): how many bases sit between
# the recognition site and the top-strand cut, and how long the 5' overhang it leaves is.
# NEB writes BsaI as GGTCTC(1/5), meaning one base of spacer and then a 4-base overhang.
# So an oligo whose 5' adaptor reads ...GGTCTC N NNNN keeps NNNN as its fused overhang.
ENZYME_CUTS: dict[str, tuple[int, int]] = {
    "BsaI": (1, 4),
    "BsmBI": (1, 4),
    "Esp3I": (1, 4),
    "BbsI": (2, 4),
    "SapI": (1, 3),
}


def cut_geometry(enzyme: str) -> tuple[int, int]:
    """``(spacer, overhang_len)`` for ``enzyme``. Raises for one we don't know how to cut,
    since guessing the geometry would put the fused overhang in the wrong place."""
    try:
        return ENZYME_CUTS[enzyme]
    except KeyError as exc:
        raise KeyError(
            f"No cut geometry known for {enzyme!r}; known: {sorted(ENZYME_CUTS)}"
        ) from exc


def contains_motif(sequence: str, motif: str, both_strands: bool = True) -> bool:
    """True if ``motif`` occurs in ``sequence``, on either strand by default.

    Plain substring matching on the uppercased sequence, not a regex, so this is for
    literal motifs. A motif on the bottom strand is found by searching the top strand for
    its reverse complement. Pass ``both_strands=False`` to look at the given strand only.

    This is the enzyme-site path, reached through ``contains_enzyme_site``. The regex
    motifs in ``spec.avoid_patterns`` do not come through here; ``check_library`` runs those
    with ``re`` on the top strand only.
    """
    seq, motif = sequence.upper(), motif.upper()
    if motif in seq:
        return True
    return both_strands and reverse_complement(motif) in seq


def recognition_site(enzyme: str) -> str:
    """The recognition sequence for ``enzyme``, or a readable error naming the ones we know.

    Every lookup goes through here. Indexing ``ENZYME_SITES`` directly used to let an
    unsupported name in ``avoid_enzymes`` get all the way past reference optimization before
    dying on a bare ``KeyError`` while stamping the first variant."""
    try:
        return ENZYME_SITES[enzyme]
    except KeyError as exc:
        raise KeyError(
            f"Unknown enzyme {enzyme!r}; known: {sorted(ENZYME_SITES)}. Names are "
            "case-sensitive."
        ) from exc


def contains_enzyme_site(sequence: str, enzyme: str) -> bool:
    """True if ``sequence`` carries a recognition site for ``enzyme`` on either strand.

    A Type IIS recognition sequence is not palindromic, so the bottom strand has to be
    searched too. An enzyme bound there cuts the sequence just the same. Raises for an
    enzyme with no entry in ``ENZYME_SITES``. Use ``count_enzyme_sites`` when the number of
    sites matters, as it does for the baseline-relative checks in ``report.py``.
    """
    return contains_motif(sequence, recognition_site(enzyme), both_strands=True)


def count_enzyme_sites(sequence: str, enzyme: str) -> int:
    """Number of recognition sites for ``enzyme`` on both strands (non-overlapping)."""
    site = recognition_site(enzyme)
    seq = sequence.upper()
    rc = reverse_complement(site)
    return seq.count(site) + (seq.count(rc) if rc != site else 0)
