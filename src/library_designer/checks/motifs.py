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


def contains_motif(sequence: str, motif: str, both_strands: bool = True) -> bool:
    seq, motif = sequence.upper(), motif.upper()
    if motif in seq:
        return True
    return both_strands and reverse_complement(motif) in seq


def contains_enzyme_site(sequence: str, enzyme: str) -> bool:
    try:
        site = ENZYME_SITES[enzyme]
    except KeyError as exc:
        raise KeyError(
            f"Unknown enzyme {enzyme!r}; known: {sorted(ENZYME_SITES)}"
        ) from exc
    return contains_motif(sequence, site, both_strands=True)


def count_enzyme_sites(sequence: str, enzyme: str) -> int:
    """Number of recognition sites for ``enzyme`` on both strands (non-overlapping)."""
    site = ENZYME_SITES[enzyme]
    seq = sequence.upper()
    rc = reverse_complement(site)
    return seq.count(site) + (seq.count(rc) if rc != site else 0)
