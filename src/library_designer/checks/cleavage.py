"""How much sequence sits outside each Type IIS recognition site.

A site flush against the end of a molecule cuts poorly. The enzyme binds duplex DNA either
side of its recognition sequence, so a site with nothing 5' of it has less to hold than it
wants and the digest is inefficient rather than wrong: the fragment is still correct when it
is released, there is just less of it. On a pooled oligo the cost lands on yield, and an
under-cut oligo contributes nothing to the assembly.

The fix is a few bases of lead-in. A tiled oligo gets it for free, since its amplification
primer sits 5' of the site. A standard library's adaptors do not, and an adaptor written as
``GGTCTC`` + spacer + overhang puts the site at base 1 of every oligo in the pool.

This is reported as an advisory rather than a failure. The design is makeable and some
protocols run it, so the call belongs to whoever is ordering.
"""
from __future__ import annotations

from ..regions import reverse_complement
from .motifs import ENZYME_SITES, recognition_site

# Bases of lead-in to ask for outside each site. Published Golden Gate protocols put a handful
# of bases there rather than none, and NEB's own adaptor designs carry several; the exact
# figure is a judgement, so this is the threshold for saying something, not a claim that 4 is
# sufficient and 3 is not.
MIN_FLANK = 4


def site_flanks(molecule: str, enzyme: str) -> tuple[int | None, int | None]:
    """``(lead, trail)`` for ``molecule``: bases 5' of the outermost forward-strand site, and
    bases 3' of the outermost reverse-strand site.

    ``None`` for an end whose site is missing, which is a different finding and reported by
    ``cut_construct`` rather than here. The outermost site is the one that matters, since it
    is the one nearest the end that has to cut.
    """
    seq = molecule.upper()
    site = recognition_site(enzyme).upper()
    site_rc = reverse_complement(site)
    fwd = seq.find(site)
    rev = seq.rfind(site_rc)
    lead = fwd if fwd != -1 else None
    trail = (len(seq) - rev - len(site_rc)) if rev != -1 else None
    return lead, trail


def design_enzyme(library) -> str | None:
    """The Type IIS enzyme this library's flanks carry.

    Taken from the design where it is stated, the tiled params or the starting vector. With
    neither, the adaptors are searched for a site the package knows, so a library with no
    vector is still checked. ``None`` when nothing carries a site.
    """
    spec = library.spec
    params = getattr(library, "tiled_params", None)
    if params is not None:
        return params.enzyme
    vec = spec.resolve_vector(params)
    if vec is not None:
        return vec.enzyme
    flanks = (spec.adaptor_5 + spec.adaptor_3).upper()
    ordered = list(spec.avoid_enzymes) + [e for e in ENZYME_SITES if e not in spec.avoid_enzymes]
    for enzyme in ordered:
        site = recognition_site(enzyme).upper()
        if site in flanks or reverse_complement(site) in flanks:
            return enzyme
    return None


def cleavage_advisories(library) -> list[str]:
    """Ends whose recognition site sits closer than ``MIN_FLANK`` to the end of the molecule.

    Reported per flank rather than per member: a standard library's adaptors are the same on
    every oligo, so naming all of them would say one thing hundreds of times. A tiled library
    is reported per tile, since each carries its own primers.
    """
    enzyme = design_enzyme(library)
    if enzyme is None:
        return []

    out: list[str] = []
    tiles = getattr(library, "tiles", None)
    if tiles is not None:
        df = library.df
        oligos = {int(t): o for t, o in zip(df.get("tile", []), df.get("oligo", []))
                  if isinstance(o, str)}
        for tile in tiles:
            oligo = oligos.get(tile.index)
            if oligo:
                out += _messages(site_flanks(oligo, enzyme), enzyme,
                                 f"tile{tile.index}'s oligo", "the forward primer",
                                 "the reverse primer")
        return out

    spec = library.spec
    if not (spec.adaptor_5 or spec.adaptor_3):
        return []
    from ..regions import assemble

    dna = next((d for d in library.df.get("variable_dna", []) if isinstance(d, str)), "")
    if not dna:
        return []
    return _messages(site_flanks(assemble(spec.adaptor_5, dna, spec.adaptor_3), enzyme),
                     enzyme, "every oligo in the pool", "adaptor_5", "adaptor_3")


def _messages(flanks, enzyme: str, what: str, five: str, three: str) -> list[str]:
    """One sentence per end that is short, naming the end, the count, and where to add bases."""
    lead, trail = flanks
    out = []
    for n, end, where in ((lead, "5'", five), (trail, "3'", three)):
        if n is not None and n < MIN_FLANK:
            out.append(
                f"the {end} {enzyme} site sits {n} base(s) from the end of {what}, so it cuts "
                f"less efficiently than one with {MIN_FLANK} or more; add lead-in bases to "
                f"{where}"
            )
    return out
