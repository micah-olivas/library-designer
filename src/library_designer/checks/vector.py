"""QC for the destination vector of a standard (untiled) library.

Three things have to line up for a one-oligo Golden Gate library to clone. Each adaptor
must carry exactly one recognition site, pointing at the insert. The fused overhang the
digested oligo keeps must be the overhang the cut vector presents, otherwise nothing
ligates (or it ligates the wrong way round). And the assembled destination vector must
carry only the two intended sites, since the whole plasmid sees the enzyme.

The first two are read out of the adaptors and the located insert (see
``layout/destination.py``); the third is counted on the emitted vector, the same check
tiling runs per tile.
"""
from __future__ import annotations


def check_vector(library) -> dict:
    """Destination-vector findings for a standard library, as lists of readable strings.

    ``adaptor_issues`` covers the adaptors and their fit to the plasmid, ``overhang_issues``
    the two fused overhangs, ``vector_extra_sites`` stray recognition sites on the emitted
    vector. Those three fail the report; ``advisories`` are informational.

    A library with no adaptors at all is not being cloned from this oligo pool directly, so
    its missing cut sites are reported as an advisory instead of a failure, and the vector's
    overhangs are drawn from the backbone flanking the insert.
    """
    from ..layout.destination import build_destination
    from ..layout.tiled import vector_site_positions

    spec = library.spec
    vec = spec.vector
    enzyme = vec.enzyme
    no_adaptors = not spec.adaptor_5 and not spec.adaptor_3

    try:
        # strict=False: an overhang collision is reported below as an overhang issue rather
        # than swallowing the rest of the checks.
        dv = build_destination(library, strict=False)
    except ValueError as exc:
        # A geometry or locating problem, reported rather than raised so check() still
        # returns a full report the user can act on.
        return {
            "adaptor_issues": [f"destination vector: {exc}"],
            "overhang_issues": [],
            "overhang_advisories": [],
            "vector_extra_sites": [],
            "advisories": [],
        }

    from .report import verbatim_advisories

    advisories: list[str] = []
    if spec.cds is not None or vec.use_vector_cds:
        advisories += verbatim_advisories(spec, library.reference)
    adaptor_issues = list(dv.issues)
    if no_adaptors:
        advisories.append(
            f"no adaptors are set, so nothing carries the {enzyme} sites that would release "
            f"the insert; the destination vector drops the whole CDS and takes both fused "
            f"overhangs ({dv.overhang_5}, {dv.overhang_3}) from the backbone"
        )
        adaptor_issues = []
    cut = dv.cut
    if cut is not None:
        # The load-bearing check: the oligo's overhangs against the vector's. They differ
        # exactly when the adaptor bases that ride along past the cut are not the backbone
        # bases flanking the insert, so name both sequences.
        for end, oligo_ovh, vector_ovh in (
            ("5'", cut.overhang_5, dv.overhang_5),
            ("3'", cut.overhang_3, dv.overhang_3),
        ):
            if oligo_ovh != vector_ovh:
                adaptor_issues.append(
                    f"the {end} fused overhang the oligo carries ({oligo_ovh}) is not the one "
                    f"the cut destination vector presents ({vector_ovh}), so the insert will "
                    f"not ligate. Match the adaptor to the plasmid flanking the insert."
                )

    # The two fused overhangs against each other, scored on how much homology they share
    # rather than on exact identity alone (see checks/overhangs.py).
    from .overhangs import overhang_findings

    overhang_issues, overhang_advisories = overhang_findings(library)

    pos = vector_site_positions(dv.sequence, enzyme, dv.topology == "circular")
    vector_extra_sites = (
        [] if len(pos) == 2
        else [f"destination vector: {len(pos)} {enzyme} site(s) (expected 2) at {pos}"]
    )
    return {
        "adaptor_issues": adaptor_issues,
        "overhang_issues": overhang_issues,
        "overhang_advisories": overhang_advisories,
        "vector_extra_sites": vector_extra_sites,
        "advisories": advisories,
    }
