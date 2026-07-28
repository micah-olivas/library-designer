"""The destination vector for a standard (untiled) library.

A tiled library needs one destination vector per tile (see ``layout/tiled.py``). A
standard library carries the whole CDS on one oligo, so it needs exactly one: the
starting plasmid with the CDS replaced by a Golden Gate drop-out.

Where that drop-out starts and ends is decided by the adaptors, not by us. The oligo's
Type IIS sites sit in the adaptors, so digesting the oligo leaves a fixed pair of fused
overhangs, and the cut vector has to present the same two. ``cut_construct`` reads those
overhangs out of the adaptors, then ``build_destination`` places the drop-out so the
vector presents them. Two conventions both work and both turn up in real designs:

- the adaptor ends at the recognition site plus its spacer, so the overhang is the first
  (last) four bases of the CDS. The vector then keeps those four bases and the drop-out
  replaces the rest of the CDS.
- the adaptor also spells the four overhang bases, drawn from the backbone flanking the
  insert. The vector then drops the whole CDS, which is what tiling does at a CDS end.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

from ..checks.motifs import ENZYME_SITES, cut_geometry
from ..regions import reverse_complement
from .vector_io import DestinationContext, assemble_vector, flanks, locating_kwargs, resolve_destination


def _find_all(hay: str, needle: str) -> list[int]:
    out, i = [], hay.find(needle)
    while i != -1:
        out.append(i)
        i = hay.find(needle, i + 1)
    return out


@dataclass
class InsertCut:
    """How the Type IIS enzyme releases the insert from an assembled construct.

    Positions are indices on the construct (``adaptor_5 + variable + adaptor_3``):
    ``cut_5`` is the first base the released fragment keeps and ``cut_3`` is one past its
    last. ``keep_5`` / ``keep_3`` are how many adaptor bases ride along on each end, which
    is what decides where the vector's drop-out has to start and end.
    """

    enzyme: str
    cut_5: int
    cut_3: int
    overhang_5: str
    overhang_3: str
    keep_5: int
    keep_3: int


def cut_construct(
    adaptor_5: str, variable: str, adaptor_3: str, enzyme: str
) -> tuple[InsertCut | None, list[str]]:
    """Work out where ``enzyme`` cuts the assembled construct, as ``(cut, issues)``.

    ``cut`` is None when the adaptors do not describe a usable pair of cuts (no site, an
    ambiguous one, or a cut that falls inside the coding region); ``issues`` says why in
    words. Findings are returned rather than raised so QC can report every one of them at
    once and the caller can still emit a whole-CDS drop-out vector.
    """
    site = ENZYME_SITES[enzyme].upper()
    site_rc = reverse_complement(site)
    spacer, overhang = cut_geometry(enzyme)
    a5, a3 = adaptor_5.upper(), adaptor_3.upper()
    var = variable.upper()
    issues: list[str] = []

    fwd5, rev5 = _find_all(a5, site), _find_all(a5, site_rc)
    fwd3, rev3 = _find_all(a3, site), _find_all(a3, site_rc)
    if rev5:
        issues.append(
            f"adaptor_5 carries a reverse-strand {enzyme} site, which cuts away from the "
            "insert instead of releasing it"
        )
    if fwd3:
        issues.append(
            f"adaptor_3 carries a forward-strand {enzyme} site, which cuts away from the "
            "insert instead of releasing it"
        )
    for label, hits, strand in (("adaptor_5", fwd5, "forward"), ("adaptor_3", rev3, "reverse")):
        if not hits:
            issues.append(
                f"{label} carries no {strand}-strand {enzyme} site, so nothing releases "
                f"that end of the insert (expected {site if strand == 'forward' else site_rc})"
            )
        elif len(hits) > 1:
            issues.append(
                f"{label} carries {len(hits)} {strand}-strand {enzyme} sites, so which one "
                "cuts is ambiguous"
            )
    if len(fwd5) != 1 or len(rev3) != 1:
        return None, issues

    # The site nearest the insert is the one that releases it, and with exactly one site
    # per adaptor that is the site we found.
    cut_5 = fwd5[0] + len(site) + spacer
    cut_3 = len(a5) + len(var) + rev3[0] - spacer
    keep_5 = len(a5) - cut_5
    keep_3 = cut_3 - (len(a5) + len(var))
    if keep_5 < 0:
        issues.append(
            f"the 5' {enzyme} cut falls {-keep_5} base(s) inside the coding region, so the "
            "oligo would lose the start of the CDS. Move the site further from the variable "
            "region in adaptor_5."
        )
    if keep_3 < 0:
        issues.append(
            f"the 3' {enzyme} cut falls {-keep_3} base(s) inside the coding region, so the "
            "oligo would lose the end of the CDS. Move the site further from the variable "
            "region in adaptor_3."
        )
    if keep_5 < 0 or keep_3 < 0:
        return None, issues

    construct = a5 + var + a3
    return (
        InsertCut(
            enzyme=enzyme,
            cut_5=cut_5,
            cut_3=cut_3,
            overhang_5=construct[cut_5:cut_5 + overhang],
            overhang_3=construct[cut_3 - overhang:cut_3],
            keep_5=keep_5,
            keep_3=keep_3,
        ),
        issues,
    )


@dataclass
class DestinationVector:
    """The one plasmid to clone for a standard library.

    ``sequence`` is the full destination vector, the starting plasmid with
    ``[start, end)`` of the CDS replaced by the drop-out. A circular plasmid is read from an
    origin upstream of the insert's promoter (see ``vector_io.choose_origin``), so the
    cassette is not split across the map's ends. ``overhang_5`` / ``overhang_3``
    are what the *cut vector* presents at the two junctions, which is what the oligo's own
    overhangs have to equal.
    """

    sequence: str
    topology: str
    start: int                  # drop-out window on the reference CDS
    end: int
    overhang_5: str
    overhang_3: str
    cut: InsertCut | None       # None when the adaptors describe no usable cut
    dest: DestinationContext
    issues: list[str] = field(default_factory=list)
    reference_len: int = 0      # the frozen CDS the drop-out window is measured against

    @property
    def length(self) -> int:
        """Size of the destination vector in bases."""
        return len(self.sequence)

    @property
    def overhangs_match(self) -> bool | None:
        """Whether the oligo's fused overhangs are the ones the cut vector presents, or None
        when the adaptors carry no usable cut to compare."""
        if self.cut is None:
            return None
        return (self.cut.overhang_5, self.cut.overhang_3) == (self.overhang_5, self.overhang_3)

    def __str__(self) -> str:
        from .vector_io import origin_in_source

        source = Path(self.dest.source_path).name or "the starting vector"
        d = self.dest
        rows = [
            ("starting plasmid", f"{len(d.full_seq)} bp {d.topology}"
                                 + (", insert on the minus strand" if d.flipped else "")),
            ("insert located", f"{d.start}-{d.end} ({len(d.located_region)} bp)"),
            ("drop-out window", f"{self.start}-{self.end} of the {self.reference_len} bp reference CDS"),
            ("vector overhangs", f"5' {self.overhang_5}   3' {self.overhang_3}"),
        ]
        if self.cut is None:
            rows.append(("oligo overhangs", "none, the adaptors carry no cut site"))
        else:
            verdict = "match" if self.overhangs_match else "MISMATCH, nothing will ligate"
            rows.append(("oligo overhangs",
                         f"5' {self.cut.overhang_5}   3' {self.cut.overhang_3}   ({verdict})"))
            rows.append(("adaptor bases kept",
                         f"{self.cut.keep_5} at the 5' end, {self.cut.keep_3} at the 3'"))
        read_from = (f", read from base {origin_in_source(d)} of {source}"
                     if d.topology == "circular" or d.flipped else "")
        rows.append(("emitted vector", f"{self.length} bp{read_from}"))

        width = max(len(k) for k, _ in rows)
        lines = [f"Destination vector for {source}"]
        lines += [f"  {k.ljust(width)}  {v}" for k, v in rows]
        lines += [f"  ! {m}" for m in self.issues]
        return "\n".join(lines)

    def __repr__(self) -> str:      # plain output in a notebook rather than a dataclass dump
        return self.__str__()

    def _repr_html_(self) -> str:
        from html import escape

        return f"<pre style='margin:0;line-height:1.4'>{escape(str(self))}</pre>"


def build_destination(library) -> DestinationVector:
    """Build the destination vector for a standard (untiled) library.

    Needs a starting vector on the spec and a codon-optimized reference. The drop-out
    window is placed so the cut vector presents exactly the overhangs the digested oligo
    carries; with no usable cut in the adaptors the whole CDS drops out and the overhangs
    come from the backbone flanking the insert. Raises when the adaptors ask for something
    the plasmid cannot present.
    """
    spec = library.spec
    if getattr(library, "tiles", None) is not None:
        raise ValueError(
            "This is a tiled library, which has one destination vector per tile. Read "
            "them off library.tiles, or export with to_vectors() / to_vector_maps()."
        )
    vec = spec.vector
    if vec is None:
        raise ValueError(
            "No starting vector to build a destination vector from. Set "
            "spec.starting_vector to the plasmid you clone into."
        )
    if getattr(library, "kind", "scan") == "sequence_set":
        raise ValueError(
            "A sequence set has no single reference CDS to drop out: its members are "
            "independent full-length sequences, so each one implies its own vector. Build "
            "the destination vector for whichever member you are cloning."
        )
    reference = library.reference
    if reference is None:
        raise ValueError("Library is not codon-optimized yet, call codon_optimize() first.")

    dest = resolve_destination(vec.path, **locating_kwargs(spec))
    overhang = cut_geometry(vec.enzyme)[1]
    cut, issues = cut_construct(spec.adaptor_5, reference, spec.adaptor_3, vec.enzyme)

    if cut is None:
        # No cut to honor, so fall back to the tiled convention at a CDS end: the whole CDS
        # drops out and the backbone supplies both overhangs.
        s, e = 0, len(reference)
        keep_5 = keep_3 = overhang
    else:
        keep_5, keep_3 = cut.keep_5, cut.keep_3
        if keep_5 > overhang or keep_3 > overhang:
            raise ValueError(
                f"the adaptors keep more bases past the {vec.enzyme} cut ({keep_5} at the 5' "
                f"end, {keep_3} at the 3') than the {overhang} bp fused overhang, so the "
                "destination vector would have to give up backbone bases the oligo "
                "re-supplies. Trim the adaptors so each ends at its overhang."
            )
        s, e = overhang - keep_5, len(reference) - (overhang - keep_3)
    if e <= s:
        raise ValueError(
            f"the reference CDS ({len(reference)} bp) is too short to carry both "
            f"{overhang} bp fused overhangs, so there is nothing left to drop out."
        )

    # The retained arms are drawn from the *reference*, not from whatever CDS the plasmid
    # holds today, so the overhangs the cut vector presents are the ones the oligos carry
    # even when the reference was codon-optimized away from the plasmid's own sequence.
    cds_region = reference[:s] + vec.vector_insert.upper() + reference[e:]
    flank_5, flank_3 = flanks(dest, keep_5, keep_3)
    return DestinationVector(
        sequence=assemble_vector(dest, cds_region),
        topology=dest.topology,
        start=s,
        end=e,
        reference_len=len(reference),
        overhang_5=flank_5 + reference[:s],
        overhang_3=reference[e:] + flank_3,
        cut=cut,
        dest=dest,
        issues=issues,
    )
