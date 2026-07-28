"""Read a destination-plasmid file and work out where the CDS insert goes.

Every library is cloned into something. A standard library drops its one oligo into one
destination vector (see ``layout/destination.py``); tiled assembly drops each amplified
sublibrary into its own vector carrying the rest of the CDS (``layout/tiled.py``). To emit
those as real plasmids rather than a bare CDS cassette, we need the backbone the user
clones into. This module reads that backbone from a ``.gb`` / ``.dna`` / ``.fasta`` file,
finds the insert locus, and reduces everything to three facts the assembler needs: the
bases just 5' of the insert, the bases just 3' of it, and whether the molecule is circular.

Reading the files needs BioPython. One dependency covers
GenBank, SnapGene ``.dna``, and FASTA, and it surfaces the topology and feature
annotations we use to locate the insert and to carry features onto the emitted maps.
BioPython reads ``.dna`` but cannot write it, so emitted maps are GenBank.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from functools import lru_cache
from pathlib import Path

from ..regions import reverse_complement

_EXT_FORMAT = {
    ".gb": "genbank", ".gbk": "genbank", ".genbank": "genbank",
    ".dna": "snapgene",
    ".fa": "fasta", ".fasta": "fasta", ".fna": "fasta", ".seq": "fasta",
}


@dataclass
class Feature:
    """One annotation read off a destination file, kept as a plain record rather than a
    BioPython object so the rest of the package does not depend on that type.

    ``[start, end)`` is the span on the file's sequence as read, and ``strand`` is 1,
    -1, or None, whatever the file gives. ``wraps_origin`` marks a location stored in
    more than one part (a GenBank ``join()``), which on a circular plasmid means the
    feature crosses base 1, so its real span is not ``[start, end)``. Locating the
    insert refuses such a feature rather than reading it as one stretch.
    """

    type: str
    start: int          # [start, end) on the file's linear sequence
    end: int
    strand: int | None
    label: str | None
    wraps_origin: bool = False   # stored as a join() crossing base 1, so the span is not [start, end)


@dataclass
class DestinationContext:
    """A located insert site in a destination plasmid, everything the assembler needs.

    ``full_seq`` is the backbone as read (uppercase). ``[start, end)`` is the located
    insert region on it (the CDS already there, a stuffer, or an annotated site).
    ``located_region`` is those bases, used verbatim as the reference when
    ``use_vector_cds`` is set. ``origin`` is where an emitted circular vector starts
    reading (see ``choose_origin``). Flanks, overhangs, and offsets are derived
    topology-aware by the helpers below, so callers never slice ``full_seq`` themselves.
    """
    full_seq: str
    topology: str          # "circular" | "linear"
    start: int
    end: int
    located_region: str
    features: list[Feature] = field(default_factory=list)
    source_path: str = ""
    origin: int = 0        # base of full_seq the emitted map starts at (circular only)
    flipped: bool = False  # the file was reverse-complemented to put the CDS on the forward strand


def _format_for(path: str) -> str:
    ext = Path(path).suffix.lower()
    fmt = _EXT_FORMAT.get(ext)
    if fmt is None:
        raise ValueError(
            f"Unsupported vector file extension {ext!r} for {path!r}; "
            f"use one of {sorted(set(_EXT_FORMAT))}."
        )
    return fmt


def _require_biopython():
    try:
        from Bio import SeqIO  # noqa: F401
    except ImportError as exc:
        raise ImportError(
            "Reading a starting vector (.gb / .dna / .fasta) needs BioPython, a base "
            "dependency. Reinstall library-designer if it is missing."
        ) from exc
    return SeqIO


def read_vector_file(path: str) -> tuple[str, str, list[Feature]]:
    """``(sequence, topology, features)`` from a GenBank / SnapGene / FASTA file.

    Topology comes from the file (GenBank ``LOCUS`` / SnapGene flag); FASTA has none,
    so it defaults to ``linear``. Features keep their type, span, strand, and label."""
    fmt = _format_for(path)
    SeqIO = _require_biopython()
    rec = SeqIO.read(str(path), fmt)
    seq = str(rec.seq).upper()
    if set(seq) - set("ACGTN"):
        bad = "".join(sorted(set(seq) - set("ACGTN")))
        raise ValueError(f"{path!r} sequence contains non-ACGTN characters: {bad!r}.")
    topo = rec.annotations.get("topology")
    topology = topo if topo in ("circular", "linear") else "linear"
    features: list[Feature] = []
    for f in rec.features:
        quals = f.qualifiers or {}
        label = None
        # regulatory_class last: it only fills in for a GenBank "regulatory" feature that
        # carries no name of its own, where the class ("promoter", "terminator") is the
        # only thing to call it.
        for key in ("label", "gene", "product", "note", "regulatory_class"):
            if quals.get(key):
                label = quals[key][0]
                break
        # A feature that crosses the file's origin is stored as a join() of two parts, and
        # BioPython reports its span as 0..len(seq), which would read as "the whole molecule".
        # Flag it so locating can refuse rather than silently swallow the backbone.
        parts = getattr(f.location, "parts", [f.location])
        features.append(
            Feature(
                type=f.type,
                start=int(f.location.start),
                end=int(f.location.end),
                strand=f.location.strand,
                label=label,
                wraps_origin=len(parts) > 1,
            )
        )
    return seq, topology, features


def _find_unique(haystack: str, needle: str) -> int:
    """Index of ``needle`` in ``haystack`` if it occurs exactly once, else raise."""
    i = haystack.find(needle)
    if i == -1:
        return -1
    if haystack.find(needle, i + 1) != -1:
        raise ValueError("sequence occurs more than once in the vector; the locus is ambiguous.")
    return i


def locate_insert(
    seq: str,
    features: list[Feature],
    *,
    search_cds: str | None = None,
    insert_label: str | None = None,
    anchors: tuple[str, str] | None = None,
) -> tuple[int, int, int]:
    """Find ``(start, end, strand)`` of the insert site on ``seq``.

    ``strand`` is ``+1`` when the CDS reads on the given (forward) strand and ``-1``
    when it reads on the reverse strand (a very common cloning orientation). The caller
    normalizes a reverse-strand hit into the CDS-sense frame.

    Precedence: an exact unique match of ``search_cds`` (forward or reverse strand),
    then an annotated feature named ``insert_label`` (or the sole CDS feature when
    ``insert_label`` is ``None`` and no CDS was given), then two unique ``anchors``
    bracketing the site. Raises with a specific reason when nothing resolves or a match
    is ambiguous.
    """
    if search_cds:
        needle = search_cds.upper()
        i = _find_unique(seq, needle)
        if i != -1:
            return i, i + len(needle), 1
        j = _find_unique(seq, reverse_complement(needle))
        if j != -1:
            return j, j + len(needle), -1
        raise ValueError(
            "the WT CDS was not found in the starting vector (either strand). Check that "
            "the vector holds this exact CDS, or locate the site with insert_label "
            "(an annotated feature) or insert_anchors."
        )

    def _span(f: Feature) -> tuple[int, int, int]:
        if f.wraps_origin:
            raise ValueError(
                f"the insert feature {f.label or f.type!r} crosses the origin of the file, so "
                "its span is not a single stretch of sequence. Rotate the plasmid in your "
                "sequence editor so the insert does not cross base 1, or bracket the site "
                "with insert_anchors."
            )
        return f.start, f.end, (f.strand or 1)

    if insert_label is not None:
        hits = [f for f in features if f.label == insert_label]
        if not hits:
            raise ValueError(
                f"no feature labelled {insert_label!r} in the vector; "
                f"labels present: {sorted({f.label for f in features if f.label})}."
            )
        if len(hits) > 1:
            raise ValueError(f"feature label {insert_label!r} is not unique in the vector.")
        return _span(hits[0])

    cds_feats = [f for f in features if f.type == "CDS"]
    if len(cds_feats) == 1:
        return _span(cds_feats[0])
    if len(cds_feats) > 1:
        raise ValueError(
            "the vector has more than one CDS feature; name the insert with "
            f"insert_label (labels present: {sorted({f.label for f in cds_feats if f.label})}) "
            "or bracket it with insert_anchors."
        )

    if anchors:
        a5, a3 = anchors[0].upper(), anchors[1].upper()
        i = _find_unique(seq, a5)
        if i == -1:
            raise ValueError(f"5' anchor {a5!r} not found in the vector.")
        start = i + len(a5)
        j = seq.find(a3, start)
        if j == -1:
            raise ValueError(f"3' anchor {a3!r} not found 3' of the 5' anchor in the vector.")
        if seq.find(a3, j + 1) != -1:
            raise ValueError(f"3' anchor {a3!r} is not unique 3' of the 5' anchor.")
        return start, j, 1

    raise ValueError(
        "could not locate the insert site: give the WT CDS (spec.cds), an annotated "
        "feature (insert_label), a single CDS feature, or insert_anchors."
    )


ORIGIN_UPSTREAM = 50   # bases to back off from the promoter when placing the map's origin


def _is_promoter(f: Feature) -> bool:
    """A promoter that could drive a forward-strand insert. Covers the SnapGene/GenBank
    ``promoter`` type and the ``regulatory`` type whose class says promoter, which
    ``read_vector_file`` surfaces as the label."""
    return (f.type == "promoter" or "promoter" in (f.label or "").lower()) and f.strand != -1


def choose_origin(seq: str, features: list[Feature], start: int, end: int, topology: str,
                  upstream: int = ORIGIN_UPSTREAM) -> int:
    """Which base an emitted circular vector should start reading at.

    A viewer draws a plasmid map from base 1, so whatever sits at the origin is split
    across the two ends of the map. Starting at the insert cuts the cassette in half right
    at the Golden Gate sites, the part you actually want to read. So the origin goes
    ``upstream`` bases before the promoter that drives the insert, the nearest one 5' of it,
    which leaves promoter and insert intact and in reading order. With no promoter
    annotated we back off the same distance from the insert itself. The result is then
    nudged to a feature boundary if it would otherwise cut through an annotation, since a
    feature spanning the origin cannot be drawn (or written to GenBank) in one piece.
    Linear molecules are never rotated.
    """
    if topology != "circular":
        return start
    n = len(seq)
    upstream_of = [f for f in features if f.end <= start or f.start >= end]   # not on the insert
    promoters = [f for f in upstream_of if _is_promoter(f)]
    # "The first one upstream": the smallest gap from a promoter's 3' end back to the
    # insert, measured around the circle so a promoter behind the file's origin counts too.
    anchor = min(promoters, key=lambda f: (start - f.end) % n).start if promoters else start
    origin = (anchor - upstream) % n
    for _ in range(8):                     # features can nest; a few passes settle it
        cut = next((f for f in upstream_of if f.start < origin < f.end), None)
        if cut is None:
            break
        origin = cut.start % n
    return start if start <= origin < end else origin


def _flip_features(features: list[Feature], n: int) -> list[Feature]:
    """Features re-expressed on the reverse-complemented sequence of length ``n``."""
    return [
        Feature(
            type=f.type,
            start=n - f.end,
            end=n - f.start,
            strand=(None if f.strand is None else -f.strand),
            label=f.label,
        )
        for f in features
    ]


def resolve_destination(
    path: str,
    *,
    topology_override: str | None = None,
    search_cds: str | None = None,
    insert_label: str | None = None,
    anchors: tuple[str, str] | None = None,
) -> DestinationContext:
    """Read ``path`` and locate the insert, returning a ``DestinationContext``.

    Cached per (path, the file's mtime and size, locating inputs) so ``build_reference``
    and ``tile_library`` share one read and one locate, keeping the reference and the
    emitted vectors consistent, while a plasmid edited on disk is read again rather than
    served stale from an earlier call in the same session.
    """
    try:
        st = Path(path).stat()
        stamp: tuple[int, int] | None = (st.st_mtime_ns, st.st_size)
    except OSError:
        stamp = None            # let the reader below raise the real error
    return _resolve_destination(
        str(path), stamp, topology_override=topology_override, search_cds=search_cds,
        insert_label=insert_label, anchors=anchors,
    )


@lru_cache(maxsize=32)
def _resolve_destination(
    path: str,
    stamp: tuple[int, int] | None,
    *,
    topology_override: str | None = None,
    search_cds: str | None = None,
    insert_label: str | None = None,
    anchors: tuple[str, str] | None = None,
) -> DestinationContext:
    seq, topology, features = read_vector_file(path)
    if topology_override:
        topology = topology_override
    start, end, strand = locate_insert(
        seq, features, search_cds=search_cds, insert_label=insert_label, anchors=anchors
    )
    if strand == -1:
        # The CDS reads on the reverse strand of the file. Re-express the whole plasmid
        # (and its features) on the CDS-sense strand so tiling, overhangs, and the emitted
        # map all live in one frame with the CDS on the forward strand.
        n = len(seq)
        seq = reverse_complement(seq)
        features = _flip_features(features, n)
        start, end = n - end, n - start
    return DestinationContext(
        full_seq=seq,
        topology=topology,
        start=start,
        end=end,
        located_region=seq[start:end],
        features=features,
        source_path=path,
        origin=choose_origin(seq, features, start, end, topology),
        flipped=(strand == -1),
    )


def locating_kwargs(spec, params=None) -> dict:
    """The insert-locating inputs for ``resolve_destination``, derived from a spec so
    that reference extraction, tiling, and the destination vector agree on the locus.

    When the reference is the vector's own CDS we cannot search for it (we are about to
    extract it), so locating leans on the label / CDS feature / anchors. Otherwise the
    known CDS (``spec.cds``) is the search key, falling back to label / anchors.

    ``params`` is the tiled-assembly params actually in use, ``spec.tiled`` when omitted.
    Pass it so a library tiled with an explicit ``tile(params)`` locates the same insert
    its layout used."""
    vec = spec.resolve_vector(params)
    if vec is None:
        raise ValueError(
            "No starting vector to locate an insert in. Set spec.starting_vector (or "
            "tiled.starting_vector) to the destination plasmid."
        )
    return dict(
        topology_override=vec.topology,
        search_cds=None if vec.use_vector_cds else spec.cds,
        insert_label=vec.insert_label,
        anchors=vec.insert_anchors,
    )


def flanks(dest: DestinationContext, n5: int, n3: int) -> tuple[str, str]:
    """The ``n5`` backbone bases immediately 5' of the insert and the ``n3`` immediately
    3' of it. Wraps the origin for a circular plasmid; either count may be 0."""
    s, e, S = dest.start, dest.end, dest.full_seq
    if dest.topology == "circular":
        arc = S[e:] + S[:s]          # the whole non-insert backbone, as one arc
        if not arc:
            raise ValueError(
                "the located insert is the entire plasmid, so there is no backbone to clone "
                "into. Check that the insert feature or anchors mark the coding region only."
            )
        if len(arc) < max(n5, n3):
            raise ValueError("backbone is shorter than the overhang length; cannot draw terminal overhangs.")
        return (arc[len(arc) - n5:] if n5 else ""), arc[:n3]
    left, right = S[:s], S[e:]
    if len(left) < n5 or len(right) < n3:
        raise ValueError(
            "the insert sits within one overhang length of a linear end; extend the "
            "flanking sequence in the starting vector so a full terminal overhang exists."
        )
    return (left[len(left) - n5:] if n5 else ""), right[:n3]


def terminal_contexts(dest: DestinationContext, overhang_len: int) -> tuple[str, str]:
    """The two backbone overhangs at the CDS ends: the ``overhang_len`` bases just 5'
    of the insert and just 3' of it. Wraps the origin for a circular plasmid."""
    return flanks(dest, overhang_len, overhang_len)


def origin_in_source(dest: DestinationContext) -> int:
    """The 1-based base of the *file* that an emitted map starts reading at.

    When the insert was found on the minus strand the whole plasmid was reverse-complemented
    so the CDS reads forward, and the emitted map runs the other way round the circle than
    the file does. This gives the coordinate on the file's own numbering either way, so a
    manifest can point back at the plasmid the user supplied."""
    i = dest.origin if dest.topology == "circular" else 0
    return (len(dest.full_seq) - i) if dest.flipped else (i + 1)


def insert_offset(dest: DestinationContext) -> int:
    """Index in the emitted vector where the CDS region begins, so callers can place
    features on what ``assemble_vector`` returns without redoing the rotation."""
    if dest.topology != "circular":
        return dest.start
    return (dest.start - dest.origin) % len(dest.full_seq)


def assemble_vector(dest: DestinationContext, cds_region: str) -> str:
    """The full destination-vector sequence with ``cds_region`` at the insert locus.

    For a circular plasmid the returned string is a linear representation of the circle
    read from ``dest.origin`` (see ``choose_origin``), which sits upstream of the insert's
    promoter so the map does not split the cassette. ``insert_offset`` says where
    ``cds_region`` landed."""
    s, e, S = dest.start, dest.end, dest.full_seq
    if dest.topology != "circular":
        return S[:s] + cds_region + S[e:]
    p = dest.origin
    if p <= s:
        return S[p:s] + cds_region + S[e:] + S[:p]
    return S[p:] + S[:s] + cds_region + S[e:p]
