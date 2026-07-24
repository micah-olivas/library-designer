"""Read a destination-plasmid file and work out where the CDS insert goes.

Tiled assembly drops each amplified sublibrary into its own destination vector that
carries the rest of the CDS (see ``layout/tiled.py``). To emit those vectors as real
plasmids rather than a bare CDS cassette, we need the backbone the user clones into.
This module reads that backbone from a ``.gb`` / ``.dna`` / ``.fasta`` file, finds the
insert locus, and reduces everything to three facts the assembler needs: the bases
just 5' of the insert, the bases just 3' of it, and whether the molecule is circular.

Reading the files needs BioPython (behind the ``tiled`` extra). One dependency covers
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
    type: str
    start: int          # [start, end) on the file's linear sequence
    end: int
    strand: int | None
    label: str | None


@dataclass
class DestinationContext:
    """A located insert site in a destination plasmid, everything the assembler needs.

    ``full_seq`` is the backbone as read (uppercase). ``[start, end)`` is the located
    insert region on it (the CDS already there, a stuffer, or an annotated site).
    ``located_region`` is those bases, used verbatim as the reference when
    ``use_vector_cds`` is set. Flanks and overhangs are derived topology-aware by the
    helpers below, so callers never slice ``full_seq`` themselves.
    """
    full_seq: str
    topology: str          # "circular" | "linear"
    start: int
    end: int
    located_region: str
    features: list[Feature] = field(default_factory=list)
    source_path: str = ""


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
            "Reading a starting vector (.gb / .dna / .fasta) needs BioPython. "
            "Install the tiled extra:  pip install 'library-designer[tiled]'  "
            "(or  uv sync --extra tiled)."
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
        for key in ("label", "gene", "product", "note"):
            if quals.get(key):
                label = quals[key][0]
                break
        features.append(
            Feature(
                type=f.type,
                start=int(f.location.start),
                end=int(f.location.end),
                strand=f.location.strand,
                label=label,
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
            "the vector holds this exact CDS, or locate the site with tiled.insert_label "
            "or tiled.insert_anchors."
        )

    if insert_label is not None:
        hits = [f for f in features if f.label == insert_label]
        if not hits:
            raise ValueError(
                f"no feature labelled {insert_label!r} in the vector; "
                f"labels present: {sorted({f.label for f in features if f.label})}."
            )
        if len(hits) > 1:
            raise ValueError(f"feature label {insert_label!r} is not unique in the vector.")
        return hits[0].start, hits[0].end, (hits[0].strand or 1)

    cds_feats = [f for f in features if f.type == "CDS"]
    if len(cds_feats) == 1:
        return cds_feats[0].start, cds_feats[0].end, (cds_feats[0].strand or 1)
    if len(cds_feats) > 1:
        raise ValueError(
            "the vector has more than one CDS feature; name the insert with "
            "tiled.insert_label or bracket it with tiled.insert_anchors."
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
        "feature (tiled.insert_label), a single CDS feature, or tiled.insert_anchors."
    )


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


@lru_cache(maxsize=32)
def resolve_destination(
    path: str,
    *,
    topology_override: str | None = None,
    search_cds: str | None = None,
    insert_label: str | None = None,
    anchors: tuple[str, str] | None = None,
) -> DestinationContext:
    """Read ``path`` and locate the insert, returning a ``DestinationContext``.

    Cached per (path, locating inputs) so ``build_reference`` and ``tile_library`` share
    one read and one locate, keeping the reference and the emitted vectors consistent.
    """
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
    )


def locating_kwargs(spec) -> dict:
    """The insert-locating inputs for ``resolve_destination``, derived from a spec so
    that reference extraction and tiling agree on the locus.

    When the reference is the vector's own CDS we cannot search for it (we are about to
    extract it), so locating leans on the label / CDS feature / anchors. Otherwise the
    known CDS (``spec.cds``) is the search key, falling back to label / anchors."""
    t = spec.tiled
    anchors = tuple(t.insert_anchors) if t.insert_anchors else None
    search_cds = None if t.use_vector_cds else spec.cds
    return dict(
        topology_override=t.topology,
        search_cds=search_cds,
        insert_label=t.insert_label,
        anchors=anchors,
    )


def terminal_contexts(dest: DestinationContext, overhang_len: int) -> tuple[str, str]:
    """The two backbone overhangs at the CDS ends: the ``overhang_len`` bases just 5'
    of the insert and just 3' of it. Wraps the origin for a circular plasmid."""
    o = overhang_len
    s, e, S = dest.start, dest.end, dest.full_seq
    if dest.topology == "circular":
        arc = S[e:] + S[:s]          # the whole non-insert backbone, as one arc
        if len(arc) < o:
            raise ValueError("backbone is shorter than the overhang length; cannot draw terminal overhangs.")
        return arc[-o:], arc[:o]
    left, right = S[:s], S[e:]
    if len(left) < o or len(right) < o:
        raise ValueError(
            "the insert sits within one overhang length of a linear end; extend the "
            "flanking sequence in the starting vector so a full terminal overhang exists."
        )
    return left[-o:], right[:o]


def assemble_vector(dest: DestinationContext, cds_region: str) -> str:
    """The full destination-vector sequence with ``cds_region`` at the insert locus.

    For a circular plasmid the returned string is a linear representation of the circle
    with its origin at the CDS start (backbone end joins the CDS start)."""
    s, e, S = dest.start, dest.end, dest.full_seq
    if dest.topology == "circular":
        return cds_region + S[e:] + S[:s]
    return S[:s] + cds_region + S[e:]
