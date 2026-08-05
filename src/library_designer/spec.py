"""A ``LibrarySpec`` plus a generator fully determines a library, so the spec doubles
as the design-specs record and can be loaded from TOML.
"""
from __future__ import annotations

import tomllib
from dataclasses import asdict, dataclass, field, replace
from html import escape as _esc
from pathlib import Path


def _preview(seq: str, keep: int = 12) -> str:
    """Middle-elide a long sequence for display: head...tail."""
    return seq if len(seq) <= 2 * keep + 1 else f"{seq[:keep]}...{seq[-keep:]}"

# Shine-Dalgarno-like motifs to avoid: ribosome-binding-site core, an 8-12 nt
# spacer, then a start codon. Regex, matching the MBO-038 amber-scan design.
DEFAULT_AVOID_PATTERNS: tuple[str, ...] = (
    "GGAGG.{8,12}[AG]TG",
    "AGGAG.{8,12}ATG",
)


@dataclass
class CodonOptimizationParams:
    """Codon-optimization knobs. Captured in the design specs so an order is self-documenting."""

    species: str = "e_coli"
    method: str = "use_best_codon"   # DNA Chisel: use_best_codon | match_codon_usage | harmonize_rca
    gc_min: float | None = None      # fraction, e.g. 0.30
    gc_max: float | None = None      # fraction, e.g. 0.68 (IDT eBlock GC ceiling)
    gc_window: int | None = None     # windowed GC in bp; None = whole sequence
    max_random_iters: int = 100_000
    # What to do when the preferred codon for a mutated residue would introduce a
    # restricted motif. True steps down that residue's usage ranking to the next codon
    # that avoids it, so the variant is still makeable at a rarer codon. False tries the
    # preferred codon only and records the variant as a failure, so no rarer codon is
    # substituted without you asking for it. A pinned literal codon (an amber TAG) has no
    # synonymous alternative, so either way it goes in verbatim and is flagged if it
    # introduces a motif. Scan libraries only. A SequenceSet has no single mutated
    # position, so DNA Chisel resolves the whole CDS instead.
    synonymous_fallback: bool = True


def optimization_line(params: dict) -> str:
    """The codon-optimization params on one line, keeping only the ones that shaped the run.

    Takes the params as a dict (``asdict`` of ``CodonOptimizationParams``, or the block
    read back from a design-specs record) so both ``LibrarySpec`` and ``LibrarySummary``
    print them the same way. An unset GC bound is left out, and so is the iteration cap
    under ``use_best_codon``, which picks each codon outright rather than searching.
    ``synonymous_fallback`` appears only when it is off, since that changes which variants
    are makeable and belongs in the record.
    """
    bits = [str(params.get("species")), str(params.get("method"))]
    gc_min, gc_max = params.get("gc_min"), params.get("gc_max")
    if gc_min is not None or gc_max is not None:
        gc = (
            f"gc {gc_min}-{gc_max}" if gc_min is not None and gc_max is not None
            else f"gc min {gc_min}" if gc_min is not None
            else f"gc max {gc_max}"
        )
        if params.get("gc_window"):
            gc += f" over {params['gc_window']} bp windows"
        bits.append(gc)
    if params.get("method") != "use_best_codon" and params.get("max_random_iters"):
        bits.append(f"{params['max_random_iters']} iters")
    # Absent from an older record, where the fallback was unconditional, so default to on.
    if not params.get("synonymous_fallback", True):
        bits.append("no synonymous fallback")
    return ", ".join(bits)


@dataclass
class StartingVectorParams:
    """The destination plasmid the library is cloned into, and how to find the insert site.

    Set ``spec.starting_vector`` to a file path for the common case; pass this dataclass
    instead when you need to name the insert site, name the cloning enzyme, or freeze the
    plasmid's own CDS. Either way the pipeline can emit the destination vector to clone
    (``to_vectors`` / ``to_vector_maps``) and check the adaptors against it, because it
    knows both the backbone and the overhangs the digested oligo will carry.

    A ``[tiled]`` block carries the same fields (``TiledAssemblyParams``); when both are
    set the tiled block wins, since that is what the layout was built with. See
    ``LibrarySpec.resolve_vector``.
    """

    path: str | None = None              # the destination plasmid (.gb / .dna / .fasta)
    use_vector_cds: bool = False         # freeze the CDS already in the plasmid as the reference (clone from it)
    insert_label: str | None = None      # locate the insert site by an annotated feature label (else the WT CDS is matched)
    insert_anchors: tuple[str, str] | None = None  # fallback: (5', 3') unique sequences bracketing the insert site
    topology: str | None = None          # override the file topology ("circular" / "linear"); FASTA defaults to linear
    enzyme: str = "BsaI"                 # Golden Gate Type IIS enzyme the adaptors carry
    vector_insert: str = "AGAGACCAAAAGGTCTCA"   # drop-out placeholder that replaces the CDS in the destination vector


AMINO_ACIDS = frozenset("ACDEFGHIKLMNPQRSTVWY")


def _clean_sequence(value, field_name: str, alphabet: frozenset[str], what: str) -> str:
    """Normalize a pasted sequence: whitespace out, uppercase, then check the alphabet.

    Pasting from a FASTA or a paper brings newlines and spaces along, and a lowercase
    sequence used to survive all the way to a bogus variant name. Anything outside the
    alphabet is refused here rather than later: DNA Chisel reverse-translates ``X``, ``B``,
    and ``Z`` to an arbitrary residue without complaining, so an unchecked ambiguity code
    designs a library for a different protein."""
    seq = "".join(str(value).split()).upper()
    bad = sorted(set(seq) - alphabet)
    if bad:
        raise ValueError(
            f"{field_name} contains {what} character(s) {bad}. Give "
            + ("one-letter residues from the 20 canonical amino acids, with no stop ('*') "
               "or ambiguity ('X', 'B', 'Z') codes."
               if alphabet is AMINO_ACIDS else "A, C, G, or T only.")
        )
    return seq


def _as_vector_params(value) -> StartingVectorParams | None:
    """Normalize ``spec.starting_vector`` (a path, a ``StartingVectorParams``, or None)."""
    if value is None:
        return None
    if isinstance(value, StartingVectorParams):
        anchors = value.insert_anchors
        if anchors is not None and not isinstance(anchors, tuple):
            return replace(value, insert_anchors=tuple(anchors))
        return value
    if isinstance(value, dict):
        data = dict(value)
        if isinstance(data.get("insert_anchors"), list):
            data["insert_anchors"] = tuple(data["insert_anchors"])
        return StartingVectorParams(**data)
    return StartingVectorParams(path=str(value))


@dataclass
class TiledAssemblyParams:
    """Layout parameters for *tiled assembly*, a long-CDS library where each single
    mutant rides on a short oligo covering only its tile window, is amplified out of
    the shared pool by a tile-specific orthogonal primer pair, and is dropped by
    Golden Gate into a destination vector carrying the rest of the WT CDS.

    Tiling is pure layout: which variants exist and their sequences are already fixed
    by the frozen-reference stamp (see ``optimize/backbone.py``). Defaults reproduce
    the glucokinase (GCK-800) design. Captured in the design specs.
    """

    oligo_budget: int = 300            # hard cap on final oligo length (bp), incl. primers + sites
    enzyme: str = "BsaI"               # Golden Gate Type IIS enzyme (recognition from checks.motifs)
    overhang_len: int = 4              # fused-overhang length left after digestion
    spacer_5: str = "A"                # base(s) between the 5' recognition site and the overhang
    spacer_3: str = "T"                # base(s) between the overhang and the 3' recognition site
    vector_insert: str = "AGAGACCAAAAGGTCTCA"   # BsaI drop-out placeholder inserted into each destination vector
    primer_set: str = "subramanian2018"  # bundled set name or a path to a primer-set CSV (see primers.py)
    primer_length: int = 20            # per-primer length assumed when sizing tiles to the budget
    tile_size: int | None = None       # override the per-tile window (bp); None means it is derived from the budget
    # Move the tile boundaries, within the budget, to the positions whose fused overhangs
    # share the least homology. The overhangs are read off the CDS at the boundaries, so
    # where a boundary falls is the only handle on them (see layout/boundaries.py). Off by
    # default, since it changes the windows and so the oligos a library emits.
    optimize_overhangs: bool = False
    # Even the pool out to one oligo length. Tiles differ in size, so the oligos do too, and
    # moving the boundaries for better overhangs widens the spread. The filler goes between
    # each primer and the recognition site beside it, outside what the enzyme releases, so it
    # is amplified with the oligo and then cut away (see layout/tiled.py).
    pad_oligos: bool = False
    pad_target: int | None = None      # length to pad to; None means the longest oligo the layout needs
    # Add one WT member per tile (``WT_Tile_<i>``), the tile window straight from the
    # reference. Each tile is amplified and assembled on its own, so without these a
    # sublibrary ships with nothing unmutated to normalize against. Set False to order
    # only the mutants.
    wt_controls: bool = True

    # Destination-vector backbone. Point at the plasmid you clone into and the per-tile
    # destination vectors are emitted as the full plasmid (see layout/vector_io.py and
    # io.to_vector_maps). Without one, only the CDS-region cassette is emitted and the
    # two vector_context overhangs below are read verbatim. These five fields mirror
    # StartingVectorParams; a value set here wins over the spec-level one, since tiling
    # is what the layout was built with.
    starting_vector: str | None = None   # path to the destination plasmid (.gb / .dna / .fasta)
    use_vector_cds: bool = False         # freeze the CDS already in the starting vector as the reference (clone from it)
    insert_label: str | None = None      # locate the insert site by an annotated feature label (else the WT CDS is matched)
    insert_anchors: tuple[str, str] | None = None  # fallback: (5', 3') unique sequences bracketing the insert site
    topology: str | None = None          # override the file topology ("circular" / "linear"); FASTA defaults to linear
    vector_context_5: str = "TGGA"     # 5'-terminal overhang when NO starting vector is given (else derived from it)
    vector_context_3: str = "GGAG"     # 3'-terminal overhang when NO starting vector is given


@dataclass
class LibrarySpec:
    """Declarative description of a variant library. A spec plus a generator fully
    determines what gets ordered.

    Only ``name`` is required. A ``SubstitutionScan`` needs ``protein_sequence`` (or
    ``uniprot`` to fetch it) and ``substitutions``. A ``SequenceSet`` is handed its
    members separately and reads the spec for the adaptors and the optimization params.
    The remaining fields shape the library or the order form, and each one carries its
    meaning in the comment on its line.

    Several fields are normalized on assignment, at construction and on later assignment
    alike. A path given to ``starting_vector`` becomes a ``StartingVectorParams``, so
    reading the field back gives the dataclass and not the string you set. A dict given
    to ``optimization`` or ``tiled`` becomes its dataclass, which is what lets a spec
    round-trip through its own design-specs JSON. ``protein_sequence`` and ``cds`` lose
    their whitespace and are uppercased, and anything outside their alphabet is refused,
    since an unchecked ambiguity code would quietly design a library for a different
    protein.

    Setting ``uniprot`` with no ``protein_sequence`` fetches the entry when the spec is
    constructed, and only then. See ``resolve_uniprot``.
    """

    name: str
    protein_sequence: str = ""               # the single WT protein for a SubstitutionScan; unused by a SequenceSet
    # A UniProt accession to take the protein from instead of pasting it. Resolved once, on
    # construction, when protein_sequence is left empty; the sequence is then stored above so
    # the spec is self-contained and the build does not depend on UniProt later. An explicit
    # protein_sequence wins and the accession is kept as provenance. See uniprot.py.
    uniprot: str | None = None
    uniprot_entry: dict | None = None        # filled in by the lookup: entry name, organism, SV
    substitutions: list[str] = field(         # amino acids ("A") and/or literal codons ("TAG"); scan only
        default_factory=list
    )

    # Optional shaping
    cds: str | None = None                   # native WT CDS to use verbatim as the reference (else codon-optimized)
    # The plasmid this library is cloned into. A path is the common case; pass a
    # StartingVectorParams to name the insert site, the cloning enzyme, or to freeze the
    # plasmid's own CDS. With one set, the standard pipeline can emit the destination
    # vector (to_vectors / to_vector_maps) and QC checks the adaptors against it.
    starting_vector: str | StartingVectorParams | None = None
    tiled: TiledAssemblyParams | None = None  # when set, use tiled-assembly layout (long CDS split across oligos)
    truncation: int = 0                      # N-terminal residues to drop
    adaptor_5: str = ""                      # 5' flanking element (lowercased on export)
    adaptor_3: str = ""                      # 3' flanking element
    avoid_patterns: list[str] = field(
        default_factory=lambda: list(DEFAULT_AVOID_PATTERNS)
    )
    avoid_enzymes: list[str] = field(default_factory=lambda: ["BsaI"])
    optimization: CodonOptimizationParams = field(default_factory=CodonOptimizationParams)
    platform: str | None = None              # synthesis platform: "pooled"/"arrayed" or a provider slug ("twist_oligo_pools"); see methods.py
    max_oligo_length: int | None = None      # hard synthesis-length cap (bp, incl. adaptors)
    # Codon optimization is stochastic, so a seed is applied on every run and recorded in
    # the design specs: the same spec gives the same library on any machine. Set to None
    # to opt out and follow the ambient RNG instead.
    seed: int | None = 0

    def __post_init__(self) -> None:
        # Resolve the accession only when there is no sequence to use, so constructing a spec
        # reaches the network exactly when you asked it to and never otherwise.
        if self.uniprot and not self.protein_sequence:
            self.resolve_uniprot()

    def resolve_uniprot(self, refresh: bool = False, **kwargs) -> "LibrarySpec":
        """Fill ``protein_sequence`` from ``spec.uniprot``, recording which entry it came
        from. Called on construction when the sequence is empty; call it yourself after
        setting ``spec.uniprot`` later, or with ``refresh=True`` to fetch the entry again.

        Called directly it overwrites ``protein_sequence`` whether or not one is set,
        unlike the lookup at construction, which runs only when the field is empty. The
        FASTA is cached on disk, so a repeat call works offline unless ``refresh=True``."""
        from .uniprot import fetch

        if not self.uniprot:
            raise ValueError("No accession to resolve. Set spec.uniprot first.")
        entry = fetch(self.uniprot, refresh=refresh, **kwargs)
        self.protein_sequence = entry.sequence
        self.uniprot_entry = entry.record()
        return self

    @classmethod
    def from_uniprot(cls, accession: str, *, refresh: bool = False, timeout: float = 30.0,
                     cache: str | Path | None = None, **kwargs) -> "LibrarySpec":
        """A spec for a UniProt entry, its protein filled in from the database.

        ``name`` defaults to the entry name (``ACYP1_HUMAN``). Everything else is as for the
        constructor::

            spec = LibrarySpec.from_uniprot("P07311", substitutions=["A"], truncation=1)

        UniProt sequences carry the initiator methionine, so pass ``truncation=1`` when your
        construct starts at the second residue."""
        from .uniprot import fetch

        entry = fetch(accession, refresh=refresh, timeout=timeout, cache=cache)
        kwargs.setdefault("name", entry.entry_name)
        spec = cls(uniprot=accession, protein_sequence=entry.sequence, **kwargs)
        spec.uniprot_entry = entry.record()
        return spec

    def __setattr__(self, name, value):
        # Normalize on the way in, at construction and on later assignment alike, so the rest
        # of the package sees one shape for each field. A bare path becomes
        # StartingVectorParams (which also lets `spec.starting_vector = "p.gb"` be followed by
        # `spec.starting_vector.use_vector_cds = True`), a nested dict becomes its dataclass
        # (so a spec round-trips through its own design-specs JSON), and a pasted sequence
        # loses its whitespace and case before anything can misread it.
        if name == "starting_vector":
            value = _as_vector_params(value)
        elif name == "optimization" and isinstance(value, dict):
            value = CodonOptimizationParams(**value)
        elif name == "tiled" and isinstance(value, dict):
            data = dict(value)
            if isinstance(data.get("insert_anchors"), list):
                data["insert_anchors"] = tuple(data["insert_anchors"])
            value = TiledAssemblyParams(**data)
        elif name == "protein_sequence" and value:
            value = _clean_sequence(value, "protein_sequence", AMINO_ACIDS, "non-amino-acid")
        elif name == "cds" and value:
            value = _clean_sequence(value, "cds", frozenset("ACGT"), "non-ACGT")
        object.__setattr__(self, name, value)

    @property
    def designed_sequence(self) -> str:
        """The protein the library encodes, which is ``protein_sequence`` with the first
        ``truncation`` residues dropped. Most libraries leave ``truncation`` at 0, where
        this is the whole protein. Read the name as "the protein being designed" rather
        than "the truncated protein"."""
        if self.truncation and self.truncation >= len(self.protein_sequence):
            raise ValueError(
                f"truncation ({self.truncation}) removes the whole {len(self.protein_sequence)} "
                "aa protein_sequence, leaving nothing to design."
            )
        return self.protein_sequence[self.truncation:]

    @property
    def truncated_sequence(self) -> str:
        """Alias for ``designed_sequence``, kept so existing scripts and notebooks still
        run. The name misleads when ``truncation`` is 0, which is most libraries, so new
        code should say ``designed_sequence``."""
        return self.designed_sequence

    def protein_description(self) -> str:
        """How to name the designed protein in a message, e.g. "protein_sequence" or
        "the truncated protein_sequence (truncation=6)". Truncation is mentioned only when
        the spec truncates, so a spec that does not never reads as if it did."""
        if not self.truncation:
            return "protein_sequence"
        return f"the truncated protein_sequence (truncation={self.truncation})"

    def resolve_vector(
        self, tiled: TiledAssemblyParams | None = None
    ) -> StartingVectorParams | None:
        """The starting vector in play, or None if the spec names no plasmid.

        Spec-level ``starting_vector`` fields merged with a tiled block's, which win
        field by field because a tiled library's vectors are built from the params the
        layout actually ran with. ``tiled`` defaults to ``spec.tiled``; pass the params
        given to ``tile()`` so an explicit ``tile(params)`` resolves what it laid out.

        Two fields are not merged. ``enzyme`` and ``vector_insert`` come from the tiled
        block whenever there is one, even where that block left them at their defaults
        and the spec-level params set them, because those two are what the tiles were
        cut and flanked for.
        """
        base = _as_vector_params(self.starting_vector) or StartingVectorParams()
        t = self.tiled if tiled is None else tiled
        if t is None:
            return base if base.path else None
        path = t.starting_vector or base.path
        if not path:
            return None
        return StartingVectorParams(
            path=path,
            use_vector_cds=t.use_vector_cds or base.use_vector_cds,
            insert_label=t.insert_label if t.insert_label is not None else base.insert_label,
            insert_anchors=tuple(t.insert_anchors) if t.insert_anchors else base.insert_anchors,
            topology=t.topology if t.topology is not None else base.topology,
            enzyme=t.enzyme,               # the enzyme the tiles were actually laid out with
            vector_insert=t.vector_insert,
        )

    @property
    def vector(self) -> StartingVectorParams | None:
        """``resolve_vector()`` with no override, the vector this spec alone describes."""
        return self.resolve_vector()

    @classmethod
    def from_toml(cls, path: str | Path) -> "LibrarySpec":
        """Load a spec from a TOML file.

        Top-level keys map to the spec's fields. An ``[optimization]`` or ``[tiled]``
        table is promoted to its dataclass here, with a TOML array of ``insert_anchors``
        coerced to the tuple the layout code expects. Everything else goes through the
        constructor, so the usual normalization applies (a ``[starting_vector]`` table
        or a bare path becomes a ``StartingVectorParams``, pasted sequences are cleaned)
        and a ``uniprot`` key with no ``protein_sequence`` fetches the entry, exactly as
        it would if you built the spec by hand. See ``examples/mbo038.toml`` for a
        standard scan and ``examples/gck_tiled.toml`` for a tiled one.
        """
        data = tomllib.loads(Path(path).read_text())
        opt = data.get("optimization")
        if isinstance(opt, dict):
            data["optimization"] = CodonOptimizationParams(**opt)
        tiled = data.get("tiled")
        if isinstance(tiled, dict):
            if isinstance(tiled.get("insert_anchors"), list):
                tiled["insert_anchors"] = tuple(tiled["insert_anchors"])
            data["tiled"] = TiledAssemblyParams(**tiled)
        return cls(**data)

    def to_dict(self) -> dict:
        """The spec as plain nested dicts, the ``"spec"`` block of the design-specs record.

        It round-trips. ``__setattr__`` turns the nested dicts back into their
        dataclasses, so ``LibrarySpec(**spec.to_dict())`` rebuilds an equal spec, and so
        does the same call after the dict has been through JSON."""
        return asdict(self)   # recurses into the nested optimization params

    # --- display --------------------------------------------------------
    def _opt_line(self) -> str:
        return optimization_line(asdict(self.optimization))

    def _uniprot_line(self) -> str:
        if not self.uniprot:
            return ""
        e = self.uniprot_entry or {}
        bits = [b for b in (e.get("entry_name"), e.get("organism")) if b]
        sv = f"SV {e['sequence_version']}" if e.get("sequence_version") else ""
        detail = ", ".join(bits + ([sv] if sv else []))
        return f"{self.uniprot}" + (f" ({detail})" if detail else "")

    def _vector_line(self) -> str:
        v = self.vector
        if v is None:
            return "not set"
        how = (
            f"feature {v.insert_label!r}" if v.insert_label
            else "anchors" if v.insert_anchors
            else "the WT CDS" if self.cds
            else "the sole CDS feature"
        )
        frozen = ", CDS frozen from the vector" if v.use_vector_cds else ""
        return f"{Path(v.path).name} ({v.enzyme}, insert located by {how}{frozen})"

    def __repr__(self) -> str:
        trunc = (
            f"  (truncation {self.truncation}, {len(self.designed_sequence)} aa scanned)"
            if self.truncation
            else ""
        )
        uniprot = f"  uniprot:       {self._uniprot_line()}\n" if self.uniprot else ""
        return (
            f"LibrarySpec {self.name!r}\n"
            f"  protein:       {len(self.protein_sequence)} aa  {_preview(self.protein_sequence)}{trunc}\n"
            f"{uniprot}"
            f"  substitutions: {', '.join(self.substitutions)}\n"
            f"  adaptors:      5' {self.adaptor_5}  |  3' {self.adaptor_3}\n"
            f"  vector:        {self._vector_line()}\n"
            f"  optimization:  {self._opt_line()}\n"
            f"  platform:      {self.platform}   max_oligo_length={self.max_oligo_length}   seed={self.seed}"
        )

    def _repr_html_(self) -> str:
        trunc = (
            f" <span style='opacity:.6'>(truncation {self.truncation}, "
            f"{len(self.designed_sequence)} aa scanned)</span>"
            if self.truncation
            else ""
        )
        rows = [
            ("name", _esc(self.name)),
            (
                "protein",
                f"{len(self.protein_sequence)} aa &nbsp; <code>{_esc(_preview(self.protein_sequence))}</code>{trunc}",
            ),
            ("substitutions", " ".join(f"<code>{_esc(s)}</code>" for s in self.substitutions)),
            *(( ("uniprot", _esc(self._uniprot_line())), ) if self.uniprot else ()),
            ("adaptors", f"5' <code>{_esc(self.adaptor_5)}</code> &nbsp; 3' <code>{_esc(self.adaptor_3)}</code>"),
            ("starting_vector", _esc(self._vector_line())),
            (
                "avoid",
                f"{', '.join(map(_esc, self.avoid_enzymes)) or 'no enzymes'} "
                f"<span style='opacity:.6'>+ {len(self.avoid_patterns)} motif pattern(s)</span>",
            ),
            ("optimization", _esc(self._opt_line())),
            ("platform", _esc(self.platform) if self.platform else "not set"),
            ("max_oligo_length", self.max_oligo_length if self.max_oligo_length is not None else "not set"),
            ("seed", self.seed),
        ]
        trs = "".join(
            "<tr>"
            f"<th style='text-align:right;padding:2px 10px;opacity:.7;white-space:nowrap'>{k}</th>"
            f"<td style='text-align:left;padding:2px 10px;font-family:var(--jp-code-font-family,monospace)'>{v}</td>"
            "</tr>"
            for k, v in rows
        )
        return (
            "<table style='border-collapse:collapse'>"
            "<caption style='text-align:left;font-weight:600;padding:2px 10px'>LibrarySpec</caption>"
            f"{trs}</table>"
        )
