"""Exporters. The construct's parts live as explicit columns (``adaptor_5``,
``variable_dna``, ``adaptor_3``), and each exporter assembles from those. The same
library serializes several ways.

- ``to_full_csv`` writes the master table (the region columns plus the assembled
  uppercase sequence, length, and GC) and may include failed rows for inspection.
- ``to_usortm`` writes the uSort-M input, ``name,sequence`` in the case-encoded
  wire format (flanking lowercase), derived only at this boundary.
- ``to_vendor`` writes the synthesis-provider order form, method-aware (pooled or arrayed).
- ``to_design_specs`` (a Library method) writes the run's spec, params, seed, and versions.
"""
from __future__ import annotations

import re
from datetime import datetime
from pathlib import Path

import pandas as pd

from .regions import assemble, usortm_sequence

# Characters uSort-M forbids in variant names (file paths, FASTA headers, delimiters).
_BAD_NAME = re.compile(r"[/|>\s]")


def run_stamp(when: datetime | None = None) -> str:
    """Local wall-clock stamp for a run directory, sortable and filename-safe."""
    return (when or datetime.now().astimezone()).strftime("%Y%m%d_%H%M%S")


def run_directory(base: str | Path = "out", name: str | None = None,
                  when: datetime | None = None, create: bool = True) -> Path:
    """``base/<name>_<stamp>``, the directory one run's files go in.

    A run writes to its own dated directory so a second run cannot overwrite the first,
    and every file on disk names the run that produced it. ``when`` pins the stamp, so a
    library can name the directory after the moment its sequences were built rather than
    the moment they were written. The directory is created unless ``create=False``. With
    no ``name`` the directory is the bare stamp."""
    stamp = run_stamp(when)
    d = Path(base) / (f"{name}_{stamp}" if name else stamp)
    if create:
        d.mkdir(parents=True, exist_ok=True)
    return d


def _assembled(library) -> pd.Series:
    """Full uppercase construct per row (NA where optimization failed)."""
    return library.df.apply(
        lambda r: assemble(r["adaptor_5"], r["variable_dna"], r["adaptor_3"])
        if isinstance(r["variable_dna"], str) else pd.NA,
        axis=1,
    )


def _usortm_assembled(library) -> pd.Series:
    """Full construct per row in uSort-M wire format (flanking lowercase)."""
    return library.df.apply(
        lambda r: usortm_sequence(r["adaptor_5"], r["variable_dna"], r["adaptor_3"]),
        axis=1,
    )


def _gc_fraction(seq: str) -> float:
    s = seq.upper()
    return (s.count("G") + s.count("C")) / len(s) if s else 0.0


def _require_optimized(library) -> None:
    if "variable_dna" not in library.df.columns:
        raise ValueError("Library is not codon-optimized yet, call codon_optimize() first.")


def _require_complete(library) -> None:
    """As _require_optimized, but also refuse to emit an *order* with failed variants."""
    _require_optimized(library)
    failed = library.df["name"][library.df["variable_dna"].isna()].astype(str).tolist()
    if failed:
        raise ValueError(
            f"{len(failed)} variant(s) failed optimization and can't be exported "
            f"(e.g. {', '.join(failed[:5])}). Inspect library.summary() or call "
            f"library.drop_failed() to exclude them."
        )


def to_full_csv(library, path: str | Path) -> None:
    """Write the master table, every column of ``lib.df`` plus the assembled construct.

    Adds ``sequence`` (``adaptor_5 + variable_dna + adaptor_3``, uppercase), ``length``,
    ``gc_content`` for the variable region alone, and ``stamp_adaptiveness``, the
    relative adaptiveness of the codon this variant carries at its own position. The
    last two are rounded to three places. ``stamp_adaptiveness`` is NA wherever there is
    no stamped position, so for the wild-type control and for every member of a sequence
    set. Those four columns and the region columns are grouped 5'->3' at the right of
    the table.

    Unlike the order-form exporters, this one keeps rows whose optimization failed, with
    a blank sequence, so a failed run can still be read.
    """
    _require_optimized(library)   # failed rows kept (blank sequence) for inspection
    df = library.df.copy()
    df["sequence"] = _assembled(library)
    df["length"] = df["sequence"].map(lambda s: len(s) if isinstance(s, str) else pd.NA)
    df["gc_content"] = df["variable_dna"].map(
        lambda v: round(_gc_fraction(v), 3) if isinstance(v, str) else pd.NA
    )
    # Adaptiveness (w) of each variant's stamped codon, the tabular twin of the QC plot.
    from .optimize.backbone import relative_adaptiveness
    w = relative_adaptiveness(library.spec.optimization.species)
    df["stamp_adaptiveness"] = df.apply(
        lambda r: round(w.get(r["variable_dna"][int(r["mut_index"]) * 3:int(r["mut_index"]) * 3 + 3], 0.0), 3)
        if isinstance(r["variable_dna"], str) and not pd.isna(r["mut_index"]) else pd.NA,
        axis=1,
    )
    # Group the region columns 5'->3' next to the assembled sequence for readability.
    region = ["adaptor_5", "variable_dna", "adaptor_3", "sequence", "length", "gc_content",
              "stamp_adaptiveness"]
    front = [c for c in df.columns if c not in region]
    df = df[front + region]
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, index=False)


def to_usortm(library, path: str | Path) -> None:
    """Write the uSort-M handoff, ``name,sequence``, one row per variant.

    The sequence is case-encoded, adaptors lowercase and the variable region uppercase,
    which is how uSort-M reads the region boundaries. Case is applied here and nowhere
    else, so ``lib.df`` stays uppercase throughout.

    Refuses a library with a failed variant, since an order should not go out
    incomplete, and refuses names carrying ``/``, ``|``, ``>``, or whitespace, which
    uSort-M forbids. A tiled library raises ``NotImplementedError``, because a tiled
    pool is many sublibraries with their own per-tile references and so has no one
    variable region to write. Use ``to_oligo_pool`` for the physical order.
    """
    if getattr(library, "tiles", None) is not None:
        raise NotImplementedError(
            "to_usortm() does not support tiled libraries. A tiled pool is many "
            "sublibraries, each with its own per-tile reference and length, so it "
            "cannot be written as one uniform name,sequence variable-region block "
            "(emitting the whole CDS per variant would produce a wrong handoff). "
            "The per-tile uSort-M convention (grouping on Tile_N) is still being "
            "settled on the uSort-M side; use to_oligo_pool() for the physical pooled "
            "order in the meantime."
        )
    _require_complete(library)
    names = library.df["name"].astype(str)
    bad = names[names.str.contains(_BAD_NAME)]
    if not bad.empty:
        raise ValueError(
            "Variant names contain characters uSort-M forbids (/ | > whitespace): "
            + ", ".join(bad.head(10))
        )
    out = pd.DataFrame({"name": names, "sequence": _usortm_assembled(library)})
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(path, index=False)

    # Record what was handed off, next to the run identity, so the downstream tool can tie
    # the CSV it is reading to the run that produced it. The CSV itself stays exactly
    # `name,sequence`: uSort-M parses it strictly, so the provenance goes in the record
    # beside it rather than as extra columns or a comment line it would have to tolerate.
    library.design_specs["handoff"] = {
        "run_id": library.run_id,
        "created": library.created,
        "variants_csv": Path(path).name,
        "n_variants": len(out),
        "sha256": _sha256(path),
        "format": "name,sequence; flanking adaptors lowercase, variable region uppercase",
    }


def _sha256(path: str | Path) -> str:
    """Digest of a written file, so a later reader can tell whether it is the one the
    record describes (hand-edited, truncated, or mixed up with another run's)."""
    import hashlib

    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def to_vendor(library, path: str | Path, method: str | None = None,
              pool_name: str | None = None) -> None:
    """Write a synthesis-provider order form. If ``method`` is omitted it is derived
    from ``spec.platform``.

    - ``pooled`` uses the Twist oligo-pool style, with one shared Pool Name and an all-uppercase insert.
    - ``arrayed`` uses one row per named construct (IDT eBlocks or gene fragments).

    So the pooled form carries no variant names, only the one shared pool name taken from
    ``pool_name`` or ``spec.name``, while the arrayed form names every row. Any other
    ``method`` raises.

    For a tiled library the molecules to order are the tile oligos, not the full coding
    sequence each variant carries, so those are what goes on the form. That includes the
    per-tile WT controls; the global ``WT`` row rides on no oligo and is left off, as it is
    in ``to_oligo_pool``.

    Refuses a library with a variant that failed optimization.
    """
    _require_complete(library)
    if method is None:
        platform = library.spec.platform
        if platform is None:
            raise ValueError("Specify method='pooled'|'arrayed' or set spec.platform.")
        from .methods import platform_type
        method = platform_type(platform)

    if getattr(library, "tiles", None) is not None:
        placed = library.df["oligo"].map(lambda o: isinstance(o, str))
        seqs = library.df["oligo"][placed]
        names = library.df["name"][placed].astype(str)
    else:
        seqs = _assembled(library)   # plain uppercase construct
        names = library.df["name"].astype(str)
    if method == "pooled":
        out = pd.DataFrame({
            "Pool Name": pool_name or library.spec.name,
            "Insert Length": seqs.str.len(),
            "Insert Sequence": seqs,
        })
    elif method == "arrayed":
        out = pd.DataFrame({
            "Name": names,
            "Insert Length": seqs.str.len(),
            "Insert Sequence": seqs,
        })
    else:
        raise ValueError(f"Unknown method {method!r} (use 'pooled' or 'arrayed').")
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(path, index=False)


# --- tiled assembly -----------------------------------------------------------

def _require_tiled(library) -> None:
    if getattr(library, "tiles", None) is None:
        raise ValueError("Library is not tiled yet, call tile() first.")


def to_oligo_pool(library, path: str | Path) -> None:
    """The single pooled synthesis order for a tiled library: ``name,sequence`` (one
    assembled oligo per placed member, mutants and the per-tile ``WT_Tile_<i>`` controls
    alike). The global ``WT`` row rides on no oligo, so it is left off. Names are validated
    as for uSort-M."""
    _require_tiled(library)
    df = library.df
    mask = df["oligo"].map(lambda o: isinstance(o, str))
    names = df["name"][mask].astype(str)
    bad = names[names.str.contains(_BAD_NAME)]
    if not bad.empty:
        raise ValueError(
            "Variant names contain characters uSort-M forbids (/ | > whitespace): "
            + ", ".join(bad.head(10))
        )
    out = pd.DataFrame({"name": names, "sequence": df["oligo"][mask]})
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(path, index=False)


def _coding_slice(library) -> dict[str, tuple[str, int | None]]:
    """Per member, the coding stretch its oligo carries and the residue index it mutates.

    A standard library's oligo holds the whole CDS. A tiled library's holds only its tile
    window, and the mutated codon is counted from the start of that window, so a tiled
    member's index is rebased. ``None`` means nothing is mutated (a wild-type control).
    """
    df = library.df
    tiles = getattr(library, "tiles", None)
    idx = df["mut_index"] if "mut_index" in df.columns else [pd.NA] * len(df)
    out: dict[str, tuple[str, int | None]] = {}
    if tiles is None:
        for name, dna, i in zip(df["name"], df["variable_dna"], idx):
            if isinstance(dna, str):
                out[str(name)] = (dna, None if pd.isna(i) else int(i))
        return out
    spans = {t.index: (t.start, t.end) for t in tiles}
    for name, dna, i, ti in zip(df["name"], df["variable_dna"], idx, df["tile"]):
        if not isinstance(dna, str) or pd.isna(ti) or int(ti) not in spans:
            continue
        s, e = spans[int(ti)]
        rebased = None if pd.isna(i) else int(i) - s // 3
        out[str(name)] = (dna[s:e], rebased)
    return out


def _oligo_record(name: str, seq: str, coding: str, mut_index: int | None,
                  enzyme: str | None, label: str):
    """One oligo as an annotated SeqRecord: its coding stretch, the mutated codon, and any
    Type IIS recognition sites, on either strand.

    The sites are annotated where they are found rather than where they were meant to be, so
    a stray site inside the coding region shows up in the map as readily as the intended ones
    in the adaptors."""
    from Bio.Seq import Seq
    from Bio.SeqFeature import SeqFeature, SimpleLocation
    from Bio.SeqRecord import SeqRecord

    from .checks.motifs import recognition_site
    from .regions import reverse_complement

    feats = []
    at = seq.find(coding) if coding else -1
    if at >= 0:
        feats.append(SeqFeature(SimpleLocation(at, at + len(coding), strand=1), type="CDS",
                                qualifiers={"label": [label], "codon_start": ["1"]}))
        if mut_index is not None and 0 <= mut_index * 3 < len(coding):
            lo = at + mut_index * 3
            feats.append(SeqFeature(
                SimpleLocation(lo, lo + 3, strand=1), type="variation",
                qualifiers={"label": [f"{name} ({seq[lo:lo + 3]})"]},
            ))
    if enzyme:
        site = recognition_site(enzyme)
        for pattern, strand in ((site, 1), (reverse_complement(site), -1)):
            if strand == -1 and pattern == site:
                continue                      # a palindromic site is not two features
            for i in _find_all(seq, pattern):
                feats.append(SeqFeature(SimpleLocation(i, i + len(pattern), strand=strand),
                                        type="protein_bind",
                                        qualifiers={"label": [enzyme]}))
    feats.sort(key=lambda f: int(f.location.start))
    rec = SeqRecord(Seq(seq), id=_file_stem(name)[:16], name=_file_stem(name)[:16],
                    description=f"{label} {name} synthesis oligo",
                    annotations={"molecule_type": "DNA", "topology": "linear"})
    rec.features = feats
    return rec


def to_oligo_files(library, directory: str | Path, fmt: str = "genbank") -> int:
    """Write one file per ordered oligo into ``directory``, and return how many were written.

    One file per library member, named for the variant (``K7stop.gb`` for an amber ``K7*``,
    since ``*`` globs in a shell). The sequence is the one the order carries, so a tiled
    library writes its assembled ``oligo`` (primers, enzyme sites, tile window, and the
    per-tile ``WT_Tile_<i>`` controls) and any other library writes the whole construct with
    its adaptors. Uppercase throughout, unlike uSort-M's ``variants.csv``, which lowercases
    the flanks to mark them.

    ``fmt`` is ``"genbank"`` (the default), ``"fasta"``, or ``"both"``. A GenBank file is
    annotated with the coding stretch the oligo carries, the mutated codon, and every Type
    IIS recognition site on either strand, so the oligo can be read in a plasmid editor
    without working out the parts by eye. FASTA is sequence only, for aligners and anything
    that wants one accession per file.

    A tiled library's global ``WT`` row rides on no oligo, so it gets no file, as in
    ``to_oligo_pool``. Two variant names that would land on the same filename are refused
    rather than silently overwriting each other."""
    if fmt not in ("genbank", "fasta", "both"):
        raise ValueError(f"Unknown fmt {fmt!r} (use 'genbank', 'fasta', or 'both').")
    _require_optimized(library)
    df = library.df
    tiled = getattr(library, "tiles", None) is not None
    seqs = df["oligo"] if tiled else _assembled(library)

    pairs = [(str(n), s.upper()) for n, s in zip(df["name"], seqs) if isinstance(s, str)]
    if not pairs:
        raise ValueError(
            "No sequences to write. Optimize the library first, and for a tiled library "
            "call tile() so each member has an oligo."
        )
    # _file_stem rewrites '*' and anything the filesystem dislikes, so two names can arrive
    # at one filename. Overwriting would drop a variant from an order.
    claimed: dict[str, str] = {}
    for name, _ in pairs:
        stem = _file_stem(name)
        if stem in claimed:
            raise ValueError(
                f"Variants {claimed[stem]!r} and {name!r} both map to the same filename "
                f"({stem}). Rename one of them."
            )
        claimed[stem] = name

    d = Path(directory)
    d.mkdir(parents=True, exist_ok=True)
    if fmt in ("fasta", "both"):
        for name, seq in pairs:
            (d / f"{_file_stem(name)}.fasta").write_text(f">{name}\n{seq}\n")
    if fmt in ("genbank", "both"):
        try:
            from Bio import SeqIO
        except ImportError as exc:
            raise ImportError(
                "Writing GenBank oligos needs BioPython, a base dependency. "
                "Reinstall library-designer if it is missing."
            ) from exc

        spec = library.spec
        vec = spec.resolve_vector(library.tiled_params)
        enzyme = (library.tiled_params.enzyme if tiled
                  else vec.enzyme if vec is not None
                  else next(iter(spec.avoid_enzymes), None))
        coding = _coding_slice(library)
        for name, seq in pairs:
            cds, mut = coding.get(name, ("", None))
            SeqIO.write(_oligo_record(name, seq, cds, mut, enzyme, spec.name),
                        str(d / f"{_file_stem(name)}.gb"), "genbank")
    return len(pairs)


def to_primer_order(library, path: str | Path,
                    scale: str = "25nm", purification: str = "STD") -> None:
    """Per-tile amplification primers as a headerless IDT bulk order
    (``name,sequence,scale,purification``). Sequences are 5'->3' as ordered."""
    _require_tiled(library)
    name = library.spec.name
    rows = []
    for t in library.tiles:
        rows.append((f"fwd_Tile_{t.index}_{name}", t.fwd, scale, purification))
        rows.append((f"rev_Tile_{t.index}_{name}", t.rev, scale, purification))
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(path, index=False, header=False)


def to_vectors(library, path: str | Path) -> None:
    """Manifest of the Golden Gate destination vectors you build to run the library.

    A tiled library gets one row per tile: each ``vector_sequence`` is the full destination
    plasmid when a starting vector is set (the backbone with that tile's window dropped
    out), or the CDS-region cassette alone when it is not. A standard library gets the one
    row it needs, the plasmid with the CDS replaced by the drop-out, plus the two fused
    overhangs the cut vector presents. For annotated plasmid maps use ``to_vector_maps``."""
    if getattr(library, "tiles", None) is None:
        from .layout.destination import build_destination

        dv = build_destination(library)
        from .layout.vector_io import origin_in_source
        out = pd.DataFrame([{
            "topology": dv.topology,
            "cds_dropout_start": dv.start,       # window on the reference CDS, not on the vector
            "cds_dropout_end": dv.end,
            "overhang_5": dv.overhang_5,
            "overhang_3": dv.overhang_3,
            "vector_length": dv.length,
            # The sequence is read from an origin upstream of the insert's promoter, so name
            # the base of the starting vector it begins at.
            "origin_in_starting_vector": origin_in_source(dv.dest),
            "insert_strand_in_starting_vector": -1 if dv.dest.flipped else 1,
            "vector_sequence": dv.sequence,
        }])
        Path(path).parent.mkdir(parents=True, exist_ok=True)
        out.to_csv(path, index=False)
        return
    rows = [
        {
            "tile": t.index, "start": t.start, "end": t.end,
            "topology": t.topology,
            "fwd_primer_id": t.fwd_id, "rev_primer_id": t.rev_id,
            "vector_length": len(t.vector),
            "vector_sequence": t.vector,
        }
        for t in library.tiles
    ]
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(path, index=False)


def _find_all(hay: str, needle: str) -> list[int]:
    out, i = [], hay.find(needle)
    while i != -1:
        out.append(i)
        i = hay.find(needle, i + 1)
    return out


def _map_targets(library) -> list[dict]:
    """The destination vectors to emit maps for, one per tile or the single one a standard
    library needs, each reduced to what the annotator below draws on."""
    tiles = getattr(library, "tiles", None)
    if tiles is None:
        from .layout.destination import build_destination

        dv = build_destination(library)
        return [{
            "file": "destination.gb", "id": "destination", "label": "CDS drop-out",
            "start": dv.start, "end": dv.end, "vector": dv.sequence, "topology": dv.topology,
            "extra": {},
        }]
    return [{
        "file": f"tile{t.index}_destination.gb", "id": f"tile{t.index}",
        "label": f"tile {t.index} (window {t.start}-{t.end})",
        "start": t.start, "end": t.end, "vector": t.vector, "topology": t.topology,
        "extra": {"tile": t.index, "fwd_primer_id": t.fwd_id, "rev_primer_id": t.rev_id},
    } for t in tiles]


def to_vector_maps(library, directory: str | Path) -> None:
    """Write the destination plasmids to clone, as annotated GenBank maps.

    Needs a starting vector (the backbone supplies the full plasmid to annotate). A tiled
    library gets one ``tile{i}_destination.gb`` per tile; a standard library gets a single
    ``destination.gb``. Both come with a ``destination_vectors.csv`` manifest. Each map
    carries the two BsaI sites, the drop-out, the two fused overhangs, the retained 5'/3'
    CDS arms, and any backbone features from the starting vector that do not overlap the
    insert.

    Two kinds of backbone feature are dropped. One overlapping the located insert goes,
    since the insert replaces what it annotated, and so does one straddling the origin
    the map is read from, which GenBank cannot write in one piece. The manifest's
    ``features_carried`` column says how many survived."""
    params = library.tiled_params
    vec = library.spec.resolve_vector(params)
    if vec is None:
        raise ValueError(
            "to_vector_maps needs a starting_vector set (the destination plasmid to clone "
            "into). Without a backbone, use to_vectors() for the CDS-cassette manifest."
        )
    if getattr(library, "tiles", None) is None and library.reference is None:
        raise ValueError("Library is not codon-optimized yet, call codon_optimize() first.")
    try:
        from Bio import SeqIO
        from Bio.Seq import Seq
        from Bio.SeqFeature import SeqFeature, SimpleLocation
        from Bio.SeqRecord import SeqRecord
    except ImportError as exc:
        raise ImportError(
            "Writing GenBank destination maps needs BioPython, a base dependency. "
            "Reinstall library-designer if it is missing."
        ) from exc

    from .checks.motifs import ENZYME_SITES, cut_geometry
    from .layout.destination import resolve_insert_locus
    from .layout.vector_io import insert_offset, origin_in_source
    from .regions import reverse_complement

    dest = resolve_insert_locus(library.spec, params)
    reference = library.reference
    enzyme = vec.enzyme
    stuffer = vec.vector_insert.upper()
    site = ENZYME_SITES[enzyme].upper()
    site_rc = reverse_complement(site)
    o = params.overhang_len if params is not None else cut_geometry(enzyme)[1]
    name = library.spec.name

    d = Path(directory)
    d.mkdir(parents=True, exist_ok=True)
    manifest = []
    for target in _map_targets(library):
        V, N = target["vector"], len(target["vector"])
        circular = target["topology"] == "circular"
        offset = insert_offset(dest)                    # where the CDS region sits in V
        s, e, L = target["start"], target["end"], len(stuffer)
        region_len = len(_cds_region(reference, s, e, stuffer))
        cds5 = offset + s                                 # end of the retained 5' arm / start of dropout
        drop_end = cds5 + L                               # end of the dropout / start of retained 3' arm

        feats: list = []

        def add(a: int, b: int, ftype: str, label: str, strand: int | None = None):
            a %= N
            b_mod = b % N if b != N else N
            if a < b_mod:                                 # skip anything that would wrap the origin
                feats.append(SeqFeature(SimpleLocation(a, b_mod, strand=strand),
                                        type=ftype, qualifiers={"label": [label]}))

        if s > 0:
            add(offset, cds5, "misc_feature", "CDS (5' arm, retained)")
        if e < len(reference):
            add(drop_end, drop_end + (len(reference) - e), "misc_feature", "CDS (3' arm, retained)")
        add(cds5, drop_end, "misc_feature", f"{enzyme} drop-out (replaced by the insert)")
        for m in _find_all(stuffer, site):
            add(cds5 + m, cds5 + m + len(site), "protein_bind", f"{enzyme} site", strand=1)
        for m in _find_all(stuffer, site_rc):
            add(cds5 + m, cds5 + m + len(site_rc), "protein_bind", f"{enzyme} site", strand=-1)
        add(cds5 - o, cds5, "misc_feature", "fused overhang 5'")
        add(drop_end, drop_end + o, "misc_feature", "fused overhang 3'")

        # Carry over backbone features that neither overlap the located insert nor cross
        # the file origin, remapping their coordinates into the emitted vector.
        carried = 0
        for f in dest.features:
            if f.end <= f.start or f.start >= f.end:
                continue
            if not (f.end <= dest.start or f.start >= dest.end):
                continue                                  # overlaps the insert, drop it
            def _remap(p: int) -> int:
                # Into the emitted frame: the two backbone arms sit either side of the CDS
                # region, which starts at `offset` (the rotation the origin implies).
                if circular:
                    if p < dest.start:
                        return (offset - (dest.start - p)) % N
                    return (offset + region_len + (p - dest.end)) % N
                return p if p < dest.start else p + (N - len(dest.full_seq))
            # Take the end from the length, not from _remap: a feature finishing on the
            # vector's last base has an exclusive end of N, which the wrap would fold to 0
            # and silently drop the feature. A feature that straddles the origin overshoots
            # N here and is skipped, which is right, GenBank cannot draw it in one piece.
            a = _remap(f.start)
            b = a + (f.end - f.start)
            if 0 <= a < b <= N:
                feats.append(SeqFeature(SimpleLocation(a, b, strand=f.strand),
                                        type=f.type, qualifiers={"label": [f.label or f.type]}))
                carried += 1

        # The map of a circular plasmid is read from an origin upstream of the insert's
        # promoter, so say which base of the starting vector that was. Without it the
        # emitted coordinates cannot be related back to the plasmid you supplied.
        origin_base = origin_in_source(dest)
        strand = " on the minus strand" if dest.flipped else ""
        read_from = (f", read from base {origin_base} of {Path(vec.path).name}{strand}"
                     if circular or dest.flipped else "")
        rec = SeqRecord(
            Seq(V), id=target["id"], name=f"{target['id']}_dest",
            description=f"{name} destination vector, {target['label']}{read_from}",
            annotations={"molecule_type": "DNA", "topology": target["topology"]},
        )
        rec.features = feats
        SeqIO.write(rec, str(d / target["file"]), "genbank")
        manifest.append({
            "vector": target["id"], **target["extra"], "file": target["file"],
            "topology": target["topology"], "vector_length": N,
            "origin_in_starting_vector": origin_base,
            "features_carried": carried,
        })
    pd.DataFrame(manifest).to_csv(d / "destination_vectors.csv", index=False)


def _cds_region(reference: str, s: int, e: int, stuffer: str) -> str:
    return reference[:s] + stuffer + reference[e:]


_UNSAFE = re.compile(r"[^A-Za-z0-9._-]")


def _file_stem(name: str) -> str:
    """A variant name as a filename. An amber stop is ``*``, which is legal on this
    filesystem but globs in a shell and is rejected on Windows, so spell it out."""
    return _UNSAFE.sub("_", name.replace("*", "stop"))


def to_assembled_vectors(library, directory: str | Path, fmt: str = "both") -> None:
    """Write the clone every variant assembles into: the simulated Golden Gate product.

    These are the plasmids to sequence against. One annotated GenBank per variant (the
    backbone features carried over, the coding sequence, and the mutated codon marked), a
    single ``all_clones.fasta`` holding the same set, the parent plasmid as
    ``parent_WT.gb`` for reference, and an ``assembled_vectors.csv`` manifest recording each
    clone's file, length, its codon change, and whether it differs from the parent only where
    it should (``as_intended``, which for the wild-type control means not at all).

    ``fmt`` is ``"genbank"``, ``"fasta"``, or ``"both"``. Needs a starting vector and
    adaptors that release the insert; a variant the reaction cannot make (see ``check()``) is
    left out, and the number skipped is raised as a warning."""
    if fmt not in ("genbank", "fasta", "both"):
        raise ValueError(f"Unknown fmt {fmt!r} (use 'genbank', 'fasta', or 'both').")
    _require_optimized(library)
    try:
        from Bio import SeqIO
        from Bio.Seq import Seq
        from Bio.SeqFeature import SeqFeature, SimpleLocation
        from Bio.SeqRecord import SeqRecord
    except ImportError as exc:
        raise ImportError(
            "Writing assembled clones needs BioPython, a base dependency. "
            "Reinstall library-designer if it is missing."
        ) from exc

    from .checks.assembly import assembled_products, parent_vector
    from .layout.destination import resolve_insert_locus
    from .layout.vector_io import insert_offset, origin_in_source

    spec = library.spec
    vec = spec.resolve_vector(library.tiled_params)
    parent = parent_vector(library)
    dest = resolve_insert_locus(spec, library.tiled_params)
    reference = library.reference
    locus_at, n = insert_offset(dest), len(parent)
    # The locus holds the reference plus whatever bases the adaptors re-supply, and the coding
    # region starts past the 5' padding. Both numbers are needed: features are remapped around
    # the whole locus, while the CDS and the mutated codon are placed from the coding start.
    locus_len = dest.pad_5 + len(reference) + dest.pad_3
    cds_at = locus_at + dest.pad_5
    source = Path(vec.path).name

    # Every clone is the same molecule apart from one codon, so the backbone annotation is
    # worked out once and reused. Coordinates are the parent's frame, which is the frame
    # assembled_products returns its products in.
    # The clone carries the reference where the plasmid carried its own insert, so a feature
    # 3' of the locus shifts by the difference in their lengths. Getting this from the located
    # region instead would misplace every downstream annotation whenever the two differ.
    def _remap(p: int) -> int:
        if p < dest.start:
            return p if dest.topology != "circular" else (locus_at - (dest.start - p)) % n
        downstream = locus_len + (p - dest.end)
        return (locus_at + downstream) % n if dest.topology == "circular" else dest.start + downstream

    shared: list = [SeqFeature(SimpleLocation(cds_at, cds_at + len(reference), strand=1),
                               type="CDS", qualifiers={"label": [spec.name], "codon_start": ["1"]})]
    for f in dest.features:
        if f.end <= f.start or not (f.end <= dest.start or f.start >= dest.end):
            continue                                     # overlaps the insert, the CDS covers it
        a = _remap(f.start)
        b = a + (f.end - f.start)
        if 0 <= a < b <= n:
            shared.append(SeqFeature(SimpleLocation(a, b, strand=f.strand), type=f.type,
                                     qualifiers={"label": [f.label or f.type]}))

    d = Path(directory)
    d.mkdir(parents=True, exist_ok=True)
    genbank, fasta = fmt in ("genbank", "both"), fmt in ("fasta", "both")

    def _record(name: str, seq: str, feats: list, note: str) -> "SeqRecord":
        rec = SeqRecord(Seq(seq), id=_file_stem(name)[:16], name=_file_stem(name)[:16],
                        description=f"{spec.name} {name} assembled into {source}: {note}",
                        annotations={"molecule_type": "DNA", "topology": dest.topology})
        rec.features = feats
        return rec

    if genbank:
        SeqIO.write(_record("WT_parent", parent, list(shared), "the parent plasmid, reference CDS"),
                    str(d / "parent_WT.gb"), "genbank")

    rows, records, written = [], [], 0
    for name, product, mut_index in assembled_products(library):
        diffs = [i for i, (a, b) in enumerate(zip(product, parent)) if a != b]
        lo, hi = (min(diffs), max(diffs) + 1) if diffs else (cds_at, cds_at)
        codon_lo = cds_at + mut_index * 3 if mut_index is not None else None
        feats = list(shared)
        if codon_lo is not None:
            wt_codon = parent[codon_lo:codon_lo + 3]
            new_codon = product[codon_lo:codon_lo + 3]
            feats.append(SeqFeature(
                SimpleLocation(codon_lo, codon_lo + 3, strand=1), type="variation",
                qualifiers={"label": [f"{name} ({wt_codon}>{new_codon})"]},
            ))
        else:
            wt_codon = new_codon = ""
        stem = _file_stem(name)
        if genbank:
            note = f"{wt_codon}>{new_codon} at codon {mut_index + 1}" if codon_lo is not None \
                else "wild-type control"
            SeqIO.write(_record(name, product, feats, note), str(d / f"{stem}.gb"), "genbank")
        if fasta:
            records.append(_record(name, product, [], ""))
        # For a variant, "as intended" means the only difference sits in its own codon. For
        # the wild-type control it means no difference at all.
        as_intended = (not diffs) if codon_lo is None else (
            bool(diffs) and lo >= codon_lo and hi <= codon_lo + 3
        )
        rows.append({
            "name": name, "file": f"{stem}.gb" if genbank else "",
            "length": len(product), "topology": dest.topology,
            "codon": (mut_index + 1) if mut_index is not None else pd.NA,
            "wt_codon": wt_codon, "clone_codon": new_codon,
            "bases_changed": len(diffs),
            "changed_at": f"{lo + 1}-{hi}" if diffs else "",
            "as_intended": as_intended,
        })
        written += 1

    if fasta:
        SeqIO.write(records, str(d / "all_clones.fasta"), "fasta")
    manifest = pd.DataFrame(rows)
    if not manifest.empty:
        manifest["codon"] = manifest["codon"].astype("Int64")   # no float codon numbers
    manifest.to_csv(d / "assembled_vectors.csv", index=False)
    skipped = int(library.df["variable_dna"].notna().sum()) - written
    if skipped > 0:
        import warnings

        warnings.warn(
            f"{skipped} variant(s) were left out of {d}: the reaction cannot make them. "
            "Run library.check() to see why.",
            RuntimeWarning, stacklevel=2,
        )
