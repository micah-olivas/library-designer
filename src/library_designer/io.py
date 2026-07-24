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
from pathlib import Path

import pandas as pd

from .regions import assemble, usortm_sequence

# Characters uSort-M forbids in variant names (file paths, FASTA headers, delimiters).
_BAD_NAME = re.compile(r"[/|>\s]")


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
    if getattr(library, "tiles", None) is not None:
        raise NotImplementedError(
            "to_usortm() does not support tiled libraries. A tiled pool is many "
            "sublibraries, each with its own per-tile reference and length, so it "
            "cannot be written as one uniform name,sequence variable-region block "
            "(emitting the whole CDS per variant would silently corrupt the handoff). "
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


def to_vendor(library, path: str | Path, method: str | None = None,
              pool_name: str | None = None) -> None:
    """Write a synthesis-provider order form. If ``method`` is omitted it is derived
    from ``spec.platform``.

    - ``pooled`` uses the Twist oligo-pool style, with one shared Pool Name and an all-uppercase insert.
    - ``arrayed`` uses one row per named construct (IDT eBlocks or gene fragments).
    """
    _require_complete(library)
    if method is None:
        platform = library.spec.platform
        if platform is None:
            raise ValueError("Specify method='pooled'|'arrayed' or set spec.platform.")
        from .methods import platform_type
        method = platform_type(platform)

    seqs = _assembled(library)   # plain uppercase construct
    if method == "pooled":
        out = pd.DataFrame({
            "Pool Name": pool_name or library.spec.name,
            "Insert Length": seqs.str.len(),
            "Insert Sequence": seqs,
        })
    elif method == "arrayed":
        out = pd.DataFrame({
            "Name": library.df["name"].astype(str),
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
    assembled oligo per placed variant). Variant names are validated as for uSort-M."""
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
    """Manifest of the per-tile Golden Gate destination vectors, one row per tile.

    Each ``vector_sequence`` is the full destination plasmid when ``tiled.starting_vector``
    is set (the backbone with that tile's window dropped out), or the CDS-region cassette
    alone when it is not. The row count is the number of separate destination vectors you
    build to run the library. For annotated plasmid maps use ``to_vector_maps``."""
    _require_tiled(library)
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


def to_vector_maps(library, directory: str | Path) -> None:
    """Write one annotated GenBank map per tile, the destination plasmids to clone.

    Requires ``tiled.starting_vector`` (the backbone supplies the full plasmid to
    annotate). Emits ``tile{i}_destination.gb`` for
    every tile plus a ``destination_vectors.csv`` manifest. Each map carries the two BsaI
    sites, the drop-out, the two fused overhangs, the retained 5'/3' CDS arms, and any
    backbone features from the starting vector that do not overlap the insert."""
    _require_tiled(library)
    params = library.tiled_params
    if not params or not params.starting_vector:
        raise ValueError(
            "to_vector_maps needs tiled.starting_vector set (the destination plasmid to "
            "clone into). Without a backbone, use to_vectors() for the CDS-cassette manifest."
        )
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

    from .checks.motifs import ENZYME_SITES
    from .layout.vector_io import locating_kwargs, resolve_destination
    from .regions import reverse_complement

    dest = resolve_destination(params.starting_vector,
                               **locating_kwargs(library.spec, params))
    reference = library.reference
    stuffer = params.vector_insert.upper()
    site = ENZYME_SITES[params.enzyme].upper()
    site_rc = reverse_complement(site)
    o = params.overhang_len
    name = library.spec.name

    d = Path(directory)
    d.mkdir(parents=True, exist_ok=True)
    manifest = []
    for t in library.tiles:
        V, N = t.vector, len(t.vector)
        circular = t.topology == "circular"
        offset = 0 if circular else dest.start          # where the CDS region sits in V
        s, e, L = t.start, t.end, len(stuffer)
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
        add(cds5, drop_end, "misc_feature", f"{params.enzyme} drop-out (replaced by tile)")
        for m in _find_all(stuffer, site):
            add(cds5 + m, cds5 + m + len(site), "protein_bind", f"{params.enzyme} site", strand=1)
        for m in _find_all(stuffer, site_rc):
            add(cds5 + m, cds5 + m + len(site_rc), "protein_bind", f"{params.enzyme} site", strand=-1)
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
                if circular:
                    return (len(_cds_region(reference, s, e, stuffer)) + (p - dest.end)) % N if p >= dest.end \
                        else (N - dest.start + p)
                return p if p < dest.start else p + (N - len(dest.full_seq))
            a, b = _remap(f.start), _remap(f.end)
            if 0 <= a < b <= N:
                feats.append(SeqFeature(SimpleLocation(a, b, strand=f.strand),
                                        type=f.type, qualifiers={"label": [f.label or f.type]}))
                carried += 1

        rec = SeqRecord(
            Seq(V), id=f"tile{t.index}", name=f"tile{t.index}_dest",
            description=f"{name} tiled destination vector, tile {t.index} (window {s}-{e})",
            annotations={"molecule_type": "DNA", "topology": t.topology},
        )
        rec.features = feats
        SeqIO.write(rec, str(d / f"tile{t.index}_destination.gb"), "genbank")
        manifest.append({
            "tile": t.index, "file": f"tile{t.index}_destination.gb",
            "topology": t.topology, "vector_length": N,
            "fwd_primer_id": t.fwd_id, "rev_primer_id": t.rev_id,
            "features_carried": carried,
        })
    pd.DataFrame(manifest).to_csv(d / "destination_vectors.csv", index=False)


def _cds_region(reference: str, s: int, e: int, stuffer: str) -> str:
    return reference[:s] + stuffer + reference[e:]
