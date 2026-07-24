"""Tiled-assembly layout.

Given a frozen WT reference and its stamped variants (see ``optimize/backbone.py``),
split the CDS into overlapping tile windows sized to the oligo budget, assign each
variant to the tile that contains its codon, and assemble the final synthesis oligo:

    fwd · SITE · spacer5 · [overhang5] · TILE · [overhang3] · spacer3 · SITE_rc · revcomp(rev)

The two ``SITE`` copies are the Golden Gate (Type IIS) recognition sequences; after
digestion the two overhang regions (drawn from the flanking WT CDS) become the fused
overhangs that ligate the tile into a per-tile destination vector carrying the rest of
the WT CDS. Primers come from a validated orthogonal set (see ``primers.py``) so each
tile's sublibrary can be selectively amplified out of the shared pool.

Two site copies are the *only* recognition sites allowed in the oligo. Because adding
the flanking primers/sites/spacers can create an unintended site spanning a junction,
assembly is screened on the *fully assembled* oligo, and pooled primers that would form
a junction site are skipped during assignment.
"""
from __future__ import annotations

from dataclasses import dataclass

import pandas as pd

from ..checks.motifs import ENZYME_SITES
from ..primers import PrimerSet, load_primer_set
from ..regions import reverse_complement
from ..spec import TiledAssemblyParams


@dataclass
class TileInfo:
    index: int
    start: int          # tile window [start, end) on the reference CDS
    end: int
    fwd_id: str
    fwd: str            # forward primer, 5'->3' as ordered (sits verbatim at the oligo 5' end)
    rev_id: str
    rev: str            # reverse primer, 5'->3' as ordered (its reverse-complement sits at the 3' end)
    vector: str         # destination vector: the backbone (or, with no starting vector, the CDS
                        # region alone) with this window replaced by the drop-out insert
    topology: str = "linear"   # "circular" for a real plasmid backbone, else "linear"


def _overhead(params: TiledAssemblyParams, rec_len: int) -> int:
    """Non-tile bases on each oligo: two primers, two sites, two spacers, two overhangs."""
    return (
        2 * params.primer_length
        + 2 * rec_len
        + len(params.spacer_5)
        + len(params.spacer_3)
        + 2 * params.overhang_len
    )


def compute_tiles(cds_len: int, params: TiledAssemblyParams) -> list[tuple[int, int]]:
    """Contiguous, gap-free tile windows covering the whole CDS.

    Boundaries fall on codon multiples so a mutated codon is never split across two
    tiles, and windows are balanced to the fewest tiles whose size fits the oligo
    budget. The 4-bp Golden Gate overhang each oligo carries is drawn from the
    flanking WT CDS, so adjacent tiles need no explicit sequence overlap.
    """
    if cds_len % 3 != 0:
        raise ValueError(f"Reference CDS length ({cds_len}) is not a multiple of 3.")
    rec_len = len(ENZYME_SITES[params.enzyme])
    max_tile = params.tile_size or (params.oligo_budget - _overhead(params, rec_len))
    max_tile_codons = max_tile // 3
    if max_tile_codons < 1:
        raise ValueError(
            f"Oligo budget ({params.oligo_budget} bp) leaves no room for a tile after "
            f"primers and sites (needs > {_overhead(params, rec_len)} bp)."
        )

    codons = cds_len // 3
    tile_n = -(-codons // max_tile_codons)          # ceil: fewest tiles that fit the budget
    per = -(-codons // tile_n)                       # ceil: balance codons across tiles
    tiles = []
    for i in range(tile_n):
        s = i * per * 3
        if s >= cds_len:
            break
        tiles.append((s, min(cds_len, (i + 1) * per * 3)))

    covered = bytearray(cds_len)
    for s, e in tiles:
        covered[s:e] = b"\x01" * (e - s)
    if not all(covered):
        raise ValueError(f"Tiling leaves CDS position {covered.index(0)} uncovered (internal error).")
    return tiles


def assign_tile(mut_index: int, tiles: list[tuple[int, int]]) -> int | None:
    """Index of the first tile whose window fully contains the codon at ``mut_index``
    (0-based residue), or None if the codon is split across every tile boundary."""
    cs, ce = mut_index * 3, mut_index * 3 + 3
    for i, (s, e) in enumerate(tiles):
        if s <= cs and ce <= e:
            return i
    return None


def tile_contexts(reference: str, start: int, end: int,
                  params: TiledAssemblyParams) -> tuple[str, str]:
    """The two Golden Gate overhang sequences for a tile, flanking WT CDS bases, or
    the configured vector context when the tile abuts a CDS end.

    QC uses this too, so the overhangs it checks are the ones the oligos carry."""
    o = params.overhang_len
    ctx5 = reference[start - o:start] if start >= o else params.vector_context_5
    ctx3 = reference[end:end + o] if end + o <= len(reference) else params.vector_context_3
    return ctx5, ctx3


def _check_primer_length(primer_set: PrimerSet, params: TiledAssemblyParams) -> None:
    """``params.primer_length`` is what sizes the tiles against the oligo budget (see
    ``_overhead``), so a set whose primers are longer than that would push every oligo
    past the budget. Refuse up front rather than emit an unorderable pool."""
    lengths = [len(s) for _, s in primer_set.primers]
    lengths += [len(s) for _, fwd, rev in primer_set.pairs for s in (fwd, rev)]
    longest = max(lengths, default=0)
    if longest > params.primer_length:
        raise ValueError(
            f"Primer set {primer_set.name!r} contains a {longest} bp primer but "
            f"tiled.primer_length is {params.primer_length}, and that is what sizes the "
            "tiles against the oligo budget. Set tiled.primer_length to the real primer "
            "length so the budget accounts for it."
        )


def assemble_oligo(reference: str, variant_cds: str, start: int, end: int,
                   fwd: str, rev: str, params: TiledAssemblyParams) -> str:
    """The full synthesis oligo for one variant/tile (all uppercase)."""
    rec = ENZYME_SITES[params.enzyme].upper()
    rec_rc = reverse_complement(rec)
    ctx5, ctx3 = tile_contexts(reference, start, end, params)
    tile = variant_cds[start:end]
    return (
        fwd + rec + params.spacer_5 + ctx5
        + tile
        + ctx3 + params.spacer_3 + rec_rc + reverse_complement(rev)
    )


def site_positions(oligo: str, enzyme: str) -> tuple[list[int], list[int]]:
    """(forward-strand hits, reverse-strand hits) of ``enzyme``'s recognition site."""
    rec = ENZYME_SITES[enzyme].upper()
    rec_rc = reverse_complement(rec)

    def _find(sub: str) -> list[int]:
        return [i for i in range(len(oligo) - len(sub) + 1) if oligo[i:i + len(sub)] == sub]

    return _find(rec), _find(rec_rc)


def extra_sites(oligo: str, fwd_len: int, rev_len: int, enzyme: str) -> bool:
    """True if the oligo carries any recognition site beyond the two intended ones
    (forward site right after the fwd primer; reverse site right before the rev primer)."""
    rec_len = len(ENZYME_SITES[enzyme])
    fwd_hits, rev_hits = site_positions(oligo, enzyme)
    want_fwd = fwd_len
    want_rev = len(oligo) - rev_len - rec_len
    return fwd_hits != [want_fwd] or rev_hits != [want_rev]


def _count_sites(s: str, enzyme: str) -> tuple[int, int]:
    fwd, rev = site_positions(s, enzyme)
    return len(fwd), len(rev)


def _end5_ok(fwd: str, ctx5: str, params: TiledAssemblyParams) -> bool:
    """The forward primer must add no site to the 5' region (fwd|site|spacer|overhang)
    beyond the one intended forward site."""
    rec = ENZYME_SITES[params.enzyme].upper()
    nf, nr = _count_sites(fwd + rec + params.spacer_5 + ctx5, params.enzyme)
    return nf == 1 and nr == 0


def _end3_ok(rev: str, ctx3: str, params: TiledAssemblyParams) -> bool:
    """The reverse primer must add no site to the 3' region (overhang|spacer|site|rev)
    beyond the one intended reverse site."""
    rec_rc = reverse_complement(ENZYME_SITES[params.enzyme].upper())
    nf, nr = _count_sites(ctx3 + params.spacer_3 + rec_rc + reverse_complement(rev), params.enzyme)
    return nf == 0 and nr == 1


def _assign_primers(tiles: list[tuple[int, int]], reference: str,
                    primer_set: PrimerSet, params: TiledAssemblyParams) -> list[tuple[str, str, str, str]]:
    """Choose (fwd_id, fwd, rev_id, rev) per tile. A pool is drawn primer-by-primer,
    skipping any that would form a junction site; a paired set is used as given.

    The 5' and 3' junctions are independent (each involves only its own primer, site,
    spacer, and overhang), so the two ends are validated separately."""
    n = len(tiles)
    if primer_set.capacity < n:
        raise ValueError(
            f"Primer set {primer_set.name!r} supplies {primer_set.capacity} tile(s) of primers "
            f"but the CDS needs {n}. Provide a larger set (tiled.primer_set=<path>)."
        )
    if primer_set.kind == "paired":
        return [(pid, fwd, pid, rev) for pid, fwd, rev in primer_set.pairs[:n]]

    pool = list(primer_set.primers)
    cursor = 0

    def draw(ok) -> tuple[str, str]:
        nonlocal cursor
        while cursor < len(pool):
            pid, seq = pool[cursor]
            cursor += 1
            if ok(seq):
                return pid, seq
        raise ValueError(
            f"Ran out of usable primers in {primer_set.name!r} while avoiding junction "
            "sites; provide a larger primer set (tiled.primer_set=<path>)."
        )

    out: list[tuple[str, str, str, str]] = []
    for start, end in tiles:
        ctx5, ctx3 = tile_contexts(reference, start, end, params)
        fid, fwd = draw(lambda s: _end5_ok(s, ctx5, params))
        rid, rev = draw(lambda s: _end3_ok(s, ctx3, params))
        out.append((fid, fwd, rid, rev))
    return out


def _vector_site_positions(vector: str, enzyme: str, circular: bool) -> list[int]:
    """Positions of every ``enzyme`` recognition site (both strands) on a destination
    vector, wrapping the origin when it is circular."""
    site = ENZYME_SITES[enzyme].upper()
    rec_rc = reverse_complement(site)
    scan = vector + vector[: len(site) - 1] if circular else vector
    hits = [i for i in range(len(scan) - len(site) + 1) if scan[i:i + len(site)] == site]
    if rec_rc != site:
        hits += [i for i in range(len(scan) - len(rec_rc) + 1) if scan[i:i + len(rec_rc)] == rec_rc]
    return sorted({i % len(vector) for i in hits})


def _vector_extra_sites(tiles: list[TileInfo], enzyme: str) -> list[str]:
    """Human-readable findings for tiles whose destination vector carries an ``enzyme``
    site beyond the two intended drop-out sites (which would misdirect assembly)."""
    out: list[str] = []
    for t in tiles:
        pos = _vector_site_positions(t.vector, enzyme, t.topology == "circular")
        if len(pos) != 2:
            out.append(f"tile{t.index}: {len(pos)} {enzyme} site(s) (expected 2) at {pos}")
    return out


def tile_library(library, params: TiledAssemblyParams) -> dict:
    """Compute tiles, assign primers, and assemble an oligo for every placeable variant.

    Returns a dict with ``tiles`` (list[TileInfo]), per-row ``tile``/``oligo`` columns,
    the loaded ``primer_set``, ``unplaced`` variant names (WT control / codons split
    across all tile boundaries), and ``vector_extra_sites`` findings.

    With ``params.starting_vector`` set, each ``TileInfo.vector`` is the full destination
    plasmid (backbone with the tile window dropped out) and the two terminal overhangs are
    drawn from the backbone flanking the CDS rather than the ``vector_context_*`` defaults."""
    from dataclasses import replace

    reference = library.reference
    if reference is None:
        raise ValueError("Library is not codon-optimized yet, call codon_optimize() first.")

    dest = None
    if params.starting_vector:
        from .vector_io import assemble_vector, locating_kwargs, resolve_destination, terminal_contexts

        dest = resolve_destination(params.starting_vector,
                                   **locating_kwargs(library.spec, params))
        term5, term3 = terminal_contexts(dest, params.overhang_len)
        params = replace(params, vector_context_5=term5, vector_context_3=term3)

    tiles_coords = compute_tiles(len(reference), params)
    primer_set = load_primer_set(params.primer_set, params.enzyme)
    _check_primer_length(primer_set, params)
    assignments = _assign_primers(tiles_coords, reference, primer_set, params)

    topology = dest.topology if dest is not None else "linear"
    tiles = []
    for i, ((s, e), (fid, fwd, rid, rev)) in enumerate(zip(tiles_coords, assignments)):
        cds_region = reference[:s] + params.vector_insert + reference[e:]
        vector = assemble_vector(dest, cds_region) if dest is not None else cds_region
        tiles.append(TileInfo(index=i, start=s, end=e, fwd_id=fid, fwd=fwd,
                              rev_id=rid, rev=rev, vector=vector, topology=topology))

    # The whole destination vector runs through the enzyme at assembly, so a stray site
    # anywhere (backbone, CDS arms, or a splice junction) breaks it. When the user chose
    # to keep the vector's own CDS (use_vector_cds) we flag this in QC; otherwise, when
    # we control the reference, a stray site can only be the user's backbone, so we stop.
    vector_extra = _vector_extra_sites(tiles, params.enzyme)
    tiled_spec = library.spec.tiled
    if vector_extra and not (tiled_spec and tiled_spec.use_vector_cds):
        raise ValueError(
            f"destination vector(s) carry {params.enzyme} site(s) beyond the two intended "
            "drop-out sites, which would misdirect Golden Gate assembly: "
            + "; ".join(vector_extra)
            + ". Domesticate the backbone (remove the extra sites) before tiling, or set "
            "tiled.use_vector_cds to keep the sequence and only flag it."
        )

    df = library.df
    tile_col: list = []
    oligo_col: list = []
    unplaced: list[str] = []
    for name, dna, idx in zip(df["name"], df["variable_dna"], df["mut_index"]):
        if pd.isna(idx) or pd.isna(dna):
            tile_col.append(pd.NA)
            oligo_col.append(pd.NA)
            if pd.isna(idx) and not pd.isna(dna):
                unplaced.append(str(name))       # WT control: no position to tile
            continue
        ti = assign_tile(int(idx), tiles_coords)
        if ti is None:
            tile_col.append(pd.NA)
            oligo_col.append(pd.NA)
            unplaced.append(str(name))
            continue
        t = tiles[ti]
        oligo_col.append(assemble_oligo(reference, str(dna), t.start, t.end, t.fwd, t.rev, params))
        tile_col.append(ti)

    return {
        "tiles": tiles,
        "tile": tile_col,
        "oligo": oligo_col,
        "primer_set": primer_set,
        "unplaced": unplaced,
        "vector_extra_sites": vector_extra,
        "params": params,     # resolved (vector-derived terminal overhangs filled in)
    }
