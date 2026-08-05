"""Tiled-assembly layout.

Given a frozen WT reference and its stamped variants (see ``optimize/backbone.py``),
split the CDS into overlapping tile windows sized to the oligo budget, assign each
variant to the tile that contains its codon, and assemble the final synthesis oligo:

    fwd · [pad5] · SITE · spacer5 · [overhang5] · TILE · [overhang3] · spacer3 · SITE_rc · [pad3] · revcomp(rev)

The two ``SITE`` copies are the Golden Gate (Type IIS) recognition sequences; after
digestion the two overhang regions (drawn from the flanking WT CDS) become the fused
overhangs that ligate the tile into a per-tile destination vector carrying the rest of
the WT CDS. Primers come from a validated orthogonal set (see ``primers.py``) so each
tile's sublibrary can be selectively amplified out of the shared pool.

Two site copies are the *only* recognition sites allowed in the oligo. Because adding
the flanking primers/sites/spacers can create an unintended site spanning a junction,
assembly is screened on the *fully assembled* oligo, and pooled primers that would form
a junction site are skipped during assignment.

``pad5`` / ``pad3`` are empty unless ``tiled.pad_oligos`` is set. Tiles differ in size, so
without them the pool holds oligos of several lengths, and moving the boundaries for better
overhangs (``layout/boundaries.py``) widens that spread. The filler brings every oligo to one
length. It sits outside both recognition sites, so it is amplified with the oligo and then
cut away, and nothing padded reaches the assembled plasmid.
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
    """One tile of a tiled library, the window plus what it takes to amplify and clone it.

    ``lib.tiles`` is a list of these, one per tile, and ``to_primer_order`` and
    ``to_vectors`` write from them. ``index`` counts from 0 at the 5' end of the CDS.
    ``start`` and ``end`` are the window ``[start, end)`` on the frozen reference, always
    on codon boundaries, so ``lib.reference[t.start:t.end]`` is the WT stretch this tile
    covers.

    ``fwd_id`` and ``rev_id`` name the pair drawn from the orthogonal primer set, and
    ``fwd`` and ``rev`` are those primers 5'->3' as you would order them. ``fwd`` sits
    verbatim at the 5' end of every oligo in this tile's sublibrary and the reverse
    complement of ``rev`` at the 3' end, so the pair pulls this tile out of the shared
    pool.

    ``vector`` is the destination vector this tile assembles into, the reference with
    this one window swapped for ``vector_insert`` and, when a starting vector is set,
    that whole CDS region put back into the backbone. ``topology`` is ``"circular"``
    when the vector is a real plasmid backbone and ``"linear"`` when there is none, and
    it is what tells QC whether to scan the vector across its origin.

    ``pad_5`` and ``pad_3`` are filler that brings this tile's oligo up to the pool's common
    length, empty unless ``tiled.pad_oligos`` is set. They sit between the primer and the
    recognition site at each end, which is outside everything the enzyme keeps, so the
    digest cuts them away and nothing padded reaches the assembled product.
    """

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
    pad_5: str = ""     # filler between the forward primer and the 5' recognition site
    pad_3: str = ""     # filler between the 3' recognition site and the reverse primer

    @property
    def lead(self) -> str:
        """Everything 5' of the recognition site: the forward primer and its pad."""
        return self.fwd + self.pad_5

    @property
    def trail(self) -> str:
        """Everything 3' of the recognition site: the pad and the reverse primer as it sits
        on the oligo (reverse-complemented)."""
        return self.pad_3 + reverse_complement(self.rev)


def _overhead(params: TiledAssemblyParams, rec_len: int) -> int:
    """Non-tile bases on each oligo: two primers, two sites, two spacers, two overhangs."""
    return (
        2 * params.primer_length
        + 2 * rec_len
        + len(params.spacer_5)
        + len(params.spacer_3)
        + 2 * params.overhang_len
    )


def max_tile_codons(params: TiledAssemblyParams) -> int:
    """The longest tile the oligo budget allows, in codons. ``tile_size`` caps it directly
    when set, otherwise it is what the budget leaves after the primers, sites, spacers, and
    overhangs."""
    rec_len = len(ENZYME_SITES[params.enzyme])
    max_tile = params.tile_size or (params.oligo_budget - _overhead(params, rec_len))
    codons = max_tile // 3
    if codons < 1:
        raise ValueError(
            f"Oligo budget ({params.oligo_budget} bp) leaves no room for a tile after "
            f"primers and sites (needs > {_overhead(params, rec_len)} bp)."
        )
    return codons


def tile_windows(reference: str, params: TiledAssemblyParams) -> list[tuple[int, int]]:
    """The tile windows a library is laid out with.

    The balanced split from ``compute_tiles``, or, with ``tiled.optimize_overhangs`` set, the
    boundary positions whose fused overhangs share the least homology (see
    ``layout/boundaries.py``). The search keeps the tile count and every constraint the
    balanced split obeys, and falls back to it whenever it cannot do better, so switching the
    flag on never costs a design anything but the boundaries moving.
    """
    baseline = compute_tiles(len(reference), params)
    if not params.optimize_overhangs:
        return baseline
    from .boundaries import search_windows

    return search_windows(reference, params, baseline, max_tile_codons(params))


def compute_tiles(cds_len: int, params: TiledAssemblyParams) -> list[tuple[int, int]]:
    """Contiguous, gap-free tile windows covering the whole CDS.

    Boundaries fall on codon multiples so a mutated codon is never split across two
    tiles, and windows are balanced to the fewest tiles whose size fits the oligo
    budget. The 4-bp Golden Gate overhang each oligo carries is drawn from the
    flanking WT CDS, so adjacent tiles need no explicit sequence overlap. That is also
    why a terminal tile has to be at least ``overhang_len`` long once there is more than
    one tile: the boundary it creates needs a full overhang of CDS on both sides.
    """
    if cds_len % 3 != 0:
        raise ValueError(f"Reference CDS length ({cds_len}) is not a multiple of 3.")
    max_codons = max_tile_codons(params)

    codons = cds_len // 3
    tile_n = -(-codons // max_codons)               # ceil: fewest tiles that fit the budget
    # Spread the codons as evenly as the count allows: the first `wide` tiles take one
    # extra codon. Repeating a ceil-sized window instead would leave the remainder in the
    # last tile, which for some CDS lengths is a stub too short to carry an overhang.
    per, wide = divmod(codons, tile_n)
    sizes = [per + 1] * wide + [per] * (tile_n - wide)

    min_codons = -(-params.overhang_len // 3)        # ceil: codons a full overhang spans
    if tile_n > 1 and min(sizes[0], sizes[-1]) < min_codons:
        raise ValueError(
            f"Splitting a {cds_len} bp CDS into {tile_n} tiles gives a terminal tile of "
            f"{min(sizes[0], sizes[-1]) * 3} bp, shorter than the {params.overhang_len} bp "
            "Golden Gate overhang, so its fused overhang cannot be drawn from the flanking "
            "CDS. Raise tiled.tile_size (or tiled.oligo_budget) to use fewer, longer tiles."
        )

    tiles, at = [], 0
    for n in sizes:
        tiles.append((at, at + n * 3))
        at += n * 3

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

    Which one applies depends on *where the boundary is*, not on whether a full overhang
    happens to fit: an internal boundary always takes CDS bases, because that is what its
    destination vector presents at the junction. ``compute_tiles`` keeps every terminal
    tile at least ``overhang_len`` long, so those bases always exist.

    QC uses this too, so the overhangs it checks are the ones the oligos carry."""
    o = params.overhang_len
    ctx5 = params.vector_context_5 if start <= 0 else reference[start - o:start]
    ctx3 = params.vector_context_3 if end >= len(reference) else reference[end:end + o]
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


# Filler for length padding: no BsaI or BsmBI site on either strand, no run longer than two,
# and an even base mix, so a slice of it can sit next to a recognition site without becoming
# part of one. `_pad` still checks, and slides along the filler when a junction goes wrong.
_FILLER = ("TAGCATCAGTTACGCATGACTAGCTTACGATCAGCATTGACGTACTAGCATGTCAGTACGT"
           "ATCGTTACGCATGCTAAGCTAGTCATGCGTAACTGATGCTAGTACGTCAATGCTAGTACG")


def _pad(n: int, before: str, after: str, enzyme: str) -> str:
    """``n`` filler bases to sit between ``before`` and ``after`` without completing a site.

    Slides a window along ``_FILLER`` until the junction it makes is clean, so the choice is
    deterministic and the same design always pads the same way. Falls back to the first
    window if nothing is clean, since ``check()`` reports the site either way and a silent
    length mismatch would be worse."""
    if n <= 0:
        return ""
    rec = ENZYME_SITES[enzyme].upper()
    rec_rc = reverse_complement(rec)
    doubled = _FILLER * (n // len(_FILLER) + 2)
    for offset in range(len(_FILLER)):
        pad = doubled[offset:offset + n]
        window = before[-(len(rec) - 1):] + pad + after[:len(rec) - 1]
        if rec not in window and rec_rc not in window:
            return pad
    return doubled[:n]


def pad_target(windows: list[tuple[int, int]], params: TiledAssemblyParams) -> int | None:
    """The one length every oligo in the pool is brought up to, or None when padding is off.

    ``tiled.pad_target`` sets it outright. Otherwise it is the longest oligo the layout
    already produces, which evens the pool out without making any oligo longer than it had to
    be. A target the budget cannot hold is refused rather than silently exceeded."""
    if not params.pad_oligos:
        return None
    overhead = _overhead(params, len(ENZYME_SITES[params.enzyme]))
    longest = overhead + max(e - s for s, e in windows)
    target = params.pad_target or longest
    if target < longest:
        raise ValueError(
            f"tiled.pad_target is {target} bp but the longest oligo this layout needs is "
            f"{longest} bp, so it cannot be padded down. Raise pad_target, or lower "
            "tiled.tile_size / tiled.oligo_budget to use shorter tiles."
        )
    if target > params.oligo_budget:
        raise ValueError(
            f"tiled.pad_target is {target} bp, past the {params.oligo_budget} bp oligo "
            "budget. Raise tiled.oligo_budget or lower pad_target."
        )
    return target


def pad_lengths(tile_len: int, params: TiledAssemblyParams,
                target: int | None) -> tuple[int, int]:
    """How many filler bases each end of one tile's oligo takes to reach ``target``.

    Split as evenly as the shortfall allows, the odd base going to the 5' end."""
    if target is None:
        return 0, 0
    short = target - (_overhead(params, len(ENZYME_SITES[params.enzyme])) + tile_len)
    if short <= 0:
        return 0, 0
    return short - short // 2, short // 2


def assemble_oligo(reference: str, variant_cds: str, start: int, end: int,
                   fwd: str, rev: str, params: TiledAssemblyParams,
                   pad_5: str = "", pad_3: str = "") -> str:
    """The full synthesis oligo for one variant/tile (all uppercase).

    ``pad_5`` / ``pad_3`` are optional filler that evens the pool out to one length. They go
    between each primer and the recognition site beside it, outside the region the enzyme
    releases."""
    rec = ENZYME_SITES[params.enzyme].upper()
    rec_rc = reverse_complement(rec)
    ctx5, ctx3 = tile_contexts(reference, start, end, params)
    tile = variant_cds[start:end]
    return (
        fwd + pad_5 + rec + params.spacer_5 + ctx5
        + tile
        + ctx3 + params.spacer_3 + rec_rc + pad_3 + reverse_complement(rev)
    )


def wt_oligo(reference: str, tile: TileInfo, params: TiledAssemblyParams) -> str:
    """The WT control oligo for one tile: that tile's window taken straight from the frozen
    reference, carrying the same primers, sites, overhangs, and padding as the tile's mutants.

    This is also the molecule the assembly simulation digests to work out what a clean
    reaction yields for the tile, so the control that ships is the one QC reasons about."""
    return assemble_oligo(reference, reference, tile.start, tile.end,
                          tile.fwd, tile.rev, params, tile.pad_5, tile.pad_3)


_MUT_COLS = ("position", "wt_residue", "mut_residue", "codon", "mut_index")


def wt_control_rows(df: pd.DataFrame, reference: str, protein: str,
                    tiles: list[TileInfo], params: TiledAssemblyParams) -> pd.DataFrame:
    """One WT control member per tile, named ``WT_Tile_<i>`` and carrying that tile's
    WT oligo.

    A tiled pool is N sublibraries that get amplified and assembled separately, so the
    library's single global WT row rides on no oligo and every sublibrary would ship with
    nothing unmutated to normalize against. Each row is a copy of the first row, so it keeps
    the adaptors and whatever else the generator wrote, with the mutation fields cleared and
    the whole reference as its coding sequence."""
    base = df.iloc[0].copy()
    base["protein"] = protein
    base["variable_dna"] = reference
    for col in _MUT_COLS:
        if col in base.index:
            base[col] = pd.NA

    rows = []
    for t in tiles:
        row = base.copy()
        row["name"] = f"WT_Tile_{t.index}"
        row["tile"] = t.index
        row["oligo"] = wt_oligo(reference, t, params)
        rows.append(row)
    return pd.DataFrame(rows).reset_index(drop=True)


def site_positions(oligo: str, enzyme: str) -> tuple[list[int], list[int]]:
    """(forward-strand hits, reverse-strand hits) of ``enzyme``'s recognition site.

    A palindromic site reads the same on both strands, so each hit appears in both lists."""
    rec = ENZYME_SITES[enzyme].upper()
    rec_rc = reverse_complement(rec)

    def _find(sub: str) -> list[int]:
        return [i for i in range(len(oligo) - len(sub) + 1) if oligo[i:i + len(sub)] == sub]

    return _find(rec), _find(rec_rc)


def extra_sites(oligo: str, lead_len: int, trail_len: int, enzyme: str) -> bool:
    """True if the oligo carries any recognition site beyond the two intended ones.

    ``lead_len`` is everything before the forward site (the primer, plus the pad when the
    pool is evened out to one length) and ``trail_len`` everything after the reverse site,
    which is what puts the two intended positions where they belong.

    For a palindromic site both intended positions show up on both strands, so compare
    the combined set of positions instead of each strand's list. No enzyme in
    ``ENZYME_SITES`` is palindromic today, but a user-supplied one can be."""
    rec = ENZYME_SITES[enzyme].upper()
    rec_len = len(rec)
    fwd_hits, rev_hits = site_positions(oligo, enzyme)
    want_fwd = lead_len
    want_rev = len(oligo) - trail_len - rec_len
    if rec == reverse_complement(rec):
        return set(fwd_hits) | set(rev_hits) != {want_fwd, want_rev}
    return fwd_hits != [want_fwd] or rev_hits != [want_rev]


def _count_sites(s: str, enzyme: str) -> tuple[int, int]:
    fwd, rev = site_positions(s, enzyme)
    return len(fwd), len(rev)


def _end5_ok(fwd: str, ctx5: str, params: TiledAssemblyParams, pad_5: str = "") -> bool:
    """The forward primer must add no site to the 5' region (fwd|pad|site|spacer|overhang)
    beyond the one intended forward site."""
    rec = ENZYME_SITES[params.enzyme].upper()
    nf, nr = _count_sites(fwd + pad_5 + rec + params.spacer_5 + ctx5, params.enzyme)
    return nf == 1 and nr == 0


def _end3_ok(rev: str, ctx3: str, params: TiledAssemblyParams, pad_3: str = "") -> bool:
    """The reverse primer must add no site to the 3' region (overhang|spacer|site|pad|rev)
    beyond the one intended reverse site."""
    rec_rc = reverse_complement(ENZYME_SITES[params.enzyme].upper())
    nf, nr = _count_sites(
        ctx3 + params.spacer_3 + rec_rc + pad_3 + reverse_complement(rev), params.enzyme
    )
    return nf == 0 and nr == 1


def _assign_primers(tiles: list[tuple[int, int]], reference: str, primer_set: PrimerSet,
                    params: TiledAssemblyParams,
                    target: int | None = None) -> list[tuple[str, str, str, str, str, str]]:
    """Choose (fwd_id, fwd, rev_id, rev, pad_5, pad_3) per tile. A pool is drawn
    primer-by-primer, skipping any that would form a junction site; a paired set is used as
    given.

    The 5' and 3' junctions are independent (each involves only its own primer, pad, site,
    spacer, and overhang), so the two ends are validated separately. A pad's length is fixed
    by the tile before any primer is picked, but its bases sit against the primer, so the pad
    is built per candidate and the junction judged with it in place."""
    n = len(tiles)
    if primer_set.capacity < n:
        raise ValueError(
            f"Primer set {primer_set.name!r} supplies {primer_set.capacity} tile(s) of primers "
            f"but the CDS needs {n}. Provide a larger set (tiled.primer_set=<path>)."
        )
    rec = ENZYME_SITES[params.enzyme].upper()
    rec_rc = reverse_complement(rec)
    pads = [pad_lengths(e - s, params, target) for s, e in tiles]

    if primer_set.kind == "paired":
        out = []
        for (pid, fwd, rev), (n5, n3) in zip(primer_set.pairs[:n], pads):
            out.append((pid, fwd, pid, rev,
                        _pad(n5, fwd, rec, params.enzyme),
                        _pad(n3, rec_rc, reverse_complement(rev), params.enzyme)))
        return out

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

    out = []
    for (start, end), (n5, n3) in zip(tiles, pads):
        ctx5, ctx3 = tile_contexts(reference, start, end, params)
        fid, fwd = draw(lambda s: _end5_ok(s, ctx5, params, _pad(n5, s, rec, params.enzyme)))
        rid, rev = draw(
            lambda s: _end3_ok(s, ctx3, params,
                               _pad(n3, rec_rc, reverse_complement(s), params.enzyme))
        )
        out.append((fid, fwd, rid, rev,
                    _pad(n5, fwd, rec, params.enzyme),
                    _pad(n3, rec_rc, reverse_complement(rev), params.enzyme)))
    return out


def vector_site_positions(vector: str, enzyme: str, circular: bool) -> list[int]:
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
        pos = vector_site_positions(t.vector, enzyme, t.topology == "circular")
        if len(pos) != 2:
            out.append(f"tile{t.index}: {len(pos)} {enzyme} site(s) (expected 2) at {pos}")
    return out


def tile_library(library, params: TiledAssemblyParams) -> dict:
    """Compute tiles, assign primers, and assemble an oligo for every placeable variant.

    Returns a dict with ``tiles`` (list[TileInfo]), per-row ``tile``/``oligo`` columns,
    the loaded ``primer_set``, ``unplaced`` variant names (the global WT control / codons
    split across all tile boundaries), ``wt_controls`` (the per-tile WT rows to append,
    None when ``params.wt_controls`` is off), and ``vector_extra_sites`` findings.

    With ``params.starting_vector`` set, each ``TileInfo.vector`` is the full destination
    plasmid (backbone with the tile window dropped out) and the two terminal overhangs are
    drawn from the backbone flanking the CDS rather than the ``vector_context_*`` defaults."""
    from dataclasses import replace

    reference = library.reference
    if reference is None:
        raise ValueError("Library is not codon-optimized yet, call codon_optimize() first.")

    dest = None
    vec = library.spec.resolve_vector(params)
    if vec is not None:
        from .vector_io import assemble_vector, locating_kwargs, resolve_destination, terminal_contexts

        dest = resolve_destination(vec.path, **locating_kwargs(library.spec, params))
        term5, term3 = terminal_contexts(dest, params.overhang_len)
        # Record the resolved path on the params too, so the layout the library carries
        # names its own backbone even when the vector came from the spec rather than the
        # tiled block.
        params = replace(params, starting_vector=vec.path,
                         vector_context_5=term5, vector_context_3=term3)

    from .boundaries import windows_cost

    tiles_coords = tile_windows(reference, params)
    # Both scores, so the run record says what the boundaries cost and what they would
    # have cost on the balanced split, whether or not the search ran.
    balanced = compute_tiles(len(reference), params)
    overhang_cost = windows_cost(reference, params, tiles_coords)
    overhang_cost_unsearched = (
        overhang_cost if tiles_coords == balanced
        else windows_cost(reference, params, balanced)
    )
    primer_set = load_primer_set(params.primer_set, params.enzyme)
    _check_primer_length(primer_set, params)
    # One length for the whole pool, when asked for. Fixed before the primers are drawn,
    # since a pad's length follows from the tile and only its bases follow from the primer.
    target = pad_target(tiles_coords, params)
    assignments = _assign_primers(tiles_coords, reference, primer_set, params, target)

    topology = dest.topology if dest is not None else "linear"
    tiles = []
    for i, ((s, e), (fid, fwd, rid, rev, p5, p3)) in enumerate(zip(tiles_coords, assignments)):
        cds_region = reference[:s] + params.vector_insert + reference[e:]
        vector = assemble_vector(dest, cds_region) if dest is not None else cds_region
        tiles.append(TileInfo(index=i, start=s, end=e, fwd_id=fid, fwd=fwd,
                              rev_id=rid, rev=rev, vector=vector, topology=topology,
                              pad_5=p5, pad_3=p3))

    # The whole destination vector runs through the enzyme at assembly, so a stray site
    # anywhere (backbone, CDS arms, or a splice junction) breaks it. When the user chose
    # to keep the vector's own CDS (use_vector_cds) we flag this in QC; otherwise, when
    # we control the reference, a stray site can only be the user's backbone, so we stop.
    vector_extra = _vector_extra_sites(tiles, params.enzyme)
    resolved_vector = library.spec.resolve_vector(params)
    if vector_extra and not (resolved_vector and resolved_vector.use_vector_cds):
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
                # The global WT control has no position, so no one tile contains it. The
                # per-tile WT controls below are what actually ship in the pool.
                unplaced.append(str(name))
            continue
        ti = assign_tile(int(idx), tiles_coords)
        if ti is None:
            tile_col.append(pd.NA)
            oligo_col.append(pd.NA)
            unplaced.append(str(name))
            continue
        t = tiles[ti]
        oligo_col.append(assemble_oligo(reference, str(dna), t.start, t.end,
                                        t.fwd, t.rev, params, t.pad_5, t.pad_3))
        tile_col.append(ti)

    wt_controls = None
    if params.wt_controls:
        wt_controls = wt_control_rows(df, reference, library.spec.designed_sequence,
                                      tiles, params)

    return {
        "tiles": tiles,
        "tile": tile_col,
        "oligo": oligo_col,
        "primer_set": primer_set,
        "unplaced": unplaced,
        "vector_extra_sites": vector_extra,
        "wt_controls": wt_controls,   # extra rows to append, one WT member per tile
        "params": params,     # resolved (vector-derived terminal overhangs filled in)
        "overhang_cost": overhang_cost,
        "overhang_cost_unsearched": overhang_cost_unsearched,
    }
