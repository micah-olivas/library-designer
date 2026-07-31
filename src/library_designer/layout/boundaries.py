"""Choosing where the tile boundaries fall, so the fused overhangs tell each other apart.

``compute_tiles`` splits a CDS into the fewest codon-aligned windows that fit the oligo
budget and balances them. Nothing in that picks the overhangs. They are whatever bases
happen to sit at the boundaries, so a tile can end up with two ends that cannot tell each
other apart and no way to see it coming.

What counts is a tile against itself. Each tile is amplified out of the pool with its own
primer pair and assembled into the vector built around its own window, so its reaction holds
its own fragments and nothing else. Homology between two different tiles' overhangs is free,
and the search spends the layout's slack only on the comparisons that can actually misfire:
a tile's two ends against each other, and each end against its own reverse complement.

There is usually slack. The budget sets a *maximum* tile, the balanced split sits under it,
and every spare codon is a boundary that could move. Shifting one boundary by a codon
changes both overhangs it creates, so a search over the boundary positions is a search over
the overhang set. This module runs that search, scoring a layout with the same functions QC
grades it with (``checks/overhangs.py``), so the two cannot disagree.

The search is a beam over the boundaries left to right. It is a candidate generator, not the
judge: every layout it surfaces is rescored whole by ``windows_cost``, and the balanced split
is scored alongside them, so the result is never worse than what the design would have had
without the search.

A boundary can also spoil an oligo. The overhang sits between the spacer and the tile, so a
boundary whose bases spell part of a recognition site next to the spacer makes an oligo no
primer can rescue. That is priced in at the same weight as a collision.
"""
from __future__ import annotations

from functools import lru_cache

from ..checks.overhangs import risk_of, self_risk, shared_bases
from ..regions import reverse_complement

BEAM = 200              # partial layouts carried forward past each boundary
MAX_CANDIDATES = 96     # boundary positions tried per step, spread over the feasible range

# What one comparison costs. A collision between the two ends of a single tile, a palindrome,
# and a boundary that spoils the oligo all make the wrong product the expected one, so they
# dominate. Below them sits a single mismatch, which T4 ligase joins at a rate you will see on
# a plate, and below that anything merely above target.
_FATAL = 1_000_000
_HIGH = 10_000
_WATCH = 1


@lru_cache(maxsize=None)
def _worst_risk(a: str, b: str) -> str:
    """The tier of the closer of the two orientations two overhangs can meet in."""
    worst = max(shared_bases(a, b), shared_bases(a, reverse_complement(b)))
    return risk_of(worst, len(a))


@lru_cache(maxsize=None)
def _palindrome(seq: str) -> bool:
    return self_risk(shared_bases(seq, reverse_complement(seq)), len(seq)) == "collision"


def _pair_cost(a: tuple[int, str], b: tuple[int, str]) -> int:
    """Cost of one pair of ends, each ``(tile index, sequence)``.

    Two ends in different tiles cost nothing. They are amplified and assembled separately,
    each into the vector built around its own window, so they never meet and homology between
    them is not a hazard to spend the layout's slack on."""
    (tile_a, seq_a), (tile_b, seq_b) = a, b
    if tile_a != tile_b:
        return 0
    risk = _worst_risk(seq_a, seq_b)
    if risk == "collision":
        return _FATAL
    if risk == "high":
        return _HIGH
    return _WATCH if risk == "watch" else 0


def _count(seq: str, sub: str) -> int:
    return sum(1 for i in range(len(seq) - len(sub) + 1) if seq[i:i + len(sub)] == sub)


@lru_cache(maxsize=None)
def _junction_cost(ctx_5: str, ctx_3: str, enzyme: str, spacer_5: str, spacer_3: str) -> int:
    """Whether the bases a boundary contributes spell a recognition site of their own.

    The overhang rides between the spacer and the tile, so a site formed there belongs to the
    boundary and no choice of primer removes it. Checked without the primers, which is what
    makes it a property of the boundary alone."""
    from ..checks.motifs import ENZYME_SITES

    rec = ENZYME_SITES[enzyme].upper()
    rec_rc = reverse_complement(rec)
    cost = 0
    head = rec + spacer_5 + ctx_5                    # the 5' end of the tile's oligo
    tail = ctx_3 + spacer_3 + rec_rc                 # the 3' end
    if (_count(head, rec), _count(head, rec_rc)) != (1, 0 if rec != rec_rc else 1):
        cost += _FATAL
    if (_count(tail, rec), _count(tail, rec_rc)) != (0 if rec != rec_rc else 1, 1):
        cost += _FATAL
    return cost


def layout_ends(reference: str, params, windows: list[tuple[int, int]]) -> list[tuple[int, str]]:
    """Every fused overhang a layout leaves, as ``(tile index, sequence)``.

    Reads them with ``tile_contexts``, the same helper the oligos and QC use, so the search
    cannot score a layout differently from the way it will be built."""
    from .tiled import tile_contexts

    ends: list[tuple[int, str]] = []
    for i, (s, e) in enumerate(windows):
        c5, c3 = tile_contexts(reference, s, e, params)
        ends.append((i, c5))
        ends.append((i, c3))
    return ends


def windows_cost(reference: str, params, windows: list[tuple[int, int]]) -> int:
    """What a whole layout costs: every pair of its overhangs, every palindrome among them,
    and every boundary that spoils an oligo.

    This is the judge. The beam only proposes."""
    ends = layout_ends(reference, params, windows)
    total = sum(_FATAL for _, seq in ends if _palindrome(seq))
    for i, a in enumerate(ends):
        for b in ends[i + 1:]:
            total += _pair_cost(a, b)
    for (_s, e), (nxt_s, _e) in zip(windows, windows[1:]):
        # The internal boundary between two tiles: the 3' overhang of the one ending here and
        # the 5' overhang of the one starting there.
        o = params.overhang_len
        total += _junction_cost(reference[nxt_s - o:nxt_s], reference[e:e + o],
                                params.enzyme, params.spacer_5, params.spacer_3)
    return total


def _feasible(j: int, n: int, codons: int, max_codons: int, term: int) -> tuple[int, int]:
    """``(lo, hi)`` codon positions the ``j``-th internal boundary can take.

    Boundary ``j`` closes tile ``j-1`` and opens tile ``j``, so ``j`` tiles sit before it and
    ``n - j`` after. Each side has to hold its tiles at no more than ``max_codons`` apiece,
    and the two terminal tiles need ``term`` codons to carry an overhang drawn from the CDS.
    """
    lo = max(term + (j - 1), codons - (n - j) * max_codons)
    hi = min(j * max_codons, codons - term - (n - j - 1))
    return lo, hi


def _candidates(lo: int, hi: int) -> list[int]:
    """Boundary positions to try, thinned to ``MAX_CANDIDATES`` when the range is wide."""
    span = hi - lo + 1
    if span <= 0:
        return []
    if span <= MAX_CANDIDATES:
        return list(range(lo, hi + 1))
    step = span / MAX_CANDIDATES
    return sorted({lo + int(i * step) for i in range(MAX_CANDIDATES)})


def search_windows(reference: str, params, baseline: list[tuple[int, int]],
                   max_codons: int) -> list[tuple[int, int]]:
    """Tile windows whose overhangs tell each other apart as well as this CDS allows.

    Same tile count and same constraints as ``baseline``, which is the balanced split. The
    boundaries move, nothing else. ``baseline`` is scored against every candidate and wins
    ties, so turning the search on can only hold a layout steady or improve it.
    """
    n = len(baseline)
    o = params.overhang_len
    codons = len(reference) // 3
    term = -(-o // 3)                       # codons a terminal tile needs to carry an overhang
    if n < 2:
        return baseline                     # one tile has no internal boundary to move

    # Seed with the 5' terminal overhang, which the vector fixes and no boundary can change.
    seed_5 = params.vector_context_5 if baseline[0][0] <= 0 else reference[:o]
    start = (_FATAL if _palindrome(seed_5) else 0, (), ((0, seed_5),))
    beam = [start]

    for j in range(1, n):
        lo, hi = _feasible(j, n, codons, max_codons, term)
        nxt = []
        for cost, chosen, ends in beam:
            prev = chosen[-1] if chosen else 0
            floor = max(lo, prev + (term if j == 1 else 1))
            ceiling = min(hi, prev + max_codons)
            for p in _candidates(floor, ceiling):
                b = 3 * p
                ctx_3 = reference[b:b + o]          # closes tile j-1
                ctx_5 = reference[b - o:b]          # opens tile j
                fresh = ((j - 1, ctx_3), (j, ctx_5))
                add = _junction_cost(ctx_5, ctx_3, params.enzyme,
                                     params.spacer_5, params.spacer_3)
                add += sum(_FATAL for _, s in fresh if _palindrome(s))
                add += _pair_cost(fresh[0], fresh[1])
                for e in ends:
                    add += _pair_cost(e, fresh[0]) + _pair_cost(e, fresh[1])
                nxt.append((cost + add, chosen + (p,), ends + fresh))
        if not nxt:
            return baseline                 # no layout satisfies the constraints; keep the split
        # Sort on the cost and then the boundaries, so an equal-cost tie always resolves the
        # same way and a given CDS always tiles identically.
        nxt.sort(key=lambda row: (row[0], row[1]))
        beam = nxt[:BEAM]

    # The beam proposes; windows_cost judges. Scoring the baseline the same way is what makes
    # the search safe to turn on.
    best, best_cost = baseline, windows_cost(reference, params, baseline)
    for _cost, chosen, _ends in beam:
        cuts = [0] + [3 * p for p in chosen] + [len(reference)]
        windows = list(zip(cuts, cuts[1:]))
        c = windows_cost(reference, params, windows)
        if c < best_cost:
            best, best_cost = windows, c
    return best
