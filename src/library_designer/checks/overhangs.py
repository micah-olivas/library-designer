"""The fused Golden Gate overhangs a design leaves, and how well they tell each other apart.

A Type IIS reaction is only directional if the overhangs it works with are distinct. Two
failure modes follow from overhangs that are too much alike.

* The cut vector's own two ends anneal to each other and it re-closes empty. That is a
  unimolecular reaction, so even a poor match gives a strong background of parent plasmid.
* The insert anneals the other way round and goes in backwards, giving a clone with the
  tile reverse-complemented.

Both are read off the same number, how many aligned bases two overhangs share. The target
here is at most 1 shared base out of 4, which is the orthogonality a designed overhang set
is chosen for. Sharing all 4 makes the wrong product the expected one, so it fails QC.
Sharing 3, a single mismatch, is the case T4 ligase actually joins at a measurable rate, so
it is called out by name. Sharing 2 is above target but not a real hazard on its own, so it
is only counted, and the full table is where it gets reviewed.

Only the exact case fails because a tiled design does not choose its overhangs. They are
read off the CDS at the tile boundaries, so the fix for a near-match is to move the boundary
with ``tiled.tile_size``, not to pick different bases.

**The unit of concern is one tile.** Each tile is amplified out of the pool with its own
primer pair and assembled into its own destination vector, so a tile's reaction contains its
own fragments and nothing else. Tile 0's overhangs and tile 3's overhangs never meet, and
they cannot be pooled either, since every tile needs the particular vector carrying the rest
of the CDS around its own window. So the only comparisons that mean anything are a tile's two
ends against each other, and each end against its own reverse complement. Cross-tile pairs
are listed when asked for, and left ungraded.

``overhang_table`` gives one row per end and ``pair_table`` one row per pair of ends, which
is what a notebook shows. ``overhang_findings`` turns the same numbers into the strings QC
prints.
"""
from __future__ import annotations

from dataclasses import dataclass

from ..regions import reverse_complement

# Aligned bases two overhangs may share and still count as orthogonal. One out of four is
# the working rule; anything above it is reported, and a full match fails.
MAX_SHARED = 1


@dataclass(frozen=True)
class OverhangEnd:
    """One fused overhang, named by the reaction it belongs to and which end it sits on.

    ``reaction`` is ``"tile0"``, ``"tile1"``, and so on for a tiled library, or ``"insert"``
    for a standard one with a single destination vector. ``end`` is ``"5'"`` or ``"3'"``, and
    ``seq`` is the overhang written on the top strand, the convention ``checks/assembly.py``
    uses, where two ends join when their sequences are equal.
    """

    reaction: str
    end: str
    seq: str

    @property
    def label(self) -> str:
        return f"{self.reaction} {self.end}"


def shared_bases(a: str, b: str) -> int:
    """How many aligned positions two overhangs match at.

    Overhangs of different lengths cannot anneal, so they score 0 rather than raising."""
    if len(a) != len(b):
        return 0
    return sum(x == y for x, y in zip(a, b))


def overhang_ends(library) -> list[OverhangEnd]:
    """Every fused overhang the design leaves, in reaction order.

    A tiled library has two per tile, read off the CDS at the tile boundaries or off the
    backbone where a tile abuts a CDS end. A standard library cloned into a starting vector
    has the one pair its destination vector presents. Anything else has no Golden Gate
    reaction to describe, so the list comes back empty.
    """
    tiles = getattr(library, "tiles", None)
    if tiles is not None:
        from ..layout.tiled import tile_contexts

        params = library.tiled_params
        ref = library.reference
        if ref is None:
            return []
        out: list[OverhangEnd] = []
        for t in tiles:
            c5, c3 = tile_contexts(ref, t.start, t.end, params)
            out.append(OverhangEnd(f"tile{t.index}", "5'", c5))
            out.append(OverhangEnd(f"tile{t.index}", "3'", c3))
        return out

    if library.spec.vector is not None and getattr(library, "reference", None):
        from ..layout.destination import build_destination

        try:
            # strict=False, since reporting the collision is the whole point of this table.
            dv = build_destination(library, strict=False)
        except ValueError:
            # A geometry problem the vector check already reports. Nothing to describe here.
            return []
        return [OverhangEnd("insert", "5'", dv.overhang_5),
                OverhangEnd("insert", "3'", dv.overhang_3)]
    return []


def self_shared(end: OverhangEnd) -> int:
    """How many bases an overhang shares with its own reverse complement.

    Sharing all of them means it is palindromic, so a second copy of the fragment anneals to
    it head to tail and the pieces concatemerize. An even-length overhang pairs its positions
    up, so this count is always even."""
    return shared_bases(end.seq, reverse_complement(end.seq))


def risk_of(shared: int, width: int) -> str:
    """The tier a shared-base count between two overhangs falls in.

    ``"collision"`` when they match at every base, so the wrong product is the expected one.
    ``"high"`` at one mismatch, which T4 ligase joins at a rate you will see on a plate.
    ``"watch"`` above ``MAX_SHARED``, worth a look but not a hazard by itself. ``"ok"``
    otherwise."""
    if shared >= width:
        return "collision"
    if shared >= width - 1:
        return "high"
    return "watch" if shared > MAX_SHARED else "ok"


def self_risk(shared: int, width: int) -> str:
    """The tier an overhang's homology with its own reverse complement falls in.

    Graded on a shorter scale than ``risk_of``, because an even-length overhang pairs its
    positions up when compared with its own reverse complement, so the count can only come
    out even. A 4 bp overhang scores 4 (palindromic, and it concatemerizes), 2 (two
    mismatches, which does not anneal), or 0. Reading 2 as "above target" the way a pair is
    read would flag more than a third of all overhangs for nothing, so only a full or
    near-full match counts here.
    """
    if shared >= width:
        return "collision"
    return "high" if shared >= width - 1 else "ok"


def _plural(n: int, word: str) -> str:
    return word if n == 1 else word + "s"


def self_note(end: OverhangEnd) -> str:
    """What to say about an overhang that anneals to itself, empty when it does not."""
    n, width = self_shared(end), len(end.seq)
    risk = self_risk(n, width)
    if risk == "collision":
        return (f"{end.label} overhang {end.seq} is palindromic, so a second copy of the "
                "fragment anneals to it and the pieces concatemerize")
    if risk == "ok":
        return ""
    return (f"{end.label} overhang {end.seq} shares {n} of {width} bases with its own reverse "
            f"complement, so a second copy of the fragment can anneal to it")


def pair_note(a: OverhangEnd, b: OverhangEnd, same: int, flipped: int) -> str:
    """What to say about a pair of overhangs, empty when nothing is wrong.

    ``same`` counts the bases they share as written, which is the annealing that closes a cut
    vector on itself. ``flipped`` counts the bases one shares with the reverse complement of
    the other, which is the annealing that puts a fragment in backwards.

    Two ends in different reactions never meet, so there is nothing to say about them.
    """
    if a.reaction != b.reaction:
        return ""
    width = len(a.seq)
    hit_same = risk_of(same, width) != "ok"
    hit_flipped = risk_of(flipped, width) != "ok"
    if not (hit_same or hit_flipped):
        return ""

    figures = ([(same, "as written")] if hit_same else []) \
        + ([(flipped, "reverse-complemented")] if hit_flipped else [])
    shares = " and ".join(
        f"{n} of {width} {_plural(width, 'base')} {how}" if i == 0 else f"{n} of {width} {how}"
        for i, (n, how) in enumerate(figures)
    )

    outcomes = []
    if hit_same:
        outcomes.append("the cut vector " + ("re-closes" if same >= width else "can re-close")
                        + " empty")
    if hit_flipped:
        outcomes.append("the fragment " + ("ligates" if flipped >= width else "can ligate")
                        + " in backwards")
    return (f"{a.reaction}: the two ends ({a.seq}, {b.seq}) share {shares}, so "
            + " and ".join(outcomes))


def overhang_table(library):
    """One row per fused overhang: which reaction and end it sits on, its sequence, how many
    bases it shares with its own reverse complement, and what that means.

    A design with no Golden Gate reaction gives an empty frame with the same columns, so a
    notebook cell does not have to branch on it."""
    import pandas as pd

    cols = ["reaction", "end", "overhang", "shared_self", "risk", "note"]
    rows = []
    for e in overhang_ends(library):
        n = self_shared(e)
        rows.append({
            "reaction": e.reaction, "end": e.end, "overhang": e.seq,
            "shared_self": n, "risk": self_risk(n, len(e.seq)), "note": self_note(e),
        })
    return pd.DataFrame(rows, columns=cols)


def pair_table(library, all_pairs: bool = False):
    """The pairs of fused overhangs that share a reaction, worst first.

    One row per tile by default, since a tile's two ends are the only pair that ever meets.
    ``shared`` counts the bases they have in common as written, the annealing that would close
    the cut vector on itself. ``shared_flipped`` counts the bases one has in common with the
    reverse complement of the other, the annealing that would put the fragment in backwards.
    ``risk`` reads off the larger of the two: ``"collision"`` at a full match, ``"high"`` at
    one mismatch, ``"watch"`` above ``MAX_SHARED``, else ``"ok"``.

    ``all_pairs=True`` adds every cross-tile pair as well. Those tiles are amplified and
    assembled separately and each needs the particular vector built around its own window, so
    they never share a tube and cannot be pooled into one. Their rows carry the counts for
    reference and a ``risk`` of ``"n/a"``, because there is no hazard to grade.
    """
    import pandas as pd

    cols = ["end_a", "overhang_a", "end_b", "overhang_b", "same_reaction",
            "shared", "shared_flipped", "risk", "note"]
    ends = overhang_ends(library)
    rows = []
    for i, a in enumerate(ends):
        for b in ends[i + 1:]:
            together = a.reaction == b.reaction
            if not (together or all_pairs):
                continue
            same = shared_bases(a.seq, b.seq)
            flipped = shared_bases(a.seq, reverse_complement(b.seq))
            rows.append({
                "end_a": a.label, "overhang_a": a.seq,
                "end_b": b.label, "overhang_b": b.seq,
                "same_reaction": together,
                "shared": same, "shared_flipped": flipped,
                "risk": risk_of(max(same, flipped), len(a.seq)) if together else "n/a",
                "note": pair_note(a, b, same, flipped),
            })
    df = pd.DataFrame(rows, columns=cols)
    if df.empty:
        return df
    # Same-reaction rows first, then by homology, so the graded pairs are always on top.
    worst = df[["shared", "shared_flipped"]].max(axis=1)
    return (df.assign(_worst=worst)
              .sort_values(["same_reaction", "_worst"], ascending=[False, False],
                           kind="stable")
              .drop(columns="_worst")
              .reset_index(drop=True))


def findings_for(ends: list[OverhangEnd]) -> tuple[list[str], list[str]]:
    """``(failures, advisories)`` for a set of ends, without needing a library.

    A failure is an exact collision between the two ends of one reaction, or a palindromic
    overhang, where the expected product is the empty vector, a backwards insert, or a
    concatemer. A single mismatch is an advisory, spelled out by name. Anything merely above
    target is collapsed into a count, so a report on a long CDS stays readable.

    Only ends that share a reaction are compared. Two tiles are amplified and assembled
    separately, each into the vector built around its own window, so their overhangs never
    meet and there is nothing to report about them.
    """
    failures: list[str] = []
    advisories: list[str] = []
    watch = 0

    def record(note: str, risk: str, fatal: bool) -> None:
        nonlocal watch
        if not note:
            return
        if fatal:
            failures.append(note)
        elif risk == "watch":
            watch += 1
        else:
            advisories.append(note)

    for e in ends:
        risk = self_risk(self_shared(e), len(e.seq))
        record(self_note(e), risk, risk == "collision")

    for i, a in enumerate(ends):
        for b in ends[i + 1:]:
            if a.reaction != b.reaction:
                continue
            width = len(a.seq)
            same = shared_bases(a.seq, b.seq)
            flipped = shared_bases(a.seq, reverse_complement(b.seq))
            risk = risk_of(max(same, flipped), width)
            record(pair_note(a, b, same, flipped), risk, risk == "collision")

    if watch:
        advisories.append(
            f"{watch} {_plural(watch, 'reaction')} {'has' if watch == 1 else 'have'} two ends "
            f"sharing more than {MAX_SHARED} {_plural(MAX_SHARED, 'base')} without being a "
            "near match; see lib.overhang_pairs()"
        )
    return failures, advisories


def overhang_findings(library) -> tuple[list[str], list[str]]:
    """``(failures, advisories)`` over every fused overhang the library leaves."""
    return findings_for(overhang_ends(library))
