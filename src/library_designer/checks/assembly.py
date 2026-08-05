"""Simulate the Golden Gate reaction a library implies, as QC.

Every other check reads the sequences we wrote down. This one puts them together: cut the
oligo with the enzyme, cut the destination vector, anneal the fused overhangs, ligate, and
look at what comes out. The digest here finds its own cut sites rather than trusting the
layout, so a spacer off by one base, an overhang drawn from the wrong place, or a tile
window that does not line up with its vector shows up as a product that is not the plasmid
we meant to build.

Two levels of checking, which keeps a large library cheap. Once per reaction the whole
plasmid is assembled from the frozen WT reference and compared against the starting vector
with that reference in place. Then, because the only thing that differs between members of
one reaction is the fragment core, every member is checked at coding-sequence level. Join its
released core to the arms the vector supplies and the result has to be that variant's
intended CDS, translating to its intended protein.

Without a destination vector there is no reaction to simulate, so nothing runs.
"""
from __future__ import annotations

from dataclasses import dataclass, field, replace

from ..regions import assemble, reverse_complement
from .motifs import ENZYME_SITES, cut_geometry
from .translation import translates_to


@dataclass
class Fragment:
    """One piece of a digest. ``left`` and ``right`` are the single-stranded 5' overhangs
    it carries, empty at an end of a linear molecule the enzyme did not cut. Both partners
    at a junction name it with the same bases, so ligating is ``a.core + b.seq``."""

    left: str
    core: str
    right: str

    @property
    def seq(self) -> str:
        return self.left + self.core + self.right


def _find_all(hay: str, needle: str) -> list[int]:
    out, i = [], hay.find(needle)
    while i != -1:
        out.append(i)
        i = hay.find(needle, i + 1)
    return out


def _arc(seq: str, a: int, b: int) -> str:
    """``seq[a:b]``, wrapping the origin when the span runs off the end."""
    return seq[a:b] if a <= b else seq[a:] + seq[:b]


def junctions(seq: str, enzyme: str, circular: bool = False) -> list[tuple[int, int]]:
    """The ``(start, end)`` span of the fused overhang each cut of ``enzyme`` leaves, in
    ascending order. A site too close to the end of a linear molecule to cut is skipped."""
    site = ENZYME_SITES[enzyme].upper()
    site_rc = reverse_complement(site)
    spacer, o = cut_geometry(enzyme)
    n = len(seq)
    scan = seq + seq[: len(site) - 1] if circular else seq   # a site spanning the origin
    starts = [i + len(site) + spacer for i in _find_all(scan, site)]
    if site_rc != site:
        starts += [j - spacer - o for j in _find_all(scan, site_rc)]
    if circular:
        return sorted({(s % n, (s + o) % n) for s in starts})
    return sorted((s, s + o) for s in starts if 0 <= s and s + o <= n)


def digest(seq: str, enzyme: str, circular: bool = False) -> list[Fragment]:
    """Cut ``seq`` with ``enzyme`` into ``Fragment`` pieces, 5'->3' along the molecule.

    A circular molecule with n cuts gives n pieces, each bounded by a junction on both
    sides. A linear one with n cuts gives n+1, the outermost two keeping the molecule's own
    ends. An uncut molecule comes back as a single piece with no overhangs."""
    j = junctions(seq, enzyme, circular)
    if not j:
        return [Fragment("", seq, "")]
    if circular:
        out = []
        for k, (a, b) in enumerate(j):
            c, d = j[(k + 1) % len(j)]
            out.append(Fragment(_arc(seq, a, b), _arc(seq, b, c), _arc(seq, c, d)))
        return out
    out = [Fragment("", seq[: j[0][0]], seq[j[0][0]:j[0][1]])]
    for k in range(1, len(j)):
        a, b = j[k - 1]
        c, d = j[k]
        out.append(Fragment(seq[a:b], seq[b:c], seq[c:d]))
    a, b = j[-1]
    out.append(Fragment(seq[a:b], seq[b:], ""))
    return out


def ligate(vector_pieces: list[Fragment], insert: Fragment) -> tuple[str | None, list[str]]:
    """Ligate ``insert`` into a digested vector, returning ``(product, issues)``.

    The overhangs decide what can join what: one vector piece has to present ``insert.left``
    on its right end and ``insert.right`` on its left. For a cut plasmid that is a single
    backbone and the product closes into a circle; for a linear molecule it is the two arms
    and the product is linear. Anything else means the reaction does not go, and the issue
    says which end failed."""
    issues: list[str] = []
    upstream = [f for f in vector_pieces if f.right == insert.left and f.right]
    downstream = [f for f in vector_pieces if f.left == insert.right and f.left]
    if not upstream:
        issues.append(
            f"nothing in the cut vector presents the insert's 5' overhang ({insert.left}), "
            f"so it cannot ligate (the vector offers "
            f"{sorted({f.right for f in vector_pieces if f.right}) or 'no overhangs'})"
        )
    if not downstream:
        issues.append(
            f"nothing in the cut vector presents the insert's 3' overhang ({insert.right}), "
            f"so it cannot ligate (the vector offers "
            f"{sorted({f.left for f in vector_pieces if f.left}) or 'no overhangs'})"
        )
    if issues:
        return None, issues
    if len(upstream) > 1 or len(downstream) > 1:
        issues.append("the insert's overhangs fit the cut vector in more than one place")
        return None, issues

    up, down = upstream[0], downstream[0]
    if up is down:                       # one backbone, so the product closes
        return up.core + insert.seq, issues
    if up.left or down.right:
        issues.append("the vector arms carry further cut sites, so the product is not one molecule")
        return None, issues
    return up.core + insert.seq + down.core, issues


@dataclass
class AssemblyResult:
    """The simulated outcome of one Golden Gate reaction.

    ``product`` is the plasmid it yields, built from the WT reference, so it is the clone you
    expect to sequence. It reads in the same frame as ``parent_vector`` and every
    ``assembled_product``, so any two can be compared directly. ``problems`` maps a
    member's name to what went wrong for it, in
    words; everything not in there rebuilt its own intended coding sequence. ``n_aligned``
    counts members whose product, aligned against the parent vector, differs only at the
    intended codon."""

    label: str
    enzyme: str
    product: str | None = None
    topology: str = "circular"
    n_members: int = 0
    n_aligned: int = 0
    issues: list[str] = field(default_factory=list)
    problems: dict[str, str] = field(default_factory=dict)

    @property
    def n_correct(self) -> int:
        """How many members rebuilt their own intended coding sequence, which is
        ``n_members`` less the ones listed in ``problems``."""
        return self.n_members - len(self.problems)

    @property
    def ok(self) -> bool:
        """True when the reaction gives a product, no member is in ``problems``, and
        ``issues`` is empty. A reaction-level issue fails it even when every member
        rebuilt correctly, and ``n_aligned`` does not enter into it."""
        return not self.issues and not self.problems and self.product is not None

    def __str__(self) -> str:
        head = f"{self.label}: {'assembles' if self.ok else 'FAILS'}"
        if self.product is not None:
            head += f", {len(self.product)} bp {self.topology} product"
        head += f", {self.n_correct}/{self.n_members} members rebuild their variant"
        if self.n_aligned:
            head += f", {self.n_aligned} aligned to the parent at the intended codon only"
        return head

    def __repr__(self) -> str:
        return f"<AssemblyResult {self}>"


def _released(construct: str, enzyme: str, label: str) -> tuple[Fragment | None, list[str]]:
    """The fragment a linear construct releases: the one piece cut on both ends."""
    pieces = digest(construct, enzyme, circular=False)
    inner = [f for f in pieces if f.left and f.right]
    if len(inner) != 1:
        return None, [
            f"{label}: digesting with {enzyme} releases {len(inner)} fragment(s) cut at both "
            f"ends (expected 1), so the insert does not come out cleanly"
        ]
    return inner[0], []


def _differing_span(a: str, b: str) -> tuple[int, int] | None:
    """``(lo, hi)`` bounding every position where two equal-length strings differ, or None
    when they are identical. Uses common prefix/suffix so a whole plasmid costs one pass."""
    from os.path import commonprefix

    lo = len(commonprefix([a, b]))
    if lo == len(a):
        return None
    suffix = len(commonprefix([a[::-1], b[::-1]]))
    return lo, len(a) - suffix


def align_to_parent(product: str, parent: str, cds_start: int, cds_len: int,
                    mut_index: int | None) -> str | None:
    """Align an assembled product against the parent vector carrying the WT reference.

    Returns None when the only difference is the intended mutation (or none at all, for the
    WT control), otherwise says what else changed. This is the end-to-end statement a
    single-mutant library rests on: the clone you get back differs from the plasmid you
    started with at one codon, and it is the codon you asked for.
    """
    if len(product) != len(parent):
        return (
            f"the product is {len(product)} bp against the parent vector's {len(parent)} bp, "
            "so it is not a clean substitution"
        )
    span = _differing_span(product, parent)
    if mut_index is None:                        # the WT control must come back identical
        if span is None:
            return None
        lo, hi = span
        return f"the WT product differs from the parent vector at {hi - lo} base(s) (from base {lo + 1})"
    codon_lo = cds_start + mut_index * 3
    codon_hi = codon_lo + 3
    if span is None:
        return f"the product is identical to the parent vector, so codon {mut_index + 1} was not changed"
    lo, hi = span
    if lo < codon_lo or hi > codon_hi:
        outside = (max(0, codon_lo - lo)) + (max(0, hi - codon_hi))
        where = f"base{'s' if hi - lo > 1 else ''} {lo + 1}-{hi} of the vector map"
        return (
            f"the product differs from the parent vector at {where}, {outside} of them outside "
            f"the intended codon {mut_index + 1}"
        )
    if not 0 <= codon_lo - cds_start < cds_len:
        return f"the intended codon {mut_index + 1} falls outside the coding sequence"
    return None


def _ligated_cds(reference: str, insert: Fragment, start: int, end: int, o: int,
                 pad_5: int = 0, pad_3: int = 0) -> str:
    """The coding sequence of the ligated product: the arms the vector keeps, joined to the
    overhang bases and core the fragment brings.

    The fragment's core covers ``[start, end)`` of the coding sequence and its overhangs sit
    just outside that, so the ``o`` bases before ``start`` come from ``insert.left``. How many
    of them are coding depends on where the cut fell: at an internal tile boundary all of
    them are, at a CDS end none are (the overhang is backbone), and an adaptor that spells
    only part of the overhang splits it between the two. ``min`` picks up all three cases, so
    only the coding part is taken from the fragment and the rest of the arm comes from
    the reference.

    ``pad_5`` / ``pad_3`` are backbone bases the adaptors re-supply, which ride inside the
    fragment's core between the overhang and the coding sequence. They are not coding, so
    they come off the core before any of this."""
    if pad_5 or pad_3:
        insert = replace(insert, core=insert.core[pad_5:len(insert.core) - pad_3 or None])
    coding_5 = min(start, o)                       # overhang bases that are coding sequence
    coding_3 = min(len(reference) - end, o)
    head = reference[: start - coding_5] + insert.left[o - coding_5:]
    tail = insert.right[:coding_3] + reference[end + coding_3:]
    return head + insert.core + tail


@dataclass
class _Plan:
    """Everything a library's Golden Gate reactions need, worked out once.

    Each entry of ``reactions`` is ``(label, vector, topology, start, end, members,
    wt_construct)``, where a member is ``(name, its CDS, the residue index it mutates, the
    molecule that gets ordered)``. ``parent`` is the starting plasmid with the frozen
    reference in place, the baseline every product is aligned against, and None when there is
    no real backbone (a tiled cassette). ``dead`` carries a result that stands in for the
    whole simulation when the library cannot be assembled at all."""

    reference: str
    enzyme: str
    overhang: int
    proteins: dict[str, str]
    reactions: list[tuple]
    parent: str | None = None
    cds_at: int | None = None
    dead: "AssemblyResult | None" = None


def _plan(library) -> _Plan | None:
    """The reactions this library implies, or None when there is nothing to assemble into."""
    import pandas as pd

    from ..layout.tiled import wt_oligo
    from ..layout.destination import resolve_insert_locus
    from ..layout.vector_io import assemble_vector, insert_offset, parent_region

    spec = library.spec
    reference = library.reference
    if reference is None:
        return None
    df = library.df
    proteins = dict(zip(df["name"].astype(str), df["protein"].astype(str)))
    tiles = getattr(library, "tiles", None)

    def _member(name, dna, idx, construct):
        return str(name), dna, (None if idx is None or pd.isna(idx) else int(idx)), construct

    idx_col = df["mut_index"] if "mut_index" in df.columns else [None] * len(df)
    if tiles is not None:
        params = library.tiled_params
        vec = spec.resolve_vector(params)
        enzyme, o = params.enzyme, params.overhang_len
        by_tile: dict[int, list] = {t.index: [] for t in tiles}
        for name, dna, ti, idx, oligo in zip(
            df["name"], df["variable_dna"], df["tile"], idx_col, df["oligo"]
        ):
            if isinstance(dna, str) and isinstance(oligo, str) and ti in by_tile:
                by_tile[int(ti)].append(_member(name, dna, idx, oligo))
        reactions = [(
            f"tile{t.index}", t.vector, t.topology, t.start, t.end, by_tile[t.index],
            wt_oligo(reference, t, params), 0, 0,      # a tile's window has no padding
        ) for t in tiles]
    else:
        vec = spec.vector
        if vec is None or not (spec.adaptor_5 or spec.adaptor_3):
            # No adaptors means nothing carries the enzyme sites, so there is no reaction to
            # simulate. checks/vector.py already reports that as an advisory.
            return None
        from ..layout.destination import build_destination

        enzyme, o = vec.enzyme, cut_geometry(vec.enzyme)[1]
        try:
            # strict=False: an overhang collision is a QC finding of its own (see
            # checks/overhangs.py), and simulating anyway says concretely what the reaction
            # does with it rather than replacing every member's result with one message.
            dv = build_destination(library, strict=False)
        except ValueError:
            # The same fault is already reported in full by checks/vector.py, which owns the
            # adaptor-against-the-plasmid checks. Restating it here printed one root cause as
            # two findings, so point at that one instead of echoing its message.
            return _Plan(reference, enzyme, o, proteins, [], dead=AssemblyResult(
                "destination", enzyme,
                issues=["no assembly was simulated, because the destination vector could "
                        "not be built; see the adaptor/destination-vector issue above"],
            ))
        a5, a3 = spec.adaptor_5, spec.adaptor_3
        members = [
            _member(n, v, i, assemble(a5, v, a3))
            for n, v, i in zip(df["name"], df["variable_dna"], idx_col) if isinstance(v, str)
        ]
        reactions = [(
            "destination", dv.sequence, dv.topology, dv.start, dv.end, members,
            assemble(a5, reference, a3), dv.dest.pad_5, dv.dest.pad_3,
        )]

    # The plasmid each reaction should yield, when a real backbone is in play. Without one the
    # "vector" is a bare coding-region cassette, so there is no plasmid to compare and the
    # per-member coding-sequence check carries the simulation on its own.
    parent = cds_at = None
    if vec is not None:
        dest = resolve_insert_locus(spec, library.tiled_params)
        # parent_region puts back any bases the locus was widened by, so the parent is the
        # real starting plasmid and not one short by the adaptors' padding.
        parent = assemble_vector(dest, parent_region(dest, reference))
        cds_at = insert_offset(dest) + dest.pad_5
    return _Plan(reference, enzyme, o, proteins, reactions, parent, cds_at)


def parent_vector(library) -> str:
    """The starting plasmid with the frozen reference CDS in place: what a WT clone should
    be, and the baseline ``assembled_product`` is aligned against."""
    plan = _plan(library)
    if plan is None or plan.parent is None:
        raise ValueError(
            "No parent plasmid to build. Set spec.starting_vector to the plasmid you clone "
            "into, and codon-optimize the library first."
        )
    return plan.parent


def assembled_products(library, names=None):
    """Yield ``(name, product, mut_index)`` for every member that assembles.

    The products come back in ``parent_vector``'s frame, so any of them can be diffed
    against it directly. The vector is digested once per reaction and the frame
    rotation taken once from the WT product, so walking a whole library is cheap. Members
    that cannot be released or ligated are skipped, since ``check()`` is where those are
    reported. ``names`` restricts the walk to a set of member names."""
    plan = _plan(library)
    if plan is None or plan.parent is None:
        raise ValueError(
            "Nothing to assemble: needs a starting vector, adaptors carrying the enzyme "
            "sites, and a codon-optimized library."
        )
    for label, vector, topology, _s, _e, members, wt_construct, _p5, _p3 in plan.reactions:
        chosen = [m for m in members if names is None or m[0] in names]
        if not chosen:
            continue
        circular = topology == "circular"
        pieces = digest(vector, plan.enzyme, circular)
        wt_insert, _ = _released(wt_construct, plan.enzyme, label)
        if wt_insert is None:
            continue
        # Take the rotation from the WT product, which is a rotation of the parent by
        # definition. A variant's product is not, so finding it in the parent would fail.
        wt_product, _ = ligate(pieces, wt_insert)
        if wt_product is None:
            continue
        rotation = max((plan.parent + plan.parent).find(wt_product), 0) if circular else 0
        for name, _cds, mut_index, construct in chosen:
            insert, _ = _released(construct, plan.enzyme, name)
            if insert is None:
                continue
            product, _ = ligate(pieces, insert)
            if product is not None:
                yield name, _rotate(product, rotation), mut_index


def assembled_product(library, name: str) -> str:
    """The plasmid member ``name`` assembles into, read in the same frame as
    ``parent_vector`` so the two can be compared.

    This is the simulated clone: its oligo digested, ligated into the cut destination
    vector, and rotated onto the parent's origin. Raises if the member has nothing to
    assemble (a variant that failed optimization, or a tiled library's global ``WT`` row,
    which rides on no oligo)."""
    for got, product, _ in assembled_products(library, names={name}):
        if got == name:
            return product
    raise ValueError(
        f"No member named {name!r} with something to assemble. Names come from "
        "library.df['name']. A tiled library's global WT row rides on no oligo; ask for "
        "a per-tile control instead (WT_Tile_0, ...) or use parent_vector() for the WT "
        "plasmid."
    )


def simulate(library) -> list[AssemblyResult]:
    """Simulate every Golden Gate reaction this library implies, one per destination vector.

    Returns an empty list when there is nothing to assemble into: no destination vector, or
    no frozen reference to assemble against. A sequence set has no shared reference, so it
    falls in the second case even after ``codon_optimize()``."""
    plan = _plan(library)
    if plan is None:
        return []
    if plan.dead is not None:
        return [plan.dead]
    reference, enzyme, o = plan.reference, plan.enzyme, plan.overhang
    parent, cds_at, proteins = plan.parent, plan.cds_at, plan.proteins
    reactions = plan.reactions

    results = []
    for label, vector, topology, start, end, members, wt_construct, pad_5, pad_3 in reactions:
        r = AssemblyResult(label, enzyme, topology=topology, n_members=len(members))
        results.append(r)

        wt_insert, issues = _released(wt_construct, enzyme, label)
        r.issues += issues
        if wt_insert is None:
            continue
        circular = topology == "circular"
        vector_frags = digest(vector, enzyme, circular)

        rotation = 0
        if parent is not None:
            product, lig_issues = ligate(vector_frags, wt_insert)
            r.issues += [f"{label}: {m}" for m in lig_issues]
            r.product = product
            if product is None:
                continue
            if circular:
                # The product reads from wherever the backbone piece began; line it up with
                # the parent so positions mean the same thing in both. Every product this
                # module hands out is in the parent's frame, so they can all be diffed
                # against each other without the caller working out a rotation.
                rotation = (parent + parent).find(product)
                if rotation < 0 or len(product) != len(parent):
                    r.issues.append(
                        f"{label}: the assembled product ({len(product)} bp) is not the starting "
                        f"vector carrying the reference CDS ({len(parent)} bp)"
                    )
                    continue
                r.product = _rotate(product, rotation)
            elif product != parent:
                r.issues.append(
                    f"{label}: the assembled product does not match the starting vector "
                    "carrying the reference CDS"
                )
                continue

        # Every member of this reaction differs from the WT only in the fragment core, so
        # this is the ligation, member by member.
        for name, variant_cds, mut_index, construct in members:
            insert, issues = _released(construct, enzyme, name)
            if insert is None:
                r.problems[name] = "does not release a clean fragment when digested"
                continue
            if (insert.left, insert.right) != (wt_insert.left, wt_insert.right):
                end5 = "5'" if insert.left != wt_insert.left else "3'"
                r.problems[name] = (
                    f"changes a base inside the {end5} fused overhang, so the fragment no longer "
                    "anneals to the cut vector"
                )
                continue
            got = _ligated_cds(reference, insert, start, end, o, pad_5, pad_3)
            protein = proteins.get(name)
            if got != variant_cds:
                r.problems[name] = "assembles to a different coding sequence than designed"
                continue
            if protein is not None and not translates_to(got, protein):
                r.problems[name] = "assembles to a coding sequence that translates to another protein"
                continue
            if parent is None:
                continue
            member_product, _ = ligate(vector_frags, insert)
            if member_product is None:
                r.problems[name] = "does not ligate into the cut vector"
                continue
            framed = _rotate(member_product, rotation)
            bad = align_to_parent(framed, parent, cds_at, len(reference), mut_index)
            if bad:
                r.problems[name] = bad
            else:
                r.n_aligned += 1
    return results


def _rotate(product: str, rotation: int) -> str:
    """A product read from ``rotation`` bases into the parent, re-read from the parent's own
    origin, so the two line up."""
    return product if not rotation else product[-rotation:] + product[:-rotation]


def check_assembly(library) -> dict:
    """Assembly-simulation findings for the QC report, members grouped by what went wrong."""
    results = simulate(library)
    issues: list[str] = []
    checked = correct = aligned = 0
    for r in results:
        issues += r.issues
        checked += r.n_members
        correct += r.n_correct
        aligned += r.n_aligned
        by_reason: dict[str, list[str]] = {}
        for name, reason in r.problems.items():
            by_reason.setdefault(reason, []).append(name)
        for reason, names in by_reason.items():
            issues.append(
                f"{r.label}: {len(names)} member(s) {reason} ({', '.join(names[:5])}"
                + (", ..." if len(names) > 5 else "") + ")"
            )
    return {
        "assembly_issues": issues,
        "assembly_checked": checked,
        "assembly_correct": correct,
        "assembly_aligned": aligned,
    }
