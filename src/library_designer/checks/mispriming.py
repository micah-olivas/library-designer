"""Whether the constant flanks a pool is amplified with prime anywhere they should not.

A pooled library is pulled apart with the constant regions its oligos carry. A tiled oligo
carries a tile-specific primer pair, and a one-oligo pool carries the shared adaptors. Either
way, a handle only amplifies what you asked for if it has one binding site. When its 3' end
also anneals inside a variable region, inside another oligo's flank, or in the destination
vector, the polymerase extends from there too and the tube holds a product the design does
not account for.

Extension starts at the 3' end, so that is what this scores. For each handle we take the
longest duplex reaching back from its 3' base that it could form where it should not, on
either strand, and grade the design on how many bases of it pair. The whole handle occurring
a second time is a second binding site, and it fails QC. Anything shorter is reported without
failing, since whether it primes depends on the annealing temperature you run at.

A duplex may carry internal mismatches, up to ``MISMATCH_BUDGET``. A primer whose last five
bases pair, then one mismatches, then fourteen more pair is a 19-of-20 duplex that primes
readily, and reading only exact runs would score that as five bases and pass it. The last
``TERMINAL_CLEAN`` bases have to pair whatever the budget is: a mismatch at or beside the 3'
base leaves nothing for the polymerase to extend from. Mismatched duplexes are advisories
only, never failures, so no design's verdict rests on the mismatch weighting.

Where we look:

- the frozen reference CDS, which is the variable region every member is built on
- a window around each member's own mutated codon, the only sequence a member has that the
  reference does not (a run that misses the codon is already in the reference, so it is
  reported there instead)
- the constant flanks of the other oligos in the pool, so a primer that anneals to another
  tile's flank is caught even though the set it came from is orthogonal on paper
- the destination-vector backbone, the part of the plasmid outside the insert locus

One limit is worth knowing: the backbone is read as the arc outside the insert locus, so a
site spelled half in the backbone and half in the insert is not looked for.
"""
from __future__ import annotations

from bisect import bisect_right
from dataclasses import dataclass

from ..regions import reverse_complement

# Shortest 3' run worth reporting. Below about 12 bases an exact run is too short to prime on
# its own at a normal annealing temperature, and short runs turn up by chance often enough to
# bury the ones that matter.
MIN_ANNEAL = 12
# A run this long anneals at the temperatures a pool is amplified at, so it is named
# individually rather than counted.
HIGH_ANNEAL = 15
# A handle shorter than this is not an amplification primer, and a run that short turns up by
# chance in any sequence, so a finding about it would mean nothing. A handle between this and
# MIN_ANNEAL is checked for occurring in full, which is the only run it has. One the design
# itself carries is reported as unscored rather than dropped in silence; one you pass in is
# refused outright, since a primer you asked about must not come back silently unexamined.
MIN_HANDLE = 8
# Internal mismatches a duplex may carry and still be reported. One covers the case exact
# matching misses worst, a long duplex broken by a single base, without turning every few kb of
# backbone into findings. Raise it to look harder; a mismatched duplex is never a failure.
MISMATCH_BUDGET = 1
# Bases at the 3' end that have to pair whatever the budget is. Taq extends from a paired 3'
# terminus, and a mismatch in the last few bases costs orders of magnitude of extension, so a
# duplex that is imperfect there is not a priming site.
TERMINAL_CLEAN = 3
# Shortest exact stretch used to find a candidate duplex. The scan is seeded on the handle's
# own k-mers, and k is chosen so that any duplex worth reporting has to contain one of them
# (see tolerant_sites), with this as the floor.
SEED_MIN = 4


@dataclass(frozen=True)
class Handle:
    """One constant region a pool is amplified with, written 5'->3' as you would order it.

    ``seq`` is the primer as ordered, which for a 3' adaptor means the reverse complement of
    what sits on the oligo: the base at its 3' end is then the adaptor base nearest the
    variable region, which is the base extension starts from. ``site`` names the region this
    handle is meant to bind, so its intended binding site is not reported as a finding.
    """

    label: str
    seq: str
    site: str = ""


@dataclass(frozen=True)
class Part:
    """One member's stretch inside a packed region: ``[start, end)`` on the packed sequence.

    Thousands of short per-member windows are searched as one string rather than one at a
    time, and this is what turns a hit position back into the member it belongs to. ``offset``
    is where the stretch begins on the frozen reference, which is what lets a duplex found in a
    member be scored against the same place on the WT it came from. ``focus`` is where the
    member's own mutated codon sits inside the stretch, the three bases that are the only reason
    this stretch is searched at all.
    """

    label: str
    start: int
    end: int
    offset: int = 0
    focus: int = 0


@dataclass(frozen=True)
class Region:
    """A stretch of sequence a handle must not prime in.

    ``kind`` is ``"variable"`` (coding sequence the library varies), ``"flank"`` (the
    constant regions of the oligos), or ``"vector"`` (the destination backbone). ``end`` is
    ``"5'`` or ``"3'"`` on a flank, which end of the oligo it sits at, and says where its
    boundary with the variable region is. ``parts`` is set when the region packs several
    members' sequences into one string, separated by ``N`` so no run can span two of them.
    """

    label: str
    seq: str
    kind: str = "variable"
    end: str = ""
    parts: tuple[Part, ...] = ()

    def part_at(self, pos: int) -> Part | None:
        """The member's stretch a hit at ``pos`` landed in, or None on an unpacked region."""
        if not self.parts:
            return None
        i = bisect_right([p.start for p in self.parts], pos) - 1
        return self.parts[max(i, 0)]

    def locate(self, pos: int) -> tuple[str, int]:
        """``(label, position)`` for a hit at ``pos``, resolved to the member it landed in."""
        part = self.part_at(pos)
        return (self.label, pos) if part is None else (part.label, pos - part.start)


@dataclass(frozen=True)
class PrimingSite:
    """One duplex a handle's 3' end could form where it should not.

    ``length`` is how many bases pair and ``aligned`` how long the duplex is, the two being
    equal when it is a clean run. ``mismatches`` gives each mismatched base as its distance
    back from the 3' base, so ``(5,)`` is the sixth base in. ``strand`` is ``"+"`` when the
    duplex reads on the region as written, so the handle anneals to the other strand and
    extends along the region, and ``"-"`` when it reads on the complement. ``where`` names the
    member or sub-region it landed in and ``position`` is the index within it. ``matched`` is
    the handle's own 3' window, the bases involved.
    """

    handle: Handle
    region: Region
    length: int
    strand: str
    where: str
    position: int
    matched: str
    aligned: int = 0
    mismatches: tuple[int, ...] = ()

    def __post_init__(self) -> None:
        if not self.aligned:                       # a clean run: the duplex is the match
            object.__setattr__(self, "aligned", self.length)

    @property
    def full(self) -> bool:
        """The whole handle pairs here, so this is a second binding site.

        Only ever true of a clean duplex: every base pairing means none mismatched."""
        return self.length == len(self.handle.seq)

    @property
    def risk(self) -> str:
        """``"collision"`` when the whole handle pairs, ``"high"`` when enough bases pair to
        anneal at a normal annealing temperature, else ``"watch"``."""
        if self.full:
            return "collision"
        return "high" if self.length >= HIGH_ANNEAL else "watch"


def _find_all(haystack: str, needle: str) -> list[int]:
    out, i = [], haystack.find(needle)
    while i != -1:
        out.append(i)
        i = haystack.find(needle, i + 1)
    return out


def anneal_sites(handle: str, region: str, min_len: int) -> list[tuple[int, str, int]]:
    """``(length, strand, position)`` for every 3'-anchored run of ``handle`` at least
    ``min_len`` bases long that occurs in ``region``, longest first.

    One site per 3' end: a 14-base run is reported once, not again as its 13- and 12-base
    suffixes. ``position`` is where the run starts in ``region``, whichever strand it reads
    on.
    """
    hits: list[tuple[int, str, int]] = []
    anchored: set[tuple[str, int]] = set()
    for length in range(len(handle), min_len - 1, -1):
        run = handle[-length:]
        for strand, needle in (("+", run), ("-", reverse_complement(run))):
            for pos in _find_all(region, needle):
                # The handle's 3' base sits at the right end of a "+" run and the left end of
                # a "-" one, so that base identifies the site across lengths.
                anchor = pos + length - 1 if strand == "+" else pos
                if (strand, anchor) in anchored:
                    continue
                anchored.add((strand, anchor))
                hits.append((length, strand, pos))
    return hits


def _extend(handle: str, template: str, anchor: int, budget: int,
            terminal_clean: int) -> tuple[int, tuple[int, ...], int] | None:
    """The duplex the handle forms with its 3' base sitting on ``template[anchor]``.

    ``(paired, mismatch offsets, aligned)``, walking back from the 3' base and stopping when
    the mismatch budget runs out or the template does. The duplex ends on a paired base, since
    a trailing mismatch adds nothing to it. None when a mismatch falls in the last
    ``terminal_clean`` bases, where there is no paired 3' terminus to extend from.
    """
    paired, aligned, mismatches = 0, 0, []
    best: tuple[int, tuple[int, ...], int] | None = None
    for d in range(len(handle)):
        j = anchor - d
        if j < 0:
            break
        if handle[len(handle) - 1 - d] == template[j]:
            paired += 1
            aligned = d + 1
            best = (paired, tuple(mismatches), aligned)
        else:
            if d < terminal_clean:
                return None
            mismatches.append(d)
            if len(mismatches) > budget:
                break
    return best


def tolerant_sites(handle: str, region: str, min_paired: int, budget: int,
                   terminal_clean: int = TERMINAL_CLEAN) -> list[tuple[int, tuple[int, ...], int, str, int]]:
    """``(paired, mismatch offsets, aligned, strand, position)`` for every 3'-anchored duplex
    of at least ``min_paired`` paired bases the handle could form in ``region``, most paired
    first.

    Mismatches are allowed inside the duplex, up to ``budget``, so this finds what exact
    matching cannot. With ``budget`` 0 it reports the same sites ``anneal_sites`` does, at more
    cost, which is why the exact scan stays the default path.

    The scan is seeded on the handle's own k-mers. A duplex of ``min_paired`` paired bases
    broken by at most ``budget`` mismatches splits into ``budget + 1`` clean stretches, so the
    longest of them is at least ``ceil(min_paired / (budget + 1))`` bases; seeding at that
    length cannot miss a duplex that qualifies. One anchor is scored once, whichever seed
    reaches it, since the walk back from a 3' base does not depend on how it was found.
    """
    n = len(handle)
    k = min(n, max(SEED_MIN, -(-min_paired // (budget + 1))))
    seeds: dict[str, list[int]] = {}
    for off in range(n - k + 1):
        seeds.setdefault(handle[off:off + k], []).append(off)

    out = []
    for strand, template in (("+", region), ("-", reverse_complement(region))):
        scored: set[int] = set()
        for seed, offsets in seeds.items():
            for hit in _find_all(template, seed):
                for off in offsets:
                    anchor = hit - off + n - 1
                    if not 0 <= anchor < len(template) or anchor in scored:
                        continue
                    scored.add(anchor)
                    got = _extend(handle, template, anchor, budget, terminal_clean)
                    if got is None or got[0] < min_paired:
                        continue
                    paired, mismatches, aligned = got
                    # Back to the region's own coordinates. On the minus strand the template is
                    # its reverse complement, so the duplex's leftmost base in the region is the
                    # handle's 3' base.
                    pos = (anchor - aligned + 1) if strand == "+" else len(region) - 1 - anchor
                    out.append((paired, mismatches, aligned, strand, pos))
    out.sort(key=lambda h: (-h[0], len(h[1]), h[4]))
    return out


# --- what to check, and against what -----------------------------------------

def handles(library) -> list[Handle]:
    """The constant regions this library is amplified with.

    A tiled library gives two per tile, the primer pair drawn from the orthogonal set. A
    standard library gives its adaptors, the constant flanks every oligo in the pool carries.
    The 3' adaptor is returned reverse-complemented, which is the primer you would order for
    it. A library with neither has nothing to amplify with, so the list comes back empty.
    """
    tiles = getattr(library, "tiles", None)
    if tiles is not None:
        out: list[Handle] = []
        for t in tiles:
            out.append(Handle(f"tile{t.index} fwd", t.fwd.upper(), f"tile{t.index} 5' flank"))
            out.append(Handle(f"tile{t.index} rev", t.rev.upper(), f"tile{t.index} 3' flank"))
        return out

    spec = library.spec
    out = []
    if spec.adaptor_5:
        out.append(Handle("adaptor_5", spec.adaptor_5.upper(), "adaptor_5"))
    if spec.adaptor_3:
        out.append(Handle("adaptor_3", reverse_complement(spec.adaptor_3.upper()), "adaptor_3"))
    return out


def extra_handles(extra) -> list[Handle]:
    """Sequences passed in to check alongside the design's own flanks, cleaned and numbered.

    Whitespace comes out and the sequence is uppercased, so a primer pasted with spaces or in
    lowercase reads the same as one typed in. Anything that cannot be scored is refused here
    rather than passed over: a primer you asked about should never come back silently
    unexamined, which is the one answer that reads like a clean result and is not one.
    """
    out = []
    for i, raw in enumerate(extra, start=1):
        seq = "".join(str(raw).split()).upper()
        bad = sorted(set(seq) - set("ACGT"))
        if bad:
            raise ValueError(
                f"extra primer {i} ({seq!r}) contains non-ACGT character(s) {bad}. Priming "
                "sites are found by exact matching, so a degenerate primer would be scored as "
                "if its ambiguous bases match nothing. Give the sequence you would order."
            )
        if len(seq) < MIN_HANDLE:
            raise ValueError(
                f"extra primer {i} ({seq!r}) is {len(seq)} base(s), too short to score as an "
                f"amplification handle (at least {MIN_HANDLE}). A run that short turns up by "
                "chance in any sequence, so a finding would say nothing about whether the "
                "primer is specific. Pass the primer you would order."
            )
        out.append(Handle(f"extra {i}", seq))
    return out


def variable_regions(library, width: int) -> list[Region]:
    """The variable sequence a handle could prime in.

    For a library built on one frozen reference that is the reference itself, plus a window
    around each member's own codon. Every member is the reference with one codon swapped, so
    a run that does not touch that codon is in the reference already and is found there; the
    windows are what carry the bases a member has and the reference does not. ``width`` is how
    far a window reaches either side of the codon, one base short of the longest handle, so
    every run that can touch the codon fits inside one.

    A ``SequenceSet`` has no shared reference (every member is a different gene), so each
    member's coding sequence goes in whole.
    """
    df = library.df
    if "variable_dna" not in df.columns:
        return []
    reference = getattr(library, "reference", None)
    names = [str(n) for n in df["name"]]
    seqs = list(df["variable_dna"])

    if not reference:
        packed = [(n, s.upper(), 0, 0) for n, s in zip(names, seqs) if isinstance(s, str)]
        return [_pack("member sequences", packed, "variable")] if packed else []

    out = [Region("reference CDS", reference.upper(), "variable")]
    if "mut_index" not in df.columns:
        return out
    import pandas as pd

    windows: dict[tuple[str, int], tuple[str, str, int, int]] = {}
    for name, dna, idx in zip(names, seqs, df["mut_index"]):
        if not isinstance(dna, str) or pd.isna(idx):
            continue
        codon = int(idx) * 3
        a, b = max(codon - width, 0), min(codon + 3 + width, len(dna))
        window = dna[a:b].upper()
        # One entry per distinct window at a given position: a scan repeats the same
        # neighbourhood at every substitution it makes there, and the whole pool is searched as
        # one string. Keyed on the position too, so two positions that happen to spell the same
        # window keep their own, each scored against its own stretch of the reference.
        windows.setdefault((window, a), (name, window, a, codon - a))
    if windows:
        out.append(_pack("variant codons", list(windows.values()), "variable"))
    return out


def _pack(label: str, parts: list[tuple[str, str, int, int]], kind: str) -> Region:
    """Several members' sequences as one searchable region, separated by ``N``.

    ``N`` is outside the ACGT a handle is written in, so no run can span two members and the
    join cannot invent a site that no oligo carries. Each part carries where its sequence starts
    on the reference and where its mutated codon sits inside it, both 0 when there is no
    reference to relate it to.
    """
    seqs, spans, at = [], [], 0
    for name, seq, offset, focus in parts:
        spans.append(Part(name, at, at + len(seq), offset, focus))
        seqs.append(seq)
        at += len(seq) + 1
    return Region(label, "N".join(seqs), kind, parts=tuple(spans))


def flank_regions(library) -> list[Region]:
    """The constant regions of the oligos in the pool, one per end per tile.

    Each is the run of an oligo outside its tile window, exactly as ``assemble_oligo`` lays
    it down: the primer, the pad, the recognition site, the spacer, and the fused overhang.
    Every handle is compared against all of them, its own included; what is passed over is the
    one full-length match a primer has on its own flank, which is where it is meant to bind.
    That is what catches a primer pulling a second tile out of the pool, or annealing back on
    the other end of its own oligo.

    A standard library's flanks are its two adaptors, so each is checked against the other.
    """
    tiles = getattr(library, "tiles", None)
    if tiles is None:
        spec = library.spec
        return [Region(label, seq.upper(), "flank", end)
                for label, seq, end in (("adaptor_5", spec.adaptor_5, "5'"),
                                        ("adaptor_3", spec.adaptor_3, "3'"))
                if seq]

    from ..layout.tiled import tile_contexts
    from .motifs import recognition_site

    params = library.tiled_params
    reference = library.reference
    if reference is None:
        return []
    rec = recognition_site(params.enzyme).upper()
    rec_rc = reverse_complement(rec)
    out = []
    for t in tiles:
        ctx5, ctx3 = tile_contexts(reference, t.start, t.end, params)
        out.append(Region(f"tile{t.index} 5' flank",
                          (t.lead + rec + params.spacer_5 + ctx5).upper(), "flank", "5'"))
        out.append(Region(f"tile{t.index} 3' flank",
                          (ctx3 + params.spacer_3 + rec_rc + t.trail).upper(), "flank", "3'"))
    return out


def vector_regions(library) -> list[Region]:
    """The destination-vector backbone, the part of the plasmid outside the insert locus.

    The coding region is left out because it is the reference, which is checked on its own,
    and every tile's vector carries the same backbone, so one region covers the whole design.
    A circular plasmid gives one region, the arc running from the 3' end of the insert round
    the origin and back to its 5' end. A linear one gives the two flanks separately, since
    joining them would make a junction the molecule does not have.

    Empty when the spec names no plasmid, or when the insert cannot be located, which the
    vector and tiling checks report on their own.
    """
    spec = library.spec
    params = getattr(library, "tiled_params", None)
    vec = spec.resolve_vector(params)
    if vec is None:
        return []
    from ..layout.vector_io import backbone, locating_kwargs, resolve_destination

    try:
        dest = resolve_destination(vec.path, **locating_kwargs(spec, params))
    except (ValueError, OSError, ImportError):
        return []
    pieces = backbone(dest)
    if len(pieces) == 1:
        return [Region("vector backbone", pieces[0], "vector")]
    return [Region(f"vector backbone {end}", seq, "vector")
            for end, seq in zip(("5'", "3'"), pieces)]


def _at_inner_edge(region: Region, aligned: int, strand: str, pos: int) -> bool:
    """True if a duplex reaches the flank's boundary with the variable region, in the
    orientation that primes into it.

    That is where an amplification primer for the pool is meant to sit: a forward primer ends
    on the last base of the 5' flank, and a reverse primer, read on the top strand, starts on
    the first base of the 3' flank.
    """
    if region.end == "5'":
        return strand == "+" and pos + aligned == len(region.seq)
    if region.end == "3'":
        return strand == "-" and pos == 0
    return False


def _new_to_the_variant(handle: str, region: Region, pos: int, aligned: int, paired: int,
                        strand: str, reference: str, rc_reference: str, budget: int) -> bool:
    """Whether a duplex found in a member's window is something its mutation brought about.

    A member is the reference with one codon swapped, so most of what a handle finds in a window
    it finds on the WT as well, and reporting that once per member that shares the neighbourhood
    would bury the report. Two things have to hold for the member's own to count. The duplex has
    to cover the mutated codon, since everything else in the window is WT sequence and the
    reference hit speaks for it. And it has to pair more bases than the same duplex pairs at the
    same place on the reference, so a mutation that spoils a WT duplex does not read as making a
    new one.

    An exact run is settled by looking the stretch up in the reference, which is quicker and
    also covers it turning up anywhere else on the WT.
    """
    part = region.part_at(pos)
    if part is None:
        return True
    local = pos - part.start
    if not (local < part.focus + 3 and part.focus < local + aligned):
        return False              # WT sequence either side of the codon; the reference has it
    chunk = region.seq[pos:pos + aligned]
    if chunk in reference or reverse_complement(chunk) in reference:
        return False
    if not budget:
        return True
    # Where this window sits on the reference, so the WT's own bases can be scored in the frame
    # the duplex was found in. The handle's 3' base is the right-hand end of a "+" duplex and the
    # left-hand end of a "-" one.
    at = part.offset + (pos - part.start)
    if strand == "+":
        anchor = at + aligned - 1
        template = reference
    else:
        anchor = len(reference) - 1 - at
        template = rc_reference
    if not 0 <= anchor < len(template):
        return True
    wt = _extend(handle, template, anchor, budget, TERMINAL_CLEAN)
    return wt is None or wt[0] < paired


def priming_sites(library, extra=(), min_anneal: int = MIN_ANNEAL,
                  mismatches: int = MISMATCH_BUDGET) -> list[PrimingSite]:
    """Every duplex a constant flank's 3' end could form where it should not, worst first.

    ``extra`` takes sequences to check alongside the design's own flanks, for a primer you
    plan to amplify with that the library does not carry (a nested primer, or a tailed primer
    whose 3' end is the adaptor). Given 5'->3' as ordered. On a design with one pair of flanks
    a supplied primer is taken to bind the flank its 3' end lands on, so the site it is for is
    not reported back at you. One too short to score, or carrying anything but A/C/G/T, is
    refused with the reason rather than passed over (see ``extra_handles``).

    ``min_anneal`` is the fewest paired bases to report. A handle shorter than that is checked
    for pairing in full, which is the only duplex it has. ``mismatches`` is how many internal
    mismatches a duplex may carry; 0 restricts this to exact runs and takes the faster path.
    """
    # The design's own flanks can be shorter than a primer (an adaptor is often just the
    # recognition site, a spacer, and the overhang), so a short one is left out and named in
    # mispriming_findings. One passed in is refused, see extra_handles.
    hs = [h for h in handles(library) if len(h.seq) >= MIN_HANDLE] + extra_handles(extra)
    if not hs:
        return []

    width = max(len(h.seq) for h in hs) - 1
    regions = variable_regions(library, width) + flank_regions(library) + vector_regions(library)
    reference = (getattr(library, "reference", None) or "").upper()
    rc_reference = reverse_complement(reference)
    # A primer you pass in has no site of its own on record, so the flank it ends in is taken to
    # be what it binds. Only when the design has a single flank pair, where that is the one
    # thing it can be: in a tiled pool a run reaching a tile's inner edge is a primer that
    # amplifies that tile, which is the finding.
    one_pair = len([r for r in regions if r.kind == "flank"]) <= 2

    out: list[PrimingSite] = []
    for h in hs:
        for r in regions:
            least = min(min_anneal, len(h.seq))
            if mismatches:
                hits = tolerant_sites(h.seq, r.seq, least, mismatches)
            else:
                hits = [(paired, (), paired, strand, pos)
                        for paired, strand, pos in anneal_sites(h.seq, r.seq, least)]
            for paired, mism, aligned, strand, pos in hits:
                if paired == len(h.seq) and r.label == h.site:
                    continue          # the handle's own binding site, which is why it is there
                if (not h.site and one_pair and r.kind == "flank"
                        and _at_inner_edge(r, aligned, strand, pos)):
                    continue          # a supplied primer, sitting where it is meant to
                if r.label == "variant codons" and not _new_to_the_variant(
                    h.seq, r, pos, aligned, paired, strand, reference, rc_reference, mismatches
                ):
                    continue          # the WT does this too, and is reported for it
                where, at = r.locate(pos)
                out.append(PrimingSite(h, r, paired, strand, where, at,
                                       h.seq[-aligned:], aligned, mism))
    out.sort(key=lambda s: (-s.length, len(s.mismatches), s.handle.label, s.where, s.position))
    return out


# --- reporting ----------------------------------------------------------------

def site_note(site: PrimingSite) -> str:
    """What to say about one priming site."""
    at = f"{site.where} at {site.position} ({site.strand} strand)"
    if site.full:
        return (f"{site.handle.label} ({site.handle.seq}) occurs in full in {at}, so it has a "
                "second binding site and cannot amplify selectively")
    if not site.mismatches:
        return (f"{site.handle.label} has its last {site.length} bases ({site.matched}) in {at}, "
                "which its 3' end can extend from")
    # Mismatches are given as how far in from the 3' base they sit, counting that base as 1,
    # which is the number that says whether the polymerase can extend.
    # Distances are 1-based (the 3' base itself is 1), and "base" agrees with the last number
    # in the list, not with how many mismatches there are.
    distances = [d + 1 for d in site.mismatches]
    where = ", ".join(str(d) for d in distances)
    n = len(distances)
    return (f"{site.handle.label} pairs {site.length} of its last {site.aligned} bases "
            f"({site.matched}) in {at}, with {n} mismatch{'' if n == 1 else 'es'} "
            f"{where} base{'' if distances == [1] else 's'} in from its 3' end, which its 3' "
            "end can extend from")


def mispriming_table(library, extra=(), min_anneal: int = MIN_ANNEAL,
                     mismatches: int = MISMATCH_BUDGET):
    """One row per duplex a constant flank could form where it should not, worst first.

    ``handle`` names the flank and ``handle_seq`` is it as ordered, so a 3' adaptor appears
    reverse-complemented. ``paired`` is how many bases pair, ``aligned`` how long the duplex is,
    and ``mismatches`` how many bases inside it do not pair, given as their distance in from the
    3' base. ``matched`` is the handle's 3' window, the bases involved. ``where`` is the member
    or sub-region the site sits in and ``position`` the index within it. ``risk`` is
    ``"collision"`` when the whole handle pairs, ``"high"`` when enough bases pair to anneal at
    a normal annealing temperature, else ``"watch"``.

    An empty frame with the same columns comes back when there is nothing to check, so a
    notebook cell does not have to branch. See ``priming_sites`` for the arguments.
    """
    import pandas as pd

    cols = ["handle", "handle_seq", "region", "kind", "where", "position", "strand",
            "paired", "aligned", "mismatches", "matched", "risk", "note"]
    rows = [{
        "handle": s.handle.label, "handle_seq": s.handle.seq, "region": s.region.label,
        "kind": s.region.kind, "where": s.where, "position": s.position,
        "strand": s.strand, "paired": s.length, "aligned": s.aligned,
        "mismatches": len(s.mismatches), "matched": s.matched, "risk": s.risk,
        "note": site_note(s),
    } for s in priming_sites(library, extra=extra, min_anneal=min_anneal,
                             mismatches=mismatches)]
    return pd.DataFrame(rows, columns=cols)


def mispriming_findings(library) -> tuple[list[str], list[str]]:
    """``(failures, advisories)`` over every constant flank the library is amplified with.

    A handle that pairs in full somewhere it should not is a failure: it has a second binding
    site, and the amplification the design relies on is not specific. Anything shorter is an
    advisory, spelled out by name when enough bases pair to anneal at a normal annealing
    temperature and collapsed into a count otherwise, so a report on a large pool stays
    readable. A duplex carrying a mismatch is only ever an advisory, whatever its length, so no
    verdict rests on how mismatches are weighted.

    A flank too short to be scored is named as well, so a report never reads as a clean check
    of something that was not checked.
    """
    failures: list[str] = []
    advisories: list[str] = []
    watch = 0
    for h in handles(library):
        if len(h.seq) < MIN_HANDLE:
            advisories.append(
                f"{h.label} is {len(h.seq)} base(s) long, too short to score as an "
                f"amplification handle (at least {MIN_HANDLE}), so it was not checked; pass "
                "the primer you amplify with to lib.mispriming(extra=[...])"
            )
    for site in priming_sites(library):
        # A mismatched duplex cannot be a full-length one, since pairing every base leaves
        # nothing to mismatch. Tested for anyway, so the verdict cannot come to rest on the
        # mismatch weighting by some later change.
        if site.risk == "collision" and not site.mismatches:
            failures.append(site_note(site))
        elif site.risk == "high":
            advisories.append(site_note(site))
        else:
            watch += 1
    if watch:
        advisories.append(
            f"{watch} shorter duplex{'' if watch == 1 else 'es'} pairing {MIN_ANNEAL}-"
            f"{HIGH_ANNEAL - 1} bases where they should not, too few to prime on their own at "
            "a normal annealing temperature; see lib.mispriming()"
        )
    return failures, advisories
