"""QC over an optimized library: optimization success, translation round-trip,
restriction sites, forbidden motifs, and, when the spec sets them, the synthesis length cap,
the GC window, and the homopolymer limit.

The last two judge the molecule that is ordered, the assembled oligo for a tiled library and
the construct with its adaptors otherwise, which is what a vendor's spec is written against.
``ordered_molecules`` builds those sequences and the GC and homopolymer screens read them.

Restriction-site and motif checks run on the **fully assembled construct**
(``adaptor_5 + variable_dna + adaptor_3``) against an *assembled-WT baseline* (the
same counts on ``adaptor_5 + reference + adaptor_3``). Adaptors legitimately carry
restriction sites, e.g. the Golden Gate BsaI site the substitution-scan adaptors
place at the variable-region boundary, so only occurrences *beyond* the baseline
count against a variant. Screening the assembled construct (not the variable region
in isolation) catches a site that a variant's edge codon spells together with
adjacent adaptor bases across the junction, which the old variable-only check
missed. Translation round-trip is still checked on the variable region alone
(adaptors are not coding).

Two extra sets of checks switch on when the spec says more about how the library is built. A tiled
library adds the per-oligo and per-tile-vector checks in ``checks/tiled.py``. A standard
library that names a starting vector adds the adaptor-against-the-plasmid checks in
``checks/vector.py``, which is the only place a construct that looks clean on its own can
be caught not fitting the backbone it is meant to clone into.

Two more run on any library that carries constant flanks. ``checks/cleavage.py`` measures how
far each Type IIS site sits from the end of the molecule: a site flush against the end cuts
poorly, which costs yield rather than correctness, so it is an advisory. ``checks/mispriming.py``
asks whether the flanks the pool is amplified with, a tile's primer pair or the shared adaptors,
anneal anywhere they should not: inside a variable region, inside another oligo's flank, or in
the destination-vector backbone. That is a PCR failure rather than a cloning one, so it is
checked on the flanks and the sequence around them and not on the assembled plasmid.
"""
from __future__ import annotations

import re
from dataclasses import asdict, dataclass, field
from html import escape as _esc

import pandas as pd

from ..regions import assemble
from .motifs import count_enzyme_sites
from .translation import translates_to


def label_width(labels, cap: int = 20) -> int:
    """How far to pad ``label:`` so a block of values lines up.

    Measured over the labels actually present, so a short report stays tight instead of being
    spaced out by a label it never printed. A label past ``cap`` (a motif pattern, say) is left
    long rather than pushing every other value across the screen.
    """
    return max((len(label) for label in labels if len(label) <= cap), default=0) + 2


def rows_to_lines(rows: list[tuple[str, str]], indent: str = "  ", cap: int = 20) -> list[str]:
    """``[(label, value), ...]`` as aligned ``label: value`` lines."""
    width = label_width([label for label, _ in rows], cap)
    return [f"{indent}{(label + ':').ljust(width)}{value}".rstrip() for label, value in rows]


def _names(items: list[str], limit: int = 4) -> str:
    """The offending members in parentheses, e.g. ``" (I7F, S8Y, and 3 more)"``, or an
    empty string when the list is empty. Named inline so a failure costs one line, not two."""
    if not items:
        return ""
    extra = len(items) - limit
    shown = ", ".join(items[:limit])
    return f" ({shown}, and {extra} more)" if extra > 0 else f" ({shown})"


@dataclass
class CheckReport:
    """What QC found, one field per check. Printing it gives the readable report.

    Most fields are lists of the variant names that failed a check, so an empty list is a
    pass. ``assembly_checked``, ``assembly_correct``, and ``assembly_aligned`` are the
    exception, being counts of members put through the simulation.

    ``off_target_edits`` is the single-WT-reference invariant, checked on the sequences rather
    than on an assembled plasmid. Every member has to match the frozen
    reference outside the codon(s) it is meant to mutate. It runs on any library with a shared
    reference, with or without a destination vector, which is what makes it different from
    ``assembly_aligned``, the same idea read off the simulated clone and only available when
    there is a plasmid to assemble into.

    Which fields can fill in depends on the library. ``oligo_extra_sites``,
    ``oligo_over_budget``, and ``unplaced`` come from a library that has been ``tile()``'d.
    ``adaptor_issues`` comes from a standard library that names a starting vector.
    ``overhang_issues`` and ``vector_extra_sites`` come from either of those two.
    ``mispriming_issues`` comes from any library that carries constant flanks to amplify
    with, tile primers or adaptors. The
    assembly counts fill in whenever there is a destination vector to assemble into, tiled
    or not. Optimization, translation, enzyme sites, motifs, and the length cap are checked
    on every optimized library.

    ``gc_out_of_range`` and ``homopolymer_hits`` fill in only when the spec sets ``gc_bounds``
    or ``max_homopolymer``, and both are read off the molecule that is ordered, the assembled
    oligo for a tiled library and the construct with its adaptors otherwise. A homopolymer is
    counted absolutely rather than against the wild-type baseline the enzyme and motif screens
    use, since a long run already in the reference is a finding too.

    ``reference_advisories`` is informational and stays out of ``passed``. It reports what a
    reference kept verbatim carries (a native BsaI site, an avoid-motif) rather than
    treating it as a failure, because the user opted into that sequence. It comes from the
    same two branches as ``overhang_issues``, so a library with neither tiles nor a starting
    vector has no advisories to report.

    ``overhang_advisories`` is informational for the same reason. It names overhang pairs
    that share more homology than a library should carry without being an outright
    collision, which is a hazard worth reading and not a verdict, since the overhangs come
    off the CDS rather than from an orthogonal set. ``lib.overhang_pairs()`` is the full
    table behind it.

    ``mispriming_advisories`` is informational too. It names a constant flank whose 3' end
    anneals somewhere it should not without the whole flank occurring there, which primes or
    does not depending on the annealing temperature you run at. ``lib.mispriming()`` is the
    full table.

    ``cleavage_advisories`` is the last of the informational ones. It names an end whose Type
    IIS recognition site sits so close to the end of the molecule that the enzyme has little
    duplex to hold, which costs yield rather than correctness, and says where to add lead-in
    bases. It fills in for any library carrying a flank with a site in it. See
    ``checks/cleavage.py``.
    """

    n_variants: int
    optimization_failed: list[str] = field(default_factory=list)
    translation_fail: list[str] = field(default_factory=list)
    # Members whose DNA differs from the frozen reference somewhere other than the codon(s)
    # they are supposed to mutate. See off_target_edits in check_library.
    off_target_edits: list[str] = field(default_factory=list)
    enzyme_hits: dict[str, list[str]] = field(default_factory=dict)
    motif_hits: dict[str, list[str]] = field(default_factory=dict)
    length_exceeded: list[str] = field(default_factory=list)
    # Members whose ordered molecule sits outside spec.gc_bounds. Empty when the gate is off.
    gc_out_of_range: list[str] = field(default_factory=list)
    # Members whose ordered molecule carries a single-base run longer than
    # spec.max_homopolymer. Empty when the gate is off.
    homopolymer_hits: list[str] = field(default_factory=list)
    # Tiled-assembly checks (empty unless the library was .tile()'d)
    oligo_extra_sites: list[str] = field(default_factory=list)   # unintended enzyme site (incl. junctions)
    oligo_over_budget: list[str] = field(default_factory=list)   # oligo longer than the budget
    overhang_issues: list[str] = field(default_factory=list)     # tiles with degenerate/palindromic overhangs
    unplaced: list[str] = field(default_factory=list)            # non-WT variants that landed in no tile
    vector_extra_sites: list[str] = field(default_factory=list)  # destination vector has enzyme sites beyond the 2 intended
    # Adaptors vs the destination vector (empty unless spec.starting_vector is set on an
    # untiled library): missing/ambiguous cut sites, or overhangs that will not ligate.
    adaptor_issues: list[str] = field(default_factory=list)
    # A constant flank the pool is amplified with (a tile primer, an adaptor) that occurs in
    # full somewhere else it would prime. See checks/mispriming.py.
    mispriming_issues: list[str] = field(default_factory=list)
    # Assembly simulation (runs whenever there is a destination vector to assemble into):
    # digest, ligate, and confirm the product is the plasmid carrying the intended variant.
    assembly_issues: list[str] = field(default_factory=list)
    assembly_checked: int = 0            # members put through the simulation
    assembly_correct: int = 0            # members whose product carries their intended CDS
    assembly_aligned: int = 0            # members whose product differs from the parent only
                                         # at the intended codon
    # Informational only, does NOT fail the report: things kept verbatim from a chosen
    # reference (SD sites, motifs, a native BsaI) that the user opted into, not recoded.
    reference_advisories: list[str] = field(default_factory=list)
    # Informational only: overhang pairs that share more homology than they should without
    # being an outright collision. See checks/overhangs.py and lib.overhang_pairs().
    overhang_advisories: list[str] = field(default_factory=list)
    # Informational: a recognition site sitting too close to the end of the molecule to cut
    # efficiently. See checks/cleavage.py.
    cleavage_advisories: list[str] = field(default_factory=list)
    # Informational only: a constant flank whose 3' end anneals somewhere it should not
    # without the whole flank occurring there. See checks/mispriming.py and lib.mispriming().
    mispriming_advisories: list[str] = field(default_factory=list)

    @property
    def passed(self) -> bool:
        """True when every failure list is empty and every simulated assembly rebuilt its
        intended variant, which is ``assembly_correct == assembly_checked``.

        With no assembly to simulate both counts are 0 and that clause holds.
        ``reference_advisories`` does not affect the result, and neither does
        ``assembly_aligned``, which is a stricter read on the same products and is reported
        for information.
        """
        return (
            not self.optimization_failed
            and not self.translation_fail
            and not self.off_target_edits
            and not any(self.enzyme_hits.values())
            and not any(self.motif_hits.values())
            and not self.length_exceeded
            and not self.gc_out_of_range
            and not self.homopolymer_hits
            and not self.oligo_extra_sites
            and not self.oligo_over_budget
            and not self.overhang_issues
            and not self.unplaced
            and not self.vector_extra_sites
            and not self.adaptor_issues
            and not self.mispriming_issues
            and not self.assembly_issues
            and self.assembly_correct == self.assembly_checked
        )

    def text(self, count: bool = True) -> str:
        """The readable report, one check per line with the values in a column.

        ``count`` prints the variant count in the header, which ``LibrarySummary`` turns off
        because its own first line already gives it. Checks that found nothing are left out,
        apart from the two that are worth seeing pass (translation and forbidden sequences),
        so a clean report stays short enough to read at a glance.
        """
        head = "QC report: " + (f"{self.n_variants} variants, " if count else "")
        checked = self.n_variants - len(self.optimization_failed)
        rows: list[tuple[str, str]] = []

        if self.optimization_failed:
            rows.append(("optimization",
                         f"{len(self.optimization_failed)} failed{_names(self.optimization_failed)}"))
        ok = checked - len(self.translation_fail)
        rows.append(("translation", f"{ok}/{checked}{_names(self.translation_fail)}"))
        if self.off_target_edits:
            rows.append(("unintended edits",
                         f"{len(self.off_target_edits)} outside their own codon"
                         f"{_names(self.off_target_edits)}"))

        # Name the sites that were hit, and collapse the clean ones into one row rather than a
        # "none" row each, which on a passing library is most of the report.
        clean: list[str] = []
        for enz, hits in self.enzyme_hits.items():
            if hits:
                rows.append((f"{enz} sites", f"{len(hits)} hit{'' if len(hits) == 1 else 's'}{_names(hits)}"))
            else:
                clean.append(enz)
        clean_motifs = 0
        for pat, hits in self.motif_hits.items():
            if hits:
                rows.append((f"motif /{pat}/", f"{len(hits)} hit{'' if len(hits) == 1 else 's'}{_names(hits)}"))
            else:
                clean_motifs += 1
        if clean or clean_motifs:
            motifs = f"{clean_motifs} motif" + ("s" if clean_motifs > 1 else "")
            rows.append(("forbidden sequences",
                         "none (checked " + ", ".join(clean + ([motifs] if clean_motifs else [])) + ")"))

        if self.length_exceeded:
            rows.append(("length",
                         f"{len(self.length_exceeded)} over max_oligo_length{_names(self.length_exceeded)}"))
        if self.gc_out_of_range:
            rows.append(("GC window",
                         f"{len(self.gc_out_of_range)} outside gc_bounds"
                         f"{_names(self.gc_out_of_range)}"))
        if self.homopolymer_hits:
            rows.append(("homopolymers",
                         f"{len(self.homopolymer_hits)} with a run over max_homopolymer"
                         f"{_names(self.homopolymer_hits)}"))
        # Checks that report member names keep them on the row. Checks that report whole
        # sentences get the count on the row and the sentences underneath, since a paragraph
        # wedged into parentheses runs off the screen and hides the column.
        for label, hits in (("oligo enzyme sites", self.oligo_extra_sites),
                            ("oligo length", self.oligo_over_budget),
                            ("unplaced", self.unplaced)):
            if hits:
                rows.append((label, f"{len(hits)}{_names(hits)}"))
        spelled_out: list[tuple[str, list[str]]] = [
            (label, hits) for label, hits in (("tile overhangs", self.overhang_issues),
                                              ("vector enzyme sites", self.vector_extra_sites),
                                              ("adaptors vs vector", self.adaptor_issues),
                                              ("mispriming", self.mispriming_issues),
                                              ("assembly issues", self.assembly_issues))
            if hits
        ]
        rows += [(label, str(len(hits))) for label, hits in spelled_out]

        if self.assembly_checked:
            rows.append(("assembly",
                         f"{self.assembly_correct}/{self.assembly_checked} rebuild their variant"))
            if self.assembly_aligned:
                rows.append(("parent alignment",
                             f"{self.assembly_aligned}/{self.assembly_checked} differ only at "
                             "the intended codon"))

        advisories = self.advisories
        if advisories:
            rows.append(("advisories", f"{len(advisories)}, not failures"))

        detail = {label: hits for label, hits in spelled_out}
        detail["advisories"] = advisories
        width = label_width([label for label, _ in rows])
        lines = [head + ("PASS" if self.passed else "FAIL")]
        for label, value in rows:
            lines.append(f"  {(label + ':').ljust(width)}{value}".rstrip())
            for msg in detail.get(label, [])[:5]:
                lines.append("    - " + msg)
            if len(detail.get(label, [])) > 5:
                lines.append(f"    ... and {len(detail[label]) - 5} more")
        return "\n".join(lines)

    def __str__(self) -> str:
        return self.text()

    def __repr__(self) -> str:   # the readable report, not a dump of every empty field
        return self.text()

    def _repr_html_(self) -> str:
        return f"<pre style='margin:0;line-height:1.4'>{_esc(str(self))}</pre>"

    # --- structured views, for reading the result in code ------------------
    #
    # Split around the two per-pattern checks, which need a key each, so ``issues`` comes
    # out in the order ``text()`` prints and ``passed`` tests.
    _CHECKS_BEFORE_PATTERNS = ("optimization_failed", "translation_fail",
                               "off_target_edits")
    _CHECKS_AFTER_PATTERNS = (
        "length_exceeded", "gc_out_of_range", "homopolymer_hits", "oligo_extra_sites",
        "oligo_over_budget",
        "overhang_issues",
        "unplaced", "vector_extra_sites", "adaptor_issues", "mispriming_issues",
        "assembly_issues",
    )

    @property
    def issues(self) -> dict[str, list[str]]:
        """Only the checks that found something, as ``{check: entries}``. Empty exactly
        when ``passed`` is True, so ``if rep.issues:`` is the branch to write in a script.

        Keys are field names, except the two per-pattern checks, which are split out one
        key per enzyme or motif (``"enzyme_hits:BsaI"``,
        ``"motif_hits:GGAGG.{8,12}[AG]TG"``), and ``"assembly_incorrect"``, which reports
        the count mismatch that has no list of its own.

        Entries are variant names for the checks that judge members one at a time
        (optimization, translation, enzyme and motif hits, length, the GC and homopolymer
        gates, the oligo checks, ``unplaced``) and sentences for the ones that judge the
        design as a whole
        (``overhang_issues``, ``vector_extra_sites``, ``adaptor_issues``,
        ``mispriming_issues``, ``assembly_issues``). Advisories are not here, they are not
        failures; see ``advisories``.
        """
        out: dict[str, list[str]] = {}
        for name in self._CHECKS_BEFORE_PATTERNS:
            if hits := getattr(self, name):
                out[name] = list(hits)
        for field_name, per_pattern in (("enzyme_hits", self.enzyme_hits),
                                       ("motif_hits", self.motif_hits)):
            for pattern, hits in per_pattern.items():
                if hits:
                    out[f"{field_name}:{pattern}"] = list(hits)
        for name in self._CHECKS_AFTER_PATTERNS:
            if hits := getattr(self, name):
                out[name] = list(hits)
        if self.assembly_correct != self.assembly_checked:
            missing = self.assembly_checked - self.assembly_correct
            out["assembly_incorrect"] = [
                f"{missing} of {self.assembly_checked} members do not rebuild their "
                "intended variant"
            ]
        return out

    @property
    def advisories(self) -> list[str]:
        """Every informational list in one, in reporting order: reference, overhang,
        mispriming, cleavage. They never affect ``passed`` or ``issues``."""
        return (list(self.reference_advisories) + list(self.overhang_advisories)
                + list(self.mispriming_advisories) + list(self.cleavage_advisories))

    def to_dict(self) -> dict:
        """The whole report as plain data, every field plus ``passed``, ``issues``, and
        ``advisories``. JSON-serializable, for recording a run or diffing two of them."""
        out = asdict(self)
        out["passed"] = self.passed
        out["issues"] = self.issues
        out["advisories"] = self.advisories
        return out


def verbatim_advisories(spec, reference: str) -> list[str]:
    """What a reference the user chose to keep verbatim (``spec.cds`` or the plasmid's own
    CDS) carries that we would otherwise have recoded away: restricted enzyme sites and
    avoid-motifs. Informational, so both the tiled and the single-vector checks report
    these rather than failing on them. The user opted into this sequence."""
    out: list[str] = []
    for e in spec.avoid_enzymes:
        n = count_enzyme_sites(reference, e)
        if n:
            out.append(
                f"reference carries {n} internal {e} site(s), kept verbatim "
                "(a Golden Gate hazard; see the destination-vector check)"
            )
    for p in spec.avoid_patterns:
        n = len(re.findall(p, reference))
        if n:
            out.append(f"reference carries {n} /{p}/ motif(s), kept verbatim (not recoded)")
    return out


def off_target_edits(library) -> list[str]:
    """Members whose DNA strays from the frozen reference outside the codon(s) they mutate.

    The single-WT-reference invariant, and the reason the backbone-and-stamp model exists: a
    single-mutant library is only interpretable if every member matches the WT
    everywhere except its own position, otherwise a phenotype cannot be pinned on the
    substitution. Optimizing each variant separately breaks this quietly, since a
    usage-matching method lets unchanged positions drift apart between members.

    Judged per member against ``library.reference``, so it needs a library with one. A
    ``SequenceSet`` has no shared WT (every member is a different gene) and gives ``[]``.
    Intended positions come from ``mut_index``; a row with none, the wild-type control and the
    per-tile ``WT_Tile_<i>`` controls, has to match the reference exactly. Written over a set
    of intended indices rather than a single one, so a multi-substitution generator would be
    covered without changing this.
    """
    reference = getattr(library, "reference", None)
    if not reference or "variable_dna" not in library.df.columns:
        return []
    df = library.df
    idx = df["mut_index"] if "mut_index" in df.columns else [pd.NA] * len(df)
    out: list[str] = []
    for name, dna, i in zip(df["name"], df["variable_dna"], idx):
        if not isinstance(dna, str):
            continue                       # optimization failed; reported on its own
        if len(dna) != len(reference):
            out.append(str(name))
            continue
        intended = set() if pd.isna(i) else set(range(int(i) * 3, int(i) * 3 + 3))
        if any(a != b and k not in intended for k, (a, b) in enumerate(zip(dna, reference))):
            out.append(str(name))
    return out


def ordered_molecules(library) -> dict[str, str]:
    """What each member is physically ordered as, by name.

    A tiled library orders the assembled oligo, primers and enzyme sites included; anything
    else orders the whole construct with its adaptors. That is the molecule a vendor's spec is
    written against, so it is the one the GC and length gates judge. Members with nothing to
    order (a failed optimization, or the global ``WT`` row of a tiled pool) are left out.
    """
    df = library.df
    if "variable_dna" not in df.columns:
        return {}
    if getattr(library, "tiles", None) is not None:
        return {str(n): o for n, o in zip(df["name"], df.get("oligo", []))
                if isinstance(o, str)}
    return {
        str(n): assemble(a5, dna, a3)
        for n, a5, dna, a3 in zip(df["name"], df["adaptor_5"], df["variable_dna"],
                                  df["adaptor_3"])
        if isinstance(dna, str)
    }


def gc_fraction(seq: str) -> float:
    """G+C as a fraction of ``seq``, 0.0 for an empty sequence."""
    up = seq.upper()
    return (up.count("G") + up.count("C")) / len(up) if up else 0.0


def gc_table(library) -> pd.DataFrame:
    """Per-member GC, one row per molecule that is ordered.

    ``ordered_gc`` is the number the ``gc_bounds`` gate judges, the whole molecule with its
    flanks. ``variable_gc`` is the coding region alone, which the gate ignores and which sits
    lower whenever the flanks are GC-rich. ``in_bounds`` is NA when no bounds are set.

    ``sublibrary`` is the residue each member mutates to, ``"WT"`` for the control, so a
    distribution can be split by scan. A sequence set is not a scan and has no mutated residue,
    so its members go in one ``"members"`` bucket, matching ``LibrarySummary.per_sublibrary``.
    Each member is still named individually in ``name``, which for a set built from a FASTA is
    the header it came in under.
    """
    ordered = ordered_molecules(library)
    df = library.df
    set_kind = getattr(library, "kind", "scan") == "sequence_set"
    sub = dict(zip(df["name"].astype(str),
                   df["mut_residue"] if "mut_residue" in df.columns else [pd.NA] * len(df)))
    variable = dict(zip(df["name"].astype(str), df.get("variable_dna", [])))
    bounds = library.spec.gc_bounds
    rows = []
    for name, molecule in ordered.items():
        gc = gc_fraction(molecule)
        coding = variable.get(name)
        rows.append({
            "name": name,
            # A sequence set carries the column but leaves it NA on every row, so reading NA as
            # the wild-type control would file every design under "WT".
            "sublibrary": ("members" if set_kind
                           else "WT" if pd.isna(sub.get(name)) else str(sub[name])),
            "ordered_gc": gc,
            "variable_gc": gc_fraction(coding) if isinstance(coding, str) else pd.NA,
            "in_bounds": (bounds[0] <= gc <= bounds[1]) if bounds else pd.NA,
        })
    return pd.DataFrame(rows, columns=["name", "sublibrary", "ordered_gc", "variable_gc",
                                       "in_bounds"])


def longest_run(seq: str) -> int:
    """The longest single-base run in ``seq``, 0 for an empty one."""
    return max((len(m.group()) for m in re.finditer(r"(.)\1*", seq.upper())), default=0)


def homopolymer_hits(library) -> list[str]:
    """Members whose ordered molecule carries a run longer than ``spec.max_homopolymer``.

    Judged on the whole molecule, so a run spelled across an adaptor-to-CDS junction counts.
    The optimizer cannot see that junction: it constrains the coding region, and the flanks are
    added afterwards, which is why this is checked again here rather than assumed.

    Counted absolutely rather than against the wild-type baseline, unlike the enzyme and motif
    screens. An adaptor's Type IIS site is there on purpose and a native CDS's motif was
    accepted by whoever supplied it, but nobody intends a long homopolymer, so one already
    present is a finding too.
    """
    limit = library.spec.max_homopolymer
    if not limit:
        return []
    return [name for name, seq in ordered_molecules(library).items()
            if longest_run(seq) > limit]


def gc_out_of_range(library) -> list[str]:
    """Members whose ordered molecule falls outside ``spec.gc_bounds``.

    Judged on the whole molecule rather than the coding region, since that is what the
    synthesiser receives and what a vendor's GC window refers to. Empty when no bounds are
    set, which is the default: the package ships no vendor registry, so the window is stated
    by the caller (Twist recommend 0.35 to 0.65 for oligo pools).
    """
    bounds = library.spec.gc_bounds
    if not bounds:
        return []
    lo, hi = bounds
    return [name for name, seq in ordered_molecules(library).items()
            if not lo <= gc_fraction(seq) <= hi]


def check_library(library) -> CheckReport:
    spec, df = library.spec, library.df
    if "variable_dna" not in df.columns:
        raise ValueError("Library is not codon-optimized yet, call codon_optimize() first.")

    opt_failed: list[str] = []
    translation_fail: list[str] = []
    length_exceeded: list[str] = []
    enzyme_hits = {e: [] for e in spec.avoid_enzymes}
    motif_hits = {p: [] for p in spec.avoid_patterns}
    compiled = {p: re.compile(p) for p in spec.avoid_patterns}
    a5, a3 = spec.adaptor_5.upper(), spec.adaptor_3.upper()
    flank = len(a5) + len(a3)

    # Enzyme sites / forbidden motifs are judged on the *assembled construct*
    # (adaptor_5 + variable + adaptor_3) relative to an *assembled-WT baseline*
    # (adaptor_5 + reference + adaptor_3). A match already in that baseline is the
    # user's intended sequence, the adaptors' own restriction sites or a native CDS's
    # motifs, not something a variant introduced, so only a net-new occurrence counts,
    # including one that spans an adaptor<->variable junction. For a codon-optimized
    # reference with adaptor-free flanks every baseline count is 0, so this reduces to
    # the plain absolute check.
    reference = getattr(library, "reference", None) or ""
    baseline = assemble(a5, reference, a3)
    ref_enzyme = {e: count_enzyme_sites(baseline, e) for e in spec.avoid_enzymes}
    ref_motif = {p: len(rx.findall(baseline)) for p, rx in compiled.items()}

    # Length gate for a one-oligo construct: the explicit max_oligo_length. It is
    # skipped for libraries whose members are whole genes fragmented downstream, a
    # tiled library (checked per-oligo against the budget in checks.tiled) or a
    # sequence set headed for OMEGA (OMEGA splits each gene into oligos itself), where
    # a full-length member is expected to exceed any single-oligo cap. No synthesis-limit
    # registry ships with the package, so set spec.max_oligo_length yourself to enable it.
    tiled = getattr(library, "tiles", None) is not None
    whole_gene = tiled or getattr(library, "kind", "scan") == "sequence_set"
    length_cap = spec.max_oligo_length

    for name, protein, dna in zip(df["name"], df["protein"], df["variable_dna"]):
        name = str(name)
        if pd.isna(dna):
            opt_failed.append(name)
            continue
        if not translates_to(dna, protein):
            translation_fail.append(name)
        construct = assemble(a5, dna, a3)
        for enz in spec.avoid_enzymes:
            if count_enzyme_sites(construct, enz) > ref_enzyme[enz]:
                enzyme_hits[enz].append(name)
        for pat, rx in compiled.items():
            if len(rx.findall(construct)) > ref_motif[pat]:
                motif_hits[pat].append(name)
        if not whole_gene and length_cap is not None and flank + len(dna) > length_cap:
            length_exceeded.append(name)

    report = CheckReport(
        len(df), opt_failed, translation_fail, off_target_edits(library), enzyme_hits,
        motif_hits, length_exceeded, gc_out_of_range(library), homopolymer_hits(library),
    )
    if getattr(library, "tiles", None) is not None:
        from .tiled import check_tiled

        t = check_tiled(library)
        report.oligo_extra_sites = t["oligo_extra_sites"]
        report.oligo_over_budget = t["oligo_over_budget"]
        report.overhang_issues = t["overhang_issues"]
        report.unplaced = t["unplaced"]
        report.vector_extra_sites = t["vector_extra_sites"]
        report.reference_advisories = t["reference_advisories"]
        report.overhang_advisories = t["overhang_advisories"]
    elif spec.vector is not None and getattr(library, "reference", None):
        # A standard library cloned into a real plasmid: check the adaptors against it.
        from .vector import check_vector

        v = check_vector(library)
        report.adaptor_issues = v["adaptor_issues"]
        report.overhang_issues = v["overhang_issues"]
        report.vector_extra_sites = v["vector_extra_sites"]
        report.reference_advisories = v["advisories"]
        report.overhang_advisories = v["overhang_advisories"]

    # The constant flanks the pool is amplified with, against the variable regions, the other
    # oligos' flanks, and the destination-vector backbone. This runs whether or not the library
    # is tiled, since the hazard is in the PCR off the pool rather than in the cloning, and it
    # comes back empty for a library that carries no flanks at all.
    from .mispriming import mispriming_findings

    report.mispriming_issues, report.mispriming_advisories = mispriming_findings(library)
    # Runs on anything carrying flanks, tiled or not, since the site's distance from the end is
    # a property of the flank rather than of how the library is cloned.
    from .cleavage import cleavage_advisories

    report.cleavage_advisories = cleavage_advisories(library)

    # Then put the construct together: digest, ligate, and align the product against the
    # parent. Needs something to assemble into, so a library with no destination vector
    # (tiled or a starting vector) is skipped.
    if getattr(library, "tiles", None) is not None or spec.vector is not None:
        from .assembly import check_assembly

        a = check_assembly(library)
        report.assembly_issues = a["assembly_issues"]
        report.assembly_checked = a["assembly_checked"]
        report.assembly_correct = a["assembly_correct"]
        report.assembly_aligned = a["assembly_aligned"]
    return report
