"""QC over an optimized library: optimization success, translation round-trip,
restriction sites, forbidden motifs, and (if set) the synthesis length cap.

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

Two extra sets of checks switch on when the design says more about how it is built. A tiled
library adds the per-oligo and per-tile-vector checks in ``checks/tiled.py``. A standard
library that names a starting vector adds the adaptor-against-the-plasmid checks in
``checks/vector.py``, which is the only place a construct that looks clean on its own can
be caught not fitting the backbone it is meant to clone into.
"""
from __future__ import annotations

import re
from dataclasses import dataclass, field

import pandas as pd

from ..regions import assemble
from .motifs import count_enzyme_sites
from .translation import translates_to


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

    Which fields can fill in depends on the design. ``oligo_extra_sites``,
    ``oligo_over_budget``, and ``unplaced`` come from a library that has been ``tile()``'d.
    ``adaptor_issues`` comes from a standard library that names a starting vector.
    ``overhang_issues`` and ``vector_extra_sites`` come from either of those two. The
    assembly counts fill in whenever there is a destination vector to assemble into, tiled
    or not. Optimization, translation, enzyme sites, motifs, and the length cap are checked
    on every optimized library.

    ``reference_advisories`` is informational and stays out of ``passed``. It reports what a
    reference kept verbatim carries (a native BsaI site, an avoid-motif) rather than
    treating it as a failure, because the user opted into that sequence. It comes from the
    same two branches as ``overhang_issues``, so a library with neither tiles nor a starting
    vector has no advisories to report.

    ``overhang_advisories`` is informational for the same reason. It names overhang pairs
    that share more homology than the design should carry without being an outright
    collision, which is a hazard worth reading and not a verdict, since the overhangs come
    off the CDS rather than from an orthogonal set. ``lib.overhang_pairs()`` is the full
    table behind it.
    """

    n_variants: int
    optimization_failed: list[str] = field(default_factory=list)
    translation_fail: list[str] = field(default_factory=list)
    enzyme_hits: dict[str, list[str]] = field(default_factory=dict)
    motif_hits: dict[str, list[str]] = field(default_factory=dict)
    length_exceeded: list[str] = field(default_factory=list)
    # Tiled-assembly checks (empty unless the library was .tile()'d)
    oligo_extra_sites: list[str] = field(default_factory=list)   # unintended enzyme site (incl. junctions)
    oligo_over_budget: list[str] = field(default_factory=list)   # oligo longer than the budget
    overhang_issues: list[str] = field(default_factory=list)     # tiles with degenerate/palindromic overhangs
    unplaced: list[str] = field(default_factory=list)            # non-WT variants that landed in no tile
    vector_extra_sites: list[str] = field(default_factory=list)  # destination vector has enzyme sites beyond the 2 intended
    # Adaptors vs the destination vector (empty unless spec.starting_vector is set on an
    # untiled library): missing/ambiguous cut sites, or overhangs that will not ligate.
    adaptor_issues: list[str] = field(default_factory=list)
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
            and not any(self.enzyme_hits.values())
            and not any(self.motif_hits.values())
            and not self.length_exceeded
            and not self.oligo_extra_sites
            and not self.oligo_over_budget
            and not self.overhang_issues
            and not self.unplaced
            and not self.vector_extra_sites
            and not self.adaptor_issues
            and not self.assembly_issues
            and self.assembly_correct == self.assembly_checked
        )

    def text(self, count: bool = True) -> str:
        """The readable report. ``count`` prints the variant count in the header, which
        ``LibrarySummary`` turns off because its own first line already gives it."""
        head = "QC report: " + (f"{self.n_variants} variants, " if count else "")
        lines = [head + ("PASS" if self.passed else "FAIL")]
        checked = self.n_variants - len(self.optimization_failed)
        if self.optimization_failed:
            lines.append(f"  optimization: {len(self.optimization_failed)} failed{_names(self.optimization_failed)}")
        ok = checked - len(self.translation_fail)
        lines.append(f"  translation round-trip: {ok}/{checked} ok{_names(self.translation_fail)}")
        # Name the sites that were hit, and collapse the clean ones into one line rather
        # than a "no hit(s)" line each, which on a passing library is most of the report.
        clean: list[str] = []
        for enz, hits in self.enzyme_hits.items():
            if hits:
                lines.append(f"  {enz} site: {len(hits)} hit(s){_names(hits)}")
            else:
                clean.append(enz)
        clean_motifs = 0
        for pat, hits in self.motif_hits.items():
            if hits:
                lines.append(f"  motif /{pat}/: {len(hits)} hit(s){_names(hits)}")
            else:
                clean_motifs += 1
        if clean or clean_motifs:
            motifs = f"{clean_motifs} motif" + ("s" if clean_motifs > 1 else "")
            what = clean + ([motifs] if clean_motifs else [])
            lines.append("  no forbidden sequences (checked " + ", ".join(what) + ")")
        if self.length_exceeded:
            lines.append(
                f"  length: {len(self.length_exceeded)} over max_oligo_length"
                f"{_names(self.length_exceeded)}"
            )
        for label, hits in (
            ("oligo extra enzyme site(s)", self.oligo_extra_sites),
            ("oligo over budget", self.oligo_over_budget),
            ("tile overhang issue(s)", self.overhang_issues),
            ("unplaced variant(s)", self.unplaced),
            ("destination vector extra enzyme site(s)", self.vector_extra_sites),
            ("adaptor/destination-vector issue(s)", self.adaptor_issues),
            ("assembly issue(s)", self.assembly_issues),
        ):
            if hits:
                lines.append(f"  {label}: {len(hits)}{_names(hits)}")
        if self.assembly_checked:
            lines.append(
                f"  assembly simulation: {self.assembly_correct}/{self.assembly_checked} "
                "members rebuild their intended variant"
            )
            if self.assembly_aligned:
                lines.append(
                    f"  aligned to the parent vector: {self.assembly_aligned}/{self.assembly_checked} "
                    "differ only at the intended codon"
                )
        advisories = self.reference_advisories + self.overhang_advisories
        if advisories:
            lines.append(f"  advisories (not failures): {len(advisories)}")
            for msg in advisories[:5]:
                lines.append("    - " + msg)
            if len(advisories) > 5:
                lines.append(f"    ... and {len(advisories) - 5} more")
        return "\n".join(lines)

    def __str__(self) -> str:
        return self.text()


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
        len(df), opt_failed, translation_fail, enzyme_hits, motif_hits, length_exceeded
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

    # Then put the design together: digest, ligate, and align the product against the
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
