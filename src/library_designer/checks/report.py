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
"""
from __future__ import annotations

import re
from dataclasses import dataclass, field

import pandas as pd

from ..regions import assemble
from .motifs import count_enzyme_sites
from .translation import translates_to


@dataclass
class CheckReport:
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
    # Informational only, does NOT fail the report: things kept verbatim from a chosen
    # reference (SD sites, motifs, a native BsaI) that the user opted into, not recoded.
    reference_advisories: list[str] = field(default_factory=list)

    @property
    def passed(self) -> bool:
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
        )

    def __str__(self) -> str:
        lines = [
            f"QC report: {self.n_variants} variants, "
            + ("PASS" if self.passed else "FAIL")
        ]
        checked = self.n_variants - len(self.optimization_failed)
        if self.optimization_failed:
            lines.append(f"  optimization: {len(self.optimization_failed)} failed")
            lines.append("    x " + ", ".join(self.optimization_failed[:5]))
        ok = checked - len(self.translation_fail)
        lines.append(f"  translation round-trip: {ok}/{checked} ok")
        if self.translation_fail:
            lines.append("    x " + ", ".join(self.translation_fail[:5]))
        for enz, hits in self.enzyme_hits.items():
            lines.append(f"  {enz} site: {len(hits) or 'no'} hit(s)")
            if hits:
                lines.append("    x " + ", ".join(hits[:5]))
        for pat, hits in self.motif_hits.items():
            lines.append(f"  motif /{pat}/: {len(hits) or 'no'} hit(s)")
            if hits:
                lines.append("    x " + ", ".join(hits[:5]))
        if self.length_exceeded:
            lines.append(f"  length: {len(self.length_exceeded)} over max_oligo_length")
            lines.append("    x " + ", ".join(self.length_exceeded[:5]))
        for label, hits in (
            ("oligo extra enzyme site(s)", self.oligo_extra_sites),
            ("oligo over budget", self.oligo_over_budget),
            ("tile overhang issue(s)", self.overhang_issues),
            ("unplaced variant(s)", self.unplaced),
            ("destination vector extra enzyme site(s)", self.vector_extra_sites),
        ):
            if hits:
                lines.append(f"  {label}: {len(hits)}")
                lines.append("    x " + ", ".join(hits[:5]))
        if self.reference_advisories:
            lines.append(f"  advisories (informational, not a failure): {len(self.reference_advisories)}")
            for msg in self.reference_advisories[:5]:
                lines.append("    - " + msg)
        return "\n".join(lines)


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
    return report
