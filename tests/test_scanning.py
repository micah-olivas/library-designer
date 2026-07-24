"""End-to-end golden path for the flagship non-tiled substitution scan.

Loads ``examples/mbo038.toml`` and runs the full pipeline (protein -> generate ->
codon_optimize -> QC/export), which exercises the DNA Chisel codon-optimization path
that the native-CDS tiled/OMEGA fixtures skip. Also regression-tests the junction
screen: QC judges the *assembled* construct (adaptor + variable + adaptor) against an
assembled-WT baseline, so a restriction site spelled across an adaptor<->variable
junction is caught even when the variable region alone is clean.
"""
from __future__ import annotations

import re
from pathlib import Path

import pandas as pd
import pytest

from library_designer import LibrarySpec, SubstitutionScan
from library_designer.checks.motifs import count_enzyme_sites
from library_designer.checks.translation import translates_to
from library_designer.regions import assemble

REPO = Path(__file__).resolve().parents[1]
MBO038_TOML = REPO / "examples" / "mbo038.toml"


@pytest.fixture(scope="module")
def mbo038():
    spec = LibrarySpec.from_toml(MBO038_TOML)
    return SubstitutionScan(spec).generate().codon_optimize()


# --- codon-optimization round trip -------------------------------------------

def test_every_variant_round_trips_to_intended_protein(mbo038):
    # The flagship path codon-optimizes a protein (no native CDS) and must not fail.
    assert mbo038.failed == {}
    for protein, dna in zip(mbo038.df["protein"], mbo038.df["variable_dna"]):
        assert isinstance(dna, str)
        assert translates_to(dna, protein)


def test_single_wt_reference_invariant(mbo038):
    # On the *codon-optimized* reference (not a native-CDS fixture): every member is
    # byte-identical to the frozen reference except within its own stamped codon.
    ref = mbo038.reference
    saw_mutant = False
    for dna, idx in zip(mbo038.df["variable_dna"], mbo038.df["mut_index"]):
        assert len(dna) == len(ref)
        diffs = [k for k in range(len(ref)) if dna[k] != ref[k]]
        if pd.isna(idx):
            assert diffs == []                        # WT control is the reference
            continue
        codon = range(int(idx) * 3, int(idx) * 3 + 3)
        assert all(k in codon for k in diffs)         # differs only within its codon
        saw_mutant = True
    assert saw_mutant


# --- QC, including the assembled-construct junction screen (fix #1) -----------

def test_qc_passes_and_no_variant_exceeds_assembled_baseline(mbo038):
    spec = mbo038.spec
    a5, a3 = spec.adaptor_5.upper(), spec.adaptor_3.upper()
    baseline = {e: count_enzyme_sites(assemble(a5, mbo038.reference, a3), e)
                for e in spec.avoid_enzymes}
    for dna in mbo038.df["variable_dna"]:
        construct = assemble(a5, dna, a3)
        for e in spec.avoid_enzymes:
            # No variant manufactures a site (junction included) beyond the intended ones.
            assert count_enzyme_sites(construct, e) <= baseline[e]

    rep = mbo038.check()
    assert rep.passed
    assert not rep.optimization_failed and not rep.translation_fail
    assert not any(rep.enzyme_hits.values())
    assert not any(rep.motif_hits.values())
    assert not rep.length_exceeded


def test_qc_flags_restriction_site_at_adaptor_junction():
    """A stamped edge codon can spell a BsaI site with adjacent adaptor bases. The
    variable region alone is clean, so the old variable-only screen missed it; QC on
    the assembled construct must catch it (this is the fix #1 regression)."""
    # Pin GGT (Gly) at every position. adaptor_3 begins "CTC...", so the last codon of
    # the variable region + adaptor_3 spells GGT|CTC = GGTCTC (BsaI) only across the
    # junction, never inside the variable region.
    spec = LibrarySpec(
        name="junction",
        protein_sequence="MAAK",
        substitutions=["GGT"],       # literal codon -> placed verbatim, residue symbol "G"
        adaptor_3="CTCAAA",
        avoid_enzymes=["BsaI"],
        avoid_patterns=[],           # isolate the enzyme check
        seed=1,
    )
    lib = SubstitutionScan(spec).generate().codon_optimize()

    row = lib.df[lib.df["name"] == "K4G"].iloc[0]
    assert isinstance(row["variable_dna"], str)
    # Clean when screened in isolation ...
    assert count_enzyme_sites(row["variable_dna"], "BsaI") == 0
    # ... but the assembled construct carries a junction BsaI site.
    construct = assemble("", row["variable_dna"], spec.adaptor_3.upper())
    assert count_enzyme_sites(construct, "BsaI") == 1

    rep = lib.check()
    assert "K4G" in rep.enzyme_hits["BsaI"]
    assert not rep.passed


# --- exporters ----------------------------------------------------------------

def test_usortm_export_contract(mbo038, tmp_path):
    p = tmp_path / "variants.csv"
    mbo038.to_usortm(p)
    out = pd.read_csv(p)
    assert list(out.columns) == ["name", "sequence"]
    assert not out["name"].astype(str).str.contains(r"[/|>\s]", regex=True).any()

    a5, a3 = mbo038.spec.adaptor_5, mbo038.spec.adaptor_3
    for seq in out["sequence"]:
        assert seq.startswith(a5.lower()) and seq.endswith(a3.lower())
        variable = seq[len(a5):len(seq) - len(a3)]
        assert variable.isupper() and set(variable) <= set("ACGT")


def test_vendor_pooled_schema(mbo038, tmp_path):
    p = tmp_path / "order.csv"
    mbo038.to_vendor(p)   # platform = twist_oligo_pools -> pooled
    out = pd.read_csv(p)
    assert list(out.columns) == ["Pool Name", "Insert Length", "Insert Sequence"]
    assert (out["Pool Name"].astype(str) == mbo038.spec.name).all()
    assert (out["Insert Length"] == out["Insert Sequence"].str.len()).all()


# --- determinism --------------------------------------------------------------

def test_codon_optimization_is_deterministic_given_seed():
    # The design-specs/reproducibility claim: same seed -> byte-identical output.
    def run():
        spec = LibrarySpec(name="det", protein_sequence="MKAILVDE",
                           substitutions=["A"], seed=3)
        return SubstitutionScan(spec).generate().codon_optimize()

    a, b = run(), run()
    assert a.reference == b.reference
    assert list(a.df["variable_dna"]) == list(b.df["variable_dna"])


# --- design-specs hygiene -----------------------------------------------------

def test_drop_failed_also_clears_the_design_specs_record():
    """The design-specs JSON is the handoff record, so it must not keep listing
    failures for variants the library no longer holds."""
    spec = LibrarySpec(name="dropped", protein_sequence="MKAILVDE", substitutions=["A"])
    lib = SubstitutionScan(spec).generate().codon_optimize()
    lib.failed = {"K2A": "synthetic failure"}
    lib.design_specs["failed"] = dict(lib.failed)

    lib.drop_failed()
    assert lib.failed == {}
    assert "failed" not in lib.design_specs


# --- spec rendering -----------------------------------------------------------

def test_spec_html_shows_readable_placeholders():
    """An unset field renders as words, not as the stray punctuation left behind when
    the placeholder glyph was stripped."""
    spec = LibrarySpec(name="x", protein_sequence="MKV", avoid_enzymes=[])
    html = spec._repr_html_()
    for label in ("platform", "max_oligo_length"):
        cell = re.search(rf">{label}</th><td[^>]*>(.*?)</td>", html).group(1)
        assert cell == "not set"
    assert "no enzymes" in html
