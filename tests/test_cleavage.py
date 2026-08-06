"""Tests for the lead-in check (checks/cleavage.py).

A Type IIS site flush against the end of a molecule cuts poorly, which costs yield rather than
correctness, so it is an advisory. The design is makeable and some protocols run it.
"""
from __future__ import annotations

import sys

import pytest

from library_designer import LibrarySpec, SubstitutionScan
from library_designer.checks.cleavage import (
    MIN_FLANK,
    cleavage_advisories,
    design_enzyme,
    site_flanks,
)
from library_designer.regions import reverse_complement as rc

sys.path.insert(0, "tests")
from test_vectors import (  # noqa: E402
    A3_BACKBONE,
    A5_BACKBONE,
    CLEAN_CDS,
    GCK_TOML,
    _plasmid,
    _protein,
    _standard_lib,
)

SITE = "GGTCTC"


def _lib(a5, a3, path=None, **kw):
    spec = LibrarySpec(name="c", protein_sequence=_protein(CLEAN_CDS), cds=CLEAN_CDS,
                       substitutions=["A"], adaptor_5=a5, adaptor_3=a3, **kw)
    return SubstitutionScan(spec).generate().codon_optimize()


def test_site_flanks_counts_from_the_outermost_site():
    lead, trail = site_flanks("GGCGC" + SITE + "ACCCC" + rc(SITE) + "TTTT", "BsaI")
    assert (lead, trail) == (5, 4)

    # Flush at both ends, which is what an adaptor written as site-first gives.
    assert site_flanks(SITE + "CAAGC" + "GGTG" + rc(SITE), "BsaI") == (0, 0)
    # A missing site is a different finding, reported elsewhere, so it is None here.
    assert site_flanks("ACGTACGT", "BsaI") == (None, None)


def test_a_terminal_site_is_reported_naming_the_end_and_the_adaptor():
    """The message has to say which end, how short it is, and where to add bases."""
    lib = _lib("ggtctccaagc", "ggtg" + "a" + rc(SITE).lower())
    notes = cleavage_advisories(lib)
    assert len(notes) == 2
    assert any("5'" in n and "adaptor_5" in n and "0 base(s)" in n for n in notes)
    assert any("3'" in n and "adaptor_3" in n for n in notes)
    assert all("BsaI" in n and str(MIN_FLANK) in n for n in notes)

    # Advisory only: it costs yield, not correctness, so the report still passes.
    rep = lib.check()
    assert rep.cleavage_advisories == notes
    assert set(notes) <= set(rep.advisories)
    assert rep.passed and "cleavage" not in rep.issues


def test_enough_lead_in_says_nothing(tmp_path):
    lib = _standard_lib(_plasmid(tmp_path, CLEAN_CDS), A5_BACKBONE, A3_BACKBONE)
    assert site_flanks(A5_BACKBONE, "BsaI")[0] >= MIN_FLANK
    assert cleavage_advisories(lib) == []
    assert lib.check().cleavage_advisories == []


def test_only_the_short_end_is_reported():
    """One end fixed and the other not should give one message, not two or none."""
    lib = _lib("gcgtc" + "ggtctcc" + "aagc", "ggtg" + "a" + rc(SITE).lower())
    (note,) = cleavage_advisories(lib)
    assert "3'" in note and "adaptor_3" in note


def test_a_tiled_oligo_gets_its_lead_in_from_the_primer():
    """Each tiled oligo starts with its amplification primer, so the site is never terminal
    and the check has nothing to say."""
    lib = SubstitutionScan(LibrarySpec.from_toml(GCK_TOML)).generate().codon_optimize()
    lib = lib.drop_failed().tile()
    assert cleavage_advisories(lib) == []
    oligo = next(o for o in lib.df["oligo"] if isinstance(o, str))
    lead, trail = site_flanks(oligo, lib.tiled_params.enzyme)
    assert lead >= MIN_FLANK and trail >= MIN_FLANK


def test_the_enzyme_is_taken_from_the_design_then_from_the_flanks(tmp_path):
    from library_designer import StartingVectorParams

    # Stated by the vector.
    lib = _standard_lib(_plasmid(tmp_path, CLEAN_CDS), A5_BACKBONE, A3_BACKBONE)
    assert design_enzyme(lib) == "BsaI"

    # No vector and no tiles: found by searching the adaptors, so the check still runs.
    bare = _lib("ggtctccaagc", "ggtg" + "a" + rc(SITE).lower())
    assert bare.spec.vector is None
    assert design_enzyme(bare) == "BsaI"
    assert cleavage_advisories(bare)

    # Nothing carrying a site at all.
    plain = _lib("", "")
    assert design_enzyme(plain) is None
    assert cleavage_advisories(plain) == []


def test_a_library_with_no_flanks_is_not_checked():
    lib = _lib("", "")
    assert cleavage_advisories(lib) == []
