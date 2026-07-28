"""Tests for looking a protein up in UniProt (uniprot.py).

The download is replaced throughout, so the suite never touches the network. What is
exercised is everything around it: the header parse, the disk cache, the accession guard,
the error messages, and how a resolved entry lands on the spec and in the design record.
"""
from __future__ import annotations

import urllib.error

import pytest

from library_designer import LibrarySpec, SubstitutionScan
from library_designer import uniprot as up

ACYP1 = (
    ">sp|P07311|ACYP1_HUMAN Acylphosphatase-1 OS=Homo sapiens OX=9606 GN=ACYP1 PE=1 SV=2\n"
    "MAEGNTLISVDYEIFGKVQGVFFRKHTQAEGKKLGLVGWVQNTDRGTVQGQLQGPISKVR\n"
    "HMQEWLETRGSPKSHIDKANFNNEKVILKLDYSDFQIVK\n"
)
SEQUENCE = (
    "MAEGNTLISVDYEIFGKVQGVFFRKHTQAEGKKLGLVGWVQNTDRGTVQGQLQGPISKVR"
    "HMQEWLETRGSPKSHIDKANFNNEKVILKLDYSDFQIVK"
)


@pytest.fixture
def offline(monkeypatch, tmp_path):
    """Point the cache at a temp directory and count downloads, serving a canned FASTA."""
    monkeypatch.setenv("LIBRARY_DESIGNER_CACHE", str(tmp_path / "cache"))
    calls = []

    def _fake(url, timeout):
        calls.append(url)
        return ACYP1

    monkeypatch.setattr(up, "_download", _fake)
    return calls


# --- parsing ------------------------------------------------------------------

def test_parse_reads_every_header_field():
    e = up.parse_entry(ACYP1, "2026-07-25T00:00:00+00:00")
    assert (e.accession, e.entry_name) == ("P07311", "ACYP1_HUMAN")
    assert e.protein_name == "Acylphosphatase-1"
    assert (e.organism, e.gene, e.sequence_version) == ("Homo sapiens", "ACYP1", 2)
    assert e.reviewed and e.sequence == SEQUENCE
    assert "P07311" in str(e) and "99 aa" in str(e)
    assert "sequence" not in e.record()          # the record is metadata; the spec holds the residues


def test_parse_handles_an_unreviewed_entry_without_a_gene():
    text = ">tr|A0A0A0A0A0|A0A0A0A0A0_ECOLI Uncharacterized protein OS=Escherichia coli OX=562 PE=4 SV=1\nMKV\n"
    e = up.parse_entry(text, "now")
    assert not e.reviewed and e.gene is None
    assert e.protein_name == "Uncharacterized protein" and e.organism == "Escherichia coli"


def test_a_non_fasta_reply_is_rejected():
    with pytest.raises(ValueError, match="did not return a FASTA"):
        up.parse_entry("<html>rate limited</html>", "now")


def test_a_sequence_of_non_residues_is_rejected():
    with pytest.raises(ValueError, match="non-residue"):
        up.parse_entry(">sp|P1|X Protein OS=E OX=1 SV=1\nMKV123\n", "now")


# --- fetching and caching -----------------------------------------------------

def test_the_second_lookup_comes_from_the_cache(offline):
    first = up.fetch("P07311")
    second = up.fetch("P07311")
    assert first.sequence == second.sequence == SEQUENCE
    assert len(offline) == 1                     # downloaded once, cached after that
    assert up.cache_dir().joinpath("P07311.fasta").is_file()

    up.fetch("P07311", refresh=True)
    assert len(offline) == 2                     # refresh goes back to the network


def test_a_malformed_accession_never_reaches_the_network(offline):
    for bad in ("not-an-accession", "MAEGNTLISVDYE", "", "P0731"):
        with pytest.raises(ValueError, match="does not look like a UniProt accession"):
            up.fetch(bad)
    assert offline == []


def test_an_isoform_accession_is_accepted(offline):
    up.fetch("P07311-2")
    assert offline == ["https://rest.uniprot.org/uniprotkb/P07311-2.fasta"]


def test_an_unknown_entry_says_so(monkeypatch, tmp_path):
    monkeypatch.setenv("LIBRARY_DESIGNER_CACHE", str(tmp_path))

    def _404(url, timeout):
        raise urllib.error.HTTPError(url, 404, "Not Found", {}, None)

    monkeypatch.setattr(up, "_download", _404)
    with pytest.raises(ValueError, match="No UniProt entry 'P07311'"):
        up.fetch("P07311")


def test_being_offline_points_at_protein_sequence(monkeypatch, tmp_path):
    monkeypatch.setenv("LIBRARY_DESIGNER_CACHE", str(tmp_path))

    def _down(url, timeout):
        raise urllib.error.URLError("no route to host")

    monkeypatch.setattr(up, "_download", _down)
    with pytest.raises(ValueError, match="Set protein_sequence on the spec to work offline"):
        up.fetch("P07311")


# --- on the spec --------------------------------------------------------------

def test_from_uniprot_names_the_spec_after_the_entry(offline):
    spec = LibrarySpec.from_uniprot("P07311", substitutions=["A"], truncation=1)
    assert spec.name == "ACYP1_HUMAN"
    assert spec.protein_sequence == SEQUENCE
    assert spec.truncated_sequence == SEQUENCE[1:]        # truncation drops the initiator Met
    assert spec.uniprot_entry["organism"] == "Homo sapiens"
    assert "P07311" in repr(spec) and "Homo sapiens" in repr(spec)


def test_the_field_resolves_on_construction(offline):
    spec = LibrarySpec(name="mine", uniprot="P07311", substitutions=["A"])
    assert spec.protein_sequence == SEQUENCE
    assert spec.uniprot_entry["sequence_version"] == 2


def test_an_explicit_protein_sequence_wins_and_keeps_the_accession(offline):
    spec = LibrarySpec(name="mine", uniprot="P07311", protein_sequence="MKV")
    assert spec.protein_sequence == "MKV"
    assert spec.uniprot == "P07311" and spec.uniprot_entry is None
    assert offline == []                                   # nothing to fetch, so nothing was


def test_the_accession_can_be_resolved_after_construction(offline):
    spec = LibrarySpec(name="mine")
    with pytest.raises(ValueError, match="No accession to resolve"):
        spec.resolve_uniprot()
    spec.uniprot = "P07311"
    spec.resolve_uniprot()
    assert spec.protein_sequence == SEQUENCE


def test_a_toml_can_name_an_accession(offline, tmp_path):
    (tmp_path / "s.toml").write_text('name = "mine"\nuniprot = "P07311"\nsubstitutions = ["A"]\n')
    spec = LibrarySpec.from_toml(tmp_path / "s.toml")
    assert spec.protein_sequence == SEQUENCE


def test_the_design_record_says_where_the_protein_came_from(offline):
    spec = LibrarySpec.from_uniprot("P07311", name="scan", substitutions=["A"], truncation=1)
    lib = SubstitutionScan(spec).generate().codon_optimize()
    recorded = lib.design_specs["spec"]

    assert recorded["uniprot"] == "P07311"
    assert recorded["uniprot_entry"]["entry_name"] == "ACYP1_HUMAN"
    # The sequence is stored too, so the design does not depend on UniProt to be rebuilt.
    assert recorded["protein_sequence"] == SEQUENCE
    assert len(lib) == 98 * 1 - sum(1 for r in SEQUENCE[1:] if r == "A") + 1
