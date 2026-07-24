"""Tests for the SequenceSet generator and independent per-member optimization.

A SequenceSet is a library of distinct full-length proteins (orthologs, designs,
multi-mutants). Unlike a scan there is no shared reference: each member is codon-
optimized on its own. These tests cover generation, that optimization is genuinely
independent, QC, and the OMEGA FASTA hand-off.
"""
from __future__ import annotations

import pandas as pd
import pytest
from dnachisel import translate

from library_designer import LibrarySpec, SequenceSet
from library_designer.integrations import omega

# Three short, distinct proteins standing in for orthologs / designs.
PROTEINS = {
    "ortholog_a": "MKAILVGADEQRSTNWYF",
    "design_01":  "MSTNRQPLVGKAEDFHIW",
    "design_02":  "MADEQRSTNKLVGAWYFH",
}


def _spec(**kw) -> LibrarySpec:
    return LibrarySpec(name="designs", **kw)


def test_generate_shape_and_kind():
    lib = SequenceSet(_spec(), PROTEINS).generate()
    assert len(lib) == 3
    assert lib.kind == "sequence_set"
    # scan-only columns exist (for downstream compatibility) but are all NA; no WT row.
    assert set(lib.df["name"]) == set(PROTEINS)
    assert lib.df["mut_residue"].isna().all()
    assert "WT" not in set(lib.df["name"])
    assert lib.design_specs["kind"] == "sequence_set"


def test_optimization_is_independent():
    lib = SequenceSet(_spec(), PROTEINS).generate().codon_optimize()
    assert not lib.failed
    assert lib.reference is None                     # no shared reference for a sequence set
    by_name = dict(zip(lib.df["name"], lib.df["variable_dna"]))
    for name, protein in PROTEINS.items():
        dna = by_name[name]
        assert isinstance(dna, str) and len(dna) == 3 * len(protein)
        assert translate(dna) == protein             # each member encodes its own protein
    # distinct proteins give distinct CDSs (not all collapsed onto one reference)
    assert len({str(v) for v in by_name.values()}) == 3


def test_check_passes():
    lib = SequenceSet(_spec(), PROTEINS).generate().codon_optimize()
    report = lib.check()
    assert report.passed
    assert not report.translation_fail
    assert not any(report.enzyme_hits.values())      # DNA Chisel avoids the BsaI site


def test_whole_gene_skips_oligo_length_gate():
    # A full-length member is longer than any single oligo; with an explicit cap set,
    # a sequence set must NOT be flagged (OMEGA fragments each gene into oligos).
    lib = SequenceSet(_spec(max_oligo_length=30), PROTEINS).generate().codon_optimize()
    assert lib.check().length_exceeded == []


def test_to_omega_fasta(tmp_path):
    lib = SequenceSet(_spec(), PROTEINS).generate().codon_optimize()
    path = tmp_path / "designs.faa"
    n = omega.write_fasta(lib, path)
    lines = path.read_text().splitlines()
    headers = [ln[1:] for ln in lines if ln.startswith(">")]
    seqs = [ln for ln in lines if not ln.startswith(">")]
    assert n == 3 and set(headers) == set(PROTEINS)
    assert all(set(s) <= set("ACGT") and len(s) % 3 == 0 for s in seqs)


def test_summary_reports_members_bucket():
    lib = SequenceSet(_spec(), PROTEINS).generate().codon_optimize()
    s = lib.summary()
    assert s.per_sublibrary == {"members": 3}
    assert s.n_optimized == 3 and s.n_failed == 0


def test_from_fasta_roundtrip(tmp_path):
    faa = tmp_path / "in.faa"
    faa.write_text(
        "".join(f">{n} some description\n{seq}\n" for n, seq in PROTEINS.items())
    )
    lib = SequenceSet.from_fasta(_spec(), faa).generate()
    assert dict(zip(lib.df["name"], lib.df["protein"])) == PROTEINS


def test_rejects_bad_name():
    with pytest.raises(ValueError, match="FASTA header|uSort-M"):
        SequenceSet(_spec(), {"has space": "MKAIL"}).generate()


def test_rejects_non_amino_acid():
    with pytest.raises(ValueError, match="non-amino-acid"):
        SequenceSet(_spec(), {"stop_in_middle": "MKA*LV"}).generate()


def test_empty_set_errors():
    with pytest.raises(ValueError, match="at least one"):
        SequenceSet(_spec(), {}).generate()
