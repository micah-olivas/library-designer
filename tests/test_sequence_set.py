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


# --- reproducibility ----------------------------------------------------------

def test_members_are_reproducible_and_ambient_rng_independent():
    """Each member gets its own seed offset, so the whole set rebuilds identically
    regardless of where the caller's RNG sits."""
    import numpy as np

    np.random.seed(1)
    first = list(SequenceSet(_spec(), PROTEINS).generate().codon_optimize().df["variable_dna"])
    np.random.seed(888)
    second = list(SequenceSet(_spec(), PROTEINS).generate().codon_optimize().df["variable_dna"])
    assert first == second

    # Distinct offsets, not one shared stream: identical proteins would still differ only
    # by their seed, so check the offsets reach the members rather than colliding.
    assert len(set(first)) == len(first)


def test_gc_table_does_not_file_every_design_under_wt():
    """``mut_residue`` is NA on every row of a set, and NA is how a scan spells its wild-type
    control, so reading the column straight would report each design as the WT."""
    lib = SequenceSet(_spec(), PROTEINS).generate().codon_optimize()
    t = lib.gc_table()
    assert set(t["sublibrary"]) == {"members"}        # one bucket, matching per_sublibrary
    assert set(t["name"]) == set(PROTEINS)            # each design still named individually
    assert "WT" not in set(t["sublibrary"])


def test_the_gc_figure_names_a_sets_bucket_members():
    import matplotlib
    matplotlib.use("Agg")

    lib = SequenceSet(_spec(gc_bounds=(0.30, 0.70)), PROTEINS).generate().codon_optimize()
    labels = [t.get_text() for t in lib.plot_gc_distribution().axes[1].get_legend().get_texts()]
    assert "members" in labels
    assert "WT control" not in labels                 # a set has no wild-type control
    assert not any(lab.startswith("to ") for lab in labels)


def test_unseeded_sequence_set_runs():
    """seed=None must reach each member as None instead of being offset into an error."""
    import numpy as np

    np.random.seed(3)
    lib = SequenceSet(_spec(seed=None), PROTEINS).generate().codon_optimize()
    assert not lib.failed
    assert lib.design_specs["seed"] is None
    for protein, dna in zip(lib.df["protein"], lib.df["variable_dna"]):
        assert translate(dna) == protein
