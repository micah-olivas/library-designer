"""Regressions from a bug-hunting pass over the package.

Each test here failed before its fix. They are grouped by what goes wrong when the fix is
missing, because that is what decides how much a defect matters: a wrong sequence that looks
plausible gets ordered, a crash only costs an afternoon.
"""
from __future__ import annotations

import json
import sys

import pytest
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import CompoundLocation, SeqFeature, SimpleLocation
from Bio.SeqRecord import SeqRecord

from library_designer import LibrarySpec, SequenceSet, SubstitutionScan
from library_designer.regions import reverse_complement as rc

sys.path.insert(0, "tests")
from test_vectors import (  # noqa: E402
    A3_BACKBONE,
    A5_BACKBONE,
    BB3,
    BB5,
    CLEAN_CDS,
    GCK_TOML,
    _plasmid,
    _protein,
    _standard_lib,
    _write_gb,
)

PROTEIN = "MKAILVVLLYTFATANADTLLILGDSLSAG"


# --- wrong sequences that looked plausible ------------------------------------

def test_a_tiled_order_form_carries_the_oligos_not_the_whole_cds(tmp_path):
    """to_vendor used to put the full 1395 bp coding sequence on a pooled order form for a
    tiled library, which no oligo-pool vendor can make and nothing in the file betrays."""
    import pandas as pd

    lib = SubstitutionScan(LibrarySpec.from_toml(GCK_TOML)).generate().codon_optimize()
    lib.tile().drop_failed()
    lib.to_vendor(tmp_path / "order.csv")

    out = pd.read_csv(tmp_path / "order.csv")
    assert set(out["Insert Sequence"]) == set(lib.df["oligo"].dropna())
    assert out["Insert Length"].max() <= lib.tiled_params.oligo_budget
    assert len(out) == int(lib.df["oligo"].notna().sum())     # the unplaced WT is left off


def test_an_overhang_split_between_adaptor_and_cds_assembles_correctly(tmp_path):
    """An adaptor spelling only part of the fused overhang is a valid design: the rest of the
    overhang comes from the coding sequence. The splice used to reconstruct nearly the whole
    reference twice and report every member broken."""
    a5 = "GGCGC" + "GGTCTC" + "A" + BB5[-2:]        # 2 of the 4 overhang bases
    a3 = BB3[:2] + "T" + rc("GGTCTC") + "GCGCC"
    lib = _standard_lib(_plasmid(tmp_path, CLEAN_CDS), a5, a3)
    dv = lib.destination_vector()
    assert (dv.cut.keep_5, dv.cut.keep_3) == (2, 2)
    assert (dv.start, dv.end) == (2, len(lib.reference) - 2)

    rep = lib.check()
    # Only the real finding survives: the codons whose bases sit inside an overhang.
    assert all("fused overhang" in m for m in rep.assembly_issues)
    assert rep.assembly_correct == rep.assembly_checked - 1


def test_a_clone_map_places_downstream_features_by_the_reference_length(tmp_path):
    """to_assembled_vectors shifted 3' backbone features by the length of the region the
    plasmid held, not the reference that replaces it, so every one was misplaced whenever the
    two differed. It also disagreed with to_vector_maps, which had it right."""
    stuffer, term = "GGGCCC" * 5, "AACCCGCTGATCGGCACGTAAGAGGTTCCA"
    seq = BB5 + stuffer + BB3 + term
    path = _write_gb(tmp_path / "stuffer.gb", seq, features=[
        ("CDS", len(BB5), len(BB5) + len(stuffer), 1, "insert"),
        ("terminator", len(BB5) + len(stuffer) + len(BB3), len(seq), 1, "term"),
    ])
    lib = _standard_lib(path, A5_BACKBONE, A3_BACKBONE, cds=None, protein=_protein(CLEAN_CDS))
    assert len(lib.reference) != len(stuffer)              # the lengths must differ to see it

    lib.to_assembled_vectors(tmp_path / "clones", fmt="genbank")
    lib.to_vector_maps(tmp_path / "maps")
    for gb in (tmp_path / "clones" / "parent_WT.gb", tmp_path / "maps" / "destination.gb"):
        rec = SeqIO.read(str(gb), "genbank")
        feat = next(f for f in rec.features if (f.qualifiers.get("label") or [""])[0] == "term")
        assert str(rec.seq)[int(feat.location.start):int(feat.location.end)] == term


def test_an_insert_feature_crossing_the_file_origin_is_refused(tmp_path):
    """Such a feature is stored as a join(), and BioPython reports its span as the whole
    molecule. Locating used to believe that and emit a 'vector' with no backbone left."""
    seq = "ACGT" * 61 + "A"
    rec = SeqRecord(Seq(seq), id="w", name="w", description="wrap",
                    annotations={"molecule_type": "DNA", "topology": "circular"})
    rec.features.append(SeqFeature(
        CompoundLocation([SimpleLocation(154, len(seq), strand=1), SimpleLocation(0, 90, strand=1)]),
        type="CDS", qualifiers={"label": ["insert"]},
    ))
    path = str(tmp_path / "wrap.gb")
    SeqIO.write(rec, path, "genbank")

    from library_designer.layout.vector_io import resolve_destination
    with pytest.raises(ValueError, match="crosses the origin of the file"):
        resolve_destination(path, insert_label="insert")


def test_a_located_insert_covering_the_whole_plasmid_is_refused(tmp_path):
    from library_designer.layout.vector_io import DestinationContext, flanks

    whole = DestinationContext(full_seq="ACGT" * 10, topology="circular", start=0, end=40,
                               located_region="ACGT" * 10)
    with pytest.raises(ValueError, match="entire plasmid"):
        flanks(whole, 0, 0)                                # keep-0 adaptors asked for nothing


def test_the_map_origin_is_reported_on_the_file_own_strand(tmp_path):
    """With the insert on the minus strand the whole plasmid is reverse-complemented, so the
    emitted map runs the other way and the origin has to be given in the file's numbering."""
    import pandas as pd

    from library_designer.layout.vector_io import read_vector_file

    path = _plasmid(tmp_path, CLEAN_CDS, reverse=True)
    lib = _standard_lib(path, A5_BACKBONE, A3_BACKBONE)
    lib.to_vectors(tmp_path / "v.csv")
    row = pd.read_csv(tmp_path / "v.csv").iloc[0]

    assert row["insert_strand_in_starting_vector"] == -1
    file_seq, _, _ = read_vector_file(path)
    base = int(row["origin_in_starting_vector"])
    assert row["vector_sequence"][:20] == rc(file_seq)[len(file_seq) - base:][:20]


# --- silently designing a different protein -----------------------------------

@pytest.mark.parametrize("residue", ["X", "B", "Z", "*"])
def test_a_non_canonical_residue_is_refused(residue):
    """DNA Chisel reverse-translates X/B/Z to an arbitrary amino acid without complaining, so
    an unchecked one designs a library for a protein nobody asked for."""
    with pytest.raises(ValueError, match="non-amino-acid"):
        LibrarySpec(name="x", protein_sequence=f"MKAILV{residue}LYTFATANAD")


def test_a_pasted_sequence_is_cleaned_of_whitespace_and_case():
    """Pasting from a FASTA brings newlines; a lowercase protein used to produce variants
    named 'a3A' and then crash with KeyError: 'm'."""
    spec = LibrarySpec(name="lc", protein_sequence="mkail vvlly\ntfat", cds="gct aaa\nctg")
    assert spec.protein_sequence == "MKAILVVLLYTFAT"
    assert spec.cds == "GCTAAACTG"

    lib = SubstitutionScan(LibrarySpec(name="lc", protein_sequence=PROTEIN.lower(),
                                       substitutions=["A"])).generate()
    assert not any(n[0].islower() for n in lib.df["name"] if n != "WT")
    assert "A3A" not in set(lib.df["name"])          # a no-op mutant is not a variant


def test_an_unknown_avoid_enzyme_says_which_ones_are_known():
    spec = LibrarySpec(name="e", protein_sequence=PROTEIN, substitutions=["A"])
    spec.avoid_enzymes = ["EcoRI"]
    with pytest.raises(KeyError, match="Unknown enzyme 'EcoRI'"):
        SubstitutionScan(spec).generate().codon_optimize()


def test_truncating_the_whole_protein_says_so():
    spec = LibrarySpec(name="t", protein_sequence=PROTEIN, substitutions=["A"],
                       truncation=len(PROTEIN))
    with pytest.raises(ValueError, match="leaving nothing to design"):
        SubstitutionScan(spec).generate()


def test_duplicate_member_names_are_refused():
    with pytest.raises(ValueError, match="Duplicate member name"):
        SequenceSet(LibrarySpec(name="s"), [("a", PROTEIN), ("a", PROTEIN[::-1])])


# --- records that misstated what was run --------------------------------------

def test_re_optimizing_clears_the_previous_run(tmp_path):
    """A second codon_optimize() rebuilds every sequence, so a tile layout and a failure list
    from the first run describe sequences that no longer exist. The oligos used to survive and
    would have been exported against the new variants."""
    lib = SubstitutionScan(LibrarySpec.from_toml(GCK_TOML)).generate().codon_optimize().tile()
    assert lib.tiles and "oligo" in lib.df.columns
    lib.design_specs["failed"] = {"Y273A": "from the last run", "K10A": "no longer true"}

    lib.codon_optimize()
    assert lib.tiles is None and "oligo" not in lib.df.columns
    assert "tiled" not in lib.design_specs
    # The record lists this run's failures, not the previous run's.
    assert lib.design_specs.get("failed", {}) == lib.failed
    assert "K10A" not in lib.design_specs.get("failed", {})


def test_the_record_names_the_parameters_actually_used(tmp_path):
    """design_specs['spec'] was snapshotted at generate() time, so editing the spec before
    optimizing left the record describing a run that never happened."""
    spec = LibrarySpec(name="a", protein_sequence=PROTEIN, substitutions=["A"], seed=0)
    lib = SubstitutionScan(spec).generate()
    spec.optimization.species = "h_sapiens"
    spec.seed = 7
    lib.codon_optimize()

    assert lib.design_specs["spec"]["optimization"]["species"] == "h_sapiens"
    assert lib.design_specs["spec"]["seed"] == lib.design_specs["seed"] == 7


def test_a_spec_round_trips_through_its_own_design_record(tmp_path):
    spec = LibrarySpec(name="r", protein_sequence=PROTEIN, substitutions=["A"],
                       starting_vector=_plasmid(tmp_path, CLEAN_CDS))
    lib = SubstitutionScan(spec).generate().codon_optimize()
    lib.to_design_specs(tmp_path / "record.json")

    rebuilt = LibrarySpec(**json.loads((tmp_path / "record.json").read_text())["spec"])
    assert rebuilt.optimization.species == spec.optimization.species
    assert rebuilt.vector.path == spec.vector.path
    assert repr(rebuilt)                                    # used to raise on the dict
    again = SubstitutionScan(rebuilt).generate().codon_optimize()
    assert again.reference == lib.reference


# --- exports that failed halfway ----------------------------------------------

def test_a_tiled_library_exports_its_own_handoff(tmp_path):
    """export_all used to write four files and then raise on to_usortm, leaving a half-filled
    run directory. A tiled pool hands off as oligos plus primers instead."""
    lib = SubstitutionScan(LibrarySpec.from_toml(GCK_TOML)).generate().codon_optimize()
    lib.tile().drop_failed()
    lib.export_all(tmp_path / "out", plots=False)

    written = {p.name for p in lib.output_dir.iterdir()}
    assert {"GCK_oligo_pool.csv", "GCK_primers.csv", "GCK_order.csv"} <= written
    assert "variants.csv" not in written                    # no uniform variable region to write


def test_an_export_that_cannot_run_leaves_no_directory(tmp_path):
    spec = LibrarySpec(name="noplat", protein_sequence=PROTEIN, substitutions=["A"])
    lib = SubstitutionScan(spec).generate().codon_optimize()
    with pytest.raises(ValueError, match="Specify method"):
        lib.export_all(tmp_path / "out", plots=False)
    assert not (tmp_path / "out").exists()

    lib2 = SubstitutionScan(spec).generate()
    with pytest.raises(ValueError, match="not codon-optimized"):
        lib2.export_all(tmp_path / "out2", plots=False)
    assert not (tmp_path / "out2").exists()


def test_a_sequence_set_exports_without_warning_about_a_plot_it_cannot_draw(tmp_path):
    spec = LibrarySpec(name="ss", platform="pooled")
    lib = SequenceSet(spec, {"a": PROTEIN, "b": PROTEIN[::-1]}).generate().codon_optimize()
    with pytest.warns(None) if False else _no_warnings():
        lib.export_all(tmp_path / "out")
    assert (lib.output_dir / "ss_full_library.csv").is_file()
    assert not list(lib.output_dir.glob("*.png"))

    with pytest.raises(ValueError, match="every member is its own gene"):
        lib.plot_codon_usage()


class _no_warnings:
    """Assert the block emits no warning, which pytest.warns(None) no longer does."""

    def __enter__(self):
        import warnings

        self._ctx = warnings.catch_warnings(record=True)
        self.caught = self._ctx.__enter__()
        import warnings as w
        w.simplefilter("always")
        return self

    def __exit__(self, *exc):
        result = self._ctx.__exit__(*exc)
        assert not self.caught, [str(x.message) for x in self.caught]
        return result
