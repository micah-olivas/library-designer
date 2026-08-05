"""Tests for the assembly simulation (checks/assembly.py).

Every other check reads the sequences we wrote down. These exercise the one that
puts them together: cut the oligo, cut the vector, ligate, and align the product against the
parent plasmid. The tampering tests are the point of it, they introduce errors that no
design-level check can see and confirm the simulation still catches them.
"""
from __future__ import annotations

import json

import pytest

from library_designer import LibrarySpec, StartingVectorParams, SubstitutionScan, TiledAssemblyParams
from library_designer.checks.assembly import align_to_parent, digest, junctions, ligate, simulate
from library_designer.regions import reverse_complement
from test_vectors import (
    A3_BACKBONE,
    A3_CDS,
    A5_BACKBONE,
    A5_CDS,
    BB3,
    BB5,
    CLEAN_CDS,
    GCK_TOML,
    _is_rotation,
    _plasmid,
    _protein,
    _standard_lib,
)

SITE = "GGTCTC"
SITE_RC = reverse_complement(SITE)


# --- the digest primitive -----------------------------------------------------

def test_digest_cuts_one_base_past_the_site_leaving_four(tmp_path):
    """BsaI is GGTCTC(1/5): one base of spacer, then the four-base overhang it leaves."""
    seq = "AAAA" + SITE + "T" + "CGCG" + "TTTTTT"
    assert junctions(seq, "BsaI") == [(11, 15)]
    left, right = digest(seq, "BsaI")
    assert (left.left, left.core, left.right) == ("", "AAAA" + SITE + "T", "CGCG")
    assert (right.left, right.core, right.right) == ("CGCG", "TTTTTT", "")
    # The overhang bases belong to both partners, so ligating restores the molecule once.
    assert left.core + right.seq == seq


def test_digest_reads_a_reverse_strand_site_backwards(tmp_path):
    seq = "AAAAAA" + "CGCG" + "T" + SITE_RC + "TTTT"
    assert junctions(seq, "BsaI") == [(6, 10)]
    left, right = digest(seq, "BsaI")
    assert left.right == right.left == "CGCG"


def test_digest_of_a_plasmid_gives_a_backbone_and_a_dropout(tmp_path):
    lib = _standard_lib(_plasmid(tmp_path, CLEAN_CDS), A5_BACKBONE, A3_BACKBONE)
    dv = lib.destination_vector()
    pieces = digest(dv.sequence, "BsaI", circular=True)
    assert len(pieces) == 2
    dropout = [p for p in pieces if SITE in p.seq or SITE_RC in p.seq]
    assert len(dropout) == 1                                  # the sites leave with the stuffer
    backbone = next(p for p in pieces if p is not dropout[0])
    assert (backbone.right, backbone.left) == (dv.overhang_5, dv.overhang_3)


def test_an_uncut_molecule_comes_back_whole():
    assert digest("ACGT" * 10, "BsaI") == [type(digest("A", "BsaI")[0])("", "ACGT" * 10, "")]


# --- the standard reaction ----------------------------------------------------

def test_the_wt_reaction_rebuilds_the_parent_plasmid(tmp_path):
    lib = _standard_lib(_plasmid(tmp_path, CLEAN_CDS), A5_BACKBONE, A3_BACKBONE)
    (r,) = lib.simulate_assembly()

    assert r.ok and r.topology == "circular"
    # The clone you would sequence is the starting plasmid, nothing added and nothing lost.
    assert _is_rotation(r.product, BB5 + CLEAN_CDS + BB3)
    assert r.n_correct == r.n_members == len(lib)
    assert r.n_aligned == r.n_members


def test_every_variant_differs_from_the_parent_at_its_own_codon_only(tmp_path):
    lib = _standard_lib(_plasmid(tmp_path, CLEAN_CDS), A5_BACKBONE, A3_BACKBONE)
    rep = lib.check()
    assert rep.passed
    assert rep.assembly_checked == len(lib)
    assert rep.assembly_correct == rep.assembly_aligned == len(lib)
    assert "differ only at the intended codon" in str(rep)


def test_a_mutation_inside_a_fused_overhang_cannot_ligate(tmp_path):
    """With the overhang drawn from the CDS ends, a variant that changes one of those bases
    presents an overhang the cut vector does not, so it cannot clone. Only a simulation
    catches this; every design-level check passes."""
    lib = _standard_lib(_plasmid(tmp_path, CLEAN_CDS), A5_CDS, A3_CDS)
    rep = lib.check()

    assert not rep.adaptor_issues and not rep.translation_fail   # the construct itself is fine
    assert not rep.passed
    blocked = [m for m in rep.assembly_issues if "fused overhang" in m]
    assert len(blocked) == 2                                    # one finding per end
    assert "K2A" in blocked[0]                                  # first codon pair
    assert "N59A" in blocked[1] and "P60A" in blocked[1]         # last codon pair
    assert rep.assembly_correct == rep.assembly_checked - 3


def test_the_alignment_catches_a_silent_off_target_change(tmp_path):
    """A synonymous change somewhere else in the CDS. Translation still round-trips and the
    assembled sequence still matches what was designed, so only the alignment against the
    parent vector can see it."""
    lib = _standard_lib(_plasmid(tmp_path, CLEAN_CDS), A5_BACKBONE, A3_BACKBONE)
    row = lib.df.index[lib.df["name"] == "D14A"][0]
    dna = lib.df.at[row, "variable_dna"]
    assert dna[24:27] == "AAT"                                  # Asn at codon 9
    tampered = dna[:24] + "AAC" + dna[27:]                      # a synonymous AAT -> AAC
    lib.df.at[row, "variable_dna"] = tampered

    rep = lib.check()
    assert not rep.translation_fail                             # protein unchanged
    assert not rep.passed
    assert any("D14A" in m and "outside the intended codon" in m for m in rep.assembly_issues)
    assert rep.assembly_aligned == len(lib) - 1


def test_the_alignment_requires_the_wt_control_to_come_back_untouched(tmp_path):
    lib = _standard_lib(_plasmid(tmp_path, CLEAN_CDS), A5_BACKBONE, A3_BACKBONE)
    row = lib.df.index[lib.df["name"] == "WT"][0]
    dna = lib.df.at[row, "variable_dna"]
    lib.df.at[row, "variable_dna"] = dna[:24] + "AAC" + dna[27:]   # synonymous, so only the
                                                                   # alignment can see it

    rep = lib.check()
    assert any("WT product differs from the parent vector" in m for m in rep.assembly_issues)


def test_align_to_parent_reports_a_length_change():
    parent = "ACGT" * 10
    assert align_to_parent(parent, parent, 0, 40, None) is None
    assert "not a clean substitution" in align_to_parent(parent + "A", parent, 0, 40, None)


def test_align_to_parent_accepts_the_intended_codon_and_nothing_else():
    parent = "AAA" + "GGGCCC" + "TTT"
    intended = "AAA" + "GGGAAA" + "TTT"          # codon 2 of a CDS starting at base 3
    assert align_to_parent(intended, parent, 3, 6, 1) is None
    assert "was not changed" in align_to_parent(parent, parent, 3, 6, 1)
    other = "AAA" + "AAACCC" + "TTT"             # codon 1 changed instead
    assert "outside the intended codon 2" in align_to_parent(other, parent, 3, 6, 1)


# --- ligation refuses what will not go ----------------------------------------

def test_ligate_says_which_end_does_not_anneal(tmp_path):
    lib = _standard_lib(_plasmid(tmp_path, CLEAN_CDS), A5_BACKBONE, A3_BACKBONE)
    dv = lib.destination_vector()
    pieces = digest(dv.sequence, "BsaI", circular=True)
    insert = digest((A5_BACKBONE + CLEAN_CDS + A3_BACKBONE).upper(), "BsaI")[1]

    product, issues = ligate(pieces, insert)
    assert product is not None and not issues

    from dataclasses import replace
    wrong = replace(insert, left="TTTT")
    product, issues = ligate(pieces, wrong)
    assert product is None
    assert any("5' overhang (TTTT)" in m for m in issues)


# --- tiled reactions ----------------------------------------------------------

@pytest.fixture(scope="module")
def tiled_lib(tmp_path_factory):
    """A tiled library over a real backbone, so every tile has a plasmid to assemble into."""
    tmp = tmp_path_factory.mktemp("tiled")
    path = _plasmid(tmp, CLEAN_CDS)
    spec = LibrarySpec(name="syn", protein_sequence=_protein(CLEAN_CDS), cds=CLEAN_CDS,
                       substitutions=["A"], starting_vector=StartingVectorParams(path=path),
                       tiled=TiledAssemblyParams(oligo_budget=150))
    return SubstitutionScan(spec).generate().codon_optimize().tile()


def test_every_tile_assembles_into_its_own_plasmid(tiled_lib):
    results = simulate(tiled_lib)
    assert len(results) == len(tiled_lib.tiles)
    for r in results:
        assert r.ok, r
        assert _is_rotation(r.product, BB5 + CLEAN_CDS + BB3)   # each tile rebuilds the parent
        assert r.n_aligned == r.n_members > 0


def test_every_tile_gets_a_wt_control_that_rebuilds_the_parent(tiled_lib):
    """The per-tile controls are what make each sublibrary's WT clone orderable, so every
    one of them has to come back as the unmutated plasmid."""
    parent = tiled_lib.parent_vector()
    controls = [f"WT_Tile_{t.index}" for t in tiled_lib.tiles]
    assert set(controls) <= set(tiled_lib.df["name"])
    for name in controls:
        assert tiled_lib.assembled_product(name) == parent


def test_tiled_qc_reports_the_simulation(tiled_lib):
    rep = tiled_lib.check()
    assert rep.passed
    # Every placed member, across every tile, is accounted for (the global WT row is not,
    # it rides on no oligo; the per-tile WT controls are).
    assert rep.assembly_checked == int(tiled_lib.df["oligo"].notna().sum())
    assert rep.assembly_correct == rep.assembly_aligned == rep.assembly_checked


def test_a_tampered_tile_is_caught_by_the_alignment(tiled_lib):
    lib = tiled_lib
    saved = dict(lib.df["oligo"])
    row = lib.df.index[lib.df["oligo"].notna()][0]
    oligo = lib.df.at[row, "oligo"]
    # Change one base of the tile body, past the 5' primer, site, spacer, and overhang.
    at = oligo.index(SITE) + len(SITE) + 1 + 4 + 6
    swapped = "A" if oligo[at] != "A" else "C"
    lib.df.at[row, "oligo"] = oligo[:at] + swapped + oligo[at + 1:]
    try:
        rep = lib.check()
        assert not rep.passed
        assert rep.assembly_correct == rep.assembly_checked - 1
    finally:
        lib.df["oligo"] = [saved[i] for i in lib.df.index]


# --- when there is nothing to simulate ----------------------------------------

def test_cassette_mode_simulates_the_coding_sequence_only():
    """No backbone, so the per-tile vector is a bare coding-region cassette. There is no
    plasmid to align against, but each oligo still has to rebuild its variant's CDS."""
    lib = SubstitutionScan(LibrarySpec.from_toml(GCK_TOML)).generate().codon_optimize().tile()
    rep = lib.check()

    assert rep.assembly_checked == int(lib.df["oligo"].notna().sum())
    assert rep.assembly_correct == rep.assembly_checked
    assert rep.assembly_aligned == 0
    assert not rep.assembly_issues


def test_no_destination_vector_means_no_simulation(tmp_path):
    spec = LibrarySpec(name="plain", protein_sequence=_protein(CLEAN_CDS), cds=CLEAN_CDS,
                       substitutions=["A"], adaptor_5=A5_CDS, adaptor_3=A3_CDS)
    lib = SubstitutionScan(spec).generate().codon_optimize()
    rep = lib.check()
    assert rep.assembly_checked == 0 and not rep.assembly_issues
    assert rep.passed
    assert lib.simulate_assembly() == []


def test_no_adaptors_means_no_reaction_to_simulate(tmp_path):
    """A pool with no cloning flanks carries no enzyme sites, so QC says so once (as an
    advisory) rather than reporting an assembly failure for every member."""
    lib = _standard_lib(_plasmid(tmp_path, CLEAN_CDS), "", "")
    rep = lib.check()
    assert rep.passed
    assert rep.assembly_checked == 0
    assert any("no adaptors are set" in a for a in rep.reference_advisories)


# --- the assembled clone, in hand ---------------------------------------------

def test_assembled_product_lines_up_with_the_parent(tmp_path):
    """The notebook-facing pair: two sequences in one frame that can be diffed directly.
    Deriving the rotation from the variant's own product would shift the frame by a base and
    make everything look different, so this pins that it comes from the WT."""
    lib = _standard_lib(_plasmid(tmp_path, CLEAN_CDS), A5_BACKBONE, A3_BACKBONE)
    parent = lib.parent_vector()
    assert _is_rotation(parent, BB5 + CLEAN_CDS + BB3)
    assert lib.assembled_product("WT") == parent          # the WT clone is the parent

    product = lib.assembled_product("D14A")
    diffs = [i for i, (a, b) in enumerate(zip(product, parent)) if a != b]
    assert len(product) == len(parent)
    assert 1 <= len(diffs) <= 3                            # one codon, nothing else
    assert max(diffs) - min(diffs) < 3


def test_assembled_product_names_what_it_cannot_assemble(tiled_lib):
    with pytest.raises(ValueError, match="No member named"):
        tiled_lib.assembled_product("not_a_variant")
    with pytest.raises(ValueError, match="rides on no oligo"):
        tiled_lib.assembled_product("WT")                  # unplaced in a tiled library
    assert tiled_lib.assembled_product("D14A") != tiled_lib.parent_vector()


# --- exporting the whole clone set --------------------------------------------

def test_to_assembled_vectors_writes_a_clone_per_variant(tmp_path):
    import pandas as pd
    from Bio import SeqIO

    lib = _standard_lib(_plasmid(tmp_path, CLEAN_CDS, extra_feature=("rep_origin", 0, 12, 1, "ori")),
                        A5_BACKBONE, A3_BACKBONE)
    d = tmp_path / "clones"
    lib.to_assembled_vectors(d)

    manifest = pd.read_csv(d / "assembled_vectors.csv")
    assert len(manifest) == len(lib)
    assert manifest["as_intended"].all()                    # every clone, WT included
    assert (manifest["length"] == len(lib.parent_vector())).all()
    assert manifest["codon"].isna().sum() == 1              # only the WT control has no codon

    # One GenBank per variant, plus the parent, plus the combined FASTA.
    assert len(list(d.glob("*.gb"))) == len(lib) + 1
    assert len(list(SeqIO.parse(str(d / "all_clones.fasta"), "fasta"))) == len(lib)

    row = manifest[manifest["name"] == "D14A"].iloc[0]
    rec = SeqIO.read(str(d / row["file"]), "genbank")
    assert str(rec.seq) == lib.assembled_product("D14A")
    labels = [(f.qualifiers.get("label") or [""])[0] for f in rec.features]
    assert any(x == "ori" for x in labels)                  # backbone feature carried over
    assert any("D14A" in x and ">" in x for x in labels)     # the codon swap is marked
    assert row["wt_codon"] != row["clone_codon"] and row["bases_changed"] >= 1


def test_an_amber_stop_gets_a_shell_safe_filename(tmp_path):
    spec = LibrarySpec(
        name="syn", protein_sequence=_protein(CLEAN_CDS), cds=CLEAN_CDS,
        substitutions=["TAG"], adaptor_5=A5_BACKBONE, adaptor_3=A3_BACKBONE,
        starting_vector=StartingVectorParams(path=_plasmid(tmp_path, CLEAN_CDS),
                                            insert_label="insert"),
    )
    lib = SubstitutionScan(spec).generate().codon_optimize()
    d = tmp_path / "clones"
    lib.to_assembled_vectors(d, fmt="genbank")
    names = {p.name for p in d.glob("*.gb")}
    assert "K2stop.gb" in names and not any("*" in n for n in names)
    assert not (d / "all_clones.fasta").exists()            # genbank only


def test_export_all_writes_the_cloning_outputs_by_default(tmp_path):
    lib = _standard_lib(_plasmid(tmp_path, CLEAN_CDS), A5_BACKBONE, A3_BACKBONE)
    lib.export_all(tmp_path / "out", method="pooled", plots=False)
    out = lib.output_dir

    assert (out / "destination_vector.csv").is_file()
    assert (out / "vector" / "destination.gb").is_file()
    assert (out / "assembled_vectors" / "parent_WT.gb").is_file()
    assert len(list((out / "assembled_vectors").glob("*.gb"))) == len(lib) + 1

    lib.export_all(tmp_path / "plain", method="pooled", plots=False, vectors=False)
    assert not (lib.output_dir / "assembled_vectors").exists()


# --- the results stay on the library ------------------------------------------

def test_simulate_assembly_leaves_its_results_on_the_library(tmp_path):
    lib = _standard_lib(_plasmid(tmp_path, CLEAN_CDS), A5_BACKBONE, A3_BACKBONE)
    assert lib.assembly is None

    results = lib.simulate_assembly()
    assert lib.assembly is results and len(results) == 1
    assert lib.assembly[0].product == lib.parent_vector()

    record = lib.design_specs["assembly"]["reactions"][0]
    assert record["label"] == "destination" and record["enzyme"] == "BsaI"
    assert record["members"] == record["rebuilt"] == record["aligned_to_parent"] == len(lib)
    assert record["product_length"] == len(lib.parent_vector())
    assert record["members_with_problems"] == [] and record["issues"] == []
    json.dumps(lib.design_specs)                          # the record stays serializable


def test_the_record_names_the_members_that_cannot_assemble(tmp_path):
    lib = _standard_lib(_plasmid(tmp_path, CLEAN_CDS), A5_CDS, A3_CDS)
    lib.simulate_assembly()
    record = lib.design_specs["assembly"]["reactions"][0]
    assert record["members_with_problems"] == ["K2A", "N59A", "P60A"]
    assert record["rebuilt"] == record["members"] - 3


@pytest.mark.parametrize("rebuild", ["codon_optimize", "drop_failed"])
def test_rebuilding_the_library_clears_a_stale_simulation(tmp_path, rebuild):
    """The results count members and carry their sequences, so anything that rebuilds those
    sequences makes them describe a library that no longer exists."""
    lib = _standard_lib(_plasmid(tmp_path, CLEAN_CDS), A5_BACKBONE, A3_BACKBONE)
    lib.simulate_assembly()
    assert lib.assembly and "assembly" in lib.design_specs

    getattr(lib, rebuild)()
    assert lib.assembly is None and "assembly" not in lib.design_specs


def test_tiling_clears_a_simulation_of_the_untiled_design(tmp_path):
    path = _plasmid(tmp_path, CLEAN_CDS)
    spec = LibrarySpec(name="syn", protein_sequence=_protein(CLEAN_CDS), cds=CLEAN_CDS,
                       substitutions=["A"], adaptor_5=A5_BACKBONE, adaptor_3=A3_BACKBONE,
                       starting_vector=StartingVectorParams(path=path, insert_label="insert"))
    lib = SubstitutionScan(spec).generate().codon_optimize()
    lib.simulate_assembly()
    assert len(lib.assembly) == 1

    lib.tile(TiledAssemblyParams(oligo_budget=150))
    assert lib.assembly is None
    assert len(lib.simulate_assembly()) == len(lib.tiles)      # now one per tile


def test_each_tile_is_simulated_against_its_own_destination_vector(tiled_lib):
    """A tiled library runs one reaction per tile, each in its own tube: that tile's oligos
    into that tile's vector. Pairing them up wrongly has to fail, or the simulation would be
    proving nothing about the layout."""
    from library_designer.checks.assembly import _plan, _released

    lib = tiled_lib
    plan = _plan(lib)
    assert len(plan.reactions) == len(lib.tiles)

    for (label, vector, _topo, start, end, members, _wt, _p5, _p3), tile in zip(
        plan.reactions, lib.tiles
    ):
        assert label == f"tile{tile.index}"
        assert vector is tile.vector                      # the tile's own drop-out vector
        assert (start, end) == (tile.start, tile.end)     # judged against its own window
        placed = {n for n, _cds, _idx, _oligo in members}
        assert placed == set(lib.df["name"][lib.df["tile"] == tile.index])

    # Every tile's vector is distinct, and an oligo ligates into its own and no other.
    assert len({t.vector for t in lib.tiles}) == len(lib.tiles)
    for i, tile in enumerate(lib.tiles):
        oligo = lib.df["oligo"][lib.df["tile"] == i].iloc[0]
        insert, _ = _released(oligo, "BsaI", "x")
        for j, other in enumerate(lib.tiles):
            product, _ = ligate(digest(other.vector, "BsaI", circular=True), insert)
            assert (product is not None) == (i == j), f"tile{i} oligo into tile{j} vector"


def test_every_tile_reaction_rebuilds_the_same_parent_plasmid(tiled_lib):
    """Each tile carries a different window, so the three reactions are different molecules,
    but they all have to converge on the one starting plasmid."""
    results = tiled_lib.simulate_assembly()
    assert len({r.product for r in results}) == 1
    assert results[0].product == tiled_lib.parent_vector()


def test_a_tiled_variant_assembles_through_its_own_tile(tiled_lib):
    lib, parent = tiled_lib, tiled_lib.parent_vector()
    for i in range(len(lib.tiles)):
        name = lib.df["name"][lib.df["tile"] == i].iloc[3]
        product = lib.assembled_product(name)
        diffs = [k for k, (a, b) in enumerate(zip(product, parent)) if a != b]
        assert diffs and max(diffs) - min(diffs) < 3       # one codon, in its own tile
