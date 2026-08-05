"""Tests for Golden Gate overhang specificity (checks/overhangs.py).

Two overhangs that are too much alike let the cut vector re-close empty or take the
fragment in backwards, so the check has to score every pair rather than only look for
exact duplicates. These cover the arithmetic, how a finding is graded, and what QC and the
destination-vector builder do with it.
"""
from __future__ import annotations

import pytest

from library_designer import LibrarySpec, StartingVectorParams, SubstitutionScan, TiledAssemblyParams
from library_designer.checks.overhangs import (
    MAX_SHARED,
    OverhangEnd,
    findings_for,
    overhang_ends,
    risk_of,
    self_shared,
    shared_bases,
)
from library_designer.regions import reverse_complement

# The same clean, BsaI-free codon block the vector tests use: 10 codons repeated, so the
# tile boundaries land on known bases.
_SAFE = ["GCT", "AAA", "CTG", "GAT", "ACC", "TTT", "CAA", "GTT", "AAT", "CCA"]
CLEAN_CDS = "".join(_SAFE * 6)                     # 60 codons, 180 bp
BB5 = "ACGTACGTTTGCAACGGATCCACAG"                  # backbone 5' of the insert
BB3 = "TGACCTAGGCATTACGTACGTACGT"                  # backbone 3' of the insert


def _protein(cds: str) -> str:
    from Bio.Seq import Seq

    return str(Seq(cds).translate())


def _write_gb(path, seq, *, topology="circular", features=()):
    from Bio import SeqIO
    from Bio.Seq import Seq
    from Bio.SeqFeature import SeqFeature, SimpleLocation
    from Bio.SeqRecord import SeqRecord

    rec = SeqRecord(Seq(seq), id="syn", name="syn", description="synthetic",
                    annotations={"molecule_type": "DNA", "topology": topology})
    for ftype, a, b, strand, label in features:
        rec.features.append(SeqFeature(SimpleLocation(a, b, strand=strand),
                                       type=ftype, qualifiers={"label": [label]}))
    SeqIO.write(rec, str(path), "genbank")
    return str(path)


def _plasmid(tmp_path, cds, bb5=BB5, bb3=BB3, name="p.gb"):
    return _write_gb(tmp_path / name, bb5 + cds + bb3,
                     features=[("CDS", len(bb5), len(bb5) + len(cds), 1, "insert")])


def _adaptors(bb5, bb3):
    """The adaptor convention where the adaptor itself spells the backbone-derived overhang,
    so the whole CDS drops out and the vector's overhangs are the flanking backbone bases."""
    return ("GGCGC" + "GGTCTC" + "A" + bb5[-4:],
            bb3[:4] + "T" + reverse_complement("GGTCTC") + "GCGCC")


def _standard_lib(path, bb5, bb3, cds=CLEAN_CDS):
    a5, a3 = _adaptors(bb5, bb3)
    spec = LibrarySpec(
        name="syn", protein_sequence=_protein(cds), cds=cds, substitutions=["A"],
        adaptor_5=a5, adaptor_3=a3,
        starting_vector=StartingVectorParams(path=path, insert_label="insert"),
    )
    return SubstitutionScan(spec).generate().codon_optimize()


def _tiled_lib(**tiled_kw):
    params = dict(oligo_budget=300, tile_size=90)
    params.update(tiled_kw)
    spec = LibrarySpec(
        name="syn_tiled", protein_sequence=_protein(CLEAN_CDS), cds=CLEAN_CDS,
        substitutions=["A"], tiled=TiledAssemblyParams(**params),
    )
    return SubstitutionScan(spec).generate().codon_optimize().tile()


def _end(reaction, end, seq):
    return OverhangEnd(reaction, end, seq)


# --- arithmetic and grading ---------------------------------------------------

def test_shared_bases_counts_aligned_matches():
    assert shared_bases("AGGT", "AGGT") == 4
    assert shared_bases("AGGT", "AGGA") == 3
    assert shared_bases("AGGT", "TCCA") == 0
    # Overhangs of different lengths cannot anneal, so they score 0 rather than raising.
    assert shared_bases("AGGT", "AGG") == 0


def test_self_shared_finds_palindromes():
    # A palindrome reads the same on both strands, so it matches its own reverse complement
    # at every base. An even-length overhang pairs its positions up, so the count is even.
    assert self_shared(_end("t", "5'", "GGCC")) == 4
    assert self_shared(_end("t", "5'", "AGGT")) in (0, 2)
    assert self_shared(_end("t", "5'", "AAAA")) == 0


def test_risk_tiers():
    assert risk_of(4, 4) == "collision"
    assert risk_of(3, 4) == "high"
    assert risk_of(2, 4) == "watch"
    assert risk_of(MAX_SHARED, 4) == "ok"
    assert risk_of(0, 4) == "ok"


# --- how a finding is graded --------------------------------------------------

def test_identical_ends_in_one_reaction_fail():
    fails, _ = findings_for([_end("tile0", "5'", "AGGT"), _end("tile0", "3'", "AGGT")])
    assert len(fails) == 1
    assert "re-closes" in fails[0] and "tile0" in fails[0]


def test_complementary_ends_in_one_reaction_fail():
    a = "AGGT"
    fails, _ = findings_for([_end("tile0", "5'", a),
                             _end("tile0", "3'", reverse_complement(a))])
    assert len(fails) == 1
    assert "backwards" in fails[0]


def test_palindromic_end_fails():
    fails, _ = findings_for([_end("tile0", "5'", "GGCC"), _end("tile0", "3'", "AAAT")])
    assert any("palindromic" in m for m in fails)


def test_one_mismatch_in_one_reaction_is_an_advisory_not_a_failure():
    fails, advisories = findings_for([_end("tile0", "5'", "AGGT"),
                                      _end("tile0", "3'", "AGGA")])
    assert fails == []
    assert any("3 of 4" in m for m in advisories)


def test_two_tiles_sharing_an_overhang_is_not_a_finding():
    """Each tile is amplified on its own and assembled into the vector built around its own
    window, so tile0's overhangs and tile1's never meet. Identical overhangs across two tiles
    are not a hazard and must not be reported as one."""
    fails, advisories = findings_for([_end("tile0", "5'", "AGGT"), _end("tile0", "3'", "TTTC"),
                                      _end("tile1", "5'", "AGGT"), _end("tile1", "3'", "CCAA")])
    assert fails == [] and advisories == []


def test_orthogonal_ends_are_clean():
    fails, advisories = findings_for([_end("tile0", "5'", "AGGT"), _end("tile0", "3'", "TCAA")])
    assert fails == [] and advisories == []


def test_pairs_above_target_are_collapsed_into_one_count():
    """Sharing 2 of 4 is above target but not a near match, so it is counted rather than
    spelled out one line at a time."""
    fails, advisories = findings_for([_end("t0", "5'", "AGGT"), _end("t0", "3'", "AGCA")])
    assert fails == []
    assert len(advisories) == 1
    assert advisories[0].startswith("1 reaction has")
    assert "overhang_pairs()" in advisories[0]


# --- the tables ---------------------------------------------------------------

def test_tiled_tables_cover_every_end_and_pair():
    lib = _tiled_lib()
    assert len(lib.tiles) == 2

    ends = lib.overhangs()
    assert len(ends) == 4                                  # two per tile
    assert list(ends["reaction"]) == ["tile0", "tile0", "tile1", "tile1"]
    assert list(ends["end"]) == ["5'", "3'", "5'", "3'"]

    # The table reads the overhangs the oligos actually carry, not a re-derivation.
    from library_designer.layout.tiled import tile_contexts
    want = []
    for t in lib.tiles:
        want += list(tile_contexts(lib.reference, t.start, t.end, lib.tiled_params))
    assert list(ends["overhang"]) == want

    # One row per tile by default: a tile's two ends are the only pair that ever meets.
    pairs = lib.overhang_pairs()
    assert len(pairs) == 2
    assert pairs["same_reaction"].all()
    assert set(pairs["risk"]) <= {"ok", "watch", "high", "collision"}


def test_pair_table_is_sorted_worst_first():
    lib = _tiled_lib()
    worst = lib.overhang_pairs()[["shared", "shared_flipped"]].max(axis=1)
    assert list(worst) == sorted(worst, reverse=True)


def test_all_pairs_adds_the_cross_tile_rows_ungraded():
    """Asked for, the cross-tile pairs are listed with their counts. They carry no verdict,
    since the two tiles are assembled in separate tubes."""
    lib = _tiled_lib()
    pairs = lib.overhang_pairs(all_pairs=True)
    assert len(pairs) == 6                                 # C(4, 2)
    cross = pairs[~pairs["same_reaction"]]
    assert len(cross) == 4
    assert set(cross["risk"]) == {"n/a"}
    assert (cross["note"] == "").all()
    # The graded rows sort ahead of the ungraded ones.
    assert pairs["same_reaction"].iloc[0] and not pairs["same_reaction"].iloc[-1]


def test_tables_are_empty_without_a_golden_gate_reaction():
    spec = LibrarySpec(name="plain", protein_sequence=_protein(CLEAN_CDS), cds=CLEAN_CDS,
                       substitutions=["A"])
    lib = SubstitutionScan(spec).generate().codon_optimize()
    assert overhang_ends(lib) == []
    assert lib.overhangs().empty and lib.overhang_pairs().empty
    # The columns are still there, so a notebook cell does not have to branch.
    assert "risk" in lib.overhangs().columns and "risk" in lib.overhang_pairs().columns


def test_figure_renders():
    fig = _tiled_lib().plot_overhangs()
    assert fig.get_axes()


def test_figure_refuses_a_design_with_no_reaction():
    spec = LibrarySpec(name="plain", protein_sequence=_protein(CLEAN_CDS), cds=CLEAN_CDS,
                       substitutions=["A"])
    lib = SubstitutionScan(spec).generate().codon_optimize()
    with pytest.raises(ValueError, match="No fused overhangs"):
        lib.plot_overhangs()


# --- QC on a tiled library ----------------------------------------------------

def test_tiled_terminal_overhang_colliding_with_a_boundary_fails_qc():
    """tile0 takes its 5' overhang from the vector context and its 3' from the CDS at the
    boundary. Setting the context to those same bases makes the two ends identical."""
    boundary = CLEAN_CDS[90:94]
    lib = _tiled_lib(vector_context_5=boundary)
    rep = lib.check()
    assert not rep.passed
    assert any("tile0" in m and "re-closes" in m for m in rep.overhang_issues)


def test_clean_tiled_library_passes_with_no_overhang_failures():
    lib = _tiled_lib()
    rep = lib.check()
    assert rep.overhang_issues == []
    # Advisories never fail the report.
    assert rep.passed or not rep.overhang_issues


def test_advisories_stay_out_of_passed():
    """A near match is reported without failing, since a tiled design reads its overhangs
    off the CDS and cannot simply pick different bases."""
    lib = _tiled_lib()
    rep = lib.check()
    rep.overhang_advisories = ["something worth reading"]
    assert rep.passed == (not rep.overhang_issues and rep.passed)
    assert "something worth reading" in rep.text()


# --- the standard (untiled) destination vector --------------------------------

def test_standard_vector_with_clean_overhangs_builds():
    import tempfile
    from pathlib import Path

    with tempfile.TemporaryDirectory() as d:
        lib = _standard_lib(_plasmid(Path(d), CLEAN_CDS), BB5, BB3)
        dv = lib.destination_vector()
        assert (dv.overhang_5, dv.overhang_3) == (BB5[-4:], BB3[:4])
        assert lib.check().overhang_issues == []
        pairs = lib.overhang_pairs()
        assert len(pairs) == 1 and pairs.loc[0, "same_reaction"]


def _colliding_lib(tmp_path, bb3):
    return _standard_lib(_plasmid(tmp_path, CLEAN_CDS, bb5=BB5, bb3=bb3), BB5, bb3)


def test_vector_whose_own_overhangs_are_identical_raises(tmp_path):
    """The flanking backbone bases decide the overhangs here, so a plasmid that presents the
    same four on both sides cannot clone at all. That is the user's vector, not something we
    can recode, so building it stops."""
    bb3 = BB5[-4:] + "CTAGGCATTACGTACGTACG"          # 3' flank starts with the 5' overhang
    with pytest.raises(ValueError, match="overhangs that collide"):
        _colliding_lib(tmp_path, bb3).destination_vector()


def test_vector_whose_overhangs_are_complementary_raises(tmp_path):
    bb3 = reverse_complement(BB5[-4:]) + "CTAGGCATTACGTACGTACG"
    with pytest.raises(ValueError, match="backwards"):
        _colliding_lib(tmp_path, bb3).destination_vector()


def test_strict_false_builds_the_colliding_vector_anyway(tmp_path):
    bb3 = BB5[-4:] + "CTAGGCATTACGTACGTACG"
    dv = _colliding_lib(tmp_path, bb3).destination_vector(strict=False)
    assert dv.overhang_5 == dv.overhang_3


def test_qc_reports_the_collision_instead_of_blowing_up(tmp_path):
    """check() has to come back with a usable report, so it builds the vector non-strictly
    and files the collision under the overhang findings."""
    bb3 = BB5[-4:] + "CTAGGCATTACGTACGTACG"
    rep = _colliding_lib(tmp_path, bb3).check()
    assert not rep.passed
    assert any("re-closes" in m for m in rep.overhang_issues)


def test_pair_table_still_works_on_a_colliding_vector(tmp_path):
    """The table is most wanted exactly when the layout is broken, so it must not go empty."""
    bb3 = BB5[-4:] + "CTAGGCATTACGTACGTACG"
    pairs = _colliding_lib(tmp_path, bb3).overhang_pairs()
    assert len(pairs) == 1
    assert pairs.loc[0, "risk"] == "collision"


# --- the two events, from first principles ------------------------------------
#
# `shared` and `shared_flipped` are two different comparisons of the same pair of
# overhangs, standing for two different reactions. These derive both from the sticky ends
# themselves, so a change to either column has to survive the annealing arithmetic.

_COMP = {"A": "T", "T": "A", "G": "C", "C": "G"}


def _anneal(a: str, b: str) -> int:
    """Watson-Crick pairs when two 5'->3' single-stranded overhangs meet antiparallel."""
    return sum(x == _COMP[y] for x, y in zip(a, b[::-1]))


def _sticky_ends(desc_5: str, desc_3: str) -> dict[str, str]:
    """The four sticky ends a Golden Gate reaction actually presents, as ssDNA read
    5'->3', from the two top-strand overhang descriptors.

    A 5'-overhang cut leaves the right-hand fragment's overhang on the top strand and the
    left-hand fragment's on the bottom strand. So the insert carries its 5' overhang on the
    top strand and its 3' overhang on the bottom, and the cut vector carries the reverse
    complement of each at the end that receives it.
    """
    return {
        "insert_5": desc_5,
        "insert_3": reverse_complement(desc_3),
        "vector_x": reverse_complement(desc_5),   # normally takes the insert's 5' end
        "vector_y": desc_3,                       # normally takes the insert's 3' end
    }


def test_the_intended_assembly_anneals_perfectly_at_both_junctions():
    e = _sticky_ends("AAGC", "GGTG")
    assert _anneal(e["insert_5"], e["vector_x"]) == 4
    assert _anneal(e["insert_3"], e["vector_y"]) == 4


def test_shared_scores_vector_re_ligation_and_shared_flipped_scores_flipping():
    """The reported row for ends (AAGC, GGTG): 0 as written, 2 reverse-complemented.

    Read left to right the bottom strand of the 5' locus (TTCG) and the top strand of the
    3' (GGTG) look like they share a base, but that reading is not the alignment they
    anneal in. Scored antiparallel, the cut vector cannot re-close on itself at all.
    """
    desc_5, desc_3 = "AAGC", "GGTG"
    e = _sticky_ends(desc_5, desc_3)

    assert _anneal(e["vector_x"], e["vector_y"]) == shared_bases(desc_5, desc_3) == 0
    assert _anneal(e["insert_3"], e["vector_x"]) == \
        shared_bases(desc_5, reverse_complement(desc_3)) == 2


def test_both_columns_hold_for_every_overhang_pair():
    from itertools import product

    words = ["".join(p) for p in product("ACGT", repeat=4)]
    for desc_5 in words:
        for desc_3 in words[::37]:          # a spread of 7 partners each, keeps it quick
            e = _sticky_ends(desc_5, desc_3)
            # Intended assembly always works, whatever the pair.
            assert _anneal(e["insert_5"], e["vector_x"]) == 4
            assert _anneal(e["insert_3"], e["vector_y"]) == 4
            # The cut vector closing on itself.
            assert _anneal(e["vector_x"], e["vector_y"]) == shared_bases(desc_5, desc_3)
            # The fragment going in backwards.
            assert _anneal(e["insert_3"], e["vector_x"]) == \
                shared_bases(desc_5, reverse_complement(desc_3))


def test_a_pair_can_be_clean_one_way_and_not_the_other():
    """Which is the reason for two columns. Identical ends re-ligate but do not flip;
    reverse-complementary ends flip but do not re-ligate."""
    same = _sticky_ends("AAGC", "AAGC")
    assert _anneal(same["vector_x"], same["vector_y"]) == 4       # re-closes empty
    assert _anneal(same["insert_3"], same["vector_x"]) == shared_bases(
        "AAGC", reverse_complement("AAGC"))                        # not a flip risk

    flip = _sticky_ends("AAGC", reverse_complement("AAGC"))
    assert _anneal(flip["insert_3"], flip["vector_x"]) == 4        # goes in backwards
    assert _anneal(flip["vector_x"], flip["vector_y"]) == shared_bases(
        "AAGC", reverse_complement("AAGC"))                        # not a re-ligation risk
