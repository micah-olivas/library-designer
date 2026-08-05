"""Tests for the mispriming check (checks/mispriming.py).

A pooled library is amplified with the constant regions on its oligos, so those regions have
to have one binding site each. These cover the 3'-anchored matching, where the check looks,
what QC does with a finding, and the opt-in primer screen in the tiler.
"""
from __future__ import annotations

import pandas as pd
import pytest
from Bio.Seq import Seq

from library_designer import LibrarySpec, SubstitutionScan, TiledAssemblyParams
from library_designer.checks.mispriming import (
    HIGH_ANNEAL,
    MIN_ANNEAL,
    MIN_HANDLE,
    TERMINAL_CLEAN,
    Handle,
    PrimingSite,
    Region,
    anneal_sites,
    tolerant_sites,
    flank_regions,
    handles,
    mispriming_findings,
    priming_sites,
    variable_regions,
    vector_regions,
)
from library_designer.primers import load_primer_set
from library_designer.regions import reverse_complement

# A 64-codon CDS with no BsaI site, no repeated 12-mer, and no Shine-Dalgarno-like motif, so a
# planted match is the only one the check can find. 192 bp splits into three tiles at
# tile_size=90.
CDS = ("CTATGGGATTGCTTGTATCTCCTAGAATGCGCCATAAGTTCCGAACTGGGAGCTTCAAGTACGCTACATAGCACCGAGTCC"
       "TGTGTAGGAAAGGCGCGCCTAGGCGTTGCTGGCAGGGCTAAATCTGACACAAATAAGATATGTATTGCGAACTATCTCCCA"
       "CTAGCCTCTATAGACATGGGACAGCTTAAA")
PROTEIN = str(Seq(CDS).translate())
BB5 = "ACGTACGTTTGCAACGGATCCACAG"       # backbone 5' of the insert
BB3 = "TGACCTAGGCATTACGTACGTACGT"       # backbone 3' of the insert


def _write_gb(path, seq, *, topology="circular", features=()):
    from Bio import SeqIO
    from Bio.SeqFeature import SeqFeature, SimpleLocation
    from Bio.SeqRecord import SeqRecord

    rec = SeqRecord(Seq(seq), id="syn", name="syn", description="synthetic",
                    annotations={"molecule_type": "DNA", "topology": topology})
    for ftype, a, b, strand, label in features:
        rec.features.append(SeqFeature(SimpleLocation(a, b, strand=strand),
                                       type=ftype, qualifiers={"label": [label]}))
    SeqIO.write(rec, str(path), "genbank")
    return str(path)


def _plasmid(tmp_path, cds=CDS, *, topology="circular", name="p.gb"):
    return _write_gb(tmp_path / name, BB5 + cds + BB3, topology=topology,
                     features=[("CDS", len(BB5), len(BB5) + len(cds), 1, "insert")])


def _pool_csv(tmp_path, planted=(), n_clean=12, name="pool.csv"):
    """A primer-set CSV: ``planted`` first, then clean primers off the bundled orthogonal
    pool, so a design draws the planted ones unless something makes it pass over them."""
    clean = [seq for _, seq in load_primer_set("subramanian2018").primers[:n_clean]]
    path = tmp_path / name
    pd.DataFrame({"primer_id": [f"p{i}" for i in range(len(planted) + len(clean))],
                  "sequence": list(planted) + clean}).to_csv(path, index=False)
    return str(path)


def _tiled(tmp_path=None, *, planted=(), cds=CDS, **tiled_kw):
    params = dict(oligo_budget=300, tile_size=90)
    if planted:
        params["primer_set"] = _pool_csv(tmp_path, planted)
    params.update(tiled_kw)
    spec = LibrarySpec(name="syn_tiled", protein_sequence=str(Seq(cds).translate()), cds=cds,
                       substitutions=["A"], tiled=TiledAssemblyParams(**params))
    return SubstitutionScan(spec).generate().codon_optimize().tile()


def _standard(adaptor_5="GGCGCGGTCTCAACAG", adaptor_3="TGACTGAGACCGCGCC", **spec_kw):
    spec = LibrarySpec(name="syn", protein_sequence=PROTEIN, cds=CDS, substitutions=["A"],
                       adaptor_5=adaptor_5, adaptor_3=adaptor_3, **spec_kw)
    return SubstitutionScan(spec).generate().codon_optimize()


# --- the 3'-anchored match ----------------------------------------------------

def test_a_full_length_match_is_found_on_either_strand():
    handle = "ACGTTGCATTGGCCAATTGA"
    target = "TTTT" + handle + "TTTT"
    assert anneal_sites(handle, target, 12) == [(len(handle), "+", 4)]
    flipped = "TTTT" + reverse_complement(handle) + "TTTT"
    assert anneal_sites(handle, flipped, 12) == [(len(handle), "-", 4)]


def test_only_the_3_end_anchors_a_match():
    """A polymerase extends from the 3' end, so that is the end that has to be annealed. The
    same 12 bases taken off the 5' end of the handle are not a priming site."""
    handle = "ACGTTGCATTGGCCAATTGA"
    assert anneal_sites(handle, "TTTT" + handle[-12:] + "TTTT", 12) == [(12, "+", 4)]
    assert anneal_sites(handle, "TTTT" + handle[:12] + "TTTT", 12) == []


def test_one_site_per_3_end():
    """A 14-base run is one site, not also its 13- and 12-base suffixes."""
    handle = "ACGTTGCATTGGCCAATTGA"
    sites = anneal_sites(handle, "TTTT" + handle[-14:] + "TTTT", 12)
    assert sites == [(14, "+", 4)]


def test_min_len_is_respected():
    handle = "ACGTTGCATTGGCCAATTGA"
    # The bases either side must not extend the run, or the target would carry a longer one
    # than the test plants.
    target = "ACCA" + handle[-11:] + "ACCA"
    assert anneal_sites(handle, target, 12) == []
    assert anneal_sites(handle, target, 11) == [(11, "+", 4)]


def test_two_separate_sites_are_both_reported():
    handle = "ACGTTGCATTGGCCAATTGA"
    target = handle[-13:] + "AAAA" + reverse_complement(handle[-12:])
    lengths = {(length, strand) for length, strand, _ in anneal_sites(handle, target, 12)}
    assert lengths == {(13, "+"), (12, "-")}


def _site(length, handle_len=20, mismatches=()):
    return PrimingSite(Handle("h", "A" * handle_len), Region("r", "A" * 40),
                       length, "+", "r", 0, "A" * length, length + len(mismatches), mismatches)


def test_risk_tiers():
    assert _site(20).risk == "collision" and _site(20).full
    assert _site(HIGH_ANNEAL).risk == "high"
    assert _site(MIN_ANNEAL).risk == "watch"
    assert not _site(HIGH_ANNEAL).full


# --- duplexes that carry a mismatch -------------------------------------------

def test_a_long_duplex_broken_by_one_base_is_found():
    """The case exact matching scores worst: the last 5 bases pair, one mismatches, then 14
    more pair. A 19-of-20 duplex primes readily, and the exact run is only 5 bases."""
    primer = "ACGTTGCATTGGCCAATTGA"
    target = "TTTT" + primer[:14] + "C" + primer[-5:] + "TTTT"     # primer[-6] swapped
    assert [(p, m) for p, m, *_ in anneal_sites(primer, target, 12)] == []
    hits = tolerant_sites(primer, target, 12, 1)
    assert len(hits) == 1
    paired, mismatches, aligned, strand, pos = hits[0]
    assert (paired, mismatches, aligned, strand) == (19, (5,), 20, "+")
    assert target[pos:pos + aligned] == primer[:14] + "C" + primer[-5:]


def test_a_mismatch_in_the_last_bases_disqualifies_the_duplex():
    """Taq extends from a paired 3' terminus. A mismatch at or beside the 3' base leaves
    nothing to extend, however much pairs behind it."""
    primer = "ACGTTGCATTGGCCAATTGA"
    for d in range(TERMINAL_CLEAN):
        swapped = primer[: len(primer) - 1 - d] + ("C" if primer[-1 - d] != "C" else "G") \
            + primer[len(primer) - d:]
        assert tolerant_sites(primer, "TTTT" + swapped + "TTTT", 12, 1) == []
    # One base further in than the clean window is allowed.
    d = TERMINAL_CLEAN
    swapped = primer[: len(primer) - 1 - d] + ("C" if primer[-1 - d] != "C" else "G") \
        + primer[len(primer) - d:]
    assert tolerant_sites(primer, "TTTT" + swapped + "TTTT", 12, 1)


def test_the_budget_caps_how_many_mismatches_a_duplex_may_carry():
    primer = "ACGTTGCATTGGCCAATTGA"
    two = primer[:8] + "C" + primer[9:14] + "C" + primer[15:]     # bases 8 and 14 swapped
    assert tolerant_sites(primer, two, 12, 1) == []
    hits = tolerant_sites(primer, two, 12, 2)
    # Counted in from the 3' base, which sits at index 19, so the two land 11 and 5 bases in.
    assert [(p, m) for p, m, *_ in hits] == [(18, (5, 11))]


def test_a_trailing_mismatch_is_not_counted_in_the_duplex():
    """A duplex ends on a paired base: a mismatch at its far end adds nothing to it."""
    primer = "ACGTTGCATTGGCCAATTGA"
    target = "TTTT" + "C" + primer[1:] + "TTTT"                     # only base 0 mismatches
    paired, mismatches, aligned, _, _ = tolerant_sites(primer, target, 12, 1)[0]
    assert (paired, mismatches, aligned) == (19, (), 19)


def test_with_no_budget_the_tolerant_scan_reproduces_the_exact_one():
    """Two code paths score the same sites, so the fast exact scan can stay the default without
    the two drifting apart."""
    import random

    rng = random.Random(4)
    for _ in range(60):
        primer = "".join(rng.choice("ACGT") for _ in range(rng.choice([9, 15, 20])))
        region = "".join(rng.choice("ACGT") for _ in range(120))
        # Plant the primer's 3' end somewhere, on one strand or the other, so there is
        # something to find rather than testing two empty lists against each other.
        run = primer[-rng.randint(6, len(primer)):]
        at = rng.randint(0, 80)
        piece = run if rng.random() < 0.5 else reverse_complement(run)
        region = region[:at] + piece + region[at + len(piece):]
        for least in (6, 12):
            exact = [(p, s, x) for p, s, x in anneal_sites(primer, region, least)]
            tolerant = [(p, s, x) for p, m, a, s, x in tolerant_sites(primer, region, least, 0)]
            assert sorted(exact) == sorted(tolerant), (primer, region, least)


# --- which flanks are checked -------------------------------------------------

def test_tiled_handles_are_the_tile_primers():
    lib = _tiled()
    assert len(lib.tiles) == 3
    hs = handles(lib)
    assert [h.label for h in hs] == [f"tile{i} {end}" for i in range(3) for end in ("fwd", "rev")]
    assert [h.seq for h in hs] == [s for t in lib.tiles for s in (t.fwd, t.rev)]


def test_a_3_adaptor_is_checked_as_the_primer_you_would_order():
    """The 3' adaptor sits on the top strand, so the primer that pairs with it is its reverse
    complement, and the 3' end of that primer is the adaptor base nearest the variable
    region."""
    lib = _standard()
    hs = {h.label: h.seq for h in handles(lib)}
    assert hs["adaptor_5"] == lib.spec.adaptor_5.upper()
    assert hs["adaptor_3"] == reverse_complement(lib.spec.adaptor_3.upper())


def test_a_library_with_no_flanks_has_nothing_to_check():
    spec = LibrarySpec(name="plain", protein_sequence=PROTEIN, cds=CDS, substitutions=["A"])
    lib = SubstitutionScan(spec).generate().codon_optimize()
    assert handles(lib) == [] and priming_sites(lib) == []
    assert lib.mispriming().empty
    # The columns are still there, so a notebook cell does not have to branch.
    assert "risk" in lib.mispriming().columns


def test_an_adaptor_too_short_to_score_is_named_rather_than_dropped_in_silence():
    """An adaptor can be shorter than a primer, being site plus spacer plus overhang. It is not
    scored, and the report has to say so, or it reads as a clean check of something that was
    never checked."""
    short = "GGTCTCA"
    assert len(short) < MIN_HANDLE
    lib = _standard(adaptor_5=short, adaptor_3="")
    assert priming_sites(lib) == []
    failures, advisories = mispriming_findings(lib)
    assert failures == []
    assert any("adaptor_5" in m and "too short to score" in m for m in advisories)
    assert any("too short to score" in m for m in lib.check().mispriming_advisories)


def test_a_supplied_primer_too_short_to_score_is_refused():
    """A sequence you passed in has to come back with an answer. Scoring a 3-base run would
    report a site in any sequence, and reporting nothing reads like a clean result."""
    lib = _standard()
    with pytest.raises(ValueError, match="too short to score"):
        lib.mispriming(extra=["GGG"])
    with pytest.raises(ValueError, match="extra primer 2"):
        lib.mispriming(extra=["ACGTTGCATTGGCCAATTGA", "AAA"])


def test_a_supplied_primer_is_cleaned_before_it_is_scored():
    """Pasted primers arrive with spaces and in either case."""
    lib = _standard()
    run = lib.reference[40:60]
    spaced = " ".join([run[:5].lower(), run[5:12], run[12:].lower()])
    assert [s.matched for s in priming_sites(lib, extra=[spaced])] == [run]


def test_a_degenerate_supplied_primer_is_refused():
    """Matching is exact, so an ambiguity code would be scored as matching nothing."""
    lib = _standard()
    with pytest.raises(ValueError, match="non-ACGT"):
        lib.mispriming(extra=["ACGTTGCANTGGCCAATTGA"])


def test_flank_regions_are_the_oligo_outside_its_tile_window():
    """The flanks have to be what the oligos carry, or the check would score a molecule the
    pool does not hold."""
    lib = _tiled()
    regions = {r.label: r.seq for r in flank_regions(lib)}
    oligos = {int(t): o for t, o in zip(lib.df["tile"], lib.df["oligo"]) if isinstance(o, str)}
    for t in lib.tiles:
        oligo = oligos[t.index]
        assert oligo.startswith(regions[f"tile{t.index} 5' flank"])
        assert oligo.endswith(regions[f"tile{t.index} 3' flank"])


def test_the_vector_backbone_is_checked_without_the_coding_region(tmp_path):
    lib = _tiled(tmp_path, starting_vector=_plasmid(tmp_path), insert_label="insert")
    regions = vector_regions(lib)
    assert [r.label for r in regions] == ["vector backbone"]      # circular: one arc
    assert regions[0].seq == BB3 + BB5                            # from the insert round to it
    assert CDS not in regions[0].seq


def test_a_linear_backbone_is_checked_as_two_pieces(tmp_path):
    """Joining the two flanks would make a junction the molecule does not have."""
    lib = _tiled(tmp_path, starting_vector=_plasmid(tmp_path, topology="linear"),
                 insert_label="insert", topology="linear")
    assert [(r.label, r.seq) for r in vector_regions(lib)] == [
        ("vector backbone 5'", BB5), ("vector backbone 3'", BB3)
    ]


def test_no_plasmid_means_no_vector_region():
    assert vector_regions(_tiled()) == []


# --- where a match can come from ----------------------------------------------

def test_a_match_in_the_reference_is_reported_against_the_reference():
    lib = _tiled()
    planted = lib.reference[30:50]
    sites = priming_sites(lib, extra=[planted])
    assert [(s.region.label, s.length, s.position) for s in sites] == [("reference CDS", 20, 30)]
    assert sites[0].risk == "collision"


def test_a_match_a_mutated_codon_creates_is_named_after_its_member():
    """A member is the reference with one codon swapped, so a run that touches that codon is
    the only sequence it has that the reference does not."""
    lib = _tiled()
    row = lib.df[lib.df["mut_index"].notna()].iloc[5]
    codon = int(row["mut_index"]) * 3
    dna = row["variable_dna"]
    planted = dna[codon - 9:codon + 12]              # spans the mutated codon
    assert planted not in lib.reference
    sites = priming_sites(lib, extra=[planted])
    assert [(s.region.label, s.where) for s in sites] == [("variant codons", row["name"])]


def test_a_run_already_in_the_reference_is_not_reported_once_per_member():
    """Every member carries the reference outside its own codon, so a run that misses the
    codon would otherwise be reported thousands of times over."""
    lib = _tiled()
    sites = priming_sites(lib, extra=[lib.reference[60:80]])
    assert {s.region.label for s in sites} == {"reference CDS"}


def test_members_packed_into_one_region_cannot_share_a_match():
    """The per-member windows are searched as one string, and the separator has to stop a run
    from spanning two members and inventing a site no oligo carries."""
    lib = _tiled()
    packed = next(r for r in variable_regions(lib, 19) if r.label == "variant codons")
    join = packed.seq.index("N")
    # The 12 bases either side of the join, with the separator taken out: sequence no member
    # holds, so it must not be found.
    across = packed.seq[join - 6:join] + packed.seq[join + 1:join + 7]
    assert len(across) == 12 and "N" not in across
    assert anneal_sites(across, packed.seq, 12) == []


def test_a_primer_that_would_pull_out_another_tile_is_caught():
    """The set is orthogonal on paper. This is the check that says so for the pool as built."""
    lib = _tiled()
    other = lib.tiles[1].fwd
    sites = priming_sites(lib, extra=[other])
    assert [(s.region.label, s.risk) for s in sites] == [("tile1 5' flank", "collision")]


def test_a_handle_is_not_reported_at_its_own_binding_site():
    """Every primer occurs in full on its own oligo, which is what it is for."""
    lib = _tiled()
    assert [s for s in priming_sites(lib) if s.risk == "collision"] == []


def test_the_bundled_tiled_example_has_no_mispriming_failures():
    lib = _tiled()
    failures, _ = mispriming_findings(lib)
    assert failures == []


# --- QC -----------------------------------------------------------------------

def test_a_primer_that_primes_on_the_cds_fails_qc(tmp_path):
    """A primer with a second binding site cannot amplify one sublibrary out of the pool, so
    this is a failure and not an advisory."""
    planted = CDS[33:53]
    lib = _tiled(tmp_path, planted=[planted])
    assert lib.tiles[0].fwd == planted
    rep = lib.check()
    assert not rep.passed
    assert any("tile0 fwd" in m and "second binding site" in m for m in rep.mispriming_issues)
    assert "mispriming_issues" in rep.issues


def test_a_primer_that_primes_on_the_backbone_fails_qc(tmp_path):
    planted = (BB3 + BB5)[15:35]
    lib = _tiled(tmp_path, planted=[planted],
                 starting_vector=_plasmid(tmp_path), insert_label="insert")
    rep = lib.check()
    assert any("vector backbone" in m for m in rep.mispriming_issues)


def _one_off(site: str, at: int = 14) -> str:
    """``site`` with the base at ``at`` swapped, so it pairs with everything but that one."""
    return site[:at] + ("C" if site[at] != "C" else "G") + site[at + 1:]


def test_a_primer_pairing_with_one_mismatch_is_reported_but_never_fails(tmp_path):
    """The case exact matching scores as five bases, end to end: 19 of 20 pair, which primes,
    so QC has to say so. It stays an advisory, since no verdict should rest on how a mismatch
    is weighted."""
    planted = _one_off(CDS[40:60])
    lib = _tiled(tmp_path, planted=[planted])
    assert lib.tiles[0].fwd == planted
    rep = lib.check()
    assert rep.mispriming_issues == []
    assert any("tile0 fwd" in m and "mismatch" in m for m in rep.mispriming_advisories)
    top = lib.mispriming().iloc[0]
    assert (top["paired"], top["aligned"], top["mismatches"], top["risk"]) == (19, 20, 1, "high")
    # Exact scoring alone sees five bases and lets it through.
    assert lib.mispriming(mismatches=0).empty


def test_screening_passes_over_a_primer_that_would_pair_with_one_mismatch(tmp_path):
    """The screen and the check are scored on the same terms, so what one would report the
    other does not draw."""
    planted = _one_off(CDS[40:60])
    path = _pool_csv(tmp_path, [planted])
    assert _tiled(tmp_path, primer_set=path).tiles[0].fwd == planted
    screened = _tiled(tmp_path, primer_set=path, screen_primers=True)
    assert screened.tiles[0].fwd != planted
    assert screened.check().mispriming_advisories == []


def test_a_partial_3_run_is_an_advisory_and_does_not_fail_the_report():
    """Whether a partial run primes depends on the annealing temperature, so it is reported
    and left to the user."""
    lib = _tiled()
    run = lib.reference[40:40 + HIGH_ANNEAL]
    sites = priming_sites(lib, extra=["TTGCA" + run])
    assert [s.risk for s in sites] == ["high"]
    rep = lib.check()
    rep.mispriming_advisories = ["something worth reading"]
    assert rep.passed
    assert "something worth reading" in rep.text()


def test_short_runs_are_collapsed_into_one_count(tmp_path):
    # 20 bases, the length the tiles are sized for, whose last 12 are lifted off the CDS.
    lib = _tiled(tmp_path, planted=["TTGCACCA" + CDS[80:80 + MIN_ANNEAL]])
    failures, advisories = mispriming_findings(lib)
    assert failures == []
    assert any("shorter duplex" in m and "lib.mispriming()" in m for m in advisories)


def test_the_table_reads_worst_first(tmp_path):
    lib = _tiled(tmp_path, planted=[CDS[33:53]])
    tab = lib.mispriming()
    assert list(tab["paired"]) == sorted(tab["paired"], reverse=True)
    assert tab.iloc[0]["risk"] == "collision"
    assert tab.iloc[0]["handle"] == "tile0 fwd"
    assert tab.iloc[0]["mismatches"] == 0
    assert list(tab["aligned"]) == [p + m for p, m in zip(tab["paired"], tab["mismatches"])]


def test_min_anneal_can_be_lowered():
    """Scored without a mismatch budget, so the threshold is the only thing moving."""
    lib = _tiled()
    run = lib.reference[70:78]                       # 8 bases, below the default threshold
    assert lib.mispriming(extra=["TTGCACCAGTTA" + run], mismatches=0).empty
    assert len(lib.mispriming(extra=["TTGCACCAGTTA" + run], min_anneal=8, mismatches=0)) == 1


def test_a_standard_pool_checks_its_adaptors_against_the_variable_region():
    """The one-oligo case: the adaptors are the constant flanks, so they are the handles. The
    3' adaptor is scored as the primer you would order for it, its reverse complement, so
    planting that primer in the CDS means planting the adaptor reverse-complemented."""
    lib = _standard(adaptor_3=reverse_complement(CDS[100:120]))
    rep = lib.check()
    assert any("adaptor_3" in m for m in rep.mispriming_issues)
    assert not rep.passed


def test_a_partial_run_on_a_3_adaptor_is_scored_from_the_base_nearest_the_insert():
    """The reverse primer's 3' end is the 5' base of the adaptor, the one next to the variable
    region, so that is the end a partial run has to reach."""
    a3 = reverse_complement(CDS[100:116]) + "GAGACCGCGCC"
    lib = _standard(adaptor_3=a3)
    # Exact runs only, so the orientation is the only thing under test.
    sites = [s for s in priming_sites(lib, mismatches=0) if s.handle.label == "adaptor_3"]
    assert len(sites) == 1
    site = sites[0]
    assert site.region.label == "reference CDS" and site.strand == "+"
    # The run is at least the 16 bases planted, and it ends where the adaptor meets the insert.
    assert site.length >= 16 and site.position + site.length == 116
    assert site.risk == "high"
    rep = lib.check()
    assert rep.mispriming_issues == []
    assert any("adaptor_3" in m for m in rep.mispriming_advisories)


def test_the_two_adaptors_are_checked_against_each_other():
    a5 = "GGCGCGGTCTCAACAGTTGCATTGGCCAATTGA"
    lib = _standard(adaptor_5=a5, adaptor_3=reverse_complement(a5) + "TTTT")
    sites = priming_sites(lib)
    assert any(s.handle.label == "adaptor_5" and s.region.label == "adaptor_3" for s in sites)


def test_a_supplied_primer_is_not_reported_at_the_flank_it_is_for():
    """The usual way to check a real amplification primer: a tail plus the adaptor at its 3'
    end. The adaptor is where it is meant to bind, so that is not a finding."""
    lib = _standard()
    tailed_fwd = "ACACTCTTTCCCTACACGACGCTCTTCCGATCT" + lib.spec.adaptor_5.upper()
    tailed_rev = "GTGACTGGAGTTCAGACGTGT" + reverse_complement(lib.spec.adaptor_3.upper())
    assert priming_sites(lib, extra=[tailed_fwd, tailed_rev]) == []
    # It is still checked everywhere else: the same primer aimed into the coding sequence is.
    into_cds = "ACACTCTTTCCCTACACGACG" + CDS[40:60]
    assert [s.region.label for s in priming_sites(lib, extra=[into_cds])] == ["reference CDS"]


def test_a_supplied_primer_on_a_tiled_pool_reports_the_tile_it_would_amplify():
    """With many flank pairs there is no single site a supplied primer can be assumed to be
    for, and which tile it amplifies is the thing worth knowing."""
    lib = _tiled()
    tailed = "ACACTCTTTCCCTACACGACGCTCTTCCGATCT" + lib.tiles[1].fwd
    assert [s.region.label for s in priming_sites(lib, extra=[tailed])] == ["tile1 5' flank"]


def test_a_clean_standard_library_passes():
    rep = _standard().check()
    assert rep.mispriming_issues == []


def test_a_sequence_set_checks_its_adaptors_against_every_member():
    """A SequenceSet has no shared reference (every member is a different gene), so each
    member's coding sequence is checked whole."""
    from library_designer import SequenceSet

    genes = {"g1": CDS, "g2": CDS[3:] + "AAA"}
    spec = LibrarySpec(name="set", adaptor_5="GGCGCGGTCTCAACAG",
                       adaptor_3="TGACTGAGACCGCGCC")
    lib = SequenceSet(spec, genes).generate().codon_optimize()
    assert lib.reference is None
    # The members are codon-optimized, so the run has to be lifted off what they end up as.
    planted = lib.df["variable_dna"].iloc[0][20:40]
    sites = priming_sites(lib, extra=[planted])
    assert {s.region.label for s in sites} == {"member sequences"}
    assert "g1" in {s.where for s in sites}


# --- the opt-in primer screen -------------------------------------------------

def test_screening_passes_over_a_primer_that_would_prime_on_the_cds(tmp_path):
    planted = CDS[33:53]
    path = _pool_csv(tmp_path, [planted])
    drawn = _tiled(tmp_path, primer_set=path).tiles[0].fwd
    assert drawn == planted                                    # off by default

    lib = _tiled(tmp_path, primer_set=path, screen_primers=True)
    assert lib.tiles[0].fwd != planted
    assert planted not in {s for t in lib.tiles for s in (t.fwd, t.rev)}
    assert lib.check().mispriming_issues == []


def test_screening_passes_over_a_primer_that_would_prime_on_the_backbone(tmp_path):
    planted = (BB3 + BB5)[15:35]
    path = _pool_csv(tmp_path, [planted])
    kw = dict(starting_vector=_plasmid(tmp_path), insert_label="insert", primer_set=path)
    assert _tiled(tmp_path, **kw).tiles[0].fwd == planted
    assert _tiled(tmp_path, screen_primers=True, **kw).tiles[0].fwd != planted


def test_screening_keeps_the_primers_of_one_pool_off_each_other(tmp_path):
    """Two primers that anneal to each other cross-amplify the pool, so the second one drawn
    has to be checked against the first."""
    first = load_primer_set("subramanian2018").primers[0][1]
    path = _pool_csv(tmp_path, [first, reverse_complement(first)])
    drawn = {s for t in _tiled(tmp_path, primer_set=path, screen_primers=True).tiles
             for s in (t.fwd, t.rev)}
    assert not (first in drawn and reverse_complement(first) in drawn)


def test_screening_does_not_change_a_design_that_was_already_clean(tmp_path):
    plain = _tiled(tmp_path)
    screened = _tiled(tmp_path, screen_primers=True)
    assert [(t.fwd, t.rev) for t in plain.tiles] == [(t.fwd, t.rev) for t in screened.tiles]
    assert list(plain.df["oligo"]) == list(screened.df["oligo"])


def test_screening_says_so_when_it_runs_out_of_primers(tmp_path):
    path = _pool_csv(tmp_path, [], n_clean=4)          # 4 primers, 2 tiles' worth, 3 needed
    with pytest.raises(ValueError, match="supplies 2 tile"):
        _tiled(tmp_path, primer_set=path, screen_primers=True)


def test_the_screen_is_recorded_in_the_design_specs(tmp_path):
    lib = _tiled(tmp_path, screen_primers=True)
    assert lib.spec.to_dict()["tiled"]["screen_primers"] is True
    assert LibrarySpec(**lib.spec.to_dict()).tiled.screen_primers is True
