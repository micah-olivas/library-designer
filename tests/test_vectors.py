"""Tests for destination-vector backbones in tiled assembly (layout/vector_io.py).

A synthetic plasmid is written to a temp GenBank file per test, so nothing binary is
committed and the whole path (read -> locate -> splice -> QC -> emit) is exercised with
BioPython.
"""
from __future__ import annotations

from pathlib import Path

import pytest

from Bio.Seq import Seq

from library_designer import LibrarySpec, SubstitutionScan, TiledAssemblyParams
from library_designer.checks.motifs import contains_enzyme_site, count_enzyme_sites
from library_designer.checks.translation import translates_to
from library_designer.layout.vector_io import (
    assemble_vector,
    insert_offset,
    locate_insert,
    read_vector_file,
    resolve_destination,
    terminal_contexts,
)
from library_designer.regions import reverse_complement


def _is_rotation(a: str, b: str) -> bool:
    """True if ``a`` is ``b`` read from a different origin, the only freedom an emitted
    circular map has."""
    return len(a) == len(b) and a in b + b

# A clean, in-frame CDS (no BsaI site, no Shine-Dalgarno motif) and two "dirty" variants.
_SAFE = ["GCT", "AAA", "CTG", "GAT", "ACC", "TTT", "CAA", "GTT", "AAT", "CCA"]
CLEAN_CDS = "".join(_SAFE * 6)                       # 60 codons, 180 bp
assert not contains_enzyme_site(CLEAN_CDS, "BsaI")
BSAI_CDS = "".join(_SAFE * 3) + "GGTCTC" + "".join(_SAFE * 3)   # carries an internal BsaI site
SD_CDS = "".join(_SAFE * 3) + "AGGAGCGCTGCTGCTATG" + "".join(_SAFE * 3)  # carries an SD-like motif

BB5 = "ACGTACGTTTGCAACGGATCCACAG"       # backbone 5' of the insert (BsaI-clean)
BB3 = "TGACCTAGGCATTACGTACGTACGT"       # backbone 3' of the insert (BsaI-clean)

GCK_TOML = Path(__file__).resolve().parents[1] / "examples" / "gck_tiled.toml"

# A synthetic destination plasmid (committed fixture): a circular backbone carrying the
# public hAcyP1 CDS on the minus strand, so the .gb read/locate/splice/emit path runs end
# to end on a real file. Built synthetically (see tests/data/README.md), not a lab construct.
SYNTH_GB = Path(__file__).resolve().parent / "data" / "synthetic_hacyp1_destination.gb"
HACYP1_CDS = (
    "GCTGAAGGAAACACCCTGATTAGCGTTGACTATGAAATCTTTGGCAAAGTGCAGGGCGTTTTCTTCCGCAAACACACCCAG"
    "GCGGAAGGCAAAAAACTGGGCCTGGTTGGCTGGGTGCAGAACACCGATCGCGGCACCGTGCAGGGCCAGCTGCAGGGGCCG"
    "ATCAGCAAAGTGCGTCATATGCAGGAATGGCTGGAAACCCGCGGCTCTCCGAAAAGCCATATTGATAAAGCGAACTTCAAC"
    "AACGAAAAAGTGATTCTGAAACTGGATTACAGCGATTTCCAGATTGTGAAA"
)


def _protein(cds: str) -> str:
    return str(Seq(cds).translate())


def _reconstitution(lib):
    """(checked, ok): excise each variant's tile between the BsaI sites, drop it back into
    the reference, and count how many rebuild the exact stamped mutant CDS."""
    p = lib.tiled_params
    ref, o = lib.reference, p.overhang_len
    rec, rec_rc = "GGTCTC", reverse_complement("GGTCTC")
    checked = ok = 0
    for _, r in lib.df.iterrows():
        if not isinstance(r["oligo"], str):
            continue
        t = lib.tiles[int(r["tile"])]
        i = r["oligo"].index(rec) + len(rec) + len(p.spacer_5)
        j = r["oligo"].rindex(rec_rc) - len(p.spacer_3)
        tile_seq = r["oligo"][i:j][o:-o]
        product = ref[: t.start] + tile_seq + ref[t.end:]
        checked += 1
        ok += product == r["variable_dna"] and translates_to(product, r["protein"])
    return checked, ok


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


def _plasmid(tmp_path, cds, *, topology="circular", reverse=False, extra_feature=None, name="p.gb"):
    """Write a plasmid BB5 + cds + BB3 with a CDS feature; ``reverse`` clones the whole
    thing on the minus strand (so the CDS reads on the reverse strand of the file)."""
    fwd = BB5 + cds + BB3
    cds_start, cds_end = len(BB5), len(BB5) + len(cds)
    feats = [("CDS", cds_start, cds_end, 1, "insert")]
    if extra_feature:
        feats.append(extra_feature)
    if reverse:
        n = len(fwd)
        fwd = reverse_complement(fwd)
        feats = [(t, n - b, n - a, -st, lab) for (t, a, b, st, lab) in feats]
    return _write_gb(tmp_path / name, fwd, topology=topology, features=feats)


# --- reading & locating -------------------------------------------------------

def test_read_reports_topology_and_features(tmp_path):
    path = _plasmid(tmp_path, CLEAN_CDS, extra_feature=("promoter", 0, 10, 1, "prom"))
    seq, topology, features = read_vector_file(path)
    assert topology == "circular"
    assert {f.label for f in features} >= {"insert", "prom"}


def test_locate_forward_and_reverse_strand(tmp_path):
    fwd = _plasmid(tmp_path, CLEAN_CDS, name="fwd.gb")
    seq, _, feats = read_vector_file(fwd)
    s, e, strand = locate_insert(seq, feats, search_cds=CLEAN_CDS)
    assert (seq[s:e], strand) == (CLEAN_CDS, 1)

    rev = _plasmid(tmp_path, CLEAN_CDS, reverse=True, name="rev.gb")
    dest = resolve_destination(rev, search_cds=CLEAN_CDS)
    # After normalization the CDS reads forward and equals the given sequence.
    assert dest.located_region == CLEAN_CDS


def test_missing_cds_raises(tmp_path):
    path = _plasmid(tmp_path, CLEAN_CDS)
    with pytest.raises(ValueError, match="not found"):
        resolve_destination(path, search_cds="ATG" + "AAA" * 20)


# --- flanks, overhangs, splice ------------------------------------------------

def test_terminal_contexts_and_assemble_circular(tmp_path):
    path = _plasmid(tmp_path, CLEAN_CDS)
    dest = resolve_destination(path, search_cds=CLEAN_CDS)
    t5, t3 = terminal_contexts(dest, 4)
    assert t5 == BB5[-4:] and t3 == BB3[:4]           # overhangs are the flanking bases
    # Dropping the CDS back in reproduces the original plasmid, read from the chosen origin.
    rebuilt = assemble_vector(dest, CLEAN_CDS)
    assert _is_rotation(rebuilt, BB5 + CLEAN_CDS + BB3)
    assert rebuilt[insert_offset(dest):][:len(CLEAN_CDS)] == CLEAN_CDS


def test_tile_params_carrying_a_starting_vector_without_a_spec_block(tmp_path):
    """A starting vector passed through ``tile(params)`` rather than ``spec.tiled``.
    Locating the insert, QC, and the map exporter all used to dereference spec.tiled."""
    path = _plasmid(tmp_path, CLEAN_CDS)
    spec = LibrarySpec(name="syn", protein_sequence=_protein(CLEAN_CDS), cds=CLEAN_CDS,
                       substitutions=["A"])
    assert spec.tiled is None
    lib = SubstitutionScan(spec).generate().codon_optimize()
    lib.tile(TiledAssemblyParams(oligo_budget=150, starting_vector=path))

    assert lib.tiled_params.starting_vector == path
    assert all(t.topology == "circular" for t in lib.tiles)
    rep = lib.check()
    assert not rep.vector_extra_sites and not rep.overhang_issues
    checked, ok = _reconstitution(lib)
    assert checked and ok == checked

    lib.to_vector_maps(tmp_path / "maps")
    assert (tmp_path / "maps" / "destination_vectors.csv").is_file()
    assert len(list((tmp_path / "maps").glob("tile*_destination.gb"))) == len(lib.tiles)


# --- end to end: freeze the vector's own CDS ----------------------------------

@pytest.fixture
def vector_lib(tmp_path):
    path = _plasmid(tmp_path, CLEAN_CDS, extra_feature=("rep_origin", 0, 12, 1, "ori"))
    params = TiledAssemblyParams(oligo_budget=150, starting_vector=path,
                                 use_vector_cds=True, insert_label="insert")
    spec = LibrarySpec(name="syn", protein_sequence=_protein(CLEAN_CDS),
                       substitutions=["A"], tiled=params)
    return SubstitutionScan(spec).generate().codon_optimize().tile(), path


def test_vector_cds_becomes_reference_and_reconstitutes(vector_lib):
    lib, _ = vector_lib
    assert lib.reference == CLEAN_CDS
    assert len(lib.tiles) >= 2 and all(t.topology == "circular" for t in lib.tiles)
    checked, ok = _reconstitution(lib)
    assert checked > 0 and ok == checked


def test_vector_cds_qc_clean(vector_lib):
    lib, _ = vector_lib
    rep = lib.check()
    assert rep.passed
    assert not rep.vector_extra_sites and not rep.reference_advisories


# --- stray backbone BsaI: raise when we own the reference, flag when the user does ---

def test_backbone_bsai_raises_without_use_vector_cds(tmp_path):
    # A backbone that carries a BsaI site, with a clean CDS given verbatim (spec.cds).
    path = _write_gb(tmp_path / "dirty.gb",
                     BB5 + CLEAN_CDS + "GGTCTCTT" + BB3,
                     features=[("CDS", len(BB5), len(BB5) + len(CLEAN_CDS), 1, "insert")])
    params = TiledAssemblyParams(oligo_budget=150, starting_vector=path)
    spec = LibrarySpec(name="syn", protein_sequence=_protein(CLEAN_CDS), substitutions=["A"],
                       cds=CLEAN_CDS, tiled=params)
    with pytest.raises(ValueError, match="BsaI"):
        SubstitutionScan(spec).generate().codon_optimize().tile()


def test_reference_bsai_flagged_not_raised_with_use_vector_cds(tmp_path):
    # The CDS the user chose to freeze itself carries a BsaI site: kept, flagged, not raised.
    path = _plasmid(tmp_path, BSAI_CDS, name="bsai.gb")
    params = TiledAssemblyParams(oligo_budget=150, starting_vector=path,
                                 use_vector_cds=True, insert_label="insert")
    spec = LibrarySpec(name="syn", protein_sequence=_protein(BSAI_CDS),
                       substitutions=["A"], tiled=params)
    lib = SubstitutionScan(spec).generate().codon_optimize().tile()   # no raise
    assert lib.reference == BSAI_CDS
    rep = lib.check()
    assert rep.vector_extra_sites                     # critical: surfaced
    assert not rep.passed
    assert any("BsaI" in a for a in rep.reference_advisories)


def test_reference_sd_motif_flagged(tmp_path):
    path = _plasmid(tmp_path, SD_CDS, name="sd.gb")
    params = TiledAssemblyParams(oligo_budget=150, starting_vector=path,
                                 use_vector_cds=True, insert_label="insert")
    spec = LibrarySpec(name="syn", protein_sequence=_protein(SD_CDS),
                       substitutions=["A"], tiled=params)
    lib = SubstitutionScan(spec).generate().codon_optimize().tile()
    rep = lib.check()
    assert any("motif" in a for a in rep.reference_advisories)   # SD flagged, informational


# --- emitting GenBank maps ----------------------------------------------------

def test_to_vector_maps_writes_annotated_gb(vector_lib, tmp_path):
    lib, _ = vector_lib
    outdir = tmp_path / "maps"
    lib.to_vector_maps(outdir)

    import pandas as pd
    manifest = pd.read_csv(outdir / "destination_vectors.csv")
    assert len(manifest) == len(lib.tiles)

    from Bio import SeqIO
    rec = SeqIO.read(str(outdir / "tile0_destination.gb"), "genbank")
    assert rec.annotations["topology"] == "circular"
    labels = [(f.qualifiers.get("label") or [""])[0] for f in rec.features]
    assert any("drop-out" in x for x in labels)
    assert any(x == "ori" for x in labels)            # backbone feature carried over


def test_to_vector_maps_requires_backbone(tmp_path):
    # A cassette-mode tiled library (no starting vector) has no plasmid to annotate.
    spec = LibrarySpec.from_toml(GCK_TOML)
    lib = SubstitutionScan(spec).generate().codon_optimize().tile()
    with pytest.raises(ValueError, match="starting_vector"):
        lib.to_vector_maps(tmp_path / "maps")


# --- committed-plasmid smoke test (a real .gb file, synthetically built) -------

need_gb = pytest.mark.skipif(not SYNTH_GB.exists(), reason=".gb fixture not present")


def _synth_lib(path):
    """The full use_vector_cds pipeline against the committed destination-plasmid file."""
    params = TiledAssemblyParams(oligo_budget=300, starting_vector=str(path),
                                 use_vector_cds=True, insert_label="ACYP1_HUMAN")
    spec = LibrarySpec(name="hacyp1", protein_sequence=_protein(HACYP1_CDS),
                       substitutions=["A"], tiled=params)
    return SubstitutionScan(spec).generate().codon_optimize().tile()


@need_gb
def test_gb_locates_reverse_strand_cds():
    # The hAcyP1 CDS is on the minus strand of the plasmid; search-by-CDS finds and
    # normalizes it so the located region reads forward and equals the given sequence.
    dest = resolve_destination(str(SYNTH_GB), search_cds=HACYP1_CDS)
    assert dest.topology == "circular"
    assert dest.located_region == HACYP1_CDS


@need_gb
def test_gb_end_to_end_and_maps(tmp_path):
    lib = _synth_lib(SYNTH_GB)

    # The plasmid's own CDS becomes the reference, and every placed mutant reconstitutes.
    assert lib.reference == HACYP1_CDS
    assert len(lib.tiles) == 2 and all(t.topology == "circular" for t in lib.tiles)
    checked, ok = _reconstitution(lib)
    assert checked > 0 and ok == checked

    rep = lib.check()
    assert rep.passed and not rep.vector_extra_sites   # this backbone is BsaI-clean

    # Emit the destination plasmids and read one back: full-length, circular, annotated,
    # with the backbone's own features carried over.
    outdir = tmp_path / "maps"
    lib.to_vector_maps(outdir)
    from Bio import SeqIO
    rec = SeqIO.read(str(outdir / "tile0_destination.gb"), "genbank")
    assert rec.annotations["topology"] == "circular"
    assert len(rec.seq) == len(lib.tiles[0].vector)
    labels = [(f.qualifiers.get("label") or [""])[0] for f in rec.features]
    assert any("drop-out" in x for x in labels)
    assert any("msGFP2" in x for x in labels)          # a carried-over backbone feature


# --- standard (untiled) destination vector ------------------------------------
#
# One oligo carries the whole CDS, so the library needs one destination vector. Where its
# drop-out starts and ends follows from the adaptors, and both conventions in the wild are
# covered: the adaptor spelling the backbone-derived overhang, and the adaptor stopping at
# the site plus spacer so the overhang is the CDS's own first/last four bases.

A5_BACKBONE = "GGCGC" + "GGTCTC" + "A" + BB5[-4:]        # ...site, spacer, backbone overhang
A3_BACKBONE = BB3[:4] + "T" + reverse_complement("GGTCTC") + "GCGCC"
A5_CDS = "ggcgcGGTCTCC"                                   # ...site, spacer, then the CDS
A3_CDS = "C" + reverse_complement("GGTCTC") + "cggcg"


def _standard_lib(path, a5, a3, *, cds=CLEAN_CDS, protein=None, **vector_kw):
    from library_designer import StartingVectorParams

    spec = LibrarySpec(
        name="syn", protein_sequence=protein or _protein(cds), cds=cds,
        substitutions=["A"], adaptor_5=a5, adaptor_3=a3,
        starting_vector=StartingVectorParams(path=path, insert_label="insert", **vector_kw),
    )
    return SubstitutionScan(spec).generate().codon_optimize()


def _overhang_blocked(rep) -> list[str]:
    """Members the simulation says cannot ligate because their codon sits in a fused
    overhang, which is the price of drawing the overhang from the CDS ends."""
    return [m for m in rep.assembly_issues if "fused overhang" in m]


def _rebuilt(lib):
    """The plasmid the library's oligos should rebuild, in the emitted vector's frame:
    the starting plasmid with the frozen reference in place."""
    return assemble_vector(lib.destination_vector().dest, lib.reference)


def _clone(lib):
    """Simulate the assembly: digest the oligo, digest the vector, ligate. Returns the
    finished plasmid, or None when the two fused overhangs do not anneal."""
    dv = lib.destination_vector()
    spec, cut = lib.spec, dv.cut
    construct = (spec.adaptor_5 + lib.reference + spec.adaptor_3).upper()
    frag = construct[cut.cut_5:cut.cut_3]
    o5, o3 = len(cut.overhang_5), len(cut.overhang_3)
    if frag[:o5] != dv.overhang_5 or frag[len(frag) - o3:] != dv.overhang_3:
        return None
    stuffer = spec.vector.vector_insert.upper()
    i = dv.sequence.index(stuffer)
    return dv.sequence[:i] + frag[o5:len(frag) - o3] + dv.sequence[i + len(stuffer):]


def test_starting_vector_accepts_a_bare_path(tmp_path):
    """The ergonomic form: spec.starting_vector = "plasmid.gb", normalized to params."""
    path = _plasmid(tmp_path, CLEAN_CDS)
    spec = LibrarySpec(name="syn", protein_sequence=_protein(CLEAN_CDS), starting_vector=path)
    assert spec.vector.path == path
    assert spec.vector.enzyme == "BsaI" and not spec.vector.use_vector_cds
    assert LibrarySpec(name="syn", protein_sequence="MKV").vector is None


def test_starting_vector_from_toml(tmp_path):
    """Both TOML forms: ``starting_vector = "p.gb"`` and a ``[starting_vector]`` table."""
    path = _plasmid(tmp_path, CLEAN_CDS)
    head = f'name = "s"\nprotein_sequence = "{_protein(CLEAN_CDS)}"\nsubstitutions = ["A"]\n'
    (tmp_path / "bare.toml").write_text(head + f'starting_vector = "{path}"\n')
    (tmp_path / "table.toml").write_text(
        head + f'\n[starting_vector]\npath = "{path}"\ninsert_label = "insert"\n'
        'use_vector_cds = true\ninsert_anchors = ["ACAG", "TGAC"]\n'
    )
    bare = LibrarySpec.from_toml(tmp_path / "bare.toml").vector
    assert (bare.path, bare.use_vector_cds) == (path, False)

    table = LibrarySpec.from_toml(tmp_path / "table.toml").vector
    assert table.insert_label == "insert" and table.use_vector_cds
    assert table.insert_anchors == ("ACAG", "TGAC")      # a TOML list, kept as a tuple


def test_backbone_overhang_convention_drops_the_whole_cds(tmp_path):
    lib = _standard_lib(_plasmid(tmp_path, CLEAN_CDS), A5_BACKBONE, A3_BACKBONE)
    dv = lib.destination_vector()

    assert (dv.start, dv.end) == (0, len(lib.reference))       # nothing of the CDS is kept
    assert (dv.overhang_5, dv.overhang_3) == (BB5[-4:], BB3[:4])
    assert dv.cut.keep_5 == dv.cut.keep_3 == 4                 # the adaptors spell the overhangs
    assert lib.check().passed
    # Ligating the digested oligo back in rebuilds the starting plasmid exactly.
    assert _clone(lib) == _rebuilt(lib)
    assert _is_rotation(_rebuilt(lib), BB5 + CLEAN_CDS + BB3)


def test_cds_overhang_convention_keeps_four_bases_of_the_cds(tmp_path):
    """The MBO-038 convention: the adaptor stops after the site and its spacer, so the
    fused overhang is the CDS's own end and the vector has to keep those four bases."""
    lib = _standard_lib(_plasmid(tmp_path, CLEAN_CDS), A5_CDS, A3_CDS)
    dv = lib.destination_vector()

    assert (dv.start, dv.end) == (4, len(lib.reference) - 4)
    assert (dv.overhang_5, dv.overhang_3) == (CLEAN_CDS[:4], CLEAN_CDS[-4:])
    assert dv.cut.keep_5 == dv.cut.keep_3 == 0
    kept = dv.sequence[insert_offset(dv.dest):]
    assert kept.startswith(CLEAN_CDS[:4] + "AGAGACC")          # kept arm, then the drop-out
    rep = lib.check()
    # The price of this convention: the three variants whose codons sit in an overhang
    # cannot ligate. Nothing else is wrong.
    assert rep.assembly_issues == _overhang_blocked(rep)
    assert rep.assembly_correct == rep.assembly_checked - 3
    assert not (rep.adaptor_issues or rep.overhang_issues or rep.vector_extra_sites)
    assert _clone(lib) == _rebuilt(lib)
    assert _is_rotation(_rebuilt(lib), BB5 + CLEAN_CDS + BB3)


def test_linear_backbone_splices_in_place(tmp_path):
    lib = _standard_lib(_plasmid(tmp_path, CLEAN_CDS, topology="linear"), A5_CDS, A3_CDS)
    dv = lib.destination_vector()
    assert dv.topology == "linear"
    assert _clone(lib) == BB5 + CLEAN_CDS + BB3


def test_codon_optimized_reference_arms_come_from_the_reference(tmp_path):
    """With no spec.cds the reference is codon-optimized, so it differs from the CDS the
    plasmid holds. The kept arms and overhangs must follow the reference, not the plasmid,
    or the oligos would not anneal."""
    lib = _standard_lib(_plasmid(tmp_path, CLEAN_CDS), A5_CDS, A3_CDS,
                        cds=None, protein=_protein(CLEAN_CDS))
    dv = lib.destination_vector()
    assert lib.reference != CLEAN_CDS
    assert dv.overhang_5 == lib.reference[:4] and dv.overhang_3 == lib.reference[-4:]
    rep = lib.check()
    assert rep.assembly_issues == _overhang_blocked(rep)       # only the terminal codons
    assert _clone(lib) == _rebuilt(lib)
    assert _is_rotation(_rebuilt(lib), BB5 + lib.reference + BB3)


def test_missing_cut_site_in_an_adaptor_is_reported(tmp_path):
    """The 3' adaptor shipped with the MBO-038 example spells the reverse of the BsaI site
    rather than its reverse complement, so nothing cuts that end."""
    lib = _standard_lib(_plasmid(tmp_path, CLEAN_CDS), A5_CDS, "CCTCTGGcggcg")
    rep = lib.check()
    assert not rep.passed
    assert any("adaptor_3" in m and "no reverse-strand BsaI site" in m for m in rep.adaptor_issues)


def test_overhang_that_will_not_ligate_is_reported(tmp_path):
    """Four adaptor bases that are not the backbone bases flanking the insert. The construct
    itself is fine, so only a check against the plasmid can catch this."""
    lib = _standard_lib(_plasmid(tmp_path, CLEAN_CDS),
                        "GGCGC" + "GGTCTC" + "A" + "TTTT", A3_BACKBONE)
    rep = lib.check()
    assert not rep.passed
    assert any("TTTT" in m and BB5[-4:] in m for m in rep.adaptor_issues)
    assert _clone(lib) is None                       # the overhangs do not anneal


def test_cut_falling_inside_the_coding_region_is_reported(tmp_path):
    # The site sits one base too close to the variable region, so the cut would shave the
    # first base off the CDS.
    lib = _standard_lib(_plasmid(tmp_path, CLEAN_CDS), "GGCGCGGTCTC", A3_CDS)
    rep = lib.check()
    assert not rep.passed
    assert any("inside the coding region" in m for m in rep.adaptor_issues)


def test_no_adaptors_is_an_advisory_not_a_failure(tmp_path):
    """A pool with no cloning flanks still gets a destination vector (whole CDS dropped
    out, overhangs from the backbone), and QC says so without failing the library."""
    lib = _standard_lib(_plasmid(tmp_path, CLEAN_CDS), "", "")
    rep = lib.check()
    assert rep.passed and not rep.adaptor_issues
    assert any("no adaptors are set" in a for a in rep.reference_advisories)
    dv = lib.destination_vector()
    assert dv.cut is None and (dv.start, dv.end) == (0, len(lib.reference))


def test_backbone_bsai_is_reported_on_a_standard_vector(tmp_path):
    path = _write_gb(tmp_path / "dirty.gb", BB5 + CLEAN_CDS + "GGTCTCTT" + BB3,
                     features=[("CDS", len(BB5), len(BB5) + len(CLEAN_CDS), 1, "insert")])
    rep = _standard_lib(path, A5_CDS, A3_CDS).check()
    assert rep.vector_extra_sites and not rep.passed


def test_use_vector_cds_without_tiling(tmp_path):
    """Freeze the plasmid's own CDS as the reference for a plain scan, no tiling involved."""
    lib = _standard_lib(_plasmid(tmp_path, CLEAN_CDS), A5_CDS, A3_CDS,
                        cds=None, protein=_protein(CLEAN_CDS), use_vector_cds=True)
    assert lib.reference == CLEAN_CDS
    rep = lib.check()
    assert rep.assembly_issues == _overhang_blocked(rep)       # only the terminal codons
    assert _clone(lib) == _rebuilt(lib)


def test_spec_level_vector_feeds_a_tiled_layout(tmp_path):
    """A [tiled] block with no backbone of its own picks up spec.starting_vector, so the
    plasmid is declared once wherever the library is laid out."""
    path = _plasmid(tmp_path, CLEAN_CDS)
    spec = LibrarySpec(name="syn", protein_sequence=_protein(CLEAN_CDS), cds=CLEAN_CDS,
                       substitutions=["A"], starting_vector=path,
                       tiled=TiledAssemblyParams(oligo_budget=150))
    lib = SubstitutionScan(spec).generate().codon_optimize().tile()

    assert lib.tiled_params.starting_vector == path      # resolved onto the layout
    assert len(lib.tiles) >= 2 and all(t.topology == "circular" for t in lib.tiles)
    checked, ok = _reconstitution(lib)
    assert checked and ok == checked


def test_standard_exports_one_vector_and_one_map(tmp_path):
    import pandas as pd

    lib = _standard_lib(_plasmid(tmp_path, CLEAN_CDS, extra_feature=("rep_origin", 0, 12, 1, "ori")),
                        A5_CDS, A3_CDS)
    dv = lib.destination_vector()

    lib.to_vectors(tmp_path / "vectors.csv")
    rows = pd.read_csv(tmp_path / "vectors.csv")
    assert len(rows) == 1
    assert rows.loc[0, "vector_sequence"] == dv.sequence
    assert (rows.loc[0, "overhang_5"], rows.loc[0, "overhang_3"]) == (dv.overhang_5, dv.overhang_3)
    assert (rows.loc[0, "cds_dropout_start"], rows.loc[0, "cds_dropout_end"]) == (dv.start, dv.end)
    assert rows.loc[0, "origin_in_starting_vector"] == dv.dest.origin + 1

    outdir = tmp_path / "maps"
    lib.to_vector_maps(outdir)
    assert not list(outdir.glob("tile*_destination.gb"))
    manifest = pd.read_csv(outdir / "destination_vectors.csv")
    assert len(manifest) == 1 and manifest.loc[0, "file"] == "destination.gb"

    from Bio import SeqIO
    rec = SeqIO.read(str(outdir / "destination.gb"), "genbank")
    assert str(rec.seq) == dv.sequence
    assert rec.annotations["topology"] == "circular"
    labels = [(f.qualifiers.get("label") or [""])[0] for f in rec.features]
    assert any("drop-out" in x for x in labels)
    assert "ori" in labels                               # backbone feature carried over


def test_to_vectors_needs_a_starting_vector(tmp_path):
    spec = LibrarySpec(name="syn", protein_sequence=_protein(CLEAN_CDS), cds=CLEAN_CDS,
                       substitutions=["A"])
    lib = SubstitutionScan(spec).generate().codon_optimize()
    with pytest.raises(ValueError, match="starting_vector"):
        lib.to_vectors(tmp_path / "vectors.csv")
    with pytest.raises(ValueError, match="starting_vector"):
        lib.to_vector_maps(tmp_path / "maps")


# --- where the emitted map starts ---------------------------------------------
#
# A viewer draws a plasmid from base 1, so whatever sits at the origin is split across the
# two ends of the map. Reading from the insert put the Golden Gate sites there, which is
# the part you want to look at, so the origin goes upstream of the promoter.

ORI = "".join("ACGTTGCA"[i % 8] for i in range(200))
PROM = "TTGACAATTAATCATCGGCTCGTATAATGTGTGGA"          # 35 bp, a trc-like promoter
LEADER = "GGCATTACGTACGTACGTACGTACGTACGT"             # 30 bp between promoter and CDS
TERM = "AACCCGCTGATCGGCACGTAAGAGGTTCC"                # 29 bp, ends at the plasmid's last base


def _cassette_plasmid(tmp_path, *, gap: int = 60, name="cassette.gb"):
    """A plasmid laid out like a real expression construct: origin of replication, a gap,
    a promoter, a leader, the insert, a terminator running to the last base. ``gap`` sets
    how much room there is upstream of the promoter."""
    spacer = "A" * gap
    seq = ORI + spacer + PROM + LEADER + CLEAN_CDS + TERM
    p5 = len(ORI) + gap
    c5 = p5 + len(PROM) + len(LEADER)
    feats = [
        ("rep_origin", 0, len(ORI), 1, "ori"),
        ("promoter", p5, p5 + len(PROM), 1, "trc promoter"),
        ("CDS", c5, c5 + len(CLEAN_CDS), 1, "insert"),
        ("terminator", c5 + len(CLEAN_CDS), len(seq), 1, "rrnB T1"),
    ]
    return _write_gb(tmp_path / name, seq, features=feats), p5, c5


def test_origin_sits_upstream_of_the_promoter(tmp_path):
    path, prom_start, cds_start = _cassette_plasmid(tmp_path)
    dest = resolve_destination(path, search_cds=CLEAN_CDS)

    assert dest.origin == prom_start - 50          # 50 bp of room before the promoter
    assert insert_offset(dest) == cds_start - dest.origin
    # The promoter and the insert follow the origin in reading order, both intact.
    rebuilt = assemble_vector(dest, CLEAN_CDS)
    assert rebuilt.index(PROM) < rebuilt.index(CLEAN_CDS)


def test_origin_is_nudged_off_a_feature_it_would_cut(tmp_path):
    """Backing off 50 bp can land inside an upstream annotation. A feature spanning the
    origin cannot be written to GenBank at all, so the origin moves to its edge instead."""
    path, prom_start, _ = _cassette_plasmid(tmp_path, gap=20, name="tight.gb")
    dest = resolve_destination(path, search_cds=CLEAN_CDS)

    assert prom_start - 50 < len(ORI)              # the plain 50 bp back-off would cut the ori
    assert dest.origin == 0                        # so it settles on that feature's edge
    assert dest.origin != dest.start               # and still not on the insert


def test_a_linear_backbone_is_never_rotated(tmp_path):
    path, _, cds_start = _cassette_plasmid(tmp_path, name="linear.gb")
    dest = resolve_destination(path, topology_override="linear", search_cds=CLEAN_CDS)
    assert dest.origin == dest.start == cds_start
    assert assemble_vector(dest, CLEAN_CDS) == dest.full_seq


def test_map_keeps_every_backbone_feature_and_records_the_rotation(tmp_path):
    """The emitted map is rotated against the file the user gave us, so the manifest says
    which base it starts at. Every backbone feature survives the remap, including one that
    ends on the plasmid's very last base."""
    import pandas as pd

    path, _, _ = _cassette_plasmid(tmp_path)
    lib = _standard_lib(path, A5_CDS, A3_CDS)
    dv = lib.destination_vector()
    lib.to_vector_maps(tmp_path / "maps")

    manifest = pd.read_csv(tmp_path / "maps" / "destination_vectors.csv")
    assert manifest.loc[0, "origin_in_starting_vector"] == dv.dest.origin + 1
    assert manifest.loc[0, "features_carried"] == 3          # ori, promoter, terminator

    from Bio import SeqIO
    rec = SeqIO.read(str(tmp_path / "maps" / "destination.gb"), "genbank")
    labels = [(f.qualifiers.get("label") or [""])[0] for f in rec.features]
    assert {"ori", "trc promoter", "rrnB T1"} <= set(labels)
    # The drop-out and its enzyme sites sit inside the map, not across its ends.
    sites = [i for i in range(len(rec.seq) - 5)
             if str(rec.seq)[i:i + 6] in ("GGTCTC", reverse_complement("GGTCTC"))]
    assert len(sites) == 2
    assert all(20 < i < len(rec.seq) - 20 for i in sites)


def test_edited_plasmid_file_is_read_again_not_served_from_cache(tmp_path):
    """resolve_destination caches per path, so an edited plasmid used to come back
    stale for the rest of the session (a notebook re-running a cell)."""
    path = Path(_plasmid(tmp_path, CLEAN_CDS))
    first = resolve_destination(str(path), search_cds=CLEAN_CDS)
    assert first.located_region == CLEAN_CDS

    other = "".join(_SAFE * 5) + "GCTAAACTG"          # a different in-frame CDS
    assert other != CLEAN_CDS
    _plasmid(tmp_path, other, name="p.gb")            # same filename, new contents
    again = resolve_destination(str(path), search_cds=other)
    assert again.located_region == other
    # and the cache still serves an unchanged file without re-reading
    assert resolve_destination(str(path), search_cds=other) is again


# --- what the destination vector says about itself ----------------------------

def test_the_destination_vector_prints_a_readable_summary(tmp_path):
    """Evaluating it in a notebook should say everything the library turns on, without the
    caller reaching into dv.dest and dv.cut to print it themselves."""
    lib = _standard_lib(_plasmid(tmp_path, CLEAN_CDS), A5_BACKBONE, A3_BACKBONE)
    dv = lib.destination_vector()
    text = str(dv)

    assert repr(dv) == text                                   # same in a cell and in print()
    html = dv._repr_html_()                                   # escaped, so 5' becomes 5&#x27;
    assert html.startswith("<pre") and "Destination vector for p.gb" in html
    assert "p.gb" in text.splitlines()[0]                     # named after the plasmid file
    assert f"{len(dv.dest.full_seq)} bp circular" in text
    assert f"{dv.dest.start}-{dv.dest.end}" in text
    assert f"{dv.start}-{dv.end} of the {len(lib.reference)} bp" in text
    assert f"5' {dv.overhang_5}   3' {dv.overhang_3}" in text
    assert "(match)" in text and dv.overhangs_match
    assert f"{dv.cut.keep_5} at the 5' end" in text
    assert f"{dv.length} bp, read from base {dv.dest.origin + 1}" in text


def test_the_summary_calls_out_overhangs_that_will_not_ligate(tmp_path):
    lib = _standard_lib(_plasmid(tmp_path, CLEAN_CDS), "GGCGC" + "GGTCTC" + "A" + "TTTT",
                        A3_BACKBONE)
    dv = lib.destination_vector()
    assert dv.overhangs_match is False
    assert "MISMATCH, nothing will ligate" in str(dv)


def test_the_summary_survives_adaptors_with_no_cut_site(tmp_path):
    """dv.cut is None here, so anything reaching for dv.cut.overhang_5 would raise."""
    dv = _standard_lib(_plasmid(tmp_path, CLEAN_CDS), "", "").destination_vector()
    assert dv.overhangs_match is None
    text = str(dv)
    assert "the adaptors carry no cut site" in text
    assert "adaptor bases kept" not in text
    assert sum(line.strip().startswith("!") for line in text.splitlines()) == 2


def test_the_summary_flags_a_minus_strand_insert(tmp_path):
    lib = _standard_lib(_plasmid(tmp_path, CLEAN_CDS, reverse=True),
                        "GGCGC" + "GGTCTC" + "A" + BB5[-4:], A3_BACKBONE)
    assert "insert on the minus strand" in str(lib.destination_vector())


# --- backbone padding: sliding the fused overhang off the junction ---------------
#
# An adaptor may spell a few bases past its fused overhang. Those bases are backbone, the
# oligo carries them, so the destination vector gives them up and the overhang window slides
# into the backbone. This is the only direction open to a saturating scan: an overhang drawn
# from coding bases breaks the variants that mutate them.

# One base of padding at each end. The pad bases are the backbone's own, so the product is
# still the parent with the CDS swapped, but the overhang window has moved by one.
A5_PAD = "GGCGC" + "GGTCTC" + "A" + BB5[-5:]              # keep_5 = 5, overhang = BB5[-5:-1]
A3_PAD = BB3[:5] + "T" + reverse_complement("GGTCTC") + "GCGCC"   # keep_3 = 5, overhang = BB3[1:5]


def test_padding_is_read_off_the_adaptors(tmp_path):
    from library_designer.layout.destination import adaptor_padding

    lib = _standard_lib(_plasmid(tmp_path, CLEAN_CDS), A5_PAD, A3_PAD)
    assert adaptor_padding(lib.spec) == (1, 1)
    # The unpadded convention asks for nothing.
    plain = _standard_lib(_plasmid(tmp_path, CLEAN_CDS, name="q.gb"), A5_BACKBONE, A3_BACKBONE)
    assert adaptor_padding(plain.spec) == (0, 0)


def test_padding_slides_the_overhang_window_into_the_backbone(tmp_path):
    padded = _standard_lib(_plasmid(tmp_path, CLEAN_CDS), A5_PAD, A3_PAD)
    plain = _standard_lib(_plasmid(tmp_path, CLEAN_CDS, name="q.gb"), A5_BACKBONE, A3_BACKBONE)

    dv, dv0 = padded.destination_vector(), plain.destination_vector()
    assert (dv0.overhang_5, dv0.overhang_3) == (BB5[-4:], BB3[:4])       # at the junction
    assert (dv.overhang_5, dv.overhang_3) == (BB5[-5:-1], BB3[1:5])      # one base along
    # And the oligo carries what the vector presents, which is the point of moving both.
    assert dv.overhangs_match
    # The vector gave up the padding bases, so its drop-out is 2 bases wider than the CDS.
    assert dv.dest.pad_5 == 1 and dv.dest.pad_3 == 1
    assert dv.dest.end - dv.dest.start == len(padded.reference) + 2


def test_a_padded_library_assembles_back_to_the_parent(tmp_path):
    """The oligo re-supplies the bases the vector dropped, so every clone is the starting
    plasmid with its own CDS in it."""
    lib = _standard_lib(_plasmid(tmp_path, CLEAN_CDS), A5_PAD, A3_PAD)
    rep = lib.check()
    assert rep.passed, rep.issues
    assert rep.assembly_correct == rep.assembly_checked == len(lib.df)
    assert rep.assembly_aligned == rep.assembly_checked
    # The simulated WT clone is the parent plasmid.
    assert lib.assembled_product("WT") == lib.parent_vector()


def test_padding_that_does_not_match_the_plasmid_is_refused(tmp_path):
    """The padding bases end up in the product, so they have to be the plasmid's own. A
    mismatch would change the backbone, so it is raised naming both sequences."""
    wrong = "A" + BB3[1:5] + "T" + reverse_complement("GGTCTC") + "GCGCC"
    assert wrong[0] != BB3[0]                      # the pad base disagrees, the overhang does not
    with pytest.raises(ValueError, match="past its fused overhang"):
        _standard_lib(_plasmid(tmp_path, CLEAN_CDS), A5_BACKBONE, wrong).destination_vector()


def test_padding_leaves_the_reference_cds_alone(tmp_path):
    """Widening the locus must not reach the reference: use_vector_cds still extracts the
    coding region, not the coding region plus padding."""
    from library_designer import StartingVectorParams

    path = _plasmid(tmp_path, CLEAN_CDS)
    spec = LibrarySpec(
        name="syn", protein_sequence=_protein(CLEAN_CDS), substitutions=["A"],
        adaptor_5=A5_PAD, adaptor_3=A3_PAD,
        starting_vector=StartingVectorParams(path=path, insert_label="insert",
                                             use_vector_cds=True),
    )
    lib = SubstitutionScan(spec).generate().codon_optimize()
    assert lib.reference == CLEAN_CDS               # not CLEAN_CDS plus a backbone base
    assert translates_to(lib.reference, _protein(CLEAN_CDS))


# --- a truncation holds residues back for the vector to supply -------------------
#
# Truncating means the library does not encode those residues, not that the construct loses
# them. The plasmid keeps their codons and supplies them, so the assembled clone still
# encodes the whole protein_sequence and the fused overhang moves inside the retained codons.

def _held_out_lib(path, truncation, a5, **kw):
    """A library that holds ``truncation`` residues back for the vector to supply. The
    truncation is a spec field, not a vector one, so this cannot go through _standard_lib."""
    from library_designer import StartingVectorParams

    spec = LibrarySpec(
        name="held", protein_sequence=_protein(CLEAN_CDS), cds=CLEAN_CDS,
        substitutions=["A"], adaptor_5=a5, adaptor_3=A3_BACKBONE, truncation=truncation,
        starting_vector=StartingVectorParams(path=path, insert_label="insert"), **kw,
    )
    return SubstitutionScan(spec).generate().codon_optimize()


def test_the_locus_is_the_designed_region_not_the_whole_cds(tmp_path):
    from library_designer.layout.destination import resolve_insert_locus

    path = _plasmid(tmp_path, CLEAN_CDS)
    full = resolve_insert_locus(_standard_lib(path, A5_BACKBONE, A3_BACKBONE).spec)
    assert full.end - full.start == len(CLEAN_CDS)

    # With two residues held back off the N terminus the vector drops six fewer bases, and the
    # six it keeps are the codons for those residues.
    lib = _held_out_lib(path, 2, A5_BACKBONE)
    d = resolve_insert_locus(lib.spec)
    assert d.end - d.start == len(CLEAN_CDS) - 6 == len(lib.reference)
    assert d.full_seq[d.start - 6:d.start] == CLEAN_CDS[:6]
    assert d.start == full.start + 6 and d.end == full.end


def test_the_vector_presents_the_overhang_from_inside_the_retained_codons(tmp_path):
    from library_designer.layout.destination import build_destination

    lib = _held_out_lib(_plasmid(tmp_path, CLEAN_CDS), 2, A5_BACKBONE)
    dv = build_destination(lib, strict=False)
    # Not the backbone bases any more: the junction moved into the CDS.
    assert dv.overhang_5 == CLEAN_CDS[2:6]
    assert dv.overhang_5 != BB5[-4:]
    # So the untruncated adaptor no longer matches, and QC says so naming both sequences.
    assert dv.overhangs_match is False
    (msg,) = lib.check().issues["adaptor_issues"]
    assert CLEAN_CDS[2:6] in msg and BB5[-4:] in msg


def test_a_held_out_library_assembles_to_the_whole_protein(tmp_path):
    """The check that matters: the clone is the starting plasmid with the designed region
    swapped, and it still encodes every residue, including the held-out ones."""
    from dnachisel import translate

    from library_designer.layout.destination import resolve_insert_locus
    from library_designer.layout.vector_io import insert_offset

    # An adaptor carrying the overhang the shortened locus presents.
    a5 = "GGCGC" + "GGTCTC" + "A" + CLEAN_CDS[2:6]
    lib = _held_out_lib(_plasmid(tmp_path, CLEAN_CDS), 2, a5)

    rep = lib.check()
    assert rep.passed, rep.issues
    assert rep.assembly_aligned == rep.assembly_checked == len(lib.df)

    parent, wt = lib.parent_vector(), lib.assembled_product("WT")
    assert parent == wt                                  # the WT clone is the parent
    d = resolve_insert_locus(lib.spec)
    at = insert_offset(d) + d.pad_5
    assert wt[at - 6:at] == CLEAN_CDS[:6]                # supplied by the vector
    assert translate(wt[at - 6:at + len(lib.reference)]) == _protein(CLEAN_CDS)


def test_holding_residues_back_off_the_c_terminus_works_the_same_way(tmp_path):
    from library_designer.layout.destination import resolve_insert_locus

    lib = _held_out_lib(_plasmid(tmp_path, CLEAN_CDS), 2, A5_BACKBONE,
                        truncation_terminus="C")
    d = resolve_insert_locus(lib.spec)
    assert d.end - d.start == len(CLEAN_CDS) - 6
    assert d.full_seq[d.end:d.end + 6] == CLEAN_CDS[-6:]   # kept at the other end
    assert lib.reference == CLEAN_CDS[:-6]


def test_the_parent_baseline_is_the_real_plasmid_not_a_shortened_one(tmp_path):
    """A regression. Locating the whole CDS while the oligo carried only the designed region
    dropped the held-out codons out of the vector with nothing putting them back, and the
    parent was rebuilt short by the same amount, so every clone matched it and QC passed."""
    path = _plasmid(tmp_path, CLEAN_CDS)
    plain = _standard_lib(path, A5_BACKBONE, A3_BACKBONE)
    a5 = "GGCGC" + "GGTCTC" + "A" + CLEAN_CDS[2:6]
    held = _held_out_lib(path, 2, a5)
    assert len(held.parent_vector()) == len(plain.parent_vector())


def test_the_destination_map_annotates_the_held_out_codons(tmp_path):
    """Opening the map should show what the vector contributes. The held-out codons get a CDS
    feature, and the fused overhang sits inside them."""
    from Bio import SeqIO

    a5 = "GGCGC" + "GGTCTC" + "A" + CLEAN_CDS[2:6]
    lib = _held_out_lib(_plasmid(tmp_path, CLEAN_CDS), 2, a5)
    lib.to_vector_maps(tmp_path / "vector")
    rec = SeqIO.read(tmp_path / "vector" / "destination.gb", "genbank")

    (held,) = [f for f in rec.features if "held out" in f.qualifiers.get("label", [""])[0]]
    assert held.type == "CDS"
    assert str(held.extract(rec.seq)) == CLEAN_CDS[:6]        # the two codons themselves
    assert "N terminus" in held.qualifiers["label"][0]

    (o5,) = [f for f in rec.features if f.qualifiers.get("label") == ["fused overhang 5'"]]
    assert str(o5.extract(rec.seq)) == CLEAN_CDS[2:6]
    # The overhang is drawn from inside the held-out codons, so it sits within that feature.
    assert int(held.location.start) <= int(o5.location.start)
    assert int(o5.location.end) <= int(held.location.end)


def test_an_untruncated_map_has_no_held_out_feature(tmp_path):
    from Bio import SeqIO

    lib = _standard_lib(_plasmid(tmp_path, CLEAN_CDS), A5_BACKBONE, A3_BACKBONE)
    lib.to_vector_maps(tmp_path / "vector")
    rec = SeqIO.read(tmp_path / "vector" / "destination.gb", "genbank")
    assert not [f for f in rec.features if "held out" in f.qualifiers.get("label", [""])[0]]


def test_a_clone_map_spans_the_whole_cds_including_the_held_out_codons(tmp_path):
    """The clone encodes every residue, so its CDS feature has to cover the held-out codons
    too. Annotating the reference alone would leave them outside the CDS on every map."""
    from Bio import SeqIO
    from dnachisel import translate

    a5 = "GGCGC" + "GGTCTC" + "A" + CLEAN_CDS[2:6]
    lib = _held_out_lib(_plasmid(tmp_path, CLEAN_CDS), 2, a5)
    lib.to_assembled_vectors(tmp_path / "clones", fmt="genbank")

    name = str(lib.df["name"].iloc[0])
    rec = SeqIO.read(tmp_path / "clones" / f"{name}.gb", "genbank")
    cds = next(f for f in rec.features
               if f.type == "CDS" and f.qualifiers.get("label") == [lib.spec.name])
    assert len(cds.extract(rec.seq)) == len(lib.reference) + 6
    assert len(translate(str(cds.extract(rec.seq)))) == len(lib.spec.protein_sequence)

    (marked,) = [f for f in rec.features if "from the vector" in f.qualifiers.get("label", [""])[0]]
    assert str(marked.extract(rec.seq)) == CLEAN_CDS[:6]
    assert int(marked.location.start) == int(cds.location.start)   # at the 5' end of the CDS
