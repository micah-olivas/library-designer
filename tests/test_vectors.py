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
    locate_insert,
    read_vector_file,
    resolve_destination,
    terminal_contexts,
)
from library_designer.regions import reverse_complement

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
    ref, o = lib.reference, lib.spec.tiled.overhang_len
    rec, rec_rc = "GGTCTC", reverse_complement("GGTCTC")
    checked = ok = 0
    for _, r in lib.df.iterrows():
        if not isinstance(r["oligo"], str):
            continue
        t = lib.tiles[int(r["tile"])]
        i = r["oligo"].index(rec) + len(rec) + len(lib.spec.tiled.spacer_5)
        j = r["oligo"].rindex(rec_rc) - len(lib.spec.tiled.spacer_3)
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
    # Dropping the CDS back in reproduces the original plasmid (rotated to CDS origin).
    rebuilt = assemble_vector(dest, CLEAN_CDS)
    assert rebuilt == CLEAN_CDS + BB3 + BB5


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
