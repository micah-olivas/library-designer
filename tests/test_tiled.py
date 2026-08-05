"""Tests for tiled-assembly layout and the frozen-reference invariant it relies on.

The glucokinase example (native human CDS, alanine scan) is the end-to-end fixture;
it runs fast because a native CDS skips codon optimization.
"""
from __future__ import annotations

from pathlib import Path

import pandas as pd
import pytest

from library_designer import LibrarySpec, SubstitutionScan, TiledAssemblyParams
from library_designer.checks.motifs import ENZYME_SITES, contains_enzyme_site
from library_designer.checks.translation import translates_to
from library_designer.layout.tiled import compute_tiles, extra_sites, tile_contexts
from library_designer.primers import load_primer_set
from library_designer.regions import reverse_complement

REPO = Path(__file__).resolve().parents[1]
GCK_TOML = REPO / "examples" / "gck_tiled.toml"


@pytest.fixture(scope="module")
def gck():
    spec = LibrarySpec.from_toml(GCK_TOML)
    return SubstitutionScan(spec).generate().codon_optimize().tile()


# --- primer sets --------------------------------------------------------------

def test_bundled_sets_load_and_screen():
    sub = load_primer_set("subramanian2018", enzyme="BsaI")
    assert sub.kind == "pool"
    assert len(sub.primers) == 164          # 165 kept, 1 dropped for a BsaI site
    assert sub.dropped == ["SUB151"]
    assert all(not contains_enzyme_site(seq, "BsaI") for _, seq in sub.primers)

    gck = load_primer_set("gck800", enzyme="BsaI")
    assert gck.kind == "paired" and len(gck.pairs) == 7


def test_unknown_primer_set_raises():
    with pytest.raises(FileNotFoundError):
        load_primer_set("does_not_exist")


# --- tiling geometry ----------------------------------------------------------

def test_compute_tiles_cover_and_codon_aligned():
    params = TiledAssemblyParams()
    for cds_len in (300, 999, 1395, 3003):
        tiles = compute_tiles(cds_len, params)
        assert tiles[0][0] == 0 and tiles[-1][1] == cds_len          # spans the CDS
        assert all(s % 3 == 0 and e % 3 == 0 for s, e in tiles)      # codon-aligned
        for (a, b), (c, d) in zip(tiles, tiles[1:]):
            assert c == b                                            # gap-free, contiguous
        rec_len = len(ENZYME_SITES[params.enzyme])
        overhead = 2 * 20 + 2 * rec_len + 2 + 2 * params.overhang_len
        assert all((e - s) + overhead <= params.oligo_budget for s, e in tiles)


def test_non_inframe_reference_rejected():
    with pytest.raises(ValueError):
        compute_tiles(100, TiledAssemblyParams())


def test_tiles_split_evenly_and_never_leave_a_stub():
    """Every tile boundary needs a full overhang of CDS on both sides, so no terminal
    tile may be shorter than overhang_len. An even split is what keeps that true."""
    params = TiledAssemblyParams()
    o = params.overhang_len
    for codons in range(2, 8000):
        tiles = compute_tiles(codons * 3, params)
        sizes = [e - s for s, e in tiles]
        assert max(sizes) - min(sizes) <= 3          # balanced to within one codon
        if len(tiles) > 1:
            assert min(sizes[0], sizes[-1]) >= o     # both ends can carry an overhang


def test_internal_boundaries_take_cds_bases_not_the_vector_context():
    """The tile before a short final tile used to be handed the vector's terminal
    overhang, which its destination vector does not present at that junction."""
    params = TiledAssemblyParams(tile_size=9, oligo_budget=1000)
    ref = "ATG" + "AAACCCGGGTTTACGTGCATCA"[:18]     # 21 bp, 7 codons
    tiles = compute_tiles(len(ref), params)
    for start, end in tiles:
        ctx5, ctx3 = tile_contexts(ref, start, end, params)
        assert ctx5 == (params.vector_context_5 if start == 0 else ref[start - 4:start])
        assert ctx3 == (params.vector_context_3 if end == len(ref) else ref[end:end + 4])
        assert len(ctx5) == 4 and len(ctx3) == 4     # a full overhang either way


def test_tile_size_too_small_to_carry_an_overhang_is_refused():
    # One codon per tile: 3 bp cannot supply the 4 bp overhang the next tile needs.
    with pytest.raises(ValueError, match="overhang"):
        compute_tiles(21, TiledAssemblyParams(tile_size=3, oligo_budget=1000))
    # An even split rescues what a ceil-sized window would have left as a 1-codon stub:
    # 7 codons in 3-codon tiles is 3 + 2 + 2, not 3 + 3 + 1.
    assert compute_tiles(21, TiledAssemblyParams(tile_size=9, oligo_budget=1000)) == [
        (0, 9), (9, 15), (15, 21)
    ]


def test_single_tile_needs_no_cds_overhang():
    """One tile spanning the whole CDS takes both overhangs from the vector, so the
    terminal-length rule does not apply to it."""
    params = TiledAssemblyParams(oligo_budget=1000)
    assert compute_tiles(3, params) == [(0, 3)]
    ctx5, ctx3 = tile_contexts("ATG", 0, 3, params)
    assert (ctx5, ctx3) == (params.vector_context_5, params.vector_context_3)


# --- end to end ---------------------------------------------------------------

def test_gck_places_all_makeable_mutants(gck):
    mutants = int(gck.df["mut_index"].notna().sum())
    placed = int(gck.df["oligo"].map(lambda o: isinstance(o, str)).sum())
    # Every makeable mutant is tiled, plus one WT control oligo per tile.
    assert placed == mutants - len(gck.failed) + len(gck.tiles)
    assert [n for n in gck._unplaced if n != "WT"] == []  # nothing unplaced but the WT control


def test_gck_only_failures_are_reported_motif_constraints(gck):
    # A handful of positions may be unmakeable because every synonymous codon would
    # introduce a forbidden motif, those are captured in .failed, not silently placed.
    rep = gck.check()
    assert set(rep.optimization_failed) == set(gck.failed)
    assert all("restricted motif" in msg for msg in gck.failed.values())
    # everything that *was* placed is clean
    assert not rep.translation_fail
    assert not any(rep.enzyme_hits.values()) and not any(rep.motif_hits.values())
    assert not rep.oligo_extra_sites and not rep.oligo_over_budget
    assert not rep.overhang_issues and not rep.unplaced


def test_gck_oligos_within_budget_and_no_extra_sites(gck):
    budget = gck.spec.tiled.oligo_budget
    for oligo, ti in zip(gck.df["oligo"], gck.df["tile"]):
        if not isinstance(oligo, str):
            continue
        t = gck.tiles[int(ti)]
        assert len(oligo) <= budget
        assert not extra_sites(oligo, len(t.fwd), len(t.rev), "BsaI")


def test_golden_gate_reconstitutes_full_mutant_cds(gck):
    """Excise ctx5+tile+ctx3 between the BsaI sites, share the 4-bp overhangs, drop
    into the tile's vector, which must equal the stamped full mutant CDS and translate to
    the target protein."""
    ref = gck.reference
    o = gck.spec.tiled.overhang_len
    rec, rec_rc = ENZYME_SITES["BsaI"], reverse_complement(ENZYME_SITES["BsaI"])
    checked = 0
    for _, r in gck.df.iterrows():
        if not isinstance(r["oligo"], str):
            continue
        t = gck.tiles[int(r["tile"])]
        i = r["oligo"].index(rec) + len(rec) + len(gck.spec.tiled.spacer_5)
        j = r["oligo"].rindex(rec_rc) - len(gck.spec.tiled.spacer_3)
        tile_seq = r["oligo"][i:j][o:-o]
        product = ref[: t.start] + tile_seq + ref[t.end :]
        assert product == r["variable_dna"]
        assert translates_to(product, r["protein"])
        checked += 1
    assert checked > 400


def test_single_wt_reference_invariant(gck):
    ref = gck.reference
    for dna, idx in zip(gck.df["variable_dna"], gck.df["mut_index"]):
        if pd.isna(dna) or pd.isna(idx):
            continue
        codon = range(int(idx) * 3, int(idx) * 3 + 3)
        diffs = [k for k in range(len(ref)) if dna[k] != ref[k]]
        assert all(k in codon for k in diffs)   # differs only within the stamped codon


def test_tile_overhangs_distinct_and_non_palindromic(gck):
    ref = gck.reference
    for t in gck.tiles:
        c5, c3 = tile_contexts(ref, t.start, t.end, gck.tiled_params)
        assert c5 != c3
        assert c5 != reverse_complement(c5) and c3 != reverse_complement(c3)


# --- the params a library was actually tiled with ------------------------------

SHORT_PROTEIN = "MKAILVGADEQTRWYFNSHCP" * 6


def test_qc_runs_when_the_params_came_from_tile_not_the_spec():
    """``tile(params)`` on a spec with no ``[tiled]`` block. QC and the summary must read
    the params off the library; reaching for ``spec.tiled`` used to raise AttributeError."""
    spec = LibrarySpec(name="p", protein_sequence=SHORT_PROTEIN, substitutions=["A"])
    assert spec.tiled is None
    lib = SubstitutionScan(spec).generate().codon_optimize()
    lib.tile(TiledAssemblyParams(oligo_budget=300))

    assert lib.tiled_params.oligo_budget == 300
    rep = lib.check()
    assert not rep.oligo_over_budget and not rep.overhang_issues
    assert lib.summary().qc is not None


def test_qc_judges_the_params_used_not_the_spec_block():
    """A stale/tight ``spec.tiled`` budget must not condemn a layout built with another."""
    spec = LibrarySpec(name="p", protein_sequence=SHORT_PROTEIN, substitutions=["A"],
                       tiled=TiledAssemblyParams(oligo_budget=200))
    lib = SubstitutionScan(spec).generate().codon_optimize()
    lib.tile(TiledAssemblyParams(oligo_budget=1000))          # laid out generously

    assert lib.tiled_params.oligo_budget == 1000
    assert int(lib.df["oligo_length"].dropna().max()) > 200   # legitimately over spec.tiled
    assert not lib.check().oligo_over_budget                  # judged at the budget used
    assert lib.design_specs["tiled"]["params"]["oligo_budget"] == 1000


# --- per-tile WT controls ------------------------------------------------------

def test_one_wt_control_oligo_per_tile(gck):
    """Each tile is amplified out of the pool and assembled on its own, so each one needs
    its own unmutated member. The global ``WT`` row stays in the table for the record and rides on
    no oligo."""
    controls = [f"WT_Tile_{t.index}" for t in gck.tiles]
    assert [n for n in gck.df["name"] if str(n).startswith("WT_Tile_")] == controls

    wt = gck.df[gck.df["name"].isin(controls)]
    assert list(wt["tile"]) == [t.index for t in gck.tiles]
    assert all(isinstance(o, str) for o in wt["oligo"])
    assert wt["mut_index"].isna().all()                       # no mutation to report
    assert (wt["variable_dna"] == gck.reference).all()        # the whole frozen reference
    assert gck.df.loc[gck.df["name"] == "WT", "oligo"].isna().all()


def test_wt_control_oligo_carries_the_reference_window(gck):
    """A control oligo is its tile's mutant oligo with the WT window in place: the same
    primers, sites, and overhangs, so it amplifies and assembles with the rest of that
    sublibrary."""
    from library_designer.layout.tiled import assemble_oligo

    for t in gck.tiles:
        oligo = gck.df.loc[gck.df["name"] == f"WT_Tile_{t.index}", "oligo"].iloc[0]
        assert oligo == assemble_oligo(gck.reference, gck.reference, t.start, t.end,
                                       t.fwd, t.rev, gck.tiled_params)
        assert oligo.startswith(t.fwd) and oligo.endswith(reverse_complement(t.rev))
        assert not extra_sites(oligo, len(t.fwd), len(t.rev), "BsaI")
        assert len(oligo) <= gck.tiled_params.oligo_budget


def test_wt_controls_reach_the_pooled_order(gck, tmp_path):
    out = tmp_path / "oligos.csv"
    gck.to_oligo_pool(out)
    pool = pd.read_csv(out)
    assert {f"WT_Tile_{t.index}" for t in gck.tiles} <= set(pool["name"])
    assert "WT" not in set(pool["name"])          # the global row has no oligo to order


def test_wt_controls_can_be_switched_off():
    spec = LibrarySpec(name="p", protein_sequence=SHORT_PROTEIN, substitutions=["A"],
                       tiled=TiledAssemblyParams(oligo_budget=300, wt_controls=False))
    lib = SubstitutionScan(spec).generate().codon_optimize().tile()
    assert not [n for n in lib.df["name"] if str(n).startswith("WT_Tile_")]
    assert lib.design_specs["tiled"]["wt_controls"] == []


def test_retiling_rebuilds_the_controls_instead_of_stacking_them():
    """A control belongs to one layout (a tile index and that tile's primers), so tiling
    again replaces the set. Left in place they would be tiled a second time as if they were
    variants, and the pool would grow a duplicate WT member per tile."""
    spec = LibrarySpec(name="p", protein_sequence=SHORT_PROTEIN, substitutions=["A"],
                       tiled=TiledAssemblyParams(oligo_budget=300))
    lib = SubstitutionScan(spec).generate().codon_optimize().tile()
    n_rows, first = len(lib.df), len(lib.tiles)

    lib.tile(TiledAssemblyParams(oligo_budget=1000))          # fewer, longer tiles
    controls = [n for n in lib.df["name"] if str(n).startswith("WT_Tile_")]
    assert len(controls) == len(lib.tiles) < first
    assert len(lib.df) == n_rows - first + len(lib.tiles)
    assert not lib.df["name"].duplicated().any()


def test_optimizing_again_drops_the_controls():
    """New sequences mean a new layout, so the controls go with the oligos they carried."""
    spec = LibrarySpec(name="p", protein_sequence=SHORT_PROTEIN, substitutions=["A"],
                       tiled=TiledAssemblyParams(oligo_budget=300))
    lib = SubstitutionScan(spec).generate().codon_optimize().tile()
    assert [n for n in lib.df["name"] if str(n).startswith("WT_Tile_")]

    lib.codon_optimize()
    assert not [n for n in lib.df["name"] if str(n).startswith("WT_Tile_")]
    assert "oligo" not in lib.df.columns
    assert lib.df["variable_dna"].notna().any()                # rows still line up


def test_primer_set_longer_than_primer_length_refused(tmp_path):
    """primer_length sizes the tiles, so a set of longer primers must be refused rather
    than silently push every oligo past the budget."""
    csv = tmp_path / "long_primers.csv"
    csv.write_text("primer_id,sequence\n" + "".join(f"P{i},{'ACGT' * 7}\n" for i in range(8)))
    spec = LibrarySpec(name="p", protein_sequence=SHORT_PROTEIN, substitutions=["A"],
                       tiled=TiledAssemblyParams(primer_set=str(csv), primer_length=20))
    lib = SubstitutionScan(spec).generate().codon_optimize()
    with pytest.raises(ValueError, match="primer_length"):
        lib.tile()


# --- native-CDS validation ----------------------------------------------------

def test_native_cds_with_internal_bsai_site_raises():
    # ATG GGT CTC AAA -> "MGLK"; the CDS carries GGTCTC (BsaI)
    spec = LibrarySpec(name="x", protein_sequence="MGLK", cds="ATGGGTCTCAAA",
                       substitutions=["A"], tiled=TiledAssemblyParams())
    assert contains_enzyme_site(spec.cds, "BsaI")
    with pytest.raises(ValueError, match="BsaI"):
        SubstitutionScan(spec).generate().codon_optimize()


def test_native_cds_translation_mismatch_raises():
    spec = LibrarySpec(name="x", protein_sequence="MGLK", cds="ATGGGCCTGAAG",  # -> MGLK? checked below
                       substitutions=["A"])
    # deliberately claim the wrong protein
    spec.protein_sequence = "MGLW"
    with pytest.raises(ValueError, match="translate"):
        SubstitutionScan(spec).generate().codon_optimize()


# --- exporters ----------------------------------------------------------------

def test_exporters_write_expected(gck, tmp_path):
    gck.to_oligo_pool(tmp_path / "oligos.csv")
    gck.to_primer_order(tmp_path / "primers.csv")
    gck.to_vectors(tmp_path / "vectors.csv")

    pool = pd.read_csv(tmp_path / "oligos.csv")
    assert list(pool.columns) == ["name", "sequence"]
    assert len(pool) == gck.df["oligo"].map(lambda o: isinstance(o, str)).sum()

    primers = pd.read_csv(tmp_path / "primers.csv", header=None)
    assert len(primers) == 2 * len(gck.tiles)      # fwd + rev per tile

    vectors = pd.read_csv(tmp_path / "vectors.csv")
    assert len(vectors) == len(gck.tiles)
    assert all(gck.spec.tiled.vector_insert in v for v in vectors["vector_sequence"])


def test_export_requires_tiling():
    spec = LibrarySpec.from_toml(GCK_TOML)
    lib = SubstitutionScan(spec).generate().codon_optimize()   # not tiled
    with pytest.raises(ValueError, match="tile"):
        lib.to_oligo_pool("/tmp/should_not_write.csv")


def test_usortm_export_rejects_tiled(gck, tmp_path):
    # A tiled pool has no single uniform variable-region block, so the uSort-M export
    # must fail loudly rather than emit one full-CDS record per variant (fix #2).
    with pytest.raises(NotImplementedError, match="tiled"):
        gck.to_usortm(tmp_path / "should_not_write.csv")
    assert not (tmp_path / "should_not_write.csv").exists()


def test_viz_figures(gck):
    import matplotlib.pyplot as plt
    for fig in (gck.plot_tiling(), gck.plot_codon_usage()):
        assert type(fig).__name__ == "Figure"
        plt.close(fig)


# --- palindromic recognition sites --------------------------------------------

def test_extra_sites_handles_a_palindromic_recognition_site():
    """A palindromic site reads the same on both strands, so both intended positions
    appear in both strand lists. Comparing per-strand lists flagged a clean oligo."""
    from library_designer.checks import motifs
    from library_designer.layout import tiled as layout

    site = "GGCGCC"                                  # its own reverse complement
    assert reverse_complement(site) == site
    motifs.ENZYME_SITES["PalI"] = site
    try:
        fwd, rev = "A" * 20, "T" * 20
        clean = (fwd + site + "A" + "TTTT" + "ATGATGATG" + "GGAG" + "T"
                 + reverse_complement(site) + reverse_complement(rev))
        assert not layout.extra_sites(clean, len(fwd), len(rev), "PalI")
        # a third copy in the middle is a real extra site and must still be caught
        dirty = clean[:40] + site + clean[40:]
        assert layout.extra_sites(dirty, len(fwd), len(rev), "PalI")
    finally:
        del motifs.ENZYME_SITES["PalI"]


def test_a_tiled_library_writes_one_fasta_per_oligo(gck, tmp_path):
    """The per-oligo FASTAs carry what the pool carries: the assembled oligo with its
    primers and sites, not the bare variable region. The global WT rides on no oligo."""
    from library_designer.io import _file_stem

    # Straight to the exporter, not export_all: this fixture keeps its one unmakeable
    # variant, and an order refuses to go out incomplete.
    d = tmp_path / "oligos"
    gck.to_oligo_files(d, fmt="fasta")

    placed = {str(n): o for n, o in zip(gck.df["name"], gck.df["oligo"])
              if isinstance(o, str)}
    assert {p.name for p in d.iterdir()} == {f"{_file_stem(n)}.fasta" for n in placed}
    assert not (d / "WT.fasta").exists()
    assert any(n.startswith("WT_Tile_") for n in placed)   # per-tile controls do ship

    for name, oligo in placed.items():
        head, seq, _ = (d / f"{_file_stem(name)}.fasta").read_text().split("\n")
        assert head == f">{name}"
        assert seq == oligo.upper()
