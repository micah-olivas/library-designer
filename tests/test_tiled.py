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
from library_designer.layout.tiled import compute_tiles, extra_sites
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


# --- end to end ---------------------------------------------------------------

def test_gck_places_all_makeable_mutants(gck):
    mutants = int(gck.df["mut_index"].notna().sum())
    placed = int(gck.df["oligo"].map(lambda o: isinstance(o, str)).sum())
    assert placed == mutants - len(gck.failed)            # every makeable mutant tiled
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
    ref, o = gck.reference, gck.spec.tiled.overhang_len
    p = gck.spec.tiled
    for t in gck.tiles:
        c5 = ref[t.start - o : t.start] if t.start >= o else p.vector_context_5
        c3 = ref[t.end : t.end + o] if t.end + o <= len(ref) else p.vector_context_3
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
