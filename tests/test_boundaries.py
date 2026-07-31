"""Tests for the tile-boundary search (layout/boundaries.py).

The overhangs of a tiled library are whatever bases sit at the tile boundaries, so moving a
boundary is the only way to change them. These cover that the search obeys every constraint
the balanced split obeys, that it never returns a worse layout than the split it started
from, and that it actually clears a collision when one is there to clear.
"""
from __future__ import annotations

import random

import pytest

from library_designer import LibrarySpec, SubstitutionScan, TiledAssemblyParams
from library_designer.checks.motifs import contains_enzyme_site
from library_designer.layout.boundaries import search_windows, windows_cost
from library_designer.layout.tiled import compute_tiles, max_tile_codons, tile_windows

# Codons with no stop and no ambiguity, spread across the bases so a shifted boundary
# actually lands on different sequence.
_CODONS = ["GCT", "AAA", "CTG", "GAT", "ACC", "TTT", "CAA", "GTT", "AAT", "CCA",
           "TCG", "ATG", "CAC", "GAG", "TGG", "CGT", "AGC", "GGC", "TAC", "ATC"]


def _cds(n_codons: int, seed: int) -> str:
    """A pseudo-random in-frame CDS with no BsaI site, so tiling is not fighting QC."""
    rng = random.Random(seed)
    for _ in range(200):
        cds = "".join(rng.choice(_CODONS) for _ in range(n_codons))
        if not contains_enzyme_site(cds, "BsaI"):
            return cds
    raise AssertionError("could not build a BsaI-free CDS")


def _protein(cds: str) -> str:
    from Bio.Seq import Seq

    return str(Seq(cds).translate())


def _params(**kw) -> TiledAssemblyParams:
    base = dict(oligo_budget=300, optimize_overhangs=True)
    base.update(kw)
    return TiledAssemblyParams(**base)


def _search(cds, params):
    baseline = compute_tiles(len(cds), params)
    return baseline, search_windows(cds, params, baseline, max_tile_codons(params))


def _check_shape(windows, cds, params):
    """Every invariant the balanced split guarantees, re-checked on a searched layout."""
    o = params.overhang_len
    assert windows[0][0] == 0 and windows[-1][1] == len(cds)          # covers the whole CDS
    for (s, e), (nxt_s, _) in zip(windows, windows[1:]):
        assert e == nxt_s                                              # contiguous, no gaps
    for s, e in windows:
        assert s % 3 == 0 and e % 3 == 0                               # codon-aligned
        assert e > s
        assert e - s <= max_tile_codons(params) * 3                    # fits the oligo budget
    if len(windows) > 1:
        # A terminal tile has to be long enough to supply the boundary's overhang.
        assert windows[0][1] - windows[0][0] >= o
        assert windows[-1][1] - windows[-1][0] >= o


# --- the flag -----------------------------------------------------------------

def test_off_by_default():
    assert TiledAssemblyParams().optimize_overhangs is False


def test_flag_off_reproduces_the_balanced_split():
    cds = _cds(300, seed=1)
    params = _params(optimize_overhangs=False)
    assert tile_windows(cds, params) == compute_tiles(len(cds), params)


def test_flag_on_is_allowed_to_move_the_boundaries():
    cds = _cds(300, seed=1)
    params = _params()
    assert tile_windows(cds, params) != compute_tiles(len(cds), params)


# --- constraints --------------------------------------------------------------

@pytest.mark.parametrize("n_codons", [120, 200, 300, 465, 700])
def test_searched_layout_keeps_every_constraint(n_codons):
    cds = _cds(n_codons, seed=n_codons)
    params = _params()
    baseline, found = _search(cds, params)
    assert len(found) == len(baseline)          # same tile count, so the same primers and vectors
    _check_shape(found, cds, params)


@pytest.mark.parametrize("seed", range(8))
def test_search_never_does_worse_than_the_split_it_started_from(seed):
    cds = _cds(280, seed=seed)
    params = _params()
    baseline, found = _search(cds, params)
    assert windows_cost(cds, params, found) <= windows_cost(cds, params, baseline)


def test_search_is_deterministic():
    cds = _cds(300, seed=3)
    params = _params()
    assert _search(cds, params)[1] == _search(cds, params)[1]


def test_a_single_tile_has_no_boundary_to_move():
    cds = _cds(30, seed=4)
    params = _params()
    baseline, found = _search(cds, params)
    assert len(baseline) == 1 and found == baseline


def test_tile_size_is_still_respected_as_a_cap():
    cds = _cds(300, seed=5)
    params = _params(tile_size=120)
    _, found = _search(cds, params)
    assert all(e - s <= 120 for s, e in found)
    _check_shape(found, cds, params)


# --- what it is for -----------------------------------------------------------

def test_the_search_clears_a_collision_the_split_walks_into():
    """The 5' terminal overhang comes from the vector, so pointing it at the bases the
    balanced split leaves at the first boundary makes tile0's two ends identical. The
    boundary can move; the vector context cannot."""
    cds = _cds(300, seed=7)
    off = _params(optimize_overhangs=False)
    first_boundary = compute_tiles(len(cds), off)[0][1]
    collide = cds[first_boundary:first_boundary + 4]

    off = _params(optimize_overhangs=False, vector_context_5=collide)
    on = _params(optimize_overhangs=True, vector_context_5=collide)

    baseline = compute_tiles(len(cds), off)
    found = tile_windows(cds, on)
    assert found != baseline                                   # it had to move something
    assert windows_cost(cds, on, found) < windows_cost(cds, off, baseline)


def test_cross_tile_homology_costs_the_search_nothing():
    """Two tiles are assembled in separate tubes, each into the vector built around its own
    window, so homology between their overhangs is not a hazard. The search must not spend
    the layout's slack on it, or it trades away the comparisons that do matter."""
    from library_designer.layout.boundaries import _pair_cost

    assert _pair_cost((0, "AGGT"), (1, "AGGT")) == 0        # identical, different tiles
    assert _pair_cost((0, "AGGT"), (0, "AGGT")) > 0         # identical, same tile


def test_a_boundary_that_would_spoil_an_oligo_is_avoided():
    """An overhang sitting next to the spacer can spell part of a recognition site, and no
    primer can remove one that the boundary itself creates. The search prices that in."""
    params = _params(spacer_5="A", enzyme="BsaI")
    # Plant "GGTCTC" so that one candidate boundary puts its tail in the 5' overhang. The
    # searched layout must not choose a boundary whose context completes the site.
    cds = _cds(200, seed=11)
    _, found = _search(cds, params)
    o = params.overhang_len
    for s, _e in found[1:]:
        head = "GGTCTC" + params.spacer_5 + cds[s - o:s]
        assert head.count("GGTCTC") == 1                       # only the intended site


# --- through the library ------------------------------------------------------

def _lib(cds, params):
    spec = LibrarySpec(name="searched", protein_sequence=_protein(cds), cds=cds,
                       substitutions=["A"], tiled=params)
    return SubstitutionScan(spec).generate().codon_optimize().tile()


def test_a_searched_library_still_assembles():
    """Moving the boundaries changes every oligo, so the whole pipeline has to still work:
    primers assign, oligos carry two sites, and every member rebuilds its variant."""
    cds = _cds(300, seed=13)
    lib = _lib(cds, _params())
    rep = lib.check()
    assert rep.oligo_extra_sites == [] and rep.oligo_over_budget == []
    assert rep.overhang_issues == []
    assert rep.assembly_correct == rep.assembly_checked and rep.assembly_checked > 0


def test_the_searched_library_has_less_overhang_homology():
    cds = _cds(300, seed=13)
    plain = _lib(cds, _params(optimize_overhangs=False)).overhang_pairs()
    found = _lib(cds, _params()).overhang_pairs()

    def bad(pairs):
        return int((pairs["risk"].isin(["collision", "high"])).sum())

    assert bad(found) <= bad(plain)


def test_the_choice_is_recorded_in_the_design_specs():
    lib = _lib(_cds(300, seed=13), _params())
    tiled = lib.design_specs["tiled"]
    assert tiled["params"]["optimize_overhangs"] is True
    assert tiled["overhangs_optimized"] is True
    assert tiled["overhang_cost"] < tiled["overhang_cost_unsearched"]
