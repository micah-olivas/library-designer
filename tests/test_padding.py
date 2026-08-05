"""Tests for evening a tiled oligo pool out to one length (tiled.pad_oligos).

Tiles differ in size, so the oligos do too, and moving the boundaries for better overhangs
widens the spread. Padding closes it by adding filler between each primer and the recognition
site beside it. The pad sits outside what the enzyme releases: it is outside what the enzyme releases, so
the pad is amplified with the oligo and cut away before ligation. The test here is that
the assembled fragment is unchanged with the padding on.
"""
from __future__ import annotations

import pytest

from library_designer import LibrarySpec, SubstitutionScan, TiledAssemblyParams
from library_designer.checks.assembly import _released
from library_designer.checks.motifs import ENZYME_SITES, count_enzyme_sites
from library_designer.layout.tiled import _FILLER, _overhead, pad_lengths, pad_target
from library_designer.regions import reverse_complement
from test_boundaries import _cds, _protein


def _params(**kw) -> TiledAssemblyParams:
    base = dict(oligo_budget=300)
    base.update(kw)
    return TiledAssemblyParams(**base)


def _lib(cds, params):
    spec = LibrarySpec(name="padded", protein_sequence=_protein(cds), cds=cds,
                       substitutions=["A"], tiled=params)
    return SubstitutionScan(spec).generate().codon_optimize().tile()


CDS = _cds(300, seed=21)


def _lengths(lib) -> set[int]:
    return {int(x) for x in lib.df["oligo_length"].dropna()}


# --- the filler itself --------------------------------------------------------

def test_filler_carries_no_recognition_site():
    """A pad next to a recognition site must not extend it, so the filler cannot contain
    one. Doubled, since pads longer than the filler wrap around it."""
    doubled = _FILLER * 2
    for enzyme in ENZYME_SITES:
        assert count_enzyme_sites(doubled, enzyme) == 0


def test_filler_has_no_long_homopolymer():
    runs, run = 1, 1
    for a, b in zip(_FILLER, _FILLER[1:]):
        run = run + 1 if a == b else 1
        runs = max(runs, run)
    assert runs <= 3


# --- lengths ------------------------------------------------------------------

def test_off_by_default():
    p = TiledAssemblyParams()
    assert p.pad_oligos is False and p.pad_target is None


def test_unpadded_pool_holds_more_than_one_length():
    lib = _lib(CDS, _params(optimize_overhangs=True))
    assert len(_lengths(lib)) > 1                     # the problem padding exists to fix


def test_padding_gives_the_pool_a_single_length():
    lib = _lib(CDS, _params(optimize_overhangs=True, pad_oligos=True))
    assert len(_lengths(lib)) == 1


def test_padding_stays_within_the_budget():
    params = _params(optimize_overhangs=True, pad_oligos=True)
    lib = _lib(CDS, params)
    assert max(_lengths(lib)) <= params.oligo_budget
    assert lib.check().oligo_over_budget == []


def test_padded_length_is_the_longest_oligo_the_layout_needed():
    """Padding evens the pool up to the tile that was already longest, so no oligo grows
    beyond what the layout required."""
    params = _params(optimize_overhangs=True, pad_oligos=True)
    lib = _lib(CDS, params)
    windows = [(t.start, t.end) for t in lib.tiles]
    overhead = _overhead(params, len(ENZYME_SITES[params.enzyme]))
    assert _lengths(lib) == {overhead + max(e - s for s, e in windows)}


def test_an_explicit_target_is_honored():
    params = _params(optimize_overhangs=True, pad_oligos=True, pad_target=300)
    assert _lengths(_lib(CDS, params)) == {300}


def test_a_target_shorter_than_the_longest_oligo_is_refused():
    params = _params(pad_oligos=True, pad_target=100)
    with pytest.raises(ValueError, match="cannot be padded down"):
        _lib(CDS, params)


def test_a_target_past_the_budget_is_refused():
    params = _params(pad_oligos=True, pad_target=400, oligo_budget=300)
    with pytest.raises(ValueError, match="past the .* oligo budget"):
        _lib(CDS, params)


def test_pad_lengths_split_evenly_with_the_odd_base_at_the_5_end():
    params = _params()
    overhead = _overhead(params, len(ENZYME_SITES[params.enzyme]))
    assert pad_lengths(100, params, overhead + 105) == (3, 2)
    assert pad_lengths(100, params, overhead + 104) == (2, 2)
    assert pad_lengths(100, params, overhead + 100) == (0, 0)
    assert pad_lengths(100, params, None) == (0, 0)


def test_no_target_when_padding_is_off():
    assert pad_target([(0, 300)], _params()) is None


# --- the pad does not reach the product ---------------------------------------

def test_the_pad_is_cut_away_and_the_fragment_is_unchanged():
    """Why the pad goes between the primer and the site. Same windows
    either way, so any difference in the released fragment would be the pad leaking through.
    """
    plain = _lib(CDS, _params(optimize_overhangs=True))
    padded = _lib(CDS, _params(optimize_overhangs=True, pad_oligos=True))
    assert [(t.start, t.end) for t in plain.tiles] == [(t.start, t.end) for t in padded.tiles]

    by_name = dict(zip(padded.df["name"], padded.df["oligo"]))
    compared = 0
    for name, oligo in zip(plain.df["name"], plain.df["oligo"]):
        other = by_name.get(name)
        if not isinstance(oligo, str) or not isinstance(other, str):
            continue
        a, issues_a = _released(oligo, "BsaI", str(name))
        b, issues_b = _released(other, "BsaI", str(name))
        assert not issues_a and not issues_b
        assert (a.left, a.core, a.right) == (b.left, b.core, b.right)
        compared += 1
    assert compared > 0


def test_a_padded_library_still_passes_qc():
    lib = _lib(CDS, _params(optimize_overhangs=True, pad_oligos=True))
    rep = lib.check()
    assert rep.oligo_extra_sites == []            # no pad completed a recognition site
    assert rep.oligo_over_budget == []
    assert rep.assembly_correct == rep.assembly_checked and rep.assembly_checked > 0


def test_every_oligo_still_carries_exactly_two_sites():
    lib = _lib(CDS, _params(optimize_overhangs=True, pad_oligos=True))
    for oligo in lib.df["oligo"].dropna():
        assert count_enzyme_sites(oligo, "BsaI") == 2


def test_the_pad_sits_between_the_primer_and_the_site():
    lib = _lib(CDS, _params(optimize_overhangs=True, pad_oligos=True))
    rec = ENZYME_SITES["BsaI"].upper()
    for t in lib.tiles:
        oligo = lib.df.loc[lib.df["tile"] == t.index, "oligo"].dropna().iloc[0]
        assert oligo.startswith(t.fwd + t.pad_5 + rec)
        assert oligo.endswith(reverse_complement(rec) + t.pad_3 + reverse_complement(t.rev))


def test_the_wt_controls_are_padded_too():
    """A tile's WT control ships in the same pool as its mutants, so it has to be the same
    length as them."""
    lib = _lib(CDS, _params(optimize_overhangs=True, pad_oligos=True))
    wt = lib.df[lib.df["name"].astype(str).str.startswith("WT_Tile_")]
    assert len(wt) == len(lib.tiles)
    assert {int(x) for x in wt["oligo_length"]} == _lengths(lib)


def test_padding_is_deterministic():
    a = _lib(CDS, _params(optimize_overhangs=True, pad_oligos=True))
    b = _lib(CDS, _params(optimize_overhangs=True, pad_oligos=True))
    assert [(t.pad_5, t.pad_3) for t in a.tiles] == [(t.pad_5, t.pad_3) for t in b.tiles]


def test_padding_off_leaves_the_oligos_byte_identical():
    """The flag has to be inert when it is off, so an existing design is untouched."""
    a = _lib(CDS, _params(optimize_overhangs=True))
    b = _lib(CDS, _params(optimize_overhangs=True, pad_oligos=False))
    assert list(a.df["oligo"].dropna()) == list(b.df["oligo"].dropna())
    assert all(t.pad_5 == "" and t.pad_3 == "" for t in a.tiles)
