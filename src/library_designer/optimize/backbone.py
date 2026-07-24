"""Backbone-and-stamp optimization for single-substitution libraries.

Codon-optimize the wild-type CDS **once** into a frozen reference, then build each
variant by overwriting only its target codon ("stamping"). Every member is then
byte-identical to the reference except at its intended position, the single-WT-
reference invariant that uSort-M mapping and clean single-mutant interpretation
depend on. (Optimizing each variant independently does not guarantee this: with a
usage-matching method, unchanged positions drift apart across variants.)

The stamped codon is chosen by a *local, usage-ranked* pick: try the target
residue's codons from most- to least-frequent and take the first that introduces
no restricted motif (BsaI site / Shine-Dalgarno). The backbone is never re-touched.
"""
from __future__ import annotations

import copy
import re
from functools import lru_cache

import pandas as pd

from ..spec import LibrarySpec
from .codon import codon_optimize

_FREQ_CACHE: dict[str, dict[str, dict[str, float]]] = {}
_RANKED_CACHE: dict[str, dict[str, list[str]]] = {}


def _frequencies(species: str) -> dict[str, dict[str, float]]:
    """Pristine {aa: {codon: freq}} for the species, isolated from DNA Chisel.

    ``python_codon_tables.get_codons_table`` is ``lru_cache``d and returns a shared
    dict that DNA Chisel mutates in place (to log-space) during optimization. We
    clear that cache and keep our own deep copy, so codon frequencies are correct
    regardless of whether optimization has already run.
    """
    if species not in _FREQ_CACHE:
        import python_codon_tables as pct

        pct.get_codons_table.cache_clear()
        _FREQ_CACHE[species] = copy.deepcopy(pct.get_codons_table(species))
    return _FREQ_CACHE[species]


def _cds_from_vector(spec: LibrarySpec) -> str:
    """Extract the CDS already sitting in the starting vector, to freeze verbatim.

    This is the "clone the destination vectors straight from my vector" path: the
    reference *is* the vector's CDS, whatever codon optimization it carries, so the
    per-tile vectors are just the starting vector with one window dropped out. The
    located region must be in-frame and encode the protein (otherwise the locus is
    wrong), but its SD sites, motifs, and even internal BsaI are kept, not recoded.
    They are surfaced as advisories by QC (see checks/tiled.py); an internal BsaI is
    a real assembly hazard and is reported at critical severity there, but it is a
    flag, not an error, because the user chose this sequence.
    """
    from ..layout.vector_io import locating_kwargs, resolve_destination

    dest = resolve_destination(spec.tiled.starting_vector, **locating_kwargs(spec))
    return dest.located_region


def build_reference(spec: LibrarySpec, seed: int | None = None) -> str:
    """The frozen WT CDS every variant is stamped onto.

    Three sources, in order: the CDS already in the starting vector when
    ``tiled.use_vector_cds`` is set (verbatim, advisory QC); an explicit ``spec.cds``
    (verbatim native sequence, e.g. a human gene whose exact codons matter); otherwise
    the protein is codon-optimized once. A verbatim CDS is validated, it must be ACGT,
    in-frame, and encode the (truncated) protein. A ``spec.cds`` is additionally checked
    to be free of the Golden Gate / avoided enzyme sites, since for that path we will
    not silently recode a reference the user chose. Errors are raised, not swallowed,
    because a bad reference poisons every variant.
    """
    from_vector = bool(spec.tiled and spec.tiled.use_vector_cds and spec.tiled.starting_vector)
    if spec.cds is None and not from_vector:
        return codon_optimize(spec.truncated_sequence, spec, seed=seed)

    from dnachisel import translate

    from ..checks.motifs import contains_enzyme_site

    cds = (_cds_from_vector(spec) if from_vector else spec.cds).upper()
    source = "the CDS located in the starting vector" if from_vector else "spec.cds"
    if set(cds) - set("ACGT"):
        raise ValueError(f"{source} must contain only A/C/G/T.")
    if len(cds) % 3 != 0:
        raise ValueError(f"{source} length ({len(cds)}) is not a multiple of 3.")
    protein = spec.truncated_sequence
    got = translate(cds)
    if got != protein:
        raise ValueError(
            f"{source} does not translate to the (truncated) protein_sequence: "
            f"{sum(a != b for a, b in zip(got, protein))} residue(s) differ "
            f"(it encodes {len(got)} aa, protein is {len(protein)} aa)."
        )
    if not from_vector:
        # spec.cds path only: a chosen native CDS with an internal Golden Gate site is
        # an error. The use_vector_cds path keeps such a site and flags it in QC instead.
        hit = next((e for e in spec.avoid_enzymes if contains_enzyme_site(cds, e)), None)
        if hit is not None:
            raise ValueError(
                f"spec.cds contains an internal {hit} site, which would break Golden Gate "
                "assembly. Recode the offending codon(s) in the native CDS yourself, or "
                f"remove {hit} from spec.avoid_enzymes if that is intended."
            )
    return cds


def ranked_codons(species: str) -> dict[str, list[str]]:
    """{amino acid: [codons, most- to least-frequent]} for the species."""
    if species not in _RANKED_CACHE:
        _RANKED_CACHE[species] = {
            aa: [c for c, _ in sorted(freqs.items(), key=lambda kv: kv[1], reverse=True)]
            for aa, freqs in _frequencies(species).items()
        }
    return _RANKED_CACHE[species]


_ADAPT_CACHE: dict[str, dict[str, float]] = {}


def relative_adaptiveness(species: str) -> dict[str, float]:
    """{codon: w}, where w = freq(codon) / max synonymous freq. 1.0 = the optimal
    codon for its amino acid; lower = a rarer synonymous choice. Same usage source
    as the optimizer, so the QC is self-consistent."""
    if species not in _ADAPT_CACHE:
        out: dict[str, float] = {}
        for codons in _frequencies(species).values():
            mx = max(codons.values()) or 1.0
            for codon, freq in codons.items():
                out[codon] = freq / mx
        _ADAPT_CACHE[species] = out
    return _ADAPT_CACHE[species]


def codon_frequency(species: str) -> dict[str, float]:
    """{codon: absolute usage frequency} for the species (pristine, isolated from
    DNA Chisel). Unlike relative adaptiveness, this varies position to position even
    under use_best_codon, so it shows the actual codon-usage landscape."""
    return {c: f for codons in _frequencies(species).values() for c, f in codons.items()}


@lru_cache(maxsize=None)
def _ref_enzyme_count(reference: str, enzyme: str) -> int:
    from ..checks.motifs import count_enzyme_sites

    return count_enzyme_sites(reference, enzyme)


@lru_cache(maxsize=None)
def _ref_pattern_count(reference: str, pattern: str) -> int:
    return len(re.findall(pattern, reference))


def _violates(cds: str, spec: LibrarySpec, reference: str | None = None) -> bool:
    """True if the stamped codon *introduces* a restricted motif or GC-bound break.

    Judged by *count* relative to ``reference``: a match already present in the
    reference is the user's native sequence, not something the stamp created, so only
    a net-new occurrence counts against the variant. (For a codon-optimized reference
    every reference count is 0, so this is the plain absolute check, the original
    behavior.) Windowed GC is approximated by the whole-sequence fraction here.
    """
    from ..checks.motifs import count_enzyme_sites

    ref = reference or ""
    for e in spec.avoid_enzymes:
        if count_enzyme_sites(cds, e) > _ref_enzyme_count(ref, e):
            return True
    for p in spec.avoid_patterns:
        if len(re.findall(p, cds)) > _ref_pattern_count(ref, p):
            return True
    o = spec.optimization
    if o.gc_min is not None or o.gc_max is not None:
        def _out_of_bounds(seq: str) -> bool:
            gc = (seq.count("G") + seq.count("C")) / len(seq) if seq else 0.0
            return (o.gc_min is not None and gc < o.gc_min) or (o.gc_max is not None and gc > o.gc_max)

        if _out_of_bounds(cds) and not (reference and _out_of_bounds(reference)):
            return True
    return False


def _stamp(reference: str, index: int, symbol: str, pinned: str | None,
           ranked: dict[str, list[str]], spec: LibrarySpec):
    """Overwrite the codon at ``index`` on the reference. Returns (cds, None) or
    (None, error). Sense residues try synonymous codons usage-first; a pinned codon
    (e.g. amber TAG) is placed verbatim, if it can't avoid a motif, it's flagged."""
    start = index * 3
    candidates = [pinned] if pinned else ranked.get(symbol)
    if not candidates:
        return None, f"no codons known for residue {symbol!r}"
    for codon in candidates:
        cds = reference[:start] + codon + reference[start + 3:]
        if not _violates(cds, spec, reference):
            return cds, None
    what = f"pinned codon {pinned}" if pinned else f"any synonymous codon for {symbol}"
    return None, f"{what} at position {index + 1} introduces a restricted motif"


def optimize_library(library, seed: int | None = None):
    """Build the frozen WT reference and stamp every variant onto it.

    Returns ``(reference_cds, variable_dna_list, failed)`` aligned to df rows."""
    spec = library.spec
    base = spec.seed if seed is None else seed
    reference = build_reference(spec, seed=base)  # WT, once (native or codon-optimized)
    ranked = ranked_codons(spec.optimization.species)

    df = library.df
    variable: list = []
    failed: dict[str, str] = {}
    for name, symbol, codon, idx in zip(
        df["name"], df["mut_residue"], df["codon"], df["mut_index"]
    ):
        if pd.isna(idx):                       # WT control uses the reference itself
            variable.append(reference)
            continue
        pinned = None if pd.isna(codon) else str(codon)
        cds, err = _stamp(reference, int(idx), str(symbol), pinned, ranked, spec)
        variable.append(cds if cds is not None else pd.NA)
        if err:
            failed[str(name)] = err
    return reference, variable, failed
