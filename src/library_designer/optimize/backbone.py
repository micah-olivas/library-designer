"""Backbone-and-stamp optimization for single-substitution libraries.

Codon-optimize the wild-type CDS **once** into a frozen reference, then build each
variant by overwriting only its target codon ("stamping").

Three rules govern how a variant's DNA is designed.

1. **Only the mutated codon changes.** No step of the design edits a base outside the
   codon a variant is meant to mutate, including to make room for a motif elsewhere. Every
   member matches the reference apart from its own codon, which is what uSort-M mapping and
   single-mutant interpretation need: an edit elsewhere means a phenotype cannot be
   attributed to the substitution. Optimizing each variant independently does not give
   this, since a usage-matching method chooses different codons at unchanged positions in
   different members. Checked at runtime by ``checks/report.off_target_edits``.
2. **A blocked codon steps down the ranking.** When the preferred codon for the target
   residue would introduce a restricted motif (a Golden Gate site, a Shine-Dalgarno
   motif), try that residue's codons from most- to least-frequent and take the first
   that does not. Only the mutated codon changes.
   ``spec.optimization.synonymous_fallback = False`` opts out of the stepping, for callers
   who would rather be told about the position than build it at a rarer codon.
3. **Exhausting the codons raises a flag.** If no codon for the residue avoids the
   motif, the variant is not built. It is recorded in ``lib.failed`` with the reason,
   its ``variable_dna`` is NA, and QC reports it under ``optimization_failed``. It is not
   placed, and the backbone is not edited to fit it, which rule 1 forbids. A pinned literal
   codon (an amber ``TAG``) has no synonymous alternative, so it reaches this case at
   once.
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
    """Unmodified {aa: {codon: freq}} for the species, kept apart from DNA Chisel's copy.

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
    reference *is* the vector's CDS, whatever codon optimization it carries, so a
    destination vector is just the starting vector with one window dropped out. The
    located region must be in-frame and encode the protein (otherwise the locus is
    wrong), but its SD sites, motifs, and even internal BsaI are kept, not recoded.
    They are surfaced as advisories by QC (see checks/tiled.py); an internal BsaI is
    a real assembly hazard and is reported at critical severity there, as a flag
    rather than an error, since the user chose the sequence.
    """
    from ..layout.vector_io import locating_kwargs, resolve_destination

    dest = resolve_destination(spec.vector.path, **locating_kwargs(spec))
    return dest.located_region


def build_reference(spec: LibrarySpec, seed: int | None = None) -> str:
    """The frozen WT CDS every variant is stamped onto.

    Three sources, in order: the CDS already in the starting vector when
    ``use_vector_cds`` is set (verbatim, advisory QC); an explicit ``spec.cds``
    (verbatim native sequence, e.g. a human gene whose exact codons matter); otherwise
    the protein is codon-optimized once. A verbatim CDS is validated, it must be ACGT,
    in-frame, and encode ``spec.designed_sequence``. A ``spec.cds`` is additionally checked
    to be free of the Golden Gate / avoided enzyme sites, since for that path we will
    not silently recode a reference the user chose. Errors are raised, not swallowed,
    because a bad reference makes every variant wrong.
    """
    vec = spec.vector
    from_vector = bool(vec and vec.use_vector_cds)
    if spec.cds is None and not from_vector:
        return codon_optimize(spec.designed_sequence, spec, seed=seed)

    from dnachisel import translate

    from ..checks.motifs import contains_enzyme_site

    cds = (_cds_from_vector(spec) if from_vector else spec.cds).upper()
    source = "the CDS located in the starting vector" if from_vector else "spec.cds"
    if set(cds) - set("ACGT"):
        raise ValueError(f"{source} must contain only A/C/G/T.")
    if len(cds) % 3 != 0:
        raise ValueError(f"{source} length ({len(cds)}) is not a multiple of 3.")
    protein = spec.designed_sequence
    got = translate(cds)
    if got != protein and spec.truncation and got == spec.protein_sequence:
        # A full-length CDS alongside a truncation. Trimming it is the only reading that makes
        # sense: the truncation is stated once and applies to the protein and the DNA that
        # encodes it, so asking for a pre-trimmed CDS would just be asking the caller to do
        # this by hand. Codons come off the same end the residues did.
        drop = spec.truncation * 3
        cds = cds[drop:] if spec.terminus == "N" else cds[:-drop]
        got = translate(cds)
    if got != protein:
        raise ValueError(
            f"{source} does not translate to {spec.protein_description()}: "
            f"{sum(a != b for a, b in zip(got, protein))} residue(s) differ "
            f"(it encodes {len(got)} aa, protein is {len(protein)} aa). Give either the "
            "designed region or the full-length CDS, which is trimmed to match."
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
    under use_best_codon, so it shows the actual per-position usage."""
    return {c: f for codons in _frequencies(species).values() for c, f in codons.items()}


@lru_cache(maxsize=None)
def _ref_enzyme_count(reference: str, enzyme: str) -> int:
    from ..checks.motifs import count_enzyme_sites

    return count_enzyme_sites(reference, enzyme)


@lru_cache(maxsize=None)
def _ref_pattern_count(reference: str, pattern: str) -> int:
    return len(re.findall(pattern, reference))


def _violates(cds: str, spec: LibrarySpec, reference: str | None = None) -> bool:
    """True if the stamped codon *introduces* a restricted motif, an over-long single-base
    run, or a GC-bound break.

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
    for p in list(spec.avoid_patterns) + spec.homopolymer_patterns:
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
    (None, error).

    A pinned codon (e.g. amber TAG) is placed verbatim, and if it can't avoid a motif it
    is flagged, since there is no synonymous alternative. A sense residue depends on
    ``spec.optimization.synonymous_fallback``: on (the default) it tries that residue's
    codons most- to least-frequent and takes the first that avoids a motif; off, it tries
    the preferred codon only and flags the variant rather than quietly using a rarer one.
    """
    start = index * 3
    candidates = [pinned] if pinned else ranked.get(symbol)
    if not candidates:
        return None, f"no codons known for residue {symbol!r}"
    fallback = spec.optimization.synonymous_fallback
    if not pinned and not fallback:
        candidates = candidates[:1]      # the preferred codon, take it or flag it
    for codon in candidates:
        cds = reference[:start] + codon + reference[start + 3:]
        if not _violates(cds, spec, reference):
            return cds, None
    if pinned:
        what, why = f"pinned codon {pinned}", ""
    elif fallback:
        what, why = f"any synonymous codon for {symbol}", ""
    else:
        what = f"the preferred codon {candidates[0]} for {symbol}"
        why = "; synonymous_fallback is off, so no rarer codon was tried"
    return None, f"{what} at position {index + 1} introduces a restricted motif{why}"


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
