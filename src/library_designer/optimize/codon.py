"""Codon optimization: a single, unified wrapper around DNA Chisel.

Reverse-translate a protein and codon-optimize it under the spec's sequence rules
(forbidden restriction sites and motifs, optional GC window). This produces the
**frozen WT reference** in ``optimize/backbone.py``; per-variant codons are never
set here. Mutations, including pinned literal codons such as an amber ``TAG`` or any
forced codon, are stamped onto that reference by ``backbone._stamp``: it places a
pinned codon verbatim and steps a sense residue down its usage ranking only to avoid
a restricted motif. So this module never pins a codon itself.

All tunable settings live on ``spec.optimization`` (species, method, GC window,
iterations) so they are configurable and recorded in the design specs.

Runs are reproducible by default. DNA Chisel's search is stochastic, so the seed on the
spec (0 unless you change it) is applied on every call, and the same spec and protein
give the same CDS on any machine. Set ``spec.seed = None`` to opt out and draw from the
ambient RNG instead.
"""
from __future__ import annotations

from contextlib import contextmanager

import numpy as np
from dnachisel import (
    CodonOptimize,
    DnaOptimizationProblem,
    EnforceGCContent,
    reverse_translate,
)

from ..spec import CodonOptimizationParams, LibrarySpec
from .constraints import sequence_rules


def _gc_constraint(params: CodonOptimizationParams):
    if params.gc_min is None and params.gc_max is None:
        return None
    kwargs = {
        "mini": params.gc_min if params.gc_min is not None else 0.0,
        "maxi": params.gc_max if params.gc_max is not None else 1.0,
    }
    if params.gc_window:
        kwargs["window"] = params.gc_window
    return EnforceGCContent(**kwargs)


@contextmanager
def _seeded(seed: int | None):
    """Run the block with DNA Chisel's RNG seeded, then hand the caller's RNG back.

    DNA Chisel draws from ``numpy.random``'s global state, so pinning a run means
    seeding that global. Restoring the previous state on the way out keeps our seed from
    reaching into the caller's own draws, which would otherwise silently inherit it.
    ``seed=None`` leaves the RNG alone, so the run follows the ambient stream.
    """
    if seed is None:
        yield
        return
    state = np.random.get_state()
    try:
        np.random.seed(seed)
        yield
    finally:
        np.random.set_state(state)


def codon_optimize(protein: str, spec: LibrarySpec, seed: int | None = None) -> str:
    """Reverse-translate ``protein`` and codon-optimize it per ``spec.optimization``,
    subject to the spec's sequence rules. Returns the optimized coding DNA (used as the
    frozen WT reference).

    Reproducible by default: ``seed`` falls back to ``spec.seed``, so two runs of the
    same spec agree whatever the ambient RNG is doing, and the caller's RNG is left as it
    was found. ``spec.seed = None`` opts out.
    """
    params = spec.optimization
    if seed is None:
        seed = spec.seed          # the spec carries the seed, so honor it here too

    with _seeded(seed):
        dna = reverse_translate(protein, table="Standard")

        constraints = sequence_rules(spec)
        gc = _gc_constraint(params)
        if gc is not None:
            constraints.append(gc)
        objectives = [CodonOptimize(species=params.species, method=params.method)]

        problem = DnaOptimizationProblem(
            sequence=dna,
            constraints=constraints,
            objectives=objectives,
            logger=None,
        )
        problem.max_random_iters = params.max_random_iters
        problem.resolve_constraints()   # satisfy hard rules (translation, avoids, GC)
        problem.optimize()              # then improve codon usage without breaking them
        return problem.sequence
