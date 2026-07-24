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
"""
from __future__ import annotations

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


def codon_optimize(protein: str, spec: LibrarySpec, seed: int | None = None) -> str:
    """Reverse-translate ``protein`` and codon-optimize it per ``spec.optimization``,
    subject to the spec's sequence rules. Deterministic given ``seed``. Returns the
    optimized coding DNA (used as the frozen WT reference).
    """
    params = spec.optimization
    if seed is not None:
        np.random.seed(seed)

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
