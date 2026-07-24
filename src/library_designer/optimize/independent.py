"""Per-member codon optimization for sequence-set libraries.

Each member of a ``SequenceSet`` is a distinct full-length protein (an ortholog, a
generative design, a deep multi-mutant), so there is no shared wild-type reference
to stamp onto. Each protein is codon-optimized on its own, under the same sequence
rules (forbidden motifs, restriction sites, GC window) as everything else, using the
single wrapper in ``optimize/codon.py``.

The backbone-and-stamp model in ``optimize/backbone.py`` is deliberately NOT used
here. That model exists to keep single-mutant variants byte-identical outside their
one mutated codon, the invariant uSort-M single-mutant mapping depends on.
Independent members share no such invariant, and forcing one would collapse the set
onto a single sequence. (This is also why the removed per-variant optimizer was wrong
for scans but is exactly right here: unchanged positions "drifting" is meaningless
when every member is a different protein.)
"""
from __future__ import annotations

import pandas as pd

from .codon import codon_optimize


def optimize_independent(library, seed: int | None = None):
    """Codon-optimize each member's protein independently.

    Returns ``(variable_dna_list, failed)`` aligned to the library's df rows. Each
    member gets its own seed offset so a run is reproducible; a protein that can't be
    optimized within the sequence rules is recorded in ``failed`` (``variable_dna``
    NA for that row) rather than aborting the whole library.
    """
    spec = library.spec
    base = spec.seed if seed is None else seed
    df = library.df

    variable: list = []
    failed: dict[str, str] = {}
    for i, (name, protein) in enumerate(zip(df["name"], df["protein"])):
        # A per-member offset keeps members from sharing one stream; seed=None (opted out
        # of reproducibility on the spec) stays None rather than becoming an offset.
        member_seed = None if base is None else base + i
        try:
            variable.append(codon_optimize(str(protein), spec, seed=member_seed))
        except Exception as exc:   # DNA Chisel can't satisfy the constraints for this member
            variable.append(pd.NA)
            failed[str(name)] = f"{type(exc).__name__}: {exc}"
    return variable, failed
