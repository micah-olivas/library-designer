"""Library summary. QC plus the structural metadata you need to review or order a
library at a glance: variant counts per sublibrary, the codon-optimization
parameters used, the adaptor regions, and the selected synthesis platform.
"""
from __future__ import annotations

from dataclasses import asdict, dataclass
from html import escape as _esc

import pandas as pd

from .checks import CheckReport, check_library
from .regions import assemble


@dataclass
class LibrarySummary:
    """A library at a glance, with its QC report attached. Printing it gives the block the
    notebooks show.

    ``per_sublibrary`` groups a scan by mutant residue, with the WT control under ``"WT"``.
    A sequence set has no such grouping, since every member is a distinct sequence, so it
    reports a single ``"members"`` bucket instead. ``adaptors`` carries the two sequences
    and their lengths, plus the construct-length range once at least one member has been
    optimized. ``qc`` is None until ``codon_optimize()`` has run, because there are no
    sequences to check before that.

    Built by ``summarize()``, which is what ``Library.summary()`` calls.
    """

    name: str
    n_variants: int
    n_optimized: int
    n_failed: int
    per_sublibrary: dict[str, int]      # mut_residue -> count (F, Y, *, ..., WT)
    optimization: dict                  # the CodonOptimizationParams used
    adaptors: dict                      # 5'/3' sequences, lengths, construct-length range
    platform: str | None
    qc: CheckReport | None
    starting_vector: str | None = None  # the destination plasmid, when the spec names one

    @property
    def ok(self) -> bool:
        """True when nothing failed optimization and QC passed.

        False before ``codon_optimize()``, since ``qc`` is None and there is no report to
        pass. The ``n_failed`` clause is belt and braces on top of ``qc.passed``, because a
        member whose optimization failed is also listed in ``qc.optimization_failed``.
        """
        return self.n_failed == 0 and self.qc is not None and self.qc.passed

    def __str__(self) -> str:
        lines = [
            f"Library '{self.name}': {self.n_variants} variants "
            f"({self.n_optimized} optimized, {self.n_failed} failed)"
        ]
        if self.platform:
            lines.append(f"  platform: {self.platform}")
        lines.append(
            "  sublibraries: "
            + ", ".join(f"{k}={v}" for k, v in self.per_sublibrary.items())
        )
        a = self.adaptors
        crange = (
            f", construct {a['construct_len_min']}-{a['construct_len_max']} bp"
            if "construct_len_min" in a
            else ""
        )
        lines.append(
            f"  adaptors: 5' {a['5_prime']!r} ({a['len_5']} bp), "
            f"3' {a['3_prime']!r} ({a['len_3']} bp){crange}"
        )
        if self.starting_vector:
            lines.append(f"  starting vector: {self.starting_vector}")
        o = self.optimization
        lines.append(
            f"  optimization: species={o.get('species')}, method={o.get('method')}, "
            f"gc=({o.get('gc_min')}, {o.get('gc_max')}), iters={o.get('max_random_iters')}"
        )
        if self.qc is not None:
            lines.append("  " + str(self.qc).replace("\n", "\n  "))
        return "\n".join(lines)

    def __repr__(self) -> str:   # clean plain output instead of the dataclass dump
        return self.__str__()

    def _repr_html_(self) -> str:
        return f"<pre style='margin:0;line-height:1.4'>{_esc(str(self))}</pre>"


def summarize(library) -> LibrarySummary:
    df, spec = library.df, library.spec
    optimized = "variable_dna" in df.columns
    n_failed = int(df["variable_dna"].isna().sum()) if optimized else 0
    n_opt = (len(df) - n_failed) if optimized else 0

    # Per-sublibrary counts group a scan by its mutant residue. A sequence set has no
    # such grouping (every member is a distinct sequence), so report a single members bucket.
    per_sub: dict[str, int] = {}
    if getattr(library, "kind", "scan") == "sequence_set":
        per_sub["members"] = len(df)
    else:
        # Counts are added into the bucket, not assigned: unmutated rows can carry either
        # NA sentinel (NaN or pd.NA, which value_counts keeps apart), and both are WT.
        for key, val in df["mut_residue"].value_counts(dropna=False).items():
            label = "WT" if pd.isna(key) else str(key)
            per_sub[label] = per_sub.get(label, 0) + int(val)

    adaptors = {
        "5_prime": spec.adaptor_5,
        "3_prime": spec.adaptor_3,
        "len_5": len(spec.adaptor_5),
        "len_3": len(spec.adaptor_3),
    }
    if optimized and n_opt:
        lengths = [
            len(assemble(a5, v, a3))
            for a5, v, a3 in zip(df["adaptor_5"], df["variable_dna"], df["adaptor_3"])
            if isinstance(v, str)
        ]
        adaptors["construct_len_min"] = min(lengths)
        adaptors["construct_len_max"] = max(lengths)

    return LibrarySummary(
        name=spec.name,
        n_variants=len(df),
        n_optimized=n_opt,
        n_failed=n_failed,
        per_sublibrary=per_sub,
        optimization=asdict(spec.optimization),
        adaptors=adaptors,
        platform=spec.platform,
        qc=check_library(library) if optimized else None,
        starting_vector=(
            vec.path if (vec := spec.resolve_vector(getattr(library, "tiled_params", None))) else None
        ),
    )
