"""Library summary. QC plus the structural metadata you need to review or order a
library at a glance: variant counts per sublibrary, the codon-optimization
parameters used, the adaptor regions, and the selected synthesis platform.
"""
from __future__ import annotations

from dataclasses import asdict, dataclass
from html import escape as _esc
from pathlib import Path

import pandas as pd

from .checks import CheckReport, check_library
from .regions import assemble
from .checks.report import rows_to_lines
from .spec import optimization_line


def _bp(lo: int, hi: int) -> str:
    """``"298 bp"`` for one length, ``"293-296 bp"`` for a range."""
    return f"{lo} bp" if lo == hi else f"{lo}-{hi} bp"


def _short(seq: str, keep: int = 12) -> str:
    """An adaptor as printed: long ones lose their middle."""
    return seq if len(seq) <= 2 * keep + 3 else f"{seq[:keep]}...{seq[-keep:]}"


@dataclass
class LibrarySummary:
    """A library at a glance, with its QC report attached. Printing it gives the block the
    notebooks show.

    ``per_sublibrary`` groups a scan by mutant residue, with the WT control under ``"WT"``.
    A sequence set has no such grouping, since every member is a distinct sequence, so it
    reports a single ``"members"`` bucket instead. ``adaptors`` carries the two sequences
    and their lengths, plus the construct-length range once at least one member has been
    optimized. ``tiles`` fills in only for a tiled library, where the oligo and not the
    construct is what gets ordered. ``qc`` is None until ``codon_optimize()`` has run,
    because there are no sequences to check before that.

    Printing keeps only the lines that say something. An empty adaptor pair, an unset GC
    window, and a sequence set's single bucket are all left out, and a tiled library shows
    its oligos in place of a construct length.

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
    tiles: dict | None = None           # tiled layout: tile count, oligo lengths, primer set
    gc: dict | None = None              # GC of the ordered molecule: min/median/max + bounds

    @property
    def ok(self) -> bool:
        """True when nothing failed optimization and QC passed.

        False before ``codon_optimize()``, since ``qc`` is None and there is no report to
        pass. The ``n_failed`` clause is redundant with ``qc.passed``, because a
        member whose optimization failed is also listed in ``qc.optimization_failed``.
        """
        return self.n_failed == 0 and self.qc is not None and self.qc.passed

    def _head(self) -> str:
        if self.qc is None:                     # codon_optimize() has not run yet
            state = "not optimized yet"
        elif self.n_failed:
            state = f"{self.n_optimized} optimized, {self.n_failed} failed"
        else:
            state = "all optimized"
        return f"Library '{self.name}': {self.n_variants} variants, {state}"

    def _size_row(self) -> tuple[str, str] | None:
        """How long the thing you order is, as ``(label, value)``: the oligo for a tiled
        library, the construct with its adaptors for everything else. None before
        optimization, when neither is known."""
        a = self.adaptors
        if self.tiles:
            t = self.tiles
            bits = [f"{t['n_tiles']} over a {t['cds_len']} bp CDS"]
            if t.get("oligo_len_min"):
                bits.append(f"oligos {_bp(t['oligo_len_min'], t['oligo_len_max'])}")
            bits.append(f"primer set {t['primer_set']}")
            return "tiles", "   ".join(bits)
        length = (
            _bp(a["construct_len_min"], a["construct_len_max"])
            if "construct_len_min" in a
            else None
        )
        if not a["len_5"] and not a["len_3"]:            # nothing to say about adaptors
            return ("construct", length) if length else None
        pair = (
            f"5' {_short(a['5_prime'])} ({a['len_5']} bp)   "
            f"3' {_short(a['3_prime'])} ({a['len_3']} bp)"
        )
        return "adaptors", f"{pair}   construct {length}" if length else pair

    def _gc_row(self) -> tuple[str, str] | None:
        """The GC spread of the ordered molecule, with the window when one is set. None before
        optimization, when there is nothing to measure."""
        if not self.gc:
            return None
        g = self.gc
        value = f"{g['min']:.3f} to {g['max']:.3f}, median {g['median']:.3f}"
        if g["bounds"]:
            lo, hi = g["bounds"]
            # The closest either end comes to its limit, which is the number worth knowing when
            # the window is wide and the pool is tight. QC names anything actually outside.
            slack = min(g["min"] - lo, hi - g["max"])
            value += (f"   window {lo:.0%}-{hi:.0%}, "
                      + (f"{slack * 100:.1f} points to spare" if slack >= 0 else "exceeded"))
        return "gc", value

    def __str__(self) -> str:
        rows: list[tuple[str, str]] = []
        if self.platform:
            rows.append(("platform", self.platform))
        # A sequence set's one bucket just repeats the variant count, so skip it there.
        if list(self.per_sublibrary) != ["members"]:
            rows.append(("sublibraries",
                         "  ".join(f"{k}={v}" for k, v in self.per_sublibrary.items())))
        if (size := self._size_row()) is not None:
            rows.append(size)
        if (gc_row := self._gc_row()) is not None:
            rows.append(gc_row)
        if self.starting_vector:
            rows.append(("vector", Path(self.starting_vector).name))
        rows.append(("optimization", optimization_line(self.optimization)))

        lines = [self._head()] + rows_to_lines(rows)
        if self.qc is not None:
            lines.append("  " + self.qc.text(count=False).replace("\n", "\n  "))
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

    # A tiled library orders oligos, not constructs, so report the layout instead of the
    # (empty) adaptor pair. Lengths come off the df, which only has them once tiled.
    tiles = None
    if getattr(library, "tiles", None):
        params = library.tiled_params
        oligo_len = df["oligo_length"].dropna() if "oligo_length" in df.columns else pd.Series(dtype=float)
        tiles = {
            "n_tiles": len(library.tiles),
            "cds_len": len(library.reference),
            "primer_set": params.primer_set,
            "oligo_len_min": int(oligo_len.min()) if len(oligo_len) else None,
            "oligo_len_max": int(oligo_len.max()) if len(oligo_len) else None,
        }

    # GC of the molecule that is ordered, which is the number a vendor's window refers to.
    # The spread is what the summary adds: QC already names the members outside the window, so
    # this line says how much room there is rather than repeating the verdict.
    gc = None
    if optimized and n_opt:
        from .checks.report import gc_table

        table = gc_table(library)
        if not table.empty:
            gc = {
                "min": float(table["ordered_gc"].min()),
                "median": float(table["ordered_gc"].median()),
                "max": float(table["ordered_gc"].max()),
                "bounds": spec.gc_bounds,
            }

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
        tiles=tiles,
        gc=gc,
    )
