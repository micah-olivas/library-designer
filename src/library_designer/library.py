"""The Library container: a pandas DataFrame of variants plus the spec and seed
needed to regenerate it (the design-specs record). Fluent methods make it read naturally in a
notebook: ``SubstitutionScan(spec).generate().codon_optimize()``.
"""
from __future__ import annotations

import json
from importlib import metadata
from pathlib import Path

import pandas as pd

from . import io as _io
from .spec import LibrarySpec, TiledAssemblyParams


def _tool_versions() -> dict:
    out: dict = {}
    for pkg in ("library-designer", "dnachisel"):
        try:
            out[pkg] = metadata.version(pkg)
        except Exception:
            pass
    return out


class Library:
    def __init__(self, df: pd.DataFrame, spec: LibrarySpec, kind: str = "scan"):
        self.df = df.reset_index(drop=True)
        self.spec = spec
        # How the members relate to each other, which sets the optimization strategy:
        #   "scan"         single-mutant variants of one WT, stamped onto a frozen reference
        #   "sequence_set" independent full-length sequences, each optimized on its own
        self.kind = kind
        self.reference: str | None = None     # frozen WT CDS all variants are stamped onto (scan only)
        self.tiles: list | None = None         # tiled assembly: per-tile primers/window/vector (TileInfo)
        self._tiled_params: TiledAssemblyParams | None = None   # what tile() actually laid out with
        self.omega = None                      # OmegaResult from assemble_with_omega(), if run
        self.failed: dict[str, str] = {}       # variant name -> optimization error
        # Design specs: the spec (incl. optimization params) fully determines the
        # library; seed, reference, and tool versions are filled in at optimize time.
        self.design_specs: dict = {"spec": spec.to_dict(), "kind": kind}

    def __len__(self) -> int:
        return len(self.df)

    def __repr__(self) -> str:
        optimized = "variable_dna" in self.df.columns
        return f"<Library {self.spec.name!r}: {len(self)} variants, optimized={optimized}>"

    # --- pipeline -------------------------------------------------------
    def codon_optimize(self, seed: int | None = None) -> "Library":
        """Codon-optimize the library, adding a ``variable_dna`` column.

        The strategy follows ``self.kind``:

        - ``"scan"`` (SubstitutionScan): optimize the WT CDS **once** into a frozen
          reference, then stamp each variant's single target codon onto it. Every
          member is then byte-identical to the reference except at its intended
          position. The stamped codon is chosen usage-first, stepping to a rarer
          synonymous codon only to avoid a restricted motif.
        - ``"sequence_set"`` (SequenceSet): each member is an independent full-length
          protein, so each is codon-optimized on its own (no shared reference).

        A member that can't satisfy the sequence rules is recorded in ``self.failed``
        with ``variable_dna`` set to NA.
        """
        base = self.spec.seed if seed is None else seed

        if self.kind == "sequence_set":
            from .optimize.independent import optimize_independent

            variable, failed = optimize_independent(self, seed=base)
            self.df["variable_dna"] = variable
            self.failed = failed
            self.design_specs.update({"seed": base, "versions": _tool_versions()})
            if failed:
                self.design_specs["failed"] = failed
            return self

        from .optimize.backbone import optimize_library

        reference, variable, failed = optimize_library(self, seed=base)
        self.df["variable_dna"] = variable
        self.reference = reference
        self.failed = failed
        self.design_specs.update(
            {"seed": base, "reference_cds": reference, "versions": _tool_versions()}
        )
        if failed:
            self.design_specs["failed"] = failed
        return self

    def tile(self, params: TiledAssemblyParams | None = None) -> "Library":
        """Lay a (codon-optimized) long-CDS library out for tiled assembly.

        Splits the frozen reference into overlapping tile windows sized to the oligo
        budget, assigns each variant to its tile, picks a tile-specific orthogonal
        primer pair, and assembles the final Golden Gate oligo (adds ``tile``,
        ``oligo``, ``oligo_length`` columns). Per-tile primers and destination vectors
        are stored on ``self.tiles``. ``params`` defaults to ``spec.tiled`` (or built-in
        defaults). Variants with no single containing tile (the WT control) get NA.
        """
        from dataclasses import asdict

        from .layout.tiled import tile_library

        if params is None:
            params = self.spec.tiled or TiledAssemblyParams()
        result = tile_library(self, params)
        # The resolved params, terminal overhangs filled in from the starting vector when
        # there is one. Everything downstream (QC, the vector exporters) reads these off
        # the library rather than spec.tiled, so an explicit tile(params) is honored.
        params = result["params"]
        self._tiled_params = params
        self.df["tile"] = result["tile"]
        self.df["oligo"] = result["oligo"]
        self.df["oligo_length"] = self.df["oligo"].map(
            lambda o: len(o) if isinstance(o, str) else pd.NA
        )
        self.tiles = result["tiles"]
        self._primer_set = result["primer_set"]
        self._unplaced = result["unplaced"]
        self._vector_extra_sites = result["vector_extra_sites"]
        self.design_specs["tiled"] = {
            "params": asdict(params),
            "primer_set": result["primer_set"].name,
            "primers_dropped": result["primer_set"].dropped,
            "n_tiles": len(result["tiles"]),
        }
        return self

    @property
    def tiled_params(self) -> TiledAssemblyParams | None:
        """The tiled-assembly params this library was laid out with, as resolved by
        ``tile()``. QC, the vector exporters, and anything else that needs the layout
        parameters read this instead of ``spec.tiled``, so a library tiled with an
        explicit ``tile(params)`` is judged against what it was actually built with.
        Falls back to ``spec.tiled`` before ``tile()`` has run."""
        return self._tiled_params if self._tiled_params is not None else self.spec.tiled

    def check(self):
        """Run QC (optimization, translation, forbidden sites, length), returning a ``CheckReport``."""
        from .checks import check_library

        return check_library(self)

    def compare_reference(self, other_dna: str, label: str = "IDT"):
        """Cross-check an externally codon-optimized WT CDS (e.g. pasted from IDT's
        Codon Optimization Tool) against this library's frozen reference, returning a
        ``ReferenceComparison``: codon agreement, GC, mean adaptiveness, and any
        BsaI/BsmBI/Shine-Dalgarno sites the external sequence would introduce that our
        design avoids. Submit the truncated protein (``spec.truncated_sequence``) to IDT."""
        from .compare import compare_reference

        return compare_reference(self, other_dna, label=label)

    def summary(self):
        """QC + structural metadata (per-sublibrary counts, params, adaptors,
        platform), returning a ``LibrarySummary``."""
        from .summary import summarize

        return summarize(self)

    def plot_codon_usage(self, metric: str = "frequency", compare: str | None = None,
                         compare_label: str = "IDT"):
        """Return a codon-usage QC Figure (renders inline in a notebook). ``metric`` is
        'frequency' (absolute usage, default) or 'adaptiveness'. Pass ``compare`` an
        external CDS (e.g. IDT's WT optimization) to overlay it as a dashed line."""
        from .viz import codon_usage_figure

        return codon_usage_figure(self, metric=metric, compare=compare, compare_label=compare_label)

    def plot_tiling(self):
        """Return a tile-layout Figure (renders inline in a notebook). Needs a tiled
        library, call ``tile()`` first."""
        from .viz import tiling_figure

        return tiling_figure(self)

    def to_qc_plots(self, path: str | Path) -> "Library":
        """Save the codon-usage QC figure to an image file (format from extension)."""
        from .viz import codon_usage_figure

        fig = codon_usage_figure(self)
        Path(path).parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(path, dpi=150, bbox_inches="tight")
        import matplotlib.pyplot as plt

        plt.close(fig)
        return self

    def drop_failed(self) -> "Library":
        """Remove variants whose optimization failed (``variable_dna`` is NA)."""
        if "variable_dna" in self.df.columns:
            self.df = self.df[self.df["variable_dna"].notna()].reset_index(drop=True)
        self.failed = {}
        return self

    # --- exports (return self so they chain) ----------------------------
    def to_full_csv(self, path: str | Path) -> "Library":
        _io.to_full_csv(self, path)
        return self

    def to_usortm(self, path: str | Path) -> "Library":
        _io.to_usortm(self, path)
        return self

    def to_vendor(self, path: str | Path, method: str | None = None,
                  pool_name: str | None = None) -> "Library":
        _io.to_vendor(self, path, method=method, pool_name=pool_name)
        return self

    def to_oligo_pool(self, path: str | Path) -> "Library":
        _io.to_oligo_pool(self, path)
        return self

    def to_primer_order(self, path: str | Path, scale: str = "25nm",
                        purification: str = "STD") -> "Library":
        _io.to_primer_order(self, path, scale=scale, purification=purification)
        return self

    def to_vectors(self, path: str | Path) -> "Library":
        _io.to_vectors(self, path)
        return self

    def to_vector_maps(self, directory: str | Path) -> "Library":
        """Write one annotated GenBank map per tile, the destination plasmids to clone.

        Each ``.gb`` is the full destination vector (backbone with that tile's window
        dropped out) with the two BsaI sites, the drop-out, the fused overhangs, the
        retained CDS arms, and the tile window annotated, plus a ``destination_vectors.csv``
        manifest. Requires a ``tiled.starting_vector``."""
        _io.to_vector_maps(self, directory)
        return self

    # --- OMEGA (optional, separately-installed external tool) ------------
    def to_omega_fasta(self, path: str | Path) -> "Library":
        """Write each member's codon-optimized CDS as a FASTA for OMEGA gene assembly
        (see ``integrations/omega.py``). Coding region only, no adaptors."""
        from .integrations.omega import write_fasta

        write_fasta(self, path)
        return self

    def assemble_with_omega(self, params, *, omega_home=None, omega_python=None,
                            primer_source: str = "subramanian2018",
                            work_dir=None, config=None):
        """Assemble this (codon-optimized) library into oligopool Golden Gate oligos
        by shelling out to a separately-installed OMEGA (https://github.com/RomeroLab/omega).

        ``params`` is an ``OmegaParams``. Point at your OMEGA checkout with
        ``omega_home`` (or the ``OMEGA_HOME`` env var) and its interpreter with
        ``omega_python`` (or ``OMEGA_PYTHON``). Returns an ``OmegaResult`` (three
        DataFrames) and records the run in the design specs. Unlike the exporters this
        returns the result, not the library, since OMEGA's output is new data.
        """
        from dataclasses import asdict

        from .integrations.omega import assemble

        result = assemble(
            self, params, omega_home=omega_home, omega_python=omega_python,
            primer_source=primer_source, work_dir=work_dir, config=config,
        )
        self.omega = result
        self.design_specs["omega"] = {"params": asdict(params), **result.design_specs()}
        return result

    def to_design_specs(self, path: str | Path) -> "Library":
        """Write the run's design specs, spec (incl. codon-optimization params),
        seed, frozen reference CDS, and tool versions, to JSON."""
        specs = dict(self.design_specs)
        specs["n_variants"] = len(self)
        Path(path).parent.mkdir(parents=True, exist_ok=True)
        Path(path).write_text(json.dumps(specs, indent=2, default=str))
        return self

    def export_all(self, output_dir: str | Path, method: str | None = None,
                   plots: bool = True) -> "Library":
        """Write the full output set into ``output_dir``: master CSV, uSort-M
        ``variants.csv``, vendor order form, the design-specs JSON
        (``<name>_design_specs.json``, the record uSort-M reads), and the
        codon-usage QC plot. The plot is best-effort; if it can't be rendered the
        data files still export (with a warning)."""
        out = Path(output_dir)
        out.mkdir(parents=True, exist_ok=True)
        name = self.spec.name
        self.to_full_csv(out / f"{name}_full_library.csv")
        self.to_usortm(out / "variants.csv")
        self.to_vendor(out / f"{name}_order.csv", method=method)
        self.to_design_specs(out / f"{name}_design_specs.json")
        if plots:
            try:
                self.to_qc_plots(out / f"{name}_codon_usage.png")
            except Exception as exc:
                import warnings
                warnings.warn(
                    f"Skipping QC plot ({exc}).",
                    RuntimeWarning, stacklevel=2,
                )
        return self
