"""The Library container: a pandas DataFrame of variants plus the spec and seed
needed to regenerate it (the design-specs record). Methods return the library, so a run chains:
``SubstitutionScan(spec).generate().codon_optimize()``.
"""
from __future__ import annotations

import json
import re
from datetime import datetime
from importlib import metadata
from pathlib import Path

import pandas as pd

from . import io as _io
from .spec import LibrarySpec, TiledAssemblyParams


_RUN_STAMP = re.compile(r"\d{8}_\d{6}")     # the dated suffix run_directory() appends


def _specs_filename(name: str) -> str:
    """The design-specs filename for a library, in one place so writers and readers agree."""
    return f"{name}_design_specs.json"


def _tool_versions() -> dict:
    out: dict = {}
    for pkg in ("library-designer", "dnachisel"):
        try:
            out[pkg] = metadata.version(pkg)
        except Exception:
            pass
    return out


class Library:
    """A generated library, made of the variant table plus the spec behind it and the
    state derived from both.

    You get one from a generator's ``generate()`` rather than building it yourself, then
    chain the pipeline methods, which return the library, so a run reads as one line::

        lib = SubstitutionScan(spec).generate().codon_optimize()

    ``lib.df`` is the variant table, one row per member. ``lib.reference`` is the frozen
    WT CDS a scan stamps onto, and ``lib.tiles`` is the per-tile layout once ``tile()``
    has run. ``lib.failed`` maps each member whose optimization could not satisfy the
    sequence rules to its error. ``lib.design_specs`` accumulates the run record that
    ``to_design_specs()`` writes.
    """

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
        self._wt_control_names: list[str] = []  # per-tile WT rows tile() appended, so it can undo them
        self.omega = None                      # OmegaResult from assemble_with_omega(), if run
        self.assembly: list | None = None      # AssemblyResults from simulate_assembly(), if run
        self.failed: dict[str, str] = {}       # variant name -> optimization error
        self._created_at: datetime | None = None   # when the sequences were built
        self.output_dir: Path | None = None    # the run directory export_all() wrote to
        # Design specs: the spec (incl. optimization params) fully determines the
        # library; seed, reference, and tool versions are filled in at optimize time.
        self.design_specs: dict = {"spec": spec.to_dict(), "kind": kind}

    def __len__(self) -> int:
        return len(self.df)

    def __repr__(self) -> str:
        optimized = "variable_dna" in self.df.columns
        where = f", out={self.output_dir}" if self.output_dir is not None else ""
        return f"<Library {self.spec.name!r}: {len(self)} variants, optimized={optimized}{where}>"

    # --- provenance -----------------------------------------------------
    @property
    def created(self) -> str | None:
        """When this library's sequences were built, ISO 8601 with the local offset, or
        None before ``codon_optimize()``. Recorded in the design specs and used to name the
        run directory, so a file on disk and the record inside it agree."""
        return self._created_at.isoformat(timespec="seconds") if self._created_at else None

    @property
    def run_id(self) -> str | None:
        """One token naming this run: ``<name>_<stamp>``, the run directory's own name.

        This is the identity a downstream tool carries forward. It is derived from the
        moment the sequences were built, so it is stable across re-exports of one library
        and different for the next run, and it survives being copied out of the directory
        because it also names the library. None before ``codon_optimize()``."""
        from .io import run_stamp

        if self._created_at is None:
            return None
        return f"{self.spec.name}_{run_stamp(self._created_at)}"

    def runs(self, base: str | Path = "out") -> list[Path]:
        """Every dated run directory this library has written under ``base``, oldest first.

        Matched on the library's own name, so a folder holding several designs only gives
        back this one's runs and nothing has to be typed out::

            lib.runs()                 # [out/hAcyP1_20260724_154541, out/hAcyP1_20260725_091203]
        """
        root = Path(base)
        if not root.is_dir():
            return []
        prefix = f"{self.spec.name}_"
        return sorted(
            p for p in root.iterdir()
            if p.is_dir() and p.name.startswith(prefix) and _RUN_STAMP.fullmatch(p.name[len(prefix):])
        )

    def latest_run(self, base: str | Path = "out") -> Path | None:
        """The most recent run directory for this library under ``base``, or None."""
        found = self.runs(base)
        return found[-1] if found else None

    def run(self, stamp: str, base: str | Path = "out") -> Path:
        """The run directory with this date stamp, e.g. ``lib.run("20260724_154541")``.

        A stamp is what actually picks a run out; the library already knows its own name, so
        that is all you give. A leading part of the stamp is enough (``"20260724"`` for a
        run that day) as long as it names exactly one."""
        found = [p for p in self.runs(base) if p.name[len(self.spec.name) + 1:].startswith(stamp)]
        if not found:
            available = [p.name[len(self.spec.name) + 1:] for p in self.runs(base)]
            raise FileNotFoundError(
                f"No run of {self.spec.name!r} under {base} stamped {stamp!r}"
                + (f"; available: {available}" if available else " (no runs there yet)")
            )
        if len(found) > 1:
            raise ValueError(
                f"{stamp!r} matches {len(found)} runs "
                f"({[p.name for p in found]}); give more of the stamp."
            )
        return found[0]

    def run_record(self, directory: str | Path) -> dict:
        """The design specs a run wrote, found by this library's name so the caller does not
        have to rebuild the filename."""
        path = Path(directory) / _specs_filename(self.spec.name)
        if not path.is_file():
            raise FileNotFoundError(
                f"No design-specs record for {self.spec.name!r} in {directory} "
                f"(expected {path.name})."
            )
        return json.loads(path.read_text())

    def matches_run(self, directory: str | Path) -> bool:
        """Whether a past run's record describes the library now in memory.

        Worth asking before adding files to a directory somebody already ordered from: the
        design is seeded, so rebuilding it gives the same sequences, but only if the spec has
        not moved on since. Compares the frozen reference, the member count, and the protein
        the run was built from."""
        rec = self.run_record(directory)
        spec = rec.get("spec", {})
        return (
            rec.get("reference_cds") == self.reference
            and rec.get("n_variants") == len(self)
            and spec.get("name") == self.spec.name
            and spec.get("protein_sequence") == self.spec.protein_sequence
        )

    def _output_target(self, sub: str) -> Path:
        """A subdirectory of this library's run directory, for exporters called with no path.
        Uses the directory ``export_all`` last wrote to, else this run's own."""
        base = self.output_dir if self.output_dir is not None else self.run_dir()
        return Path(base) / sub

    def run_dir(self, base: str | Path = "out") -> Path:
        """Create and return this run's own output directory, ``base/<name>_<stamp>``.

        The stamp is the moment the sequences were built, so exporting the same library
        twice refreshes one directory while a fresh run gets a fresh one. Pass the result
        to the exporters that take a path::

            out = lib.run_dir("out")
            lib.to_oligo_pool(out / "oligos.csv")

        Call it after ``codon_optimize()``. Before then there is no build time to stamp
        with, so you get a directory named for the moment you asked, which a later export
        will not write to.
        """
        from .io import run_directory

        return run_directory(base, self.spec.name, when=self._created_at)

    # --- pipeline -------------------------------------------------------
    def codon_optimize(self, seed: int | None = None) -> "Library":
        """Codon-optimize the library, adding a ``variable_dna`` column.

        The strategy follows ``self.kind``:

        - ``"scan"`` (SubstitutionScan): optimize the WT CDS **once** into a frozen
          reference, then stamp each variant's single target codon onto it. Every
          member then matches the reference except at its intended
          position. The stamped codon is chosen usage-first, stepping to a rarer
          synonymous codon only to avoid a restricted motif. Set
          ``spec.optimization.synonymous_fallback = False`` to refuse that step and
          have such a position reported as a failure instead.
        - ``"sequence_set"`` (SequenceSet): each member is an independent full-length
          protein, so each is codon-optimized on its own (no shared reference).

        A member that can't satisfy the sequence rules is recorded in ``self.failed``
        with ``variable_dna`` set to NA.

        Reproducible by default: ``seed`` falls back to ``spec.seed`` (0 unless changed)
        and is recorded in the design specs, so re-running a spec rebuilds the same
        library. The caller's own ``numpy.random`` state is left untouched. Set
        ``spec.seed = None`` for an unseeded run.

        Optimizing again rebuilds every sequence, so anything derived from the last run is
        dropped first: a tile layout and its ``tile``/``oligo``/``oligo_length`` columns,
        the failure record, a simulated assembly, and the design-specs entries for all
        three. Call ``tile()`` again afterwards. The spec is re-snapshotted at the same
        moment, so params edited since ``generate()`` are the ones recorded.
        """
        base = self.spec.seed if seed is None else seed
        self._created_at = datetime.now().astimezone()
        self._reset_derived()
        # Re-snapshot the spec: it may have been edited since generate(), and the record has
        # to name the parameters this run actually used.
        self.design_specs["spec"] = self.spec.to_dict()

        if self.kind == "sequence_set":
            from .optimize.independent import optimize_independent

            variable, failed = optimize_independent(self, seed=base)
            self.df["variable_dna"] = variable
            self.failed = failed
            self.design_specs.update(
                {"seed": base, "run_id": self.run_id, "created": self.created,
                 "versions": _tool_versions()}
            )
            if failed:
                self.design_specs["failed"] = failed
            return self

        from .optimize.backbone import optimize_library

        reference, variable, failed = optimize_library(self, seed=base)
        self.df["variable_dna"] = variable
        self.reference = reference
        self.failed = failed
        self.design_specs.update(
            {"seed": base, "run_id": self.run_id, "created": self.created,
             "reference_cds": reference, "versions": _tool_versions()}
        )
        if failed:
            self.design_specs["failed"] = failed
        return self

    def _reset_derived(self) -> None:
        """Drop everything computed from a previous optimization.

        Optimizing again rebuilds every sequence, so a tile layout, its oligos, and a list of
        failures from the last run all describe sequences that no longer exist. Leaving them
        in place used to export oligos that did not contain their own variants, and a
        design-specs record listing failures for variants that had since succeeded."""
        self.failed = {}
        self.assembly = None
        for key in ("failed", "tiled", "assembly", "handoff"):
            self.design_specs.pop(key, None)
        if self.tiles is not None:
            self.tiles = None
            self._tiled_params = None
            self._drop_wt_controls()
            self.df = self.df.drop(columns=["tile", "oligo", "oligo_length"], errors="ignore")
            for attr in ("_primer_set", "_unplaced", "_vector_extra_sites"):
                if hasattr(self, attr):
                    delattr(self, attr)

    def _drop_wt_controls(self) -> None:
        """Remove the per-tile WT control rows a previous ``tile()`` appended.

        Each one belongs to a specific layout, a tile index and that tile's primers, so a
        new layout rebuilds them. Left in place they would also be tiled a second time as if
        they were variants, and the pool would grow a duplicate WT member per tile."""
        if self._wt_control_names:
            keep = ~self.df["name"].isin(self._wt_control_names)
            self.df = self.df[keep].reset_index(drop=True)
        self._wt_control_names = []

    def tile(self, params: TiledAssemblyParams | None = None) -> "Library":
        """Lay a (codon-optimized) long-CDS library out for tiled assembly.

        Splits the frozen reference into overlapping tile windows sized to the oligo
        budget, assigns each variant to its tile, picks a tile-specific orthogonal
        primer pair, and assembles the final Golden Gate oligo (adds ``tile``,
        ``oligo``, ``oligo_length`` columns). Per-tile primers and destination vectors
        are stored on ``self.tiles``. ``params`` defaults to ``spec.tiled`` (or built-in
        defaults). Variants with no single containing tile (the global WT control) get NA.

        One WT member per tile is appended, ``WT_Tile_0`` and so on, each carrying its
        tile's window straight from the reference. Every tile is amplified out of the pool
        and assembled on its own, so a sublibrary without one has nothing unmutated to
        normalize against. The global ``WT`` row stays in the table for the record and rides
        on no oligo. Set ``tiled.wt_controls = False`` to order only the mutants.

        Needs a codon-optimized library, since the tiles are cut from the frozen
        reference. The params it resolves, with the terminal overhangs filled in from the
        starting vector when there is one, are kept on ``lib.tiled_params``, which is what
        QC and the vector exporters read.
        """
        from dataclasses import asdict

        from .layout.tiled import tile_library

        if params is None:
            params = self.spec.tiled or TiledAssemblyParams()
        self.assembly = None                  # a new layout means new oligos to simulate
        self.design_specs.pop("assembly", None)
        self._drop_wt_controls()              # the last layout's controls, before laying out again
        result = tile_library(self, params)
        # The resolved params, terminal overhangs filled in from the starting vector when
        # there is one. Everything downstream (QC, the vector exporters) reads these off
        # the library rather than spec.tiled, so an explicit tile(params) is honored.
        params = result["params"]
        self._tiled_params = params
        self.df["tile"] = result["tile"]
        self.df["oligo"] = result["oligo"]
        # One WT member per tile, appended as rows of their own: each sublibrary is
        # amplified and assembled separately, so each needs its own unmutated oligo.
        wt = result["wt_controls"]
        if wt is not None and len(wt):
            self.df = pd.concat([self.df, wt], ignore_index=True)
            self._wt_control_names = [str(n) for n in wt["name"]]
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
            "wt_controls": self._wt_control_names,
            # Where the boundaries came from, and what the overhangs cost either way, so the
            # record says whether the search was on and what it bought.
            "overhangs_optimized": params.optimize_overhangs,
            "overhang_cost": result["overhang_cost"],
            "overhang_cost_unsearched": result["overhang_cost_unsearched"],
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

    def destination_vector(self, strict: bool = True):
        """The one destination vector a standard (untiled) library clones into, as a
        ``DestinationVector``: the starting plasmid with the CDS replaced by the Golden
        Gate drop-out, plus the two fused overhangs the cut vector presents.

        Needs ``spec.starting_vector`` and a codon-optimized reference. The drop-out is
        placed so those overhangs are the ones the digested oligo carries, which is read
        out of the adaptors. A tiled library has one vector per tile instead, on
        ``self.tiles``.

        Raises when the plasmid and the adaptors between them give the cut vector two
        overhangs that collide, since that vector re-closes empty or takes the insert
        backwards and there is no correct product to clone. Pass ``strict=False`` to build
        it anyway and inspect it, or read ``lib.overhang_pairs()`` for the full picture."""
        from .layout.destination import build_destination

        return build_destination(self, strict=strict)

    def simulate_assembly(self):
        """Simulate the Golden Gate reaction(s) this design implies, one ``AssemblyResult``
        per destination vector.

        Digests the oligos and the vector, anneals the fused overhangs, ligates, and aligns
        each product against the parent vector. ``result.product`` is the plasmid the WT
        reaction yields, the clone you expect to sequence. QC runs this too (see
        ``check()``); call it directly when you want the assembled sequences in hand.
        Returns an empty list when there is nothing to assemble into, either no destination
        vector or no frozen reference to assemble against. A sequence set has no shared
        reference, so it falls in the second case even after ``codon_optimize()``.

        The results are left on ``self.assembly`` and summarized in the design specs, so a
        notebook can look at them later without re-running. Each call simulates afresh, so
        the attribute always describes the library as it stands; anything that rebuilds the
        sequences (``codon_optimize``, ``tile``, ``drop_failed``) clears it."""
        from .checks.assembly import simulate

        self.assembly = simulate(self)
        if self.assembly:
            self.design_specs["assembly"] = {
                "simulated": self.created,
                "reactions": [
                    {
                        "label": r.label,
                        "enzyme": r.enzyme,
                        "product_length": len(r.product) if r.product else None,
                        "members": r.n_members,
                        "rebuilt": r.n_correct,
                        "aligned_to_parent": r.n_aligned,
                        "issues": list(r.issues),
                        # Names only: the reasons are aggregated into `issues` by check(),
                        # and a per-member dict would swamp the record on a broken design.
                        "members_with_problems": sorted(r.problems),
                    }
                    for r in self.assembly
                ],
            }
        return self.assembly

    def parent_vector(self) -> str:
        """The starting plasmid with the frozen reference CDS in place: what a wild-type
        clone should be, and the baseline ``assembled_product`` is compared against."""
        from .checks.assembly import parent_vector

        return parent_vector(self)

    def assembled_product(self, name: str) -> str:
        """The plasmid one variant assembles into, in the same frame as ``parent_vector()``.

        The simulated clone for that member: its oligo digested, ligated into the cut
        destination vector, and rotated onto the parent's origin, so the two sequences line
        up and can be diffed directly."""
        from .checks.assembly import assembled_product

        return assembled_product(self, name)

    def overhangs(self):
        """The fused Golden Gate overhangs this design leaves, one row per end.

        Two per tile for a tiled library, or the destination vector's pair for a standard
        one. ``shared_self`` is how many of its bases each overhang shares with its own
        reverse complement, so a full count means it is palindromic and a copy of the
        fragment can anneal to it. ``lib.overhang_pairs()`` is the pairwise view."""
        from .checks.overhangs import overhang_table

        return overhang_table(self)

    def overhang_pairs(self, all_pairs: bool = False):
        """The pairs of fused overhangs that share a reaction, worst first, as a DataFrame.

        This is the table to read before ordering. One row per tile, since a tile's two ends
        are the only overhangs that ever meet: it is amplified out of the pool on its own and
        assembled into the vector built around its own window. ``shared`` counts the bases the
        pair has in common as written, which is how easily the cut vector re-closes on itself,
        and ``shared_flipped`` counts them against the reverse complement, which is how easily
        the fragment goes in backwards. A pair should share at most one base of either;
        ``risk`` grades that and ``note`` says what the consequence is.

        ``all_pairs=True`` also lists the cross-tile pairs. Those never share a tube, so their
        rows carry the counts for reference and a ``risk`` of ``"n/a"``."""
        from .checks.overhangs import pair_table

        return pair_table(self, all_pairs=all_pairs)

    def plot_overhangs(self, dpi: int | None = None):
        """Return a Figure of overhang homology, every end against every other (renders
        inline in a notebook). Needs a design with a Golden Gate reaction, so ``tile()`` or a
        starting vector first. ``dpi`` raises the resolution for a printable copy."""
        from .viz import DEFAULT_DPI, overhang_figure

        return overhang_figure(self, dpi=dpi or DEFAULT_DPI)

    def check(self, fmt: str = "report"):
        """Run QC (optimization, translation, forbidden sites, length), returning a ``CheckReport``.

        Those four run for every library, and the rest depend on what it is. A tiled library gets
        the per-oligo and per-tile-vector checks, a standard library with a starting vector
        gets its adaptors checked against the plasmid, and either one gets the assembly
        simulated end to end, digested, ligated, and aligned against the parent. Raises if
        the library has not been codon-optimized.

        ``fmt`` picks what you get back. ``"report"`` (the default) is the ``CheckReport``
        itself, which prints as the readable report and carries ``passed``, ``issues``, and
        ``advisories``. ``"text"`` is that report as a plain string, for a log file or an
        email. ``"dict"`` is plain data, every field plus ``passed`` and the two parsed
        views, ready for JSON or a DataFrame.
        """
        if fmt not in ("report", "text", "dict"):
            raise ValueError(f"Unknown fmt {fmt!r} (use 'report', 'text', or 'dict').")
        from .checks import check_library

        report = check_library(self)
        if fmt == "text":
            return report.text()
        return report.to_dict() if fmt == "dict" else report

    def compare_reference(self, other_dna: str, label: str = "IDT"):
        """Cross-check an externally codon-optimized WT CDS (e.g. pasted from IDT's
        Codon Optimization Tool) against this library's frozen reference, returning a
        ``ReferenceComparison``: codon agreement, GC, mean adaptiveness, and any
        BsaI/BsmBI/Shine-Dalgarno sites the external sequence would introduce that our
        design avoids. Submit the protein the reference encodes
        (``spec.designed_sequence``) to IDT."""
        from .compare import compare_reference

        return compare_reference(self, other_dna, label=label)

    def summary(self):
        """QC + structural metadata (per-sublibrary counts, params, adaptors,
        platform), returning a ``LibrarySummary``."""
        from .summary import summarize

        return summarize(self)

    def plot_codon_usage(self, metric: str = "frequency", compare: str | None = None,
                         compare_label: str = "IDT", dpi: int | None = None):
        """Return a codon-usage QC Figure (renders inline in a notebook). ``metric`` is
        'frequency' (absolute usage, default) or 'adaptiveness'. Pass ``compare`` an
        external CDS (e.g. IDT's WT optimization) to overlay it as a dashed line. Needs the
        one shared reference a scan has, so it does not apply to a ``SequenceSet``. ``dpi``
        raises the resolution for a printable copy."""
        from .viz import DEFAULT_DPI, codon_usage_figure

        return codon_usage_figure(self, metric=metric, compare=compare,
                                  compare_label=compare_label, dpi=dpi or DEFAULT_DPI)

    def plot_codon_matrix(self, log: bool = True, reference_only: bool = False,
                          dpi: int | None = None):
        """Return a codon-map Figure (renders inline in a notebook): which codon sits at
        every position of the CDS.

        Codons run down the y axis grouped by amino acid, each group banded and ordered
        most- to least-used in the host, so a codon drawn low in its own band is one the
        design had to compromise on. A cell counts the members carrying that codon there, so
        the frozen reference reads as a path and each stamped substitution as a mark off it.
        ``reference_only=True`` plots the reference alone; ``log=False`` uses a linear colour
        scale, which suits a ``SequenceSet`` whose members spread more evenly. Cells are a
        couple of pixels wide, so raise ``dpi`` for a long CDS or a printable copy."""
        from .viz import DEFAULT_DPI, codon_matrix_figure

        return codon_matrix_figure(self, log=log, reference_only=reference_only,
                                   dpi=dpi or DEFAULT_DPI)

    def plot_tiling(self, dpi: int | None = None):
        """Return a tile-layout Figure (renders inline in a notebook). Needs a tiled
        library, call ``tile()`` first. ``dpi`` raises the resolution for a printable copy."""
        from .viz import DEFAULT_DPI, tiling_figure

        return tiling_figure(self, dpi=dpi or DEFAULT_DPI)

    def _save_figure(self, fig, path: str | Path, dpi: int | None = None) -> None:
        """Write a figure and drop it, so a long export does not accumulate open figures.

        With no ``dpi`` the file inherits the figure's own, so one number governs how a plot
        looks inline and on disk instead of the two drifting apart."""
        import matplotlib.pyplot as plt

        Path(path).parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(path, dpi=dpi or fig.dpi, bbox_inches="tight")
        plt.close(fig)

    def to_qc_plots(self, path: str | Path, dpi: int | None = None) -> "Library":
        """Save the codon-usage QC figure to an image file (format from extension). The
        codon map is a separate figure; see ``plot_codon_matrix``. ``export_all`` writes
        both. ``dpi`` sets the resolution of the file."""
        from .viz import DEFAULT_DPI, codon_usage_figure

        self._save_figure(codon_usage_figure(self, dpi=dpi or DEFAULT_DPI), path)
        return self

    def drop_failed(self) -> "Library":
        """Remove variants whose optimization failed (``variable_dna`` is NA).

        The design specs drop them too, so the record cannot list failures for variants
        the library no longer contains."""
        if "variable_dna" in self.df.columns:
            self.df = self.df[self.df["variable_dna"].notna()].reset_index(drop=True)
        self.failed = {}
        self.assembly = None                  # it counted members that are no longer here
        for key in ("failed", "assembly"):
            self.design_specs.pop(key, None)
        return self

    # --- exports (return self so they chain) ----------------------------
    def to_full_csv(self, path: str | Path) -> "Library":
        """Write the master table, every column of ``lib.df`` plus the assembled construct.

        Adds ``sequence``, its ``length``, the ``gc_content`` of the variable region, and
        ``stamp_adaptiveness``, the relative adaptiveness of the codon this variant carries
        at its own position. Rows whose optimization failed are kept with a blank sequence,
        so a failed run can still be read."""
        _io.to_full_csv(self, path)
        return self

    def to_usortm(self, path: str | Path) -> "Library":
        """Write the uSort-M handoff, ``name,sequence`` with one row per variant.

        The sequence is case-encoded, adaptors lowercase and the variable region uppercase,
        which is how uSort-M reads the region boundaries. Refuses a library with a failed
        variant, and refuses names carrying ``/``, ``|``, ``>``, or whitespace. A tiled
        library raises ``NotImplementedError``, since a tiled pool has no one variable
        region to write; use ``to_oligo_pool`` for the physical order."""
        _io.to_usortm(self, path)
        return self

    def to_vendor(self, path: str | Path, method: str | None = None,
                  pool_name: str | None = None) -> "Library":
        """Write a synthesis-provider order form, pooled or arrayed. ``method`` defaults to
        the shape ``spec.platform`` implies. The pooled form carries one shared pool name
        from ``pool_name`` or ``spec.name`` and no variant names; the arrayed form names
        every row. For a tiled library the molecules ordered are the tile oligos, not each
        variant's full construct. Refuses a library with a failed variant."""
        _io.to_vendor(self, path, method=method, pool_name=pool_name)
        return self

    def to_oligo_pool(self, path: str | Path) -> "Library":
        """Write the single pooled synthesis order for a tiled library, ``name,sequence``
        with one assembled oligo per placed variant. Names are validated as for uSort-M.
        Needs a tiled library, so call ``tile()`` first."""
        _io.to_oligo_pool(self, path)
        return self

    def to_oligo_files(self, directory: str | Path, fmt: str = "genbank") -> "Library":
        """Write one file per ordered oligo into ``directory``, named for the variant.

        Each file holds the sequence the order carries, so a tiled library writes its
        assembled oligo and anything else writes the whole construct with adaptors.
        ``fmt`` is ``"genbank"`` (annotated with the coding stretch, the mutated codon, and
        the Type IIS sites), ``"fasta"`` (sequence only), or ``"both"``. ``export_all``
        already writes these to ``oligos/``, so call this yourself to put them somewhere
        else or in the other format."""
        _io.to_oligo_files(self, directory, fmt=fmt)
        return self

    def to_primer_order(self, path: str | Path, scale: str = "25nm",
                        purification: str = "STD") -> "Library":
        """Write the per-tile amplification primers as a headerless IDT bulk order
        (``name,sequence,scale,purification``). Two rows per tile, forward and reverse,
        5'->3' as ordered. Needs a tiled library."""
        _io.to_primer_order(self, path, scale=scale, purification=purification)
        return self

    def to_vectors(self, path: str | Path) -> "Library":
        """Write a manifest of the Golden Gate destination vectors you build to run the
        library, sequences included. A tiled library gets one row per tile, a standard
        library the single row it needs. For annotated plasmid maps use
        ``to_vector_maps``."""
        _io.to_vectors(self, path)
        return self

    def to_vector_maps(self, directory: str | Path | None = None) -> "Library":
        """Write the destination plasmids to clone as annotated GenBank maps.

        Each ``.gb`` is a full destination vector (the backbone with a window dropped out)
        with the two BsaI sites, the drop-out, the fused overhangs, and the retained CDS
        arms annotated, plus a ``destination_vectors.csv`` manifest. A tiled library gets
        one map per tile, a standard library a single ``destination.gb``. Requires a
        ``starting_vector``. With no ``directory`` they go in this library's run
        directory."""
        if directory is None:
            directory = self._output_target("vectors" if self.tiles is not None else "vector")
        _io.to_vector_maps(self, directory)
        return self

    def to_assembled_vectors(self, directory: str | Path | None = None,
                             fmt: str = "both") -> "Library":
        """Write the clone every variant assembles into, the plasmids to sequence against.

        One annotated GenBank per variant with its mutated codon marked, an
        ``all_clones.fasta`` holding the same set, the parent plasmid as ``parent_WT.gb``,
        and an ``assembled_vectors.csv`` manifest saying how each clone differs from the
        parent. ``fmt`` is ``"genbank"``, ``"fasta"``, or ``"both"``. Needs a starting
        vector; ``export_all`` calls this on its own when one is set.

        With no ``directory`` the clones go in ``assembled_vectors/`` inside this library's
        run directory, so back-filling a past run reads
        ``lib.to_assembled_vectors(lib.latest_run() / "assembled_vectors")`` and a fresh one
        needs no path at all."""
        _io.to_assembled_vectors(self, directory or self._output_target("assembled_vectors"),
                                 fmt=fmt)
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
                   plots: bool = True, timestamp: bool = True, vectors: bool = True,
                   oligos: bool = True, oligo_fmt: str = "genbank",
                   dpi: int | None = None) -> "Library":
        """Write the full output set: master CSV, uSort-M ``variants.csv``, vendor order
        form, the design-specs JSON (``<name>_design_specs.json``, the record uSort-M
        reads), and two QC plots: codon usage along the CDS and the codon map (which codon
        sits at each position, grouped by amino acid). Plots are best-effort; if one can't be
        rendered the data files still export (with a warning). ``dpi`` sets their
        resolution.

        ``oligos/`` holds the same sequences one file per member, annotated GenBank by
        default (``oligo_fmt="fasta"`` or ``"both"`` for the other forms). Pass
        ``oligos=False`` to skip it when a library of thousands would make that directory
        unwieldy.

        With a starting vector set, the cloning outputs come too: ``destination_vector.csv``
        and an annotated ``vector/`` map for the plasmid to build, and ``assembled_vectors/``
        holding the clone every variant assembles into, one GenBank each plus a combined
        FASTA. Pass ``vectors=False`` to skip those, e.g. for a library of thousands where
        one file per clone is more than you want.

        Files go in a dated run directory under ``output_dir``
        (``out/<name>_20260724_143210``), so a second run cannot overwrite the first and
        every file says which run produced it. The stamp is the moment the sequences were
        built, so re-exporting one library refreshes its directory instead of making a new one.
        Pass ``timestamp=False`` to write straight into ``output_dir``, e.g. when a build
        script wants a fixed path. The directory used is left on ``self.output_dir``.

        Refuses to write anything if a variant failed optimization, since an order should
        not go out incomplete; inspect ``lib.failed`` or call ``drop_failed()``. An untiled
        library also needs a synthesis method, from ``method=`` or ``spec.platform``. Both
        are checked before the directory is created, so a library that is not ready creates
        nothing."""
        # Checked before the directory is made, so nothing is created for a library that
        # is not ready.
        _io._require_complete(self)
        tiled = self.tiles is not None
        if not tiled and method is None and self.spec.platform is None:
            raise ValueError(
                "Specify method='pooled'|'arrayed' or set spec.platform before exporting."
            )
        out = self.run_dir(output_dir) if timestamp else Path(output_dir)
        out.mkdir(parents=True, exist_ok=True)
        self.output_dir = out
        self.design_specs["output_dir"] = str(out)
        name = self.spec.name
        self.to_full_csv(out / f"{name}_full_library.csv")
        if tiled:
            # A tiled pool hands off as oligos plus their amplification primers. There is no
            # uniform variable region to write, which is why to_usortm refuses it.
            self.to_oligo_pool(out / f"{name}_oligo_pool.csv")
            self.to_primer_order(out / f"{name}_primers.csv")
        else:
            self.to_usortm(out / "variants.csv")
        self.to_vendor(out / f"{name}_order.csv", method=method)
        if oligos:
            self.to_oligo_files(out / "oligos", fmt=oligo_fmt)
        self.to_design_specs(out / _specs_filename(name))
        if vectors and tiled and self.spec.resolve_vector(self.tiled_params) is not None:
            self.to_vectors(out / "destination_vectors.csv")
            self.to_vector_maps(out / "vectors")
        elif vectors and self.spec.vector is not None and not tiled:
            self.to_vectors(out / "destination_vector.csv")
            self.to_vector_maps(out / "vector")
            self.to_assembled_vectors(out / "assembled_vectors")
        if plots:
            from .viz import DEFAULT_DPI, codon_matrix_figure, codon_usage_figure

            # The codon-usage plot needs one shared reference to plot against, which a sequence
            # set has not got, so skip it there rather than warn about it on every export. The
            # codon map counts across members, so it works either way.
            wanted = [("codon_matrix", codon_matrix_figure)]
            if self.kind != "sequence_set":
                wanted.insert(0, ("codon_usage", codon_usage_figure))
            for label, build in wanted:
                try:
                    self._save_figure(build(self, dpi=dpi or DEFAULT_DPI),
                                      out / f"{name}_{label}.png")
                except Exception as exc:
                    import warnings
                    warnings.warn(
                        f"Skipping the {label} plot ({exc}).", RuntimeWarning, stacklevel=2,
                    )
        return self
