"""The design, as data. A ``LibrarySpec`` plus a generator fully determines a
library, so the spec doubles as the design-specs record and can be loaded from TOML.
"""
from __future__ import annotations

import tomllib
from dataclasses import asdict, dataclass, field
from html import escape as _esc
from pathlib import Path


def _preview(seq: str, keep: int = 12) -> str:
    """Middle-elide a long sequence for display: head...tail."""
    return seq if len(seq) <= 2 * keep + 1 else f"{seq[:keep]}...{seq[-keep:]}"

# Shine-Dalgarno-like motifs to avoid: ribosome-binding-site core, an 8-12 nt
# spacer, then a start codon. Regex, matching the MBO-038 amber-scan design.
DEFAULT_AVOID_PATTERNS: tuple[str, ...] = (
    "GGAGG.{8,12}[AG]TG",
    "AGGAG.{8,12}ATG",
)


@dataclass
class CodonOptimizationParams:
    """Codon-optimization knobs. Captured in the design specs so an order is self-documenting."""

    species: str = "e_coli"
    method: str = "use_best_codon"   # DNA Chisel: use_best_codon | match_codon_usage | harmonize_rca
    gc_min: float | None = None      # fraction, e.g. 0.30
    gc_max: float | None = None      # fraction, e.g. 0.68 (IDT eBlock GC ceiling)
    gc_window: int | None = None     # windowed GC in bp; None = whole sequence
    max_random_iters: int = 100_000


@dataclass
class TiledAssemblyParams:
    """Layout parameters for *tiled assembly*, a long-CDS library where each single
    mutant rides on a short oligo covering only its tile window, is amplified out of
    the shared pool by a tile-specific orthogonal primer pair, and is dropped by
    Golden Gate into a destination vector carrying the rest of the WT CDS.

    Tiling is pure layout: which variants exist and their sequences are already fixed
    by the frozen-reference stamp (see ``optimize/backbone.py``). Defaults reproduce
    the glucokinase (GCK-800) design. Captured in the design specs.
    """

    oligo_budget: int = 300            # hard cap on final oligo length (bp), incl. primers + sites
    enzyme: str = "BsaI"               # Golden Gate Type IIS enzyme (recognition from checks.motifs)
    overhang_len: int = 4              # fused-overhang length left after digestion
    spacer_5: str = "A"                # base(s) between the 5' recognition site and the overhang
    spacer_3: str = "T"                # base(s) between the overhang and the 3' recognition site
    vector_insert: str = "AGAGACCAAAAGGTCTCA"   # BsaI drop-out placeholder inserted into each destination vector
    primer_set: str = "subramanian2018"  # bundled set name or a path to a primer-set CSV (see primers.py)
    primer_length: int = 20            # expected primer length (informational / QC)
    tile_size: int | None = None       # override the per-tile window (bp); None means it is derived from the budget

    # Destination-vector backbone. Point at the plasmid you clone into and the per-tile
    # destination vectors are emitted as the full plasmid (see layout/vector_io.py and
    # io.to_vector_maps). Without one, only the CDS-region cassette is emitted and the
    # two vector_context overhangs below are read verbatim.
    starting_vector: str | None = None   # path to the destination plasmid (.gb / .dna / .fasta)
    use_vector_cds: bool = False         # freeze the CDS already in the starting vector as the reference (clone from it)
    insert_label: str | None = None      # locate the insert site by an annotated feature label (else the WT CDS is matched)
    insert_anchors: tuple[str, str] | None = None  # fallback: (5', 3') unique sequences bracketing the insert site
    topology: str | None = None          # override the file topology ("circular" / "linear"); FASTA defaults to linear
    vector_context_5: str = "TGGA"     # 5'-terminal overhang when NO starting vector is given (else derived from it)
    vector_context_3: str = "GGAG"     # 3'-terminal overhang when NO starting vector is given


@dataclass
class LibrarySpec:
    """Declarative description of a variant library."""

    name: str
    protein_sequence: str = ""               # the single WT protein for a SubstitutionScan; unused by a SequenceSet
    substitutions: list[str] = field(         # amino acids ("A") and/or literal codons ("TAG"); scan only
        default_factory=list
    )

    # Optional shaping
    cds: str | None = None                   # native WT CDS to use verbatim as the reference (else codon-optimized)
    tiled: TiledAssemblyParams | None = None  # when set, use tiled-assembly layout (long CDS split across oligos)
    truncation: int = 0                      # N-terminal residues to drop
    adaptor_5: str = ""                      # 5' flanking element (lowercased on export)
    adaptor_3: str = ""                      # 3' flanking element
    avoid_patterns: list[str] = field(
        default_factory=lambda: list(DEFAULT_AVOID_PATTERNS)
    )
    avoid_enzymes: list[str] = field(default_factory=lambda: ["BsaI"])
    optimization: CodonOptimizationParams = field(default_factory=CodonOptimizationParams)
    platform: str | None = None              # synthesis platform: "pooled"/"arrayed" or a provider slug ("twist_oligo_pools"); see methods.py
    max_oligo_length: int | None = None      # hard synthesis-length cap (bp, incl. adaptors)
    seed: int = 0                            # RNG seed for reproducible optimization

    @property
    def truncated_sequence(self) -> str:
        """The protein sequence after removing the N-terminal ``truncation`` residues."""
        return self.protein_sequence[self.truncation:]

    @classmethod
    def from_toml(cls, path: str | Path) -> "LibrarySpec":
        data = tomllib.loads(Path(path).read_text())
        opt = data.get("optimization")
        if isinstance(opt, dict):
            data["optimization"] = CodonOptimizationParams(**opt)
        tiled = data.get("tiled")
        if isinstance(tiled, dict):
            if isinstance(tiled.get("insert_anchors"), list):
                tiled["insert_anchors"] = tuple(tiled["insert_anchors"])
            data["tiled"] = TiledAssemblyParams(**tiled)
        return cls(**data)

    def to_dict(self) -> dict:
        return asdict(self)   # recurses into the nested optimization params

    # --- display --------------------------------------------------------
    def _opt_line(self) -> str:
        o = self.optimization
        gc = ""
        if o.gc_min is not None or o.gc_max is not None:
            gc = f", gc=({o.gc_min}, {o.gc_max})"
            if o.gc_window:
                gc += f" window={o.gc_window}"
        return f"{o.species}, {o.method}, iters={o.max_random_iters}{gc}"

    def __repr__(self) -> str:
        trunc = (
            f"  (truncation {self.truncation}, {len(self.truncated_sequence)} aa scanned)"
            if self.truncation
            else ""
        )
        return (
            f"LibrarySpec {self.name!r}\n"
            f"  protein:       {len(self.protein_sequence)} aa  {_preview(self.protein_sequence)}{trunc}\n"
            f"  substitutions: {', '.join(self.substitutions)}\n"
            f"  adaptors:      5' {self.adaptor_5}  |  3' {self.adaptor_3}\n"
            f"  optimization:  {self._opt_line()}\n"
            f"  platform:      {self.platform}   max_oligo_length={self.max_oligo_length}   seed={self.seed}"
        )

    def _repr_html_(self) -> str:
        trunc = (
            f" <span style='opacity:.6'>(truncation {self.truncation}, "
            f"{len(self.truncated_sequence)} aa scanned)</span>"
            if self.truncation
            else ""
        )
        rows = [
            ("name", _esc(self.name)),
            (
                "protein",
                f"{len(self.protein_sequence)} aa &nbsp; <code>{_esc(_preview(self.protein_sequence))}</code>{trunc}",
            ),
            ("substitutions", " ".join(f"<code>{_esc(s)}</code>" for s in self.substitutions)),
            ("adaptors", f"5' <code>{_esc(self.adaptor_5)}</code> &nbsp; 3' <code>{_esc(self.adaptor_3)}</code>"),
            (
                "avoid",
                f"{', '.join(map(_esc, self.avoid_enzymes)) or ', '} "
                f"<span style='opacity:.6'>+ {len(self.avoid_patterns)} motif pattern(s)</span>",
            ),
            ("optimization", _esc(self._opt_line())),
            ("platform", _esc(self.platform) if self.platform else ", "),
            ("max_oligo_length", self.max_oligo_length if self.max_oligo_length is not None else ", "),
            ("seed", self.seed),
        ]
        trs = "".join(
            "<tr>"
            f"<th style='text-align:right;padding:2px 10px;opacity:.7;white-space:nowrap'>{k}</th>"
            f"<td style='text-align:left;padding:2px 10px;font-family:var(--jp-code-font-family,monospace)'>{v}</td>"
            "</tr>"
            for k, v in rows
        )
        return (
            "<table style='border-collapse:collapse'>"
            "<caption style='text-align:left;font-weight:600;padding:2px 10px'>LibrarySpec</caption>"
            f"{trs}</table>"
        )
