"""library_designer, design DNA variant libraries for pooled/arrayed synthesis.

Notebook-first Python API. Typical flow::

    from library_designer import LibrarySpec, SubstitutionScan

    spec = LibrarySpec(name="hAcyP1", protein_sequence="AEG...IVK",
                       substitutions=["F", "Y", "M", "A", "TAG"], truncation=6,
                       adaptor_5="ggcgcGGTCTCC", adaptor_3="CCTCTGGcggcg")
    lib = SubstitutionScan(spec).generate().codon_optimize()
    print(lib.check())
    lib.to_usortm("variants.csv")

Name the plasmid you clone into with ``starting_vector=`` and QC also checks the adaptors
against it, while ``to_vectors`` / ``to_vector_maps`` emit the destination vector to build::

    spec.starting_vector = "my_plasmid.gb"
    lib.destination_vector()      # the plasmid with the CDS replaced by the drop-out

For a CDS longer than one oligo, add ``tiled=TiledAssemblyParams()`` (or a ``[tiled]``
TOML block) and chain ``.tile()`` to split it across oligos with per-tile orthogonal
primers and Golden Gate assembly::

    lib = SubstitutionScan(spec).generate().codon_optimize().tile()
    lib.to_oligo_pool("oligos.csv"); lib.to_primer_order("primers.csv"); lib.to_vectors("vectors.csv")

For a library of independent full-length sequences (orthologs, generative designs,
deep multi-mutants) rather than single mutants, use ``SequenceSet``. Each member is
codon-optimized on its own, and the result feeds a whole-gene assembler like OMEGA::

    spec = LibrarySpec(name="my_designs", platform="pooled")
    lib = SequenceSet(spec, proteins={"design_01": "MKAIL...", ...}).generate().codon_optimize()
    lib.assemble_with_omega(OmegaParams(njunctions=40), omega_home="~/repos/omega")
"""
from __future__ import annotations

from .spec import (
    CodonOptimizationParams,
    LibrarySpec,
    StartingVectorParams,
    TiledAssemblyParams,
)
from .library import Library
from .generators import SequenceSet, SubstitutionScan
from .integrations.omega import OmegaParams

__version__ = "0.1.0"
__all__ = [
    "LibrarySpec",
    "CodonOptimizationParams",
    "StartingVectorParams",
    "TiledAssemblyParams",
    "Library",
    "SubstitutionScan",
    "SequenceSet",
    "OmegaParams",
    "__version__",
]
