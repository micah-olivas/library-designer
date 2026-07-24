# Library Designer

Design DNA libraries for commercial synthesis.

## Install

This project uses [uv](https://docs.astral.sh/uv/) for environment management (the Python version, the venv, and dependencies). To set up, clone the repo and run 

```bash
uv sync
```

or install the package into an existing environment:

```bash
pip install git+https://github.com/micah-olivas/library-designer
```

Launch the notebooks with the project environment active (no kernel registration needed):

```bash
uv run --with jupyterlab jupyter lab
```

## Usage

```python
from library_designer import LibrarySpec, CodonOptimizationParams, SubstitutionScan

spec = LibrarySpec(
    name="my_library",
    protein_sequence="MKAILV...GADEQ",           # your protein, one-letter, no stop codon
    substitutions=["A", "G", "P"],               # amino acids, or literal codons (e.g. "TTT")
    adaptor_5="ggcgcGGTCTCC",                    # 5' flank added to every oligo (here a BsaI cloning site)
    adaptor_3="CCTCTGGcggcg",                    # 3' flank added to every oligo

    optimization=CodonOptimizationParams(species="e_coli", method="use_best_codon", gc_max=0.68),
)
# or load a saved spec: spec = LibrarySpec.from_toml("my_spec.toml")

lib = SubstitutionScan(spec).generate().codon_optimize()
print(lib.summary())

lib.to_full_csv("out/my_library_full.csv")                  # master table with length and gc_content
lib.to_usortm("out/variants.csv")                           # name,sequence for a downstream pooling/QC plan
lib.to_vendor("out/order.csv")                              # order form, method taken from spec.platform
lib.to_design_specs("out/my_library_design_specs.json")     # design record: spec, params, seed, versions
```

By default, `codon_optimize()` optimizes the WT CDS once into a frozen reference (`lib.reference`), then stamps each variant's single codon onto it. This ensures that each library member matches the reference sequence except at its intended position. A variant whose codon can't avoid a forbidden motif is recorded in `lib.failed` and reported by QC instead of aborting the run.

`lib.summary()` returns a `LibrarySummary` with variant counts per sublibrary, the codon-optimization parameters, the adaptor regions and construct length range, the synthesis platform, and the QC report.

Substitutions accept amino acids, where the optimizer picks the codon, or literal codons (e.g. `"TAG"`), which are placed verbatim and protected. Codon-optimization parameters (species, method, GC window, iterations) live on `CodonOptimizationParams` and are recorded in `to_design_specs()`, so every order documents how it was generated.

## Tiled assembly

When the defined CDS is longer than the maximum pooled oligo/fragment length, tiled assembly splits the CDS into acceptable-length tiles and designs one pool in which each tile is a sublibrary flanked by its own orthogonal primer pair. Each sublibrary may then be amplified out of the pool and assembled by Golden Gate into a per-tile destination vector that carries the rest of the WT CDS.

Add a `[tiled]` block (or a `TiledAssemblyParams`) and chain `.tile()`:

```python
from library_designer import LibrarySpec, SubstitutionScan

spec = LibrarySpec.from_toml("my_tiled_spec.toml")   # a CDS with a [tiled] block
lib = SubstitutionScan(spec).generate().codon_optimize().tile()
print(lib.summary())                             # per-tile primers and QC, including junction sites

lib.to_oligo_pool("out/oligos.csv")              # the pooled synthesis order (name, sequence)
lib.to_primer_order("out/primers.csv")           # per-tile amplification primers (IDT bulk format)
lib.to_vectors("out/vectors.csv")                # per-tile Golden Gate destination vectors
```

The WT reference is either codon-optimized (the default) or a native CDS supplied with `cds=`. Designed sequences are screened for enzyme sites to ensure orthogonality. Orthogonal primers come from a validated set ([Subramanian et al. 2018](https://doi.org/10.1093/synbio/ysx008)) by default, or may be provided with `tiled.primer_set=<path>`.

To generate destination vectors for downstream tiled assembly, point `tiled.starting_vector` at the starting plasmid (`.gb`, `.dna`, or `.fasta`) containing the WT CDS. The parent vector is screened for stray Golden Gate assembly sites to prevent off-target cleavage. Set `use_vector_cds=true` to freeze the CDS already in the plasmid and clone the vectors straight from it, which flags undesirable motifs in the CDS instead of recoding them.

```python
spec.tiled.starting_vector = "my_destination.gb"   # .gb / .dna / .fasta
lib = SubstitutionScan(spec).generate().codon_optimize().tile()
lib.to_vector_maps("out/vectors/")                 # one annotated GenBank per tile plus a manifest
```

## Explicitly-defined Libraries

To create a library of explictly-defined sequences, use `SequenceSet` to load a sequence dictionary or multi-sequence fasta file input. By default,
each sequence is codon-optimized independently:

```python
from library_designer import LibrarySpec, SequenceSet

spec = LibrarySpec(name="my_designs", platform="pooled")
lib = SequenceSet(spec, proteins={
    "ortholog_ec": "MKAILV...",
    "design_01":   "MKGILV...",
}).generate().codon_optimize()
# or load them from a FASTA of proteins (test set located in examples/):
# lib = SequenceSet.from_fasta(spec, "examples/example_designs.faa").generate().codon_optimize()

print(lib.summary())
lib.to_full_csv("out/designs_full.csv")
```
For these libraries, library-designer works with [OMEGA](https://github.com/RomeroLab/omega) (Freschlin et al. 2026) to perform pooled Golden Gate assembly
from an oligo pool starting point. The OMEGA step runs separately as a subprocess aside from library-designer; it is never imported or bundled to abide by the terms of its GPL-3.0 license. Point at your checkout with `omega_home` (or the `OMEGA_HOME` env var). See the worked example at `notebooks/tutorials/03-omega-sequence-set.ipynb`.

```python
from library_designer import OmegaParams
result = lib.assemble_with_omega(OmegaParams(njunctions=40), omega_home="~/repos/omega")
result.oligos.to_csv("out/oligo_order.csv", index=False)   # the pooled oligo order
```

## Outputs

The same library serializes several ways.

| Method | For | Notes |
|---|---|---|
| `to_full_csv` | you / records | all metadata plus the assembled sequence |
| `to_usortm` | uSort-M | `name,sequence`, case-encoded regions, name-validated |
| `to_vendor` | synthesis provider | method-aware (pooled or arrayed) |
| `to_design_specs` | you / uSort-M | run record: spec, params, seed, reference, versions (`<name>_design_specs.json`) |
| `to_oligo_pool` | synthesis provider, tiled | one pool of assembled tile oligos |
| `to_primer_order` | synthesis provider, tiled | per-tile amplification primers (IDT) |
| `to_vectors` | cloning, tiled | per-tile destination-vector sequences (manifest CSV) |
| `to_vector_maps` | cloning, tiled | annotated GenBank plasmid per tile (needs a starting vector) |
