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
    substitutions=["F", "Y", "M", "A", "TAG"],   # amino acids, or literal codons ("TAG" is an amber stop)
    adaptor_5="ggcgcGGTCTCC",                     # 5' flank added to every oligo (here a BsaI cloning site)
    adaptor_3="CCTCTGGcggcg",
    optimization=CodonOptimizationParams(species="e_coli", method="use_best_codon", gc_max=0.68),
)
# or load a saved spec: spec = LibrarySpec.from_toml("my_spec.toml")

lib = SubstitutionScan(spec).generate().codon_optimize()
print(lib.summary())

lib.to_full_csv("out/my_library_full.csv")          # master table with length and gc_content
lib.to_usortm("out/variants.csv")                   # name,sequence for a downstream pooling/QC plan
lib.to_vendor("out/order.csv")                       # order form, method taken from spec.platform
lib.to_design_specs("out/my_library_design_specs.json")  # design record: spec, params, seed, versions
```

By default, `codon_optimize()` optimizes the WT CDS once into a frozen reference (`lib.reference`), then stamps each variant's single codon onto it. This ensures that each library member matches the reference sequence except at its intended position. A variant whose codon can't avoid a forbidden motif is recorded in `lib.failed` and reported by QC instead of aborting the run.

`lib.summary()` returns a `LibrarySummary` with variant counts per sublibrary, the codon-optimization parameters, the adaptor regions and construct length range, the synthesis platform, and the QC report.

Substitutions accept amino acids, where the optimizer picks the codon, or literal codons (e.g. `"TAG"`), which are placed verbatim and protected. Codon-optimization parameters (species, method, GC window, iterations) live on `CodonOptimizationParams` and are recorded in `to_design_specs()`, so every order documents how it was generated.

## Libraries of diverse full-length sequences

Some libraries are not single mutants of one protein but a set of distinct full-length sequences: orthologs, generative designs, or deep multi-mutants. Use `SequenceSet` for these. Each member is codon-optimized on its own, there is no shared reference, which is what a whole-gene assembler like OMEGA expects.

```python
from library_designer import LibrarySpec, SequenceSet

spec = LibrarySpec(name="my_designs", platform="pooled")
lib = SequenceSet(spec, proteins={
    "ortholog_ec": "MKAILV...",
    "design_01":   "MKGILV...",
}).generate().codon_optimize()
# or load them from a FASTA of proteins (a ready-to-run set ships in examples/):
# lib = SequenceSet.from_fasta(spec, "examples/example_designs.faa").generate().codon_optimize()

print(lib.summary())
lib.to_full_csv("out/designs_full.csv")
```

To fragment each gene into a pooled Golden Gate oligo order, route the library to [OMEGA](https://github.com/RomeroLab/omega). OMEGA is a separate GPL-3.0 tool that `library_designer` shells out to; it is never imported or bundled. Point at your checkout with `omega_home` (or the `OMEGA_HOME` env var). The worked example is `notebooks/tutorials/03-omega-sequence-set.ipynb`.

```python
from library_designer import OmegaParams
result = lib.assemble_with_omega(OmegaParams(njunctions=40), omega_home="~/repos/omega")
result.oligos.to_csv("out/oligo_order.csv", index=False)   # the pooled oligo order
```

## Tiled assembly (CDS longer than one oligo)

When the CDS is longer than a synthesis oligo, such as a full-length gene, you can't put each single mutant on its own oligo. Tiled assembly splits the CDS into oligo-sized tiles and orders one pool in which each tile is a sublibrary flanked by its own orthogonal primer pair. You amplify each sublibrary out of the pool and drop it by Golden Gate into a per-tile destination vector that carries the rest of the WT CDS.

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

The WT reference is either codon-optimized (the default) or a native CDS you supply with `cds=`, used verbatim. Supply a native CDS when the exact codons matter, for example a human gene whose codons you want to preserve. Tiles are codon-aligned, so a mutated codon never splits across a boundary. Every assembled oligo is screened for the two intended enzyme sites, and any extra site fails QC, including one created at a primer or context junction. Pooled primers that would form a junction site are skipped during assignment. Orthogonal primers come from a validated set (Subramanian 2018 by default, 164 usable after BsaI screening; see `data/primers/PROVENANCE.md`), and you can supply your own with `tiled.primer_set=<path>`.

Without a backbone, `to_vectors` emits each destination vector as the bare coding-region cassette. To get the plasmids you actually clone into, point `tiled.starting_vector` at your backbone (`.gb`, `.dna`, or `.fasta`). The WT CDS is located inside it on either strand, each tile's window is swapped for the BsaI drop-out, and `to_vector_maps` writes one annotated GenBank plasmid per tile with the enzyme sites, fused overhangs, retained CDS arms, and the backbone's own features carried over. The whole plasmid is screened for stray BsaI. Set `use_vector_cds=true` to freeze the CDS already in the plasmid and clone the vectors straight from it, which flags Shine-Dalgarno sites and motifs in that CDS instead of recoding them. Reading plasmid files and writing GenBank needs the `tiled` extra (`uv sync --extra tiled`).

```python
spec.tiled.starting_vector = "my_destination.gb"   # .gb / .dna / .fasta
lib = SubstitutionScan(spec).generate().codon_optimize().tile()
lib.to_vector_maps("out/vectors/")                 # one annotated GenBank per tile plus a manifest
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
| `to_vector_maps` | cloning, tiled | annotated GenBank plasmid per tile (needs a starting vector + the `tiled` extra) |
