# Library Designer

Design DNA libraries for commercial synthesis.

## Install

You need git and [uv](https://docs.astral.sh/uv/). uv handles the rest (the Python
version, the virtual environment, and the dependencies), so you don't have to install
Python 3.12 yourself. On macOS or Linux:

```bash
curl -LsSf https://astral.sh/uv/install.sh | sh
```

Clone the repo and build the environment:

```bash
git clone https://github.com/micah-olivas/library-designer.git
cd library-designer
uv sync
```

`uv sync` creates `.venv/` inside the repo, installs the dependencies, and installs
library-designer itself in editable mode, so edits to `src/` take effect without a
reinstall. Confirm it works:

```bash
uv run pytest -q                                          # the test suite, a few seconds
uv run python -c "from library_designer import LibrarySpec; print('ok')"
```

Run anything inside that environment by prefixing it with `uv run`, or activate the venv
once with `source .venv/bin/activate` and drop the prefix.

### Notebooks

The tutorials in `notebooks/tutorials/` run end to end on the bundled data in `examples/`,
so they are the quickest way to see a whole workflow. JupyterLab comes with the project, so
launching it needs no extra flags:

```bash
uv run jupyter lab                                  # from the repo
uv run --project ~/repos/library-designer jupyter lab    # from anywhere else
```

Working in a notebook that lives outside the repo, as an experiment folder usually does,
set the project once in your shell profile and drop that flag too:

```bash
export UV_PROJECT=~/repos/library-designer
uv run jupyter lab
```

`~/repos/library-designer/.venv/bin/jupyter lab` works from anywhere as well, which makes a
short alias easy. `uv sync --no-group notebooks` skips JupyterLab if you only want the
library.

Open `notebooks/tutorials/01-standard-scan.ipynb` first. `notebooks/templates/` holds
copy-and-edit skeletons of the same workflows; put your own working copies in
`notebooks/personal/`, which is git-ignored. See `notebooks/README.md` for what each
notebook covers.

### Installing into an existing environment

If you want the package alone, without the notebooks and example data, install it into any
Python 3.12+ environment:

```bash
pip install git+https://github.com/micah-olivas/library-designer.git
```

### OMEGA

Only the `SequenceSet` assembly step needs [OMEGA](https://github.com/RomeroLab/omega),
and library-designer calls it as a separate process rather than importing it (see
[Explicitly-defined Libraries](#explicitly-defined-libraries)). Install it on its own,
following its instructions, then point library-designer at the checkout:

```bash
git clone https://github.com/RomeroLab/omega.git ~/repos/omega
export OMEGA_HOME=~/repos/omega     # or pass omega_home="~/repos/omega" at the call site
```

The runner expects the CLI at `$OMEGA_HOME/code/omega.py` and calls it with `python`. If
OMEGA lives in its own conda env, name that interpreter with `OMEGA_PYTHON` (or
`omega_python=`).

## Usage

```python
from library_designer import LibrarySpec, CodonOptimizationParams, SubstitutionScan

spec = LibrarySpec(
    name="my_library",
    protein_sequence="MKAILV...GADEQ",           # your protein, one-letter, no stop codon
    substitutions=["A", "G", "P"],               # amino acids, or literal codons (e.g. "TTT")
    adaptor_5="ggtctccaagc",                     # 5' flank on every oligo: BsaI site, spacer, overhang
    adaptor_3="ggtgggagacc",                     # 3' flank: overhang, spacer, site reverse-complemented

    optimization=CodonOptimizationParams(species="e_coli", method="use_best_codon", gc_max=0.68),
)
# or load a saved spec: spec = LibrarySpec.from_toml("my_spec.toml")

lib = SubstitutionScan(spec).generate().codon_optimize()
print(lib.summary())

out = lib.run_dir("out")                                    # out/my_library_20260724_143210
lib.to_full_csv(out / "my_library_full.csv")                # master table with length and gc_content
lib.to_usortm(out / "variants.csv")                         # name,sequence for a downstream pooling/QC plan
lib.to_vendor(out / "order.csv")                            # order form, method taken from spec.platform
lib.to_design_specs(out / "my_library_design_specs.json")   # run record: spec, params, seed, versions

lib.export_all("out")                                       # the same set, named for you, plus the QC plots
```

Each run writes into its own dated directory, so a re-run cannot overwrite an order you
already sent and any file found later can be traced back to the run that produced it. The
stamp is the moment the sequences were built, so re-exporting one library refreshes its
directory rather than making a new one; the same value is recorded as `created` in the design
specs, and the path is left on `lib.output_dir`. Pass `export_all(..., timestamp=False)`
to write straight into the directory you name, e.g. from a build script that wants a fixed
path.

By default, `codon_optimize()` optimizes the WT CDS once into a frozen reference (`lib.reference`), then stamps each variant's single codon onto it. Each member then matches the reference except at its intended position. A variant whose codon can't avoid a forbidden motif is recorded in `lib.failed` and reported by QC instead of aborting the run.

When the preferred codon for a mutated residue would spell a forbidden motif, the stamp steps down that residue's usage ranking and takes the next codon that avoids it, so the variant is still makeable at a rarer codon. Set `CodonOptimizationParams(synonymous_fallback=False)` to be told instead. The position then lands in `lib.failed` with the codon named, and no rarer codon is substituted. A pinned literal codon has no synonymous alternative, so it is placed verbatim and flagged either way.

`lib.summary()` returns a `LibrarySummary` with variant counts per sublibrary, the codon-optimization parameters, the adaptor regions and construct length range, the synthesis platform, and the QC report. Printing it keeps only the lines that say something, so an empty adaptor pair, an unset GC bound, and a check that found nothing do not crowd out the ones that need reading. A tiled library shows its tiles and oligo lengths in place of a construct length, since the oligo is what you order.

A single-mutant library is only interpretable if every
member matches the WT outside its own codon, otherwise a phenotype cannot be
pinned on the substitution. `rep.off_target_edits` names any member that strays, checked on
the sequences themselves, so it runs whether or not there is a plasmid to clone into. `assembly_aligned` checks the same thing on the simulated clone, and needs a destination
vector.

`lib.check()` returns a `CheckReport` that shows the same readable block whether you print
it or just echo it in a notebook. To act on the result in code, read `rep.passed` for the
verdict, `rep.issues` for a dict of only the checks that found something, and
`rep.advisories` for the informational notes that never fail a report. `rep.issues` is
empty exactly when `rep.passed` is true, so `if rep.issues:` is the branch to write.

```python
rep = lib.check()
if rep.issues:
    for check, entries in rep.issues.items():
        print(check, len(entries), entries[:3])   # e.g. enzyme_hits:BsaI 1 ['K4G']

lib.check(fmt="text")     # the report as a plain string, for a log
lib.check(fmt="dict")     # every field plus passed/issues/advisories, ready for JSON
```

`export_all` writes two QC figures. `<name>_codon_usage.png` plots host codon-usage
frequency along the CDS, the WT as a line and each substitution's codon as a point.
`<name>_codon_matrix.png` is the codon map: codons down the y axis grouped by amino acid,
each group banded and ordered most- to least-used in the host, CDS position across the x
axis, and each cell counting the members carrying that codon there. The frozen reference
reads as a dark path and the stamped substitutions as pale marks off it, so a codon drawn
low in its own band is one the design had to compromise on. `lib.plot_codon_matrix()`
returns it for a notebook, with `reference_only=True` to plot the WT alone.

Every plotting call takes `dpi`, and one number governs both the inline figure and the file
written from it. The default is 200 rather than matplotlib's 100, because the codon map draws
each cell a couple of pixels wide and a long CDS goes soft at the lower setting. Raise it for
a figure headed into a manuscript: `lib.plot_codon_matrix(dpi=600)`,
`lib.to_qc_plots(path, dpi=600)`, or `lib.export_all("out", dpi=600)` for both QC images at
once.

Substitutions accept amino acids, where the optimizer picks the codon, or literal codons (e.g. `"TAG"`), which are placed verbatim and protected. Codon-optimization parameters (species, method, GC window, iterations) live on `CodonOptimizationParams` and are recorded in `to_design_specs()`, so every order documents how it was generated.

Codon optimization searches stochastically, so it is seeded on every run and the seed goes
into the run record. Re-running a spec rebuilds the same library, and your own `numpy.random` state is left as
it was. Change `spec.seed` for a different
draw, or set it to `None` to follow the ambient RNG.

## Protein from UniProt

Give an accession instead of pasting residues:

```python
spec = LibrarySpec.from_uniprot("P07311", substitutions=["A"], truncation=1)
spec.name          # ACYP1_HUMAN, the entry name, unless you pass your own
```

Or name it on a spec you are building anyway, including from TOML (`uniprot = "P07311"`):

```python
spec = LibrarySpec(name="my_library", uniprot="P07311", substitutions=["A"], truncation=1)
```

The lookup runs once, at construction, and only when `protein_sequence` is left empty, so a
spec reaches the network exactly when you asked it to. The sequence is then stored on the
spec, and the entry it came from (name, organism, gene, sequence version, when it was
fetched) is recorded alongside it in `to_design_specs()`. That matters because UniProt
entries are revised: the record holds the residues that were actually used, so the library
stays reproducible after the entry changes. The FASTA is cached under
`~/.cache/library-designer/uniprot/`, so later runs work offline; `spec.resolve_uniprot(refresh=True)`
fetches it again, and `LIBRARY_DESIGNER_CACHE` moves the cache.

UniProt sequences carry the initiator methionine. Pass `truncation=1` when your construct
starts at the second residue, as the example above does. An accession that does not match
the UniProt format is refused before any request, so a pasted protein sequence cannot
become an HTTP call.

## Starting vector

Point `starting_vector` at the plasmid you clone into (`.gb`, `.dna`, or `.fasta`) and the
library is checked and exported against the real backbone:

```python
spec = LibrarySpec(
    name="my_library",
    protein_sequence="MKAILV...GADEQ",
    substitutions=["A"],
    adaptor_5="ggtctccaagc",                     # BsaI site, spacer, then the fused overhang
    adaptor_3="ggtgggagacc",
    starting_vector="my_plasmid.gb",             # or StartingVectorParams(...) for the details
)
lib = SubstitutionScan(spec).generate().codon_optimize()
print(lib.check())                               # now also checks the adaptors against the plasmid

dv = lib.destination_vector()                    # the plasmid with the CDS dropped out

out = lib.run_dir("out")
lib.to_vectors(out / "vector.csv")               # its sequence, overhangs, and drop-out window
lib.to_vector_maps(out / "vector")               # an annotated GenBank map plus a manifest
```

The insert site is found by matching a `cds=` you supplied, else an annotated feature
(`insert_label=`), else the plasmid's sole CDS feature, else two bracketing sequences
(`insert_anchors=`). A CDS on the minus strand is normalized, so either cloning
orientation reads the same way here.

QC gains three checks once a plasmid is named. Each adaptor has to carry one Type IIS site
pointing at the insert. The fused overhang the digested oligo keeps has to be the one the
cut vector presents, which is what catches an adaptor that looks fine on its own but will
not ligate. And the assembled destination vector has to carry only the two intended sites,
since the whole plasmid sees the enzyme.

An emitted circular map is read from an origin ~50 bp upstream of the promoter that drives
the insert, the nearest one 5' of it, so the cassette and its Golden Gate sites sit in the
middle of the map instead of being split across its two ends. With no promoter annotated
the origin backs off the same distance from the insert, and it shifts to a feature boundary
rather than cut through an annotation, since a feature spanning the origin cannot be
written to GenBank. The base it started at is recorded in the manifest as
`origin_in_starting_vector`, which is how you relate the emitted coordinates back to the
plasmid you supplied.

### Overhang specificity

A Golden Gate reaction is only directional if its two fused overhangs are distinct. When they do not, the cut vector's own ends anneal and it re-closes empty,
or the fragment goes in the other way round and you clone the tile backwards. QC counts how
many of the four bases the two share, in both orientations, and checks each against its own
reverse complement to catch a palindrome. Aim for at most one shared base. A full match
fails the report; a single mismatch is an advisory, because a tiled design reads its
overhangs off the coding sequence at the tile boundaries rather than picking them from an
orthogonal set.

The unit of concern is one tile. Each is amplified out of the pool with its own primer pair
and assembled into the vector built around its own window, so a tile's reaction holds its
own fragments and nothing else. Tile 0's overhangs and tile 3's never meet, and they cannot
be pooled either, since every tile needs the particular vector that drops out its own
window. So the comparisons that mean anything are a tile's two ends against each other and
each end against its own reverse complement.

```python
lib.overhangs()                     # one row per fused overhang
lib.overhang_pairs()                # one row per tile, worst first
lib.overhang_pairs(all_pairs=True)  # cross-tile rows too, listed but ungraded
lib.plot_overhangs()                # the same as a matrix
```

On an untiled library the pair comes from your backbone and your adaptors, so no recoding
can fix a collision. `destination_vector()` raises rather than hand back a plasmid that
cannot clone; pass `strict=False` to build it and look anyway.

### Moving the boundaries

A tiled design does not pick its overhangs, so the boundary is the only handle on them. The
budget caps a tile and the balanced split sits under that cap, which leaves spare codons, and
every one of them is a boundary that could move. `tiled.optimize_overhangs` searches those
positions and keeps the layout whose overhangs share the least sequence. The balanced split
is scored against every candidate and wins ties, so the search can only hold a layout steady
or improve it. It scores only what can misfire, a tile against itself, so the slack is not
spent on cross-tile homology that never meets. On the bundled glucokinase example the
balanced split leaves tile4 with two ends one mismatch from complementary, which would let
that fragment ligate in backwards. The search clears it and leaves five of the six reactions
at or under target.

Boundaries that move make the tiles uneven, and the oligos with them. `tiled.pad_oligos`
evens the pool back out to one length with filler between each primer and the recognition
site beside it. That position is what makes it safe: the pad is outside both sites, so it is
amplified with the oligo and cut away before anything ligates, and the fragment that reaches
the vector is unchanged. `pad_target` sets the length outright,
otherwise the pool levels up to the longest oligo the layout already needed.

```python
tiled = TiledAssemblyParams(oligo_budget=300, optimize_overhangs=True, pad_oligos=True)
```

Both default to off, so an existing library keeps the boundaries and the oligos it had.
Turning the search on rewrites every oligo, so re-order the pool rather than mixing layouts.
`lib.design_specs["tiled"]` records which layout was used and what the overhangs cost either
way.

### Assembly simulation

With a destination vector in play, QC stops reading the sequences and starts putting them
together. It digests the oligo you would order and the vector you would build, anneals the
fused overhangs, ligates, and aligns the product against the parent plasmid:

```
QC report: 480 variants, PASS
  translation round-trip: 480/480 ok
  no forbidden sequences (checked BsaI, 2 motifs)
  assembly simulation: 480/480 members rebuild their intended variant
  aligned to the parent vector: 480/480 differ only at the intended codon
```

The digest finds its own cut sites rather than trusting the layout, so a spacer off by one
base, an overhang drawn from the wrong place, or a tile window that does not line up with
its vector shows up as a product that is not the plasmid you meant to build. The alignment
is the end-to-end statement a single-mutant library rests on: the clone you get back differs
from the plasmid you started with at one codon, and it is the codon you asked for. A
synonymous change anywhere else in the coding sequence passes every other check and fails
this one.

`lib.simulate_assembly()` returns the same simulation with the sequences in hand, one
`AssemblyResult` per destination vector. `result.product` is the assembled plasmid, the clone
you expect to sequence.

The clones themselves come out with everything else `export_all` writes:

```python
lib.export_all("out")           # -> out/<name>_<date>_<time>/assembled_vectors/
```

`assembled_vectors/` holds one annotated GenBank per variant, with its mutated codon marked
as a `variation` feature and the backbone's own annotations carried over, plus
`all_clones.fasta` with the same set for aligning sequencing reads against,
`parent_WT.gb` as the baseline, and a manifest naming each clone's codon change and whether
it differs from the parent only where it should. `lib.assembled_product(name)` and
`lib.parent_vector()` return any single pair as strings, in one frame, ready to diff.

One finding worth knowing about in advance. When the fused overhang is drawn from the CDS's
own ends (an adaptor that stops after the recognition site and its spacer), a variant that
mutates one of those first or last four bases presents an overhang the cut vector does not,
so it cannot clone, and the simulation names those variants. Adaptors that spell the four
overhang bases themselves, as the ones above do, avoid the problem: the whole coding
sequence rides on the fragment and no mutable base sits in an overhang.

Where the drop-out begins depends on the adaptors, and both conventions work. An adaptor
that ends at the recognition site plus its spacer leaves the CDS's own first four bases as
the overhang, so the vector keeps those four bases. An adaptor that spells four more bases,
drawn from the backbone flanking the insert, lets the whole CDS drop out. Set
`use_vector_cds=True` to freeze the CDS already in the plasmid as the reference and clone
from it, which flags motifs in that CDS instead of recoding them.

## Tiled assembly

When the defined CDS is longer than the maximum pooled oligo/fragment length, tiled assembly splits the CDS into acceptable-length tiles and designs one pool in which each tile is a sublibrary flanked by its own orthogonal primer pair. Each sublibrary may then be amplified out of the pool and assembled by Golden Gate into a per-tile destination vector that carries the rest of the WT CDS.

Add a `[tiled]` block (or a `TiledAssemblyParams`) and chain `.tile()`:

```python
from library_designer import LibrarySpec, SubstitutionScan

spec = LibrarySpec.from_toml("my_tiled_spec.toml")   # a CDS with a [tiled] block
lib = SubstitutionScan(spec).generate().codon_optimize().tile()
print(lib.summary())                             # tiles, oligo lengths, and QC, including junction sites

out = lib.run_dir("out")
lib.to_oligo_pool(out / "oligos.csv")            # the pooled synthesis order (name, sequence)
lib.to_primer_order(out / "primers.csv")         # per-tile amplification primers (IDT bulk format)
lib.to_vectors(out / "vectors.csv")              # per-tile Golden Gate destination vectors
```

The WT reference is either codon-optimized (the default) or a native CDS supplied with `cds=`. Designed sequences are screened for enzyme sites to ensure orthogonality. Orthogonal primers come from a validated set ([Subramanian et al. 2018](https://doi.org/10.1093/synbio/ysx008)) by default, or may be provided with `tiled.primer_set=<path>`.

Each tile also gets its own WT member, `WT_Tile_0` and so on, which is that tile's window taken straight from the reference and flanked by the same primers as its mutants. Since a sublibrary is amplified and assembled on its own, this is the unmutated clone you sequence and normalize that tile against. These rows are in the pooled order alongside the mutants. The single `WT` row remains in `lib.df` for the record and carries no oligo of its own. Set `tiled.wt_controls = false` to order only the mutants.

A tiled library needs one destination vector per tile, each carrying the rest of the WT CDS. Name the starting plasmid as above with `starting_vector`, or in the `[tiled]` block with `tiled.starting_vector` when the layout has its own backbone. The parent vector is screened for stray Golden Gate assembly sites to prevent off-target cleavage. Set `use_vector_cds=true` to freeze the CDS already in the plasmid and clone the vectors straight from it, which flags undesirable motifs in the CDS instead of recoding them.

```python
spec.starting_vector = "my_destination.gb"         # .gb / .dna / .fasta
lib = SubstitutionScan(spec).generate().codon_optimize().tile()
lib.to_vector_maps(lib.run_dir("out") / "vectors")                # one annotated GenBank per tile plus a manifest
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
lib.to_full_csv(lib.run_dir("out") / "designs_full.csv")
```
For these libraries, library-designer works with [OMEGA](https://github.com/RomeroLab/omega) (Freschlin et al. 2026) to perform pooled Golden Gate assembly
from an oligo pool starting point. The OMEGA step runs separately as a subprocess aside from library-designer; it is never imported or bundled to abide by the terms of its GPL-3.0 license. Point at your checkout with `omega_home` (or the `OMEGA_HOME` env var). See the worked example at `notebooks/tutorials/03-omega-sequence-set.ipynb`.

```python
from library_designer import OmegaParams
result = lib.assemble_with_omega(OmegaParams(njunctions=40), omega_home="~/repos/omega")
result.oligos.to_csv(lib.run_dir("out") / "oligo_order.csv", index=False)   # the pooled oligo order
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
| `to_vectors` | cloning | the destination vector to build (one row per tile if tiled) |
| `to_vector_maps` | cloning | annotated GenBank of that vector (needs a starting vector) |
| `to_assembled_vectors` | sequence verification | the clone every variant assembles into: one annotated GenBank each, a combined FASTA, the parent, and a manifest (standard libraries only) |
| `to_oligo_files` | plasmid editors, aligners | one file per ordered oligo, annotated GenBank (or FASTA) |

`export_all` always writes `to_full_csv`, `to_vendor`, `to_design_specs`, and
`to_oligo_files`. A standard library also gets `to_usortm`; a tiled one gets
`to_oligo_pool` and `to_primer_order` in its place, because a tiled pool has no single
variable region to write. With a starting vector set the cloning outputs come too:
`to_vectors` and `to_vector_maps` for either kind, and `to_assembled_vectors` for a
standard library. Pass `vectors=False` to leave the cloning outputs out, or
`oligos=False` to skip the per-oligo directory.

The per-oligo files land in `oligos/`, one per member. They are annotated GenBank by
default, marking the coding stretch the oligo carries, the mutated codon, and every Type
IIS site on either strand, so an oligo opens in a plasmid editor already labelled. Pass
`oligo_fmt="fasta"` for sequence only, or `"both"`. A tiled library writes its assembled
oligo with the primers and enzyme sites on it; anything else writes the whole construct,
adaptors included. An amber `K7*` becomes `K7stop.gb`, since `*` globs in a shell, while
the record inside keeps the real name.

### The uSort-M handoff

uSort-M reads two files from a run directory: `variants.csv` and the design-specs JSON
beside it. The CSV is exactly `name,sequence` with flanking adaptors lowercase, and stays
that way, because uSort-M parses it strictly. The run's identity therefore travels in the
record rather than as extra columns:

```json
"run_id":  "hAcyP1_GCE_scan_summer26_20260724_154541",
"created": "2026-07-24T15:45:41-07:00",
"handoff": {
  "run_id": "hAcyP1_GCE_scan_summer26_20260724_154541",
  "variants_csv": "variants.csv",
  "n_variants": 384,
  "sha256": "318b1cfa38d7…",
  "format": "name,sequence; flanking adaptors lowercase, variable region uppercase"
}
```

`run_id` is the run directory's own name, also on `lib.run_id`. It is derived from the
moment the sequences were built, so it is stable across re-exports of one library, distinct
for the next run, and still meaningful once a file has been copied somewhere else. A
downstream tool that echoes it into its own outputs gives you one string to trace a plate,
a pool, or a sequencing result back to the run that produced it, and
`LibrarySpec.from_toml` plus the recorded seed rebuilds that run exactly.

The `sha256` is of `variants.csv` as written. Checking it before reading tells a consumer
whether the CSV in front of it is the one the record describes, which is what catches a
hand-edited file or two runs' outputs mixed in one folder.
