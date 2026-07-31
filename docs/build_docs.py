#!/usr/bin/env python3
"""Generate a self-contained HTML reference for library_designer.

Introspects the live package (signatures, docstrings, and the ``#`` comments on
dataclass fields) and writes a single offline-friendly ``docs/index.html``. The guide's
tables and plots are rendered from the bundled examples by ``_examples.py`` and embedded
as base64 SVG or inline HTML, so the page stays one file with no external assets. If an
artifact fails to render the page still builds without it.

Usage:
    uv run python docs/build_docs.py
"""
from __future__ import annotations

import ast
import dataclasses as dc
import html
import inspect
import io
import re
import sys
import tempfile
import textwrap
import tokenize
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import library_designer as ld  # noqa: E402
from library_designer.checks.assembly import AssemblyResult  # noqa: E402
from library_designer.checks.report import CheckReport  # noqa: E402
from library_designer.compare import ReferenceComparison  # noqa: E402
from library_designer.integrations.omega import OmegaResult  # noqa: E402
from library_designer.layout.destination import DestinationVector  # noqa: E402
from library_designer.layout.tiled import TileInfo  # noqa: E402
from library_designer.summary import LibrarySummary  # noqa: E402
from library_designer.uniprot import UniProtEntry  # noqa: E402

# --------------------------------------------------------------------------- #
# Reference layout. Unlike plater, whose __all__ is free functions, this package
# exports classes, so a section is a list of classes and each one contributes a class
# card plus a card per public method and property. The last three sections hold return
# types that are not exported but are what you actually inspect in a notebook.
# --------------------------------------------------------------------------- #
REFERENCE = [
    ("design", "The design", [
        ld.LibrarySpec, ld.CodonOptimizationParams,
        ld.StartingVectorParams, ld.TiledAssemblyParams,
    ]),
    ("generators", "Generators", [ld.SubstitutionScan, ld.SequenceSet]),
    ("library", "The library", [ld.Library]),
    ("qc", "QC results", [CheckReport, LibrarySummary, AssemblyResult, ReferenceComparison]),
    ("layout", "Vectors and tiling", [DestinationVector, TileInfo]),
    ("omega", "OMEGA", [ld.OmegaParams, OmegaResult]),
    ("uniprot", "Protein lookup", [UniProtEntry]),
]

# --------------------------------------------------------------------------- #
# The guide. Each entry is (anchor, title, intro HTML, code, [(kind, key, caption)]).
# 'kind' is 'img' for a base64 SVG or 'html' for a ready-made fragment. Keys index the
# dict from _examples.render_all().
# --------------------------------------------------------------------------- #
GUIDE = [
    (
        "install", "Install",
        "You need git and <a href='https://docs.astral.sh/uv/'>uv</a>. "
        "<code>uv sync</code> builds the environment and installs library-designer in "
        "editable mode.",
        """git clone https://github.com/micah-olivas/library-designer.git
cd library-designer
uv sync

uv run pytest -q          # the test suite, a few seconds""",
        [],
    ),
    (
        "spec", "The design, as data",
        "A <code>LibrarySpec</code> plus a generator fully determines a library, so the "
        "spec doubles as the design record. Write it in Python or load it from TOML. "
        "Give it a UniProt accession instead of pasting residues and the lookup runs "
        "once, at construction.",
        """from library_designer import LibrarySpec, CodonOptimizationParams

spec = LibrarySpec(
    name="hAcyP1",
    protein_sequence="AEGNTLISVDYE...",   # or uniprot="P07311"
    substitutions=["F", "Y", "A", "TAG"], # residues and/or codons
    truncation=6,                         # drop N-terminal residues
    adaptor_5="ggtctccaagc",              # site, spacer, overhang
    adaptor_3="ggtgggagacc",
    optimization=CodonOptimizationParams(species="e_coli"),
)

spec = LibrarySpec.from_toml("mbo038.toml")   # or load one""",
        [("html", "spec-table",
          "A spec renders as a table in a notebook. This is "
          "<code>examples/mbo038.toml</code>, the hAcyP1 amber and missense scan.")],
    ),
    (
        "generate", "Generate and optimize",
        "<code>generate()</code> builds the variant table, one row per substitution per "
        "position, plus a WT control. <code>codon_optimize()</code> then optimizes the WT "
        "CDS once into a frozen reference and stamps each variant's single codon onto it, "
        "so every member matches the reference except where you asked it to differ.",
        """from library_designer import SubstitutionScan

lib = SubstitutionScan(spec).generate().codon_optimize()
len(lib)          # 451 variants
lib.reference     # the frozen WT CDS members stamp onto
lib.failed        # members whose codon hit a restricted motif""",
        [("html", "variant-table",
          "The variant table <code>lib.df</code>, one row per member. "
          "<code>variable_dna</code> appears once <code>codon_optimize()</code> has run; "
          "<code>mut_index</code> is the 0-based codon it stamps.")],
    ),
    (
        "qc", "QC",
        "<code>check()</code> runs the checks the design supports and returns a "
        "<code>CheckReport</code>. <code>summary()</code> wraps that in the structural "
        "metadata you need to review an order. Name a starting vector and this same "
        "report grows the assembly lines in <span class='io-ref'>Out 7a</span>.",
        """print(lib.summary())     # counts, params, adaptors, platform, and QC in one block
report = lib.check()     # the CheckReport on its own
report.passed            # True when every check is clean""",
        [("html", "summary-text", "<code>lib.summary()</code>."),
         ("html", "qc-text",
          "<code>lib.check()</code>. Sites and motifs are judged on the assembled "
          "construct against an assembled-WT baseline, so the adaptors' own BsaI site "
          "does not count against every variant.")],
    ),
    (
        "codon-usage", "Codon usage",
        "The QC plot compares the library's codon usage against the host table it was "
        "optimized for. Pass <code>compare=</code> an externally optimized CDS (IDT's "
        "tool, say) to overlay it.",
        """fig = lib.plot_codon_usage()      # or metric='adaptiveness'
lib.to_qc_plots("out/usage.png")

cmp = lib.compare_reference(idt_cds)   # agreement, GC, new sites""",
        [("img", "codon-usage",
          "Codon usage across the library against the <em>E. coli</em> table it was "
          "optimized for.")],
    ),
    (
        "vector", "Starting vector",
        "Point <code>starting_vector</code> at the plasmid you clone into and the design "
        "is checked and exported against the real backbone. The insert site is found by "
        "a supplied <code>cds=</code>, an annotated feature, the plasmid's sole CDS "
        "feature, or two bracketing anchors.",
        """spec.starting_vector = "my_plasmid.gb"        # .gb / .dna / .fasta
lib = SubstitutionScan(spec).generate().codon_optimize()

dv = lib.destination_vector()                 # the plasmid with the CDS dropped out
lib.to_vectors("out/vector.csv")              # sequence, overhangs, drop-out window
lib.to_vector_maps("out/vector")              # annotated GenBank plus a manifest""",
        [("html", "vector-table",
          "<code>to_vectors()</code>. The drop-out window is given on the reference CDS, "
          "and <code>origin_in_starting_vector</code> relates the emitted coordinates "
          "back to the plasmid you supplied."),
         ("html", "vector-features",
          "The features on the emitted <code>destination.gb</code>. The two BsaI sites, "
          "the drop-out, and both fused overhangs are annotated, and backbone features "
          "that do not overlap the insert are carried over.")],
    ),
    (
        "assembly", "Assembly simulation",
        "With a destination vector in play, QC stops reading the sequences and starts "
        "putting them together. It digests the oligo and the vector built in "
        "<span class='io-ref in'>In 6</span>, anneals the fused overhangs, ligates, and "
        "aligns the product against the parent plasmid. The digest finds its own cut "
        "sites rather than trusting the layout, so a spacer off by one base shows up as "
        "a product that is not the plasmid you meant to build.",
        """results = lib.simulate_assembly()       # one AssemblyResult per vector
results[0].product                      # the plasmid the WT reaction yields

lib.parent_vector()                     # what a wild-type clone should be
lib.assembled_product("K57A")           # what one variant assembles into
lib.to_assembled_vectors("out/clones")  # a GenBank per clone, plus a manifest""",
        [("html", "assembly-text",
          "The assembly lines of the QC report. The alignment is the end-to-end claim a "
          "single-mutant library rests on."),
         ("html", "clone-diff",
          "One clone against the parent plasmid. A synonymous change anywhere else in "
          "the coding sequence would pass every other check and fail this one.")],
    ),
    (
        "tiled", "Tiled assembly",
        "When the CDS is longer than one oligo, tiled assembly splits it into tiles and "
        "designs one pool in which each tile is a sublibrary flanked by its own "
        "orthogonal primer pair. Each sublibrary is amplified out of the pool and "
        "assembled by Golden Gate into a per-tile destination vector carrying the rest "
        "of the WT CDS, so you get one vector per tile rather than the single one in "
        "<span class='io-ref'>Out 6a</span>.",
        """from library_designer import TiledAssemblyParams

spec.tiled = TiledAssemblyParams(oligo_budget=300, enzyme="BsaI")
lib = SubstitutionScan(spec).generate().codon_optimize().tile()

lib.to_oligo_pool("out/oligos.csv")     # the pooled synthesis order
lib.to_primer_order("out/primers.csv")  # per-tile primers, IDT bulk format
lib.to_vectors("out/vectors.csv")       # one destination vector per tile""",
        [("img", "tiling",
          "The tile layout for the bundled glucokinase example. Each variant rides on "
          "the tile whose window contains its codon."),
         ("html", "oligo-pool", "<code>to_oligo_pool()</code>, the single pooled order."),
         ("html", "primer-order",
          "<code>to_primer_order()</code>, the per-tile amplification primers as a "
          "headerless IDT bulk order.")],
    ),
    (
        "overhangs", "Overhang specificity",
        "A Golden Gate reaction is only directional if the two fused overhangs tell each "
        "other apart. When they do not, the cut vector's own ends anneal and it re-closes "
        "empty, or the fragment goes in the other way round. QC counts how many of the four "
        "bases the two share, in both orientations, and checks each against its own reverse "
        "complement so a palindrome is caught too. Aim for at most one shared base. A full "
        "match fails the report; a single mismatch is an advisory, because a tiled design "
        "reads its overhangs off the CDS rather than picking them from an orthogonal set. "
        "The unit is one tile: each is amplified out of the pool on its own and assembled "
        "into the vector built around its own window, so two tiles' overhangs never meet and "
        "homology between them is not a hazard.",
        """lib.overhangs()                   # one row per fused overhang
lib.overhang_pairs()              # one row per tile: its two ends, and what they would do
lib.overhang_pairs(all_pairs=True)  # cross-tile rows too, listed but ungraded
lib.plot_overhangs()              # the same as a matrix

# On an untiled library the pair comes from your backbone and adaptors, so a
# collision cannot be recoded away and building the vector refuses outright:
lib.destination_vector()                # raises on colliding overhangs
lib.destination_vector(strict=False)    # build it anyway and look""",
        [("img", "overhang-matrix",
          "Every overhang against every other, worse of the two orientations. Only the cells "
          "that can misfire are graded: the boxed pairs, which are a tile's own two ends, and "
          "the diagonal, an overhang against its own reverse complement, where a hot cell is "
          "a palindrome. Cross-tile cells are shown faded, since those tiles are assembled in "
          "separate tubes."),
         ("html", "overhang-pairs",
          "<code>overhang_pairs()</code>, one row per tile, worst first. This is the table to "
          "read before ordering.")],
    ),
    (
        "boundaries", "Moving the boundaries",
        "A tiled design does not pick its overhangs, so the tile boundary is the only handle "
        "on them. The budget caps a tile and the balanced split sits under that cap, which "
        "leaves spare codons, and every one of them is a boundary that could move. "
        "<code>optimize_overhangs</code> searches those positions and keeps the layout whose "
        "overhangs share the least sequence. The balanced split is scored against every "
        "candidate and wins ties, so the search can only hold a design steady or improve it. "
        "Boundaries that move make the tiles uneven and the oligos with them, so "
        "<code>pad_oligos</code> evens the pool back out to one length. The filler goes "
        "between each primer and the recognition site beside it, outside what the enzyme "
        "releases, so it is amplified with the oligo and cut away before anything ligates. "
        "Both default to off, so an existing design keeps the boundaries and oligos it had.",
        """spec.tiled = TiledAssemblyParams(
    oligo_budget=300,
    optimize_overhangs=True,   # move the boundaries to the least similar overhangs
    pad_oligos=True,           # one oligo length across the pool
    # pad_target=300,          # else the longest oligo the layout already needed
)
lib = SubstitutionScan(spec).generate().codon_optimize().tile()
lib.design_specs["tiled"]["overhang_cost"]   # and ..._unsearched, for the record""",
        [("html", "boundary-search",
          "The bundled glucokinase example. The balanced split leaves tile4 with two ends one "
          "mismatch from complementary, so that fragment can ligate in backwards. The search "
          "clears it and leaves five of the six reactions at or under target.")],
    ),
    (
        "sequence-set", "Libraries of explicit sequences",
        "For a library of independent full-length sequences (orthologs, generative "
        "designs, deep multi-mutants) use <code>SequenceSet</code>. Each member is "
        "codon-optimized on its own, with no shared reference, which is the input a "
        "whole-gene assembler expects.",
        """from library_designer import SequenceSet

lib = SequenceSet.from_fasta(spec, "designs.faa") \\
                 .generate().codon_optimize()

# Golden Gate from an oligo pool, via a separate OMEGA:
from library_designer import OmegaParams
result = lib.assemble_with_omega(OmegaParams(njunctions=40),
                                 omega_home="~/repos/omega")""",
        [("html", "sequence-set",
          "Each member optimized independently. OMEGA runs as a separate process, never "
          "imported, to keep its GPL-3.0 licence off this package.")],
    ),
    (
        "outputs", "Outputs",
        "<code>export_all</code> always writes the master CSV, the vendor order form, and "
        "the design-specs JSON. A standard library also gets <code>variants.csv</code> for "
        "uSort-M; a tiled one gets the oligo pool and primer order in its place. With a "
        "starting vector set the cloning outputs come too. Files go in a dated run "
        "directory, so a second run cannot overwrite the first.",
        """lib.export_all("out")     # -> out/hAcyP1_20260725_141530/
lib.output_dir            # the directory it used

out = lib.run_dir("out")  # or drive the exporters yourself
lib.to_full_csv(out / "full.csv")
lib.to_design_specs(out / "design_specs.json")""",
        [],
    ),
]


# --------------------------------------------------------------------------- #
# Docstring and field-comment extraction
# --------------------------------------------------------------------------- #
def esc(s: str) -> str:
    return html.escape(s, quote=True)


def real_doc(obj) -> str | None:
    """The written docstring, or None.

    ``@dataclass`` sets ``__doc__`` to the generated signature when a class has none, and
    ``inspect.getdoc`` returns it happily, so a class with no description would otherwise
    render a signature dump where its summary goes. Filter that out.
    """
    doc = inspect.getdoc(obj)
    if not doc:
        return None
    name = getattr(obj, "__name__", "")
    if name and doc.startswith(f"{name}(") and doc.rstrip().endswith(")"):
        return None
    return doc


_comment_cache: dict[str, dict[int, str]] = {}


def _comments_by_line(path: str) -> dict[int, str]:
    """{line number: comment text} for every ``#`` comment in a source file."""
    if path not in _comment_cache:
        out: dict[int, str] = {}
        with open(path, "rb") as fh:
            for tok in tokenize.tokenize(fh.readline):
                if tok.type == tokenize.COMMENT:
                    out[tok.start[0]] = tok.string.lstrip("#").strip()
        _comment_cache[path] = out
    return _comment_cache[path]


def field_docs(cls) -> dict[str, str]:
    """Field descriptions for a dataclass, lifted from its ``#`` comments.

    This package documents dataclass fields with a trailing comment on the field line,
    and occasionally with a comment block just above it. Neither is visible to
    ``inspect``, so read them out of the source instead. A one-line block of three words
    or fewer is treated as a section heading (``# Optional shaping``) rather than a
    description, and skipped.
    """
    try:
        path = inspect.getsourcefile(cls)
        src = Path(path).read_text()
    except (TypeError, OSError):
        return {}
    comments = _comments_by_line(path)

    tree = ast.parse(src)
    node = next((n for n in ast.walk(tree)
                 if isinstance(n, ast.ClassDef) and n.name == cls.__name__), None)
    if node is None:
        return {}

    fields = [(n.target.id, n.lineno, n.end_lineno or n.lineno)
              for n in node.body if isinstance(n, ast.AnnAssign) and isinstance(n.target, ast.Name)]
    field_lines = {ln for _, start, end in fields for ln in range(start, end + 1)}

    out: dict[str, str] = {}
    for name, start, end in fields:
        trailing = next((comments[ln] for ln in range(start, end + 1) if ln in comments), None)
        if trailing:
            out[name] = trailing
            continue
        block: list[str] = []
        ln = start - 1
        while ln in comments and ln not in field_lines:
            block.insert(0, comments[ln])
            ln -= 1
        if block and not (len(block) == 1 and len(block[0].split()) <= 3):
            out[name] = " ".join(block)
    return out


# --------------------------------------------------------------------------- #
# Prose rendering. This package writes plain prose with reST inline literals, ``::``
# literal blocks, and the occasional bullet list. There are no NumPy sections anywhere,
# so there is no section parser here.
# --------------------------------------------------------------------------- #
_INLINE_DOUBLE = re.compile(r"``(.+?)``")
_INLINE_SINGLE = re.compile(r"(?<![`\w])`([^`]+?)`(?!`)")
_BULLET = re.compile(r"^\s*[-*]\s")


def inline_code(escaped: str) -> str:
    s = _INLINE_DOUBLE.sub(r"<code>\1</code>", escaped)
    return _INLINE_SINGLE.sub(r"<code>\1</code>", s)


def _blocks(text: str):
    """Split prose into ('para' | 'code' | 'bullets', payload).

    A paragraph ending in ``::`` introduces the indented literal block that follows, the
    way reST does; the trailing colons are dropped from the prose.
    """
    text = textwrap.dedent(text).strip("\n")
    paras = re.split(r"\n\s*\n", text)
    i = 0
    while i < len(paras):
        para = paras[i]
        if not para.strip():
            i += 1
            continue
        lines = para.split("\n")
        if any(_BULLET.match(ln) for ln in lines):
            yield "bullets", para
            i += 1
            continue
        if para.rstrip().endswith("::") and i + 1 < len(paras):
            lead = para.rstrip()[:-2].rstrip()
            if lead:
                yield "para", lead
            yield "code", textwrap.dedent(paras[i + 1])
            i += 2
            continue
        yield "para", para
        i += 1


def _bullets_html(para: str) -> str:
    items, cur = [], []
    for ln in para.split("\n"):
        if _BULLET.match(ln):
            if cur:
                items.append(" ".join(cur))
            cur = [_BULLET.sub("", ln).strip()]
        elif ln.strip() and cur:
            cur.append(ln.strip())
    if cur:
        items.append(" ".join(cur))
    lis = "".join(f"<li>{inline_code(esc(it))}</li>" for it in items)
    return f"<ul class=\"doc-list\">{lis}</ul>"


def flow(text: str) -> str:
    """Render a prose block (everything after the summary line) as HTML."""
    out = []
    for kind, payload in _blocks(text):
        if kind == "bullets":
            out.append(_bullets_html(payload))
        elif kind == "code":
            out.append(f'<pre class="doc-code"><code>{highlight_code(payload)}</code></pre>')
        else:
            joined = " ".join(ln.strip() for ln in payload.split("\n") if ln.strip())
            out.append(f'<p class="pdesc">{inline_code(esc(joined))}</p>')
    return "\n".join(out)


def render_docstring(doc: str | None) -> str:
    if not doc or not doc.strip():
        return '<p class="muted">No description available.</p>'
    lines = doc.split("\n")
    i, summary = 0, []
    while i < len(lines) and lines[i].strip():
        summary.append(lines[i].strip())
        i += 1
    while i < len(lines) and not lines[i].strip():
        i += 1
    parts = [f'<p class="fn-summary">{inline_code(esc(" ".join(summary)))}</p>']
    rest = "\n".join(lines[i:])
    if rest.strip():
        parts.append(flow(rest))
    return "\n".join(parts)


def summary_of(doc: str | None) -> str:
    if not doc:
        return ""
    out = []
    for ln in doc.split("\n"):
        if not ln.strip():
            break
        out.append(ln.strip())
    return " ".join(out)


# --------------------------------------------------------------------------- #
# Build-time syntax highlighting
# --------------------------------------------------------------------------- #
_PY_KEYWORDS = {
    "import", "from", "as", "for", "in", "if", "else", "elif", "return", "def", "None",
    "True", "False", "and", "or", "not", "with", "while", "lambda", "class", "print",
}
_NAMESPACE = {"lib", "spec", "ld", "library_designer", "dv", "report", "result", "results"}

_HL_RE = re.compile(
    r"(?P<comment>\#[^\n]*)"
    r"|(?P<string>'[^'\n]*'|\"[^\"\n]*\")"
    r"|(?P<number>\b\d+\.?\d*\b)"
    r"|(?P<name>[A-Za-z_][A-Za-z0-9_]*)"
)


def highlight_code(code: str) -> str:
    out, pos = [], 0
    for m in _HL_RE.finditer(code):
        out.append(esc(code[pos:m.start()]))
        kind, text, cls = m.lastgroup, m.group(), None
        if kind == "comment":
            cls = "c"
        elif kind == "string":
            cls = "s"
        elif kind == "number":
            cls = "n"
        elif kind == "name":
            if text in _PY_KEYWORDS:
                cls = "k"
            elif text in _NAMESPACE:
                cls = "b"
            elif code[m.end():m.end() + 1] == "(":
                cls = "f"
        out.append(f'<span class="tok-{cls}">{esc(text)}</span>' if cls else esc(text))
        pos = m.end()
    out.append(esc(code[pos:]))
    return "".join(out)


# --------------------------------------------------------------------------- #
# Introspection
# --------------------------------------------------------------------------- #
def members_of(cls) -> list[dict]:
    """Public methods and properties of a class, in source order where we can get it."""
    out = []
    for name, obj in vars(cls).items():
        if name.startswith("_"):
            continue
        if isinstance(obj, property):
            doc = real_doc(obj.fget)
            out.append({"name": name, "kind": "property", "sig": "", "doc": doc,
                        "summary": summary_of(doc)})
        elif isinstance(obj, (classmethod, staticmethod)):
            fn = obj.__func__
            out.append({"name": name, "kind": type(obj).__name__, "sig": _sig(fn),
                        "doc": real_doc(fn), "summary": summary_of(real_doc(fn))})
        elif inspect.isfunction(obj):
            out.append({"name": name, "kind": "method", "sig": _sig(obj),
                        "doc": real_doc(obj), "summary": summary_of(real_doc(obj))})
    return out


def _sig(fn) -> str:
    try:
        s = inspect.signature(fn)
    except (ValueError, TypeError):
        return "(...)"
    params = [p for n, p in s.parameters.items() if n not in ("self", "cls")]
    return "(" + ", ".join(str(p) for p in params) + ")"


def class_signature(cls) -> str:
    """A constructor signature, but only when it is short enough to read.

    A dataclass with 17 fields produces a 586-character line that tells you nothing; those
    get the field table below instead.
    """
    try:
        s = str(inspect.signature(cls))
    except (ValueError, TypeError):
        return ""
    return s if len(s) <= 90 else ""


# --------------------------------------------------------------------------- #
# HTML assembly
# --------------------------------------------------------------------------- #
def fields_html(cls) -> str:
    if not dc.is_dataclass(cls):
        return ""
    docs = field_docs(cls)
    rows = []
    for f in dc.fields(cls):
        typ = f.type if isinstance(f.type, str) else getattr(f.type, "__name__", str(f.type))
        if f.default is not dc.MISSING:
            default = repr(f.default)
        elif f.default_factory is not dc.MISSING:      # type: ignore[misc]
            try:
                default = repr(f.default_factory())    # type: ignore[misc]
            except Exception:
                default = "..."
        else:
            default = ""
        desc = docs.get(f.name, "")
        rows.append(
            f'<dt><code class="pname">{esc(f.name)}</code>'
            f'<span class="ptype">{esc(str(typ))}</span>'
            + (f'<span class="pdefault">= {esc(default)}</span>' if default else "")
            + "</dt>"
            + (f'<dd><p class="pdesc">{inline_code(esc(desc))}</p></dd>'
               if desc else '<dd class="muted-dd"></dd>')
        )
    return '<h4 class="doc-section">Fields</h4><dl class="params">' + "".join(rows) + "</dl>"


def class_card(cls) -> str:
    name = cls.__name__
    sig = class_signature(cls)
    parts = [f'<section class="fn-card cls-card" id="{esc(name)}" '
             f'data-name="{esc(name)}" data-summary="{esc(summary_of(real_doc(cls)))}">']
    parts.append(
        f'<code class="fn-sig"><span class="cls-kw">class</span> '
        f'<span class="fn-name">{esc(name)}</span>'
        f'<span class="fn-args">{esc(sig)}</span></code>'
    )
    parts.append(render_docstring(real_doc(cls)))
    parts.append(fields_html(cls))
    parts.append("</section>")
    return "".join(parts)


def member_card(cls, m: dict) -> str:
    anchor = f"{cls.__name__}.{m['name']}"
    tag = ('<span class="badge">property</span>' if m["kind"] == "property"
           else '<span class="badge">classmethod</span>' if m["kind"] == "classmethod"
           else "")
    return (
        f'<section class="fn-card mem-card" id="{esc(anchor)}" '
        f'data-name="{esc(anchor)}" data-summary="{esc(m["summary"])}">'
        f'<code class="fn-sig"><span class="fn-owner">{esc(cls.__name__)}.</span>'
        f'<span class="fn-name">{esc(m["name"])}</span>'
        f'<span class="fn-args">{esc(m["sig"])}</span></code>{tag}'
        f'{render_docstring(m["doc"])}'
        "</section>"
    )


CSS = """
:root {
  --bg:#F8FAFC; --surface:#FFFFFF; --text:#1E293B; --muted:#64748B; --border:#E2E8F0;
  --code-bg:#EAEFF3; --accent:#2563EB; --accent-soft:#EFF4FE; --num:#334155;
  --out-bg:#F1F5F9; --sidebar-w:264px; --measure:820px;
  --out-fg:#0F766E; --out-bd:#5EEAD4; --out-soft:#F0FDFA;
  color-scheme: light;
}
[data-theme="dark"] {
  --bg:#0D1117; --surface:#171E29; --text:#E6EDF3; --muted:#8B949E; --border:#30363D;
  --code-bg:#0F141C; --accent:#58A6FF; --accent-soft:#16304D; --num:#C9D1D9;
  --out-bg:#0F141C;
  --out-fg:#5EEAD4; --out-bd:#155E56; --out-soft:#0C2925;
  color-scheme: dark;
}
/* Embedded plots are transparent SVGs with dark axes, so back them with white in dark
   mode to keep them legible while the card around them stays dark. */
[data-theme="dark"] .guide-fig img { background:#FFF; border-radius:6px; padding:6px; }
* { box-sizing:border-box; }
html { scroll-behavior:smooth; }
body {
  margin:0; background:var(--bg); color:var(--text);
  font-family:-apple-system,BlinkMacSystemFont,"Segoe UI",Roboto,Helvetica,Arial,sans-serif;
  font-size:16px; line-height:1.6; -webkit-font-smoothing:antialiased;
  transition:background-color 200ms,color 200ms;
}
code, pre, .mono { font-family:ui-monospace,"SF Mono","JetBrains Mono",Menlo,Consolas,monospace; }

.sidebar {
  position:fixed; top:0; left:0; width:var(--sidebar-w); height:100vh; overflow-y:auto;
  border-right:1px solid var(--border); background:var(--surface); padding:20px 16px 48px;
}
main { margin-left:var(--sidebar-w); padding:48px 40px 120px; max-width:calc(var(--measure) + 80px); }
.wrap { max-width:var(--measure); }

.brand-row { display:flex; align-items:center; justify-content:space-between; gap:8px; margin:0 0 4px; }
.brand { font-weight:700; font-size:17px; letter-spacing:-0.01em; }
.theme-toggle {
  flex:none; border:1px solid var(--border); background:var(--bg); color:var(--muted);
  font:inherit; font-size:11px; font-weight:600; text-transform:uppercase; letter-spacing:.04em;
  padding:4px 9px; border-radius:6px; cursor:pointer;
  transition:background-color 180ms,color 180ms,border-color 180ms;
}
.theme-toggle:hover { color:var(--text); border-color:var(--accent); }
/* Keyboard focus. The sidebar is the page's main navigation and is all links, so they
   need a visible ring; :focus-visible keeps it off mouse clicks. */
a:focus-visible, .theme-toggle:focus-visible, .nav-fn:focus-visible, .nav-sub:focus-visible {
  outline:2px solid var(--accent); outline-offset:2px; border-radius:4px;
}
.search input {
  width:100%; padding:8px 10px; font-size:14px; color:var(--text); border:1px solid var(--border);
  border-radius:8px; background:var(--bg); margin:12px 0 18px; outline:none;
}
.search input:focus-visible { border-color:var(--accent); box-shadow:0 0 0 3px var(--accent-soft); }
.nav-h {
  display:block; font-size:12px; font-weight:600; text-transform:uppercase; letter-spacing:.06em;
  color:var(--muted); margin:18px 0 6px; text-decoration:none;
}
.nav-fn, .nav-sub {
  display:block; padding:3px 8px; margin:1px 0; border-radius:6px; font-size:13.5px;
  color:var(--text); text-decoration:none; white-space:nowrap; overflow:hidden;
  text-overflow:ellipsis; transition:background-color 180ms,color 180ms;
}
.nav-fn { font-family:ui-monospace,"SF Mono",Menlo,Consolas,monospace; }
.nav-fn.child { padding-left:20px; font-size:12.5px; color:var(--muted); }
.nav-fn:hover, .nav-sub:hover { background:var(--code-bg); }
.nav-fn.active, .nav-sub.active { background:var(--accent-soft); color:var(--accent); }

h1 { font-size:34px; font-weight:700; letter-spacing:-0.02em; margin:0 0 6px; }
.lede { color:var(--muted); font-size:17px; margin:0 0 8px; }
.version { color:var(--muted); font-size:14px; }
h2 {
  font-size:13px; font-weight:700; text-transform:uppercase; letter-spacing:.08em;
  color:var(--muted); margin:56px 0 20px; padding-bottom:8px; border-bottom:1px solid var(--border);
}
h3 { font-size:22px; font-weight:650; letter-spacing:-0.01em; margin:40px 0 16px; }
a { color:var(--accent); }

.guide-block { margin:0 0 56px; }
.guide-block h3 { margin-top:28px; }
.guide-block > p { margin:0 0 18px; }
/* The code+media row breaks out past the prose measure (which stays narrow to read) so
   both columns get real room. It shares the left edge with the text and only extends
   right into the empty margin. Stacks on narrow screens.
   align-items is stretch, not flex-start: a flex-start column shrink-wraps to its
   content's intrinsic width, and a wide table inside an overflow-x:auto box still
   reports its full unwrapped width, which dragged the whole page to 1221px. */
.guide-cols { display:flex; gap:32px; align-items:stretch; width:min(1180px, calc(100vw - 320px)); }
.guide-cols > .guide-code { flex:0 1 auto; min-width:0; max-width:52%; }
.guide-cols > .guide-figs { flex:1 1 0; min-width:0; }
.guide-code pre { margin:0; font-size:12.5px; line-height:1.5; }
.guide-stack { width:min(1180px, calc(100vw - 320px)); }
/* Wrap rather than squeeze: three artifacts in one row would clip the widest table. */
.guide-figs--row { display:flex; flex-wrap:wrap; gap:28px; align-items:stretch; }
.guide-figs--row > .io { flex:1 1 340px; min-width:0; max-width:100%; margin:26px 0 0; }
@media (max-width:980px) {
  .guide-cols { width:auto; flex-direction:column; gap:26px; }
  .guide-cols > .guide-code, .guide-cols > .guide-figs { max-width:100%; width:100%; }
  .guide-stack { width:auto; }
  .guide-figs--row { flex-direction:column; }
}
/* The badge has to escape the card to sit on its corner, so it hangs off this wrapper
   while the card itself keeps clipping its own wide content. */
.io { position:relative; min-width:0; max-width:100%; }
.guide-fig {
  margin:0; padding:16px; background:var(--surface); border:1px solid var(--border);
  border-radius:10px; overflow:hidden; max-width:100%; height:100%;
}

/* In/Out badges. Each guide section is one notebook-style cell: the code block is In n,
   and the artifacts it produces are Out n (lettered when there is more than one), so the
   prose can point at a specific block. */
.io-badge {
  position:absolute; top:-10px; left:-10px; z-index:2;
  font-family:ui-monospace,"SF Mono",Menlo,Consolas,monospace;
  font-size:10.5px; font-weight:700; letter-spacing:.04em; white-space:nowrap;
  padding:2px 8px; border-radius:999px; border:1px solid var(--border);
  background:var(--surface); color:var(--muted);
}
.io-in > .io-badge { color:var(--accent); border-color:var(--accent); background:var(--accent-soft); }
.io-out > .io-badge { color:var(--out-fg); border-color:var(--out-bd); background:var(--out-soft); }
/* Inline reference to a numbered block, e.g. "the manifest in Out 6b". */
.io-ref {
  font-family:ui-monospace,"SF Mono",Menlo,Consolas,monospace; font-size:.82em;
  font-weight:700; padding:1px 6px; border-radius:999px; white-space:nowrap;
  border:1px solid var(--out-bd); color:var(--out-fg); background:var(--out-soft);
}
.io-ref.in { border-color:var(--accent); color:var(--accent); background:var(--accent-soft); }
.io-legend {
  margin:0 0 28px; padding:12px 16px; border:1px solid var(--border); border-radius:10px;
  background:var(--surface); color:var(--muted); font-size:14.5px; line-height:1.6;
}
.guide-figs > .io { margin:0 0 26px; }
.guide-figs > .io:last-child { margin-bottom:0; }
/* Stacked layout: the code block sits directly above the first card, and the badge hangs
   10px above it, so leave clearance rather than letting them collide. */
.guide-stack > .guide-code { margin-bottom:4px; }
.guide-fig img { display:block; width:100%; height:auto; }
.guide-fig figcaption { margin:12px 2px 0; color:var(--muted); font-size:13.5px; line-height:1.5; }

.df-preview { overflow-x:auto; }
table.df { width:auto; border-collapse:collapse; font-size:12.5px; }
table.df thead th {
  text-align:left; font-weight:600; color:var(--muted); font-size:10.5px; text-transform:uppercase;
  letter-spacing:.04em; padding:0 12px 9px; border-bottom:1px solid var(--border); white-space:nowrap;
}
table.df tbody td {
  padding:7px 12px; border-bottom:1px solid var(--border); white-space:nowrap; vertical-align:middle;
}
table.df tbody tr:last-child td { border-bottom:none; }
table.df tbody tr:hover td { background:var(--accent-soft); }
table.df th.num, table.df td.num { text-align:right; }
table.df td.num, table.df td:last-child {
  font-family:ui-monospace,"SF Mono",Menlo,Consolas,monospace;
  font-variant-numeric:tabular-nums; color:var(--num);
}
table.df td.txt:first-child { font-weight:600; }

pre {
  background:var(--code-bg); border:1px solid var(--border); border-radius:8px;
  padding:14px 16px; overflow-x:auto; font-size:13.5px; line-height:1.55; margin:0 0 8px;
}
pre.code { background:transparent; }
pre.out-text {
  background:var(--out-bg); font-size:12.5px; line-height:1.5; white-space:pre;
  overflow-x:auto; margin:0;
}
pre.doc-code { background:transparent; font-size:12.5px; margin:8px 0 12px; }
p code, li code, figcaption code, dd code {
  background:var(--code-bg); padding:1px 5px; border-radius:4px; font-size:.9em;
}
pre code { background:none; padding:0; border-radius:0; font-size:inherit; }

/* Light-mode token colors are darkened from the usual palette so each one clears WCAG AA
   (4.5:1) against --code-bg. The obvious picks (#6B7280, #B45309, #2563EB, #DB2777) all
   land between 3.9 and 4.5, which is under the bar for 13.5px text. */
.tok-c { color:#5A6270; font-style:italic; }
.tok-s { color:#047857; }
.tok-n { color:#9A4708; }
.tok-k { color:#7C3AED; }
.tok-f { color:#1D4ED8; }
.tok-b { color:#BE185D; font-weight:600; }
[data-theme="dark"] .tok-c { color:#8B949E; }
[data-theme="dark"] .tok-s { color:#7EE787; }
[data-theme="dark"] .tok-n { color:#FFA657; }
[data-theme="dark"] .tok-k { color:#D2A8FF; }
[data-theme="dark"] .tok-f { color:#79C0FF; }
[data-theme="dark"] .tok-b { color:#FF7B9C; }

.fn-card {
  background:var(--surface); border:1px solid var(--border); border-radius:10px;
  padding:18px 20px; margin:0 0 14px; scroll-margin-top:24px;
}
.cls-card { border-left:3px solid var(--accent); }
.mem-card { margin-left:18px; }
.fn-sig { display:block; font-size:14px; overflow-x:auto; margin:0 0 10px; padding-bottom:2px; }
.fn-sig .fn-name { font-weight:700; color:var(--text); }
.fn-sig .fn-owner { color:var(--muted); }
.fn-sig .fn-args { color:var(--muted); }
.cls-kw { color:var(--accent); font-weight:600; }
.badge {
  display:inline-block; font-size:10px; font-weight:700; text-transform:uppercase;
  letter-spacing:.05em; color:var(--accent); background:var(--accent-soft);
  padding:2px 7px; border-radius:5px; margin:0 0 8px;
}
.fn-summary { margin:0 0 4px; font-size:15.5px; }
.doc-section {
  font-size:11px; font-weight:700; text-transform:uppercase; letter-spacing:.07em;
  color:var(--muted); margin:20px 0 10px;
}
.doc-list { margin:6px 0; padding-left:20px; font-size:14.5px; }
.doc-list li { margin:4px 0; }

dl.params { margin:4px 0 0; }
dl.params dt { margin:0 0 2px; }
dl.params dd { margin:2px 0 14px; padding:0 0 0 16px; border-left:2px solid var(--border); }
dl.params dd.muted-dd, dl.params dd:empty { border-left-color:transparent; margin-bottom:8px; }
.pname { background:var(--code-bg); padding:1.5px 7px; border-radius:5px; font-weight:700; font-size:13.5px; }
.ptype {
  font-family:ui-monospace,"SF Mono",Menlo,Consolas,monospace; color:var(--muted);
  font-size:12.5px; margin-left:6px;
}
.pdefault { font-family:ui-monospace,Menlo,monospace; color:var(--muted); font-size:12.5px; margin-left:6px; opacity:.75; }
.pdesc { margin:6px 0; font-size:14.5px; line-height:1.55; }
.muted { color:var(--muted); }
.hidden { display:none !important; }

@media (max-width:860px) {
  /* The nav is ~90 links. Left at its natural height it puts 2700px of sidebar in front
     of the first line of content, so cap it and let it scroll in place. */
  .sidebar {
    position:static; width:auto; height:auto; max-height:42vh; overflow-y:auto;
    border-right:none; border-bottom:1px solid var(--border);
    padding:16px 16px 20px;
  }
  main { margin-left:0; padding:28px 20px 80px; }
  .mem-card { margin-left:0; }
  .guide-block { margin-bottom:44px; }
}
@media (prefers-reduced-motion:reduce) {
  html { scroll-behavior:auto; }
  * { transition:none !important; }
}
"""

JS = """
const q = document.getElementById('q');
const cards = Array.from(document.querySelectorAll('.fn-card'));
const navFns = Array.from(document.querySelectorAll('.nav-fn'));
const refGroups = Array.from(document.querySelectorAll('.ref-group'));
const navGroups = Array.from(document.querySelectorAll('.nav-group'));
const guideEls = Array.from(document.querySelectorAll('[data-guide]'));

function applyFilter() {
  const term = q.value.trim().toLowerCase();
  const searching = term.length > 0;
  guideEls.forEach(el => el.classList.toggle('hidden', searching));
  cards.forEach(c => {
    const hay = (c.dataset.name + ' ' + c.dataset.summary).toLowerCase();
    c.classList.toggle('hidden', searching && !hay.includes(term));
  });
  navFns.forEach(a => {
    a.classList.toggle('hidden', searching && !a.dataset.name.toLowerCase().includes(term));
  });
  refGroups.forEach(g => {
    g.classList.toggle('hidden', !g.querySelectorAll('.fn-card:not(.hidden)').length);
  });
  navGroups.forEach(g => {
    g.classList.toggle('hidden', !g.querySelectorAll('.nav-fn:not(.hidden)').length);
  });
}
q.addEventListener('input', applyFilter);

const linkByName = {};
navFns.forEach(a => { linkByName[a.dataset.name] = a; });
const obs = new IntersectionObserver((entries) => {
  entries.forEach(e => {
    const link = linkByName[e.target.dataset.name];
    if (!link) return;
    if (e.isIntersecting) {
      navFns.forEach(a => a.classList.remove('active'));
      link.classList.add('active');
    }
  });
}, { rootMargin: '-10% 0px -80% 0px', threshold: 0 });
cards.forEach(c => obs.observe(c));

// Color theme: Auto (dark during sundown) / Light / Dark, cycled by the sidebar button.
const themeBtn = document.getElementById('theme-toggle');
const THEME_KEY = 'library-designer-theme';
const MODES = ['auto', 'light', 'dark'];
function isSundown() { const h = new Date().getHours(); return h >= 18 || h < 6; }
function themeMode() {
  const m = localStorage.getItem(THEME_KEY);
  return MODES.includes(m) ? m : 'auto';
}
function renderTheme() {
  const mode = themeMode();
  document.documentElement.setAttribute(
    'data-theme', mode === 'auto' ? (isSundown() ? 'dark' : 'light') : mode);
  themeBtn.textContent = mode.charAt(0).toUpperCase() + mode.slice(1);
}
themeBtn.addEventListener('click', () => {
  localStorage.setItem(THEME_KEY, MODES[(MODES.indexOf(themeMode()) + 1) % MODES.length]);
  renderTheme();
});
renderTheme();
setInterval(() => { if (themeMode() === 'auto') renderTheme(); }, 60000);
"""

HEAD_THEME_JS = """
(function () {
  try {
    var m = localStorage.getItem('library-designer-theme');
    var h = new Date().getHours();
    document.documentElement.setAttribute('data-theme',
      (m === 'light' || m === 'dark') ? m : ((h >= 18 || h < 6) ? 'dark' : 'light'));
  } catch (e) {}
})();
"""


def build_html(figures: dict[str, str]) -> str:
    try:
        from importlib.metadata import version
        ver = version("library-designer")
    except Exception:
        ver = getattr(ld, "__version__", "")
    # The package docstring opens "library_designer, design DNA ..."; drop the name and
    # raise the first letter only, so acronyms downstream keep their case.
    lede = (ld.__doc__ or "").strip().split("\n")[0]
    if "," in lede:
        lede = lede.split(",", 1)[1].strip()
        lede = lede[:1].upper() + lede[1:]

    # ---- sidebar ----
    side = ['<aside class="sidebar">']
    side.append('<div class="brand-row"><div class="brand">library-designer</div>'
                '<button id="theme-toggle" class="theme-toggle" type="button" '
                'title="Color theme, Auto follows sundown" aria-label="Color theme">Auto</button></div>')
    side.append('<div class="search"><input id="q" type="search" '
                'placeholder="Search the reference…" aria-label="Search"></div>')
    side.append('<a class="nav-h" href="#guide" data-guide>Guide</a>')
    for gid, title, *_ in GUIDE:
        side.append(f'<a class="nav-sub" href="#{gid}" data-guide>{esc(title)}</a>')
    side.append('<a class="nav-h" href="#reference">Reference</a>')
    for slug, title, classes in REFERENCE:
        side.append(f'<div class="nav-group" data-section="{slug}">')
        side.append(f'<div class="nav-h">{esc(title)}</div>')
        for cls in classes:
            side.append(f'<a class="nav-fn" data-name="{esc(cls.__name__)}" '
                        f'href="#{esc(cls.__name__)}">{esc(cls.__name__)}</a>')
            for m in members_of(cls):
                anchor = f"{cls.__name__}.{m['name']}"
                side.append(f'<a class="nav-fn child" data-name="{esc(anchor)}" '
                            f'href="#{esc(anchor)}">{esc(m["name"])}</a>')
        side.append("</div>")
    side.append("</aside>")

    # ---- main ----
    m = ['<main><div class="wrap">', "<header>", "<h1>library-designer</h1>"]
    if lede:
        m.append(f'<p class="lede">{esc(lede)}</p>')
    if ver:
        m.append(f'<p class="version">v{esc(ver)}</p>')
    m.append("</header>")

    m.append('<h2 id="guide" data-guide>Guide</h2>')
    m.append(
        '<p class="io-legend" data-guide>Each section below reads as one notebook cell. '
        '<span class="io-ref in">In 1</span> is the code you run and '
        '<span class="io-ref">Out 1</span> is what it produces, lettered when one cell '
        'produces more than one thing. Every output on this page was generated by running '
        'the code above it against the bundled examples.</p>'
    )
    for n, (gid, title, intro, code, media) in enumerate(GUIDE, 1):
        m.append(f'<section class="guide-block" id="{gid}" data-guide>')
        m.append(f"<h3>{esc(title)}</h3>")
        m.append(f"<p>{intro}</p>")
        code_html = (f'<div class="io io-in"><span class="io-badge">In {n}</span>'
                     f'<pre class="code"><code>{highlight_code(code)}</code></pre></div>')

        # Number the artifacts that actually rendered, so a figure that failed to build
        # cannot leave a gap in the lettering or a dangling reference in the prose.
        present = [(k, key, cap) for k, key, cap in media if figures.get(key)]
        figs = []
        for i, (kind, key, caption) in enumerate(present):
            content = figures[key]
            label = f"{n}{chr(ord('a') + i)}" if len(present) > 1 else str(n)
            body = (f'<img src="data:image/svg+xml;base64,{content}" alt="{esc(title)}">'
                    if kind == "img" else content)
            figs.append(f'<div class="io io-out" id="out-{label}">'
                        f'<span class="io-badge">Out {label}</span>'
                        f'<figure class="guide-fig">{body}'
                        f"<figcaption>{caption}</figcaption></figure></div>")

        if len(figs) > 1:
            m.append('<div class="guide-stack">')
            m.append(f'<div class="guide-code">{code_html}</div>')
            m.append(f'<div class="guide-figs guide-figs--row">{"".join(figs)}</div>')
            m.append("</div>")
        elif figs:
            m.append('<div class="guide-cols">')
            m.append(f'<div class="guide-code">{code_html}</div>')
            m.append(f'<div class="guide-figs">{"".join(figs)}</div>')
            m.append("</div>")
        else:
            m.append(code_html)
        m.append("</section>")

    m.append('<h2 id="reference">Reference</h2>')
    n_cards = 0
    for slug, title, classes in REFERENCE:
        m.append(f'<div class="ref-group" data-section="{slug}"><h3>{esc(title)}</h3>')
        for cls in classes:
            m.append(class_card(cls))
            n_cards += 1
            for mem in members_of(cls):
                m.append(member_card(cls, mem))
                n_cards += 1
        m.append("</div>")
    m.append("</div></main>")

    page = f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>library-designer reference</title>
<script>{HEAD_THEME_JS}</script>
<style>{CSS}</style>
</head>
<body>
{''.join(side)}
{''.join(m)}
<script>{JS}</script>
</body>
</html>
"""
    return page, n_cards


def main() -> int:
    print("Rendering guide artifacts from the bundled examples...")
    try:
        from _examples import render_all
        with tempfile.TemporaryDirectory() as d:
            figures = render_all(Path(d))
    except Exception as exc:                       # noqa: BLE001
        print(f"  (no guide artifacts: {type(exc).__name__}: {exc})")
        figures = {}

    page, n_cards = build_html(figures)
    out = Path(__file__).resolve().parent / "index.html"
    out.write_text(page, encoding="utf-8")
    kb = len(page.encode()) / 1024
    print(f"Wrote {out}  ({n_cards} cards, {len(figures)} artifacts, {kb:.0f} KB)")

    # An artifact that failed to render takes its Out badge with it, which would leave any
    # prose pointing at it dangling. Catch that here rather than shipping a broken pointer.
    badges = set(re.findall(r'<span class="io-badge">((?:In|Out) [^<]+)</span>', page))
    refs = {r for r in re.findall(r"<span class='io-ref(?: in)?'>([^<]+)</span>", page)}
    dangling = sorted(refs - badges)
    if dangling:
        print(f"  WARNING dangling reference(s) in the prose: {', '.join(dangling)}")
    else:
        print(f"  {len(badges)} numbered blocks, every prose reference resolves")

    undocumented = [
        f"{cls.__name__}.{mem['name']}"
        for _, _, classes in REFERENCE for cls in classes
        for mem in members_of(cls) if not mem["doc"]
    ] + [cls.__name__ for _, _, classes in REFERENCE for cls in classes if not real_doc(cls)]
    if undocumented:
        print(f"  {len(undocumented)} undocumented: {', '.join(sorted(undocumented))}")
    else:
        print("  every symbol in the reference has a docstring")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
