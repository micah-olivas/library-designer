"""Build the artifacts the guide shows beside its code blocks.

Everything here runs the real pipeline on the bundled examples in ``examples/``, so the
docs build doubles as a smoke test: if a guide figure stops rendering, the workflow it
documents is broken. Nothing is faked and nothing is cached.

Each artifact is produced by its own function and called through ``_safe``, so one
failure costs one figure rather than the whole page. Plots come back as base64 SVG
(vector, so they stay crisp at any zoom) and tables as ready-made HTML.
"""
from __future__ import annotations

import base64
import io as _io
import textwrap
from html import escape
from pathlib import Path

import matplotlib

matplotlib.use("Agg")   # headless: the build must not need a display
# Matplotlib salts the ids it gives SVG clip paths and glyphs with a fresh uuid per figure
# unless this is set, so the same figure came out different on every build. Any fixed string
# does; this one names where it is used.
matplotlib.rcParams["svg.hashsalt"] = "library-designer-docs"

import matplotlib.pyplot as plt   # noqa: E402

ROOT = Path(__file__).resolve().parent.parent
EXAMPLES = ROOT / "examples"

# A synthetic backbone for the destination-vector guide sections. The four bases each arm
# contributes are the fused overhangs the mbo038 adaptors spell (AAGC and GGTG), so the
# digested oligo and the cut vector present the same sticky ends and the whole CDS drops
# out. Both arms are BsaI-clean, checked below.
BB5 = "ACGTACGTTTGCAACGGATCCACAAGC"     # ...ends with the 5' fused overhang
BB3 = "GGTGCCTAGGCATTACGTACGTACGT"      # starts with the 3' fused overhang


def _svg(fig) -> str:
    """A matplotlib figure as a base64 SVG, so it embeds in the page with no side files.

    ``Date`` is dropped from the SVG metadata. Matplotlib writes the build time into every
    figure otherwise, so each of these blobs changed on every build and the generated page
    never matched the committed one.
    """
    buf = _io.BytesIO()
    fig.savefig(buf, format="svg", bbox_inches="tight", transparent=True,
                metadata={"Date": None})
    plt.close(fig)
    return base64.b64encode(buf.getvalue()).decode("ascii")


def _text_block(text: str, limit: int = 24) -> str:
    """A monospace panel, for report output that is already laid out as text."""
    lines = text.splitlines()
    clipped = lines[:limit] + ([f"... ({len(lines) - limit} more lines)"]
                               if len(lines) > limit else [])
    return f'<pre class="out-text">{escape(chr(10).join(clipped))}</pre>'


def df_preview_html(df, numeric: set[str] | None = None, maxlen: int = 28) -> str:
    """Render a DataFrame as a compact HTML table.

    Text columns left-align, numeric columns right-align in tabular-figure monospace.
    Long sequence cells are middle-elided so one wide column cannot push the table past
    the page.
    """
    numeric = numeric or set()

    def fmt(v) -> str:
        s = "" if v is None else str(v)
        if s in ("<NA>", "nan", "None"):
            return ""
        if len(s) > maxlen:
            keep = (maxlen - 3) // 2
            return f"{s[:keep]}...{s[-keep:]}"
        return s

    cols = list(df.columns)
    head = "".join(
        f'<th class="{"num" if c in numeric else "txt"}">{escape(str(c))}</th>' for c in cols
    )
    body = "".join(
        "<tr>" + "".join(
            f'<td class="{"num" if c in numeric else "txt"}">{escape(fmt(row[c]))}</td>'
            for c in cols
        ) + "</tr>"
        for _, row in df.iterrows()
    )
    return ('<div class="df-preview"><table class="df">'
            f"<thead><tr>{head}</tr></thead><tbody>{body}</tbody></table></div>")


def _ruler(lo: int, hi: int, column, step: int = 10) -> tuple[str, str]:
    """A coordinate line and a tick line for bases ``lo`` to ``hi - 1``, both 0-based.

    ``column(i)`` says which screen column base ``i`` prints at, so the ruler stays honest
    about which character is which base even when the sequence line has brackets in it.
    Ticks land on the 1-based coordinates that are multiples of ``step``. A number ends at
    its own tick, or starts there when ending there would run it off the left of the line.
    """
    nums: list[str] = []
    ticks: list[str] = []

    def place(row: list[str], text: str, at: int) -> None:
        row += " " * (at + len(text) - len(row))
        row[at:at + len(text)] = text

    for i in range(lo, hi):
        pos = i + 1
        if pos % step:
            continue
        col = column(i)
        place(ticks, "|", col)
        label = str(pos)
        start = col - len(label) + 1
        place(nums, label, start if start >= column(lo) else col)
    return "".join(nums), "".join(ticks)


# --------------------------------------------------------------------------- #
# The three worked libraries. Built once each and reused across the artifacts.
# --------------------------------------------------------------------------- #
_CACHE: dict = {}


def scan_library():
    """The mbo038 hAcyP1 amber + missense scan, codon-optimized."""
    if "scan" not in _CACHE:
        from library_designer import LibrarySpec, SubstitutionScan

        spec = LibrarySpec.from_toml(EXAMPLES / "mbo038.toml")
        _CACHE["scan"] = SubstitutionScan(spec).generate().codon_optimize()
    return _CACHE["scan"]


def tiled_library():
    """The GCK tiled-assembly example, codon-optimized and tiled."""
    if "tiled" not in _CACHE:
        from library_designer import LibrarySpec, SubstitutionScan

        spec = LibrarySpec.from_toml(EXAMPLES / "gck_tiled.toml")
        _CACHE["tiled"] = SubstitutionScan(spec).generate().codon_optimize().tile()
    return _CACHE["tiled"]


def tiled_optimized_library():
    """The same GCK example with the boundary search and oligo padding turned on, which is
    what the overhang section compares against the plain layout."""
    if "tiled_opt" not in _CACHE:
        from dataclasses import replace

        from library_designer import LibrarySpec, SubstitutionScan

        spec = LibrarySpec.from_toml(EXAMPLES / "gck_tiled.toml")
        spec.tiled = replace(spec.tiled, optimize_overhangs=True, pad_oligos=True)
        _CACHE["tiled_opt"] = SubstitutionScan(spec).generate().codon_optimize().tile()
    return _CACHE["tiled_opt"]


def vector_library(tmp: Path):
    """The scan again, this time cloned into a synthetic plasmid built around its own
    reference, so the destination vector and the assembly simulation have a real backbone
    to work against."""
    if "vector" not in _CACHE:
        from library_designer import LibrarySpec, StartingVectorParams, SubstitutionScan
        from library_designer.checks.motifs import contains_enzyme_site

        base = scan_library()
        reference = base.reference
        assert not contains_enzyme_site(BB5, "BsaI") and not contains_enzyme_site(BB3, "BsaI")

        plasmid = _write_genbank(tmp / "synthetic_backbone.gb", BB5 + reference + BB3,
                                 cds=(len(BB5), len(BB5) + len(reference)))
        spec = LibrarySpec.from_toml(EXAMPLES / "mbo038.toml")
        spec.cds = reference                      # clone the reference we already built
        spec.starting_vector = StartingVectorParams(path=str(plasmid), insert_label="insert")
        _CACHE["vector"] = SubstitutionScan(spec).generate().codon_optimize()
    return _CACHE["vector"]


def _write_genbank(path: Path, seq: str, cds: tuple[int, int]) -> Path:
    from Bio import SeqIO
    from Bio.Seq import Seq
    from Bio.SeqFeature import SeqFeature, SimpleLocation
    from Bio.SeqRecord import SeqRecord

    rec = SeqRecord(Seq(seq), id="synthetic", name="synthetic",
                    description="synthetic destination backbone for the docs build",
                    annotations={"molecule_type": "DNA", "topology": "circular"})
    rec.features = [
        SeqFeature(SimpleLocation(cds[0], cds[1], strand=1), type="CDS",
                   qualifiers={"label": ["insert"]}),
        SeqFeature(SimpleLocation(0, 12, strand=1), type="promoter",
                   qualifiers={"label": ["T7 promoter"]}),
    ]
    path.parent.mkdir(parents=True, exist_ok=True)
    SeqIO.write(rec, str(path), "genbank")
    return path


# --------------------------------------------------------------------------- #
# Artifacts
# --------------------------------------------------------------------------- #
def art_spec_table() -> str:
    return scan_library().spec._repr_html_()


def art_variant_table() -> str:
    lib = scan_library()
    # The adaptor columns are the same on every row, so leave them out of the preview and
    # keep the table narrow enough to read without scrolling.
    cols = ["name", "position", "wt_residue", "mut_residue", "codon", "mut_index",
            "variable_dna"]
    return df_preview_html(lib.df[cols].head(6), numeric={"position", "mut_index"})


def art_summary_text() -> str:
    return _text_block(str(scan_library().summary()), limit=16)


def art_qc_text() -> str:
    return _text_block(str(scan_library().check()), limit=12)


def art_codon_usage() -> str:
    return _svg(scan_library().plot_codon_usage())


def art_tiling() -> str:
    return _svg(tiled_library().plot_tiling())


def art_oligo_pool(tmp: Path) -> str:
    import pandas as pd

    lib = tiled_library()
    p = tmp / "oligos.csv"
    lib.to_oligo_pool(p)
    return df_preview_html(pd.read_csv(p).head(5), maxlen=28)


def art_primer_order(tmp: Path) -> str:
    import pandas as pd

    lib = tiled_library()
    p = tmp / "primers.csv"
    lib.to_primer_order(p)
    df = pd.read_csv(p, header=None, names=["name", "sequence", "scale", "purification"])
    # scale and purification are the same on every row, so leave them out of the preview.
    return df_preview_html(df[["name", "sequence"]].head(6), maxlen=26)


def art_vector_table(tmp: Path) -> str:
    import pandas as pd

    lib = vector_library(tmp)
    p = tmp / "vector.csv"
    lib.to_vectors(p)
    df = pd.read_csv(p)
    keep = [c for c in df.columns if c != "vector_sequence"]
    return df_preview_html(df[keep], numeric={
        "cds_dropout_start", "cds_dropout_end", "vector_length",
        "origin_in_starting_vector", "insert_strand_in_starting_vector"})


def art_vector_features(tmp: Path) -> str:
    """The feature table of the emitted destination map, which is what you open in a
    plasmid viewer."""
    from Bio import SeqIO

    lib = vector_library(tmp)
    d = tmp / "vectormaps"
    lib.to_vector_maps(d)
    rec = SeqIO.read(str(d / "destination.gb"), "genbank")
    rows = [
        (f.qualifiers.get("label", [f.type])[0], f.type,
         f"{int(f.location.start) + 1}-{int(f.location.end)}",
         "+" if f.location.strand != -1 else "-")
        for f in rec.features
    ]
    head = "".join(f"<th>{escape(h)}</th>" for h in ("label", "type", "span", "strand"))
    body = "".join("<tr>" + "".join(f"<td>{escape(c)}</td>" for c in r) + "</tr>"
                   for r in rows)
    return ('<div class="df-preview"><table class="df">'
            f"<thead><tr>{head}</tr></thead><tbody>{body}</tbody></table></div>")


def art_assembly_text(tmp: Path) -> str:
    lib = vector_library(tmp)
    report = lib.check()
    lines = [ln for ln in str(report).splitlines()
             if "assembly" in ln or "aligned" in ln or ln.startswith("QC report")]
    return _text_block("\n".join(lines) or str(report), limit=10)


def art_clone_diff(tmp: Path) -> str:
    """The one codon that separates a clone from the parent, which is the whole claim a
    single-mutant library rests on."""
    lib = vector_library(tmp)
    parent = lib.parent_vector()
    # A variant from the middle of the CDS, so the diff has context on both sides.
    placed = lib.df[lib.df["mut_index"].notna()]
    name = str(placed["name"].iloc[len(placed) // 2])
    clone = lib.assembled_product(name)
    diffs = [i for i, (a, b) in enumerate(zip(parent, clone)) if a != b]
    lo, hi = min(diffs), max(diffs) + 1
    pad = 21
    a, b = max(0, lo - pad), min(len(parent), hi + pad)
    marker = "".join("^" if a + i in diffs else " " for i in range(b - a))
    gutter = len("parent  ")

    def column(i: int) -> int:
        """Where base ``i`` lands once the gutter and the brackets are counted in."""
        return gutter + (i - a) + (1 if i >= lo else 0) + (1 if i >= hi else 0)

    nums, ticks = _ruler(a, b, column)
    rows = [nums, ticks] if ticks.strip() else []   # a plasmid too short to hold a tick
    rows += [
        f"parent  {parent[a:lo]}[{parent[lo:hi]}]{parent[hi:b]}",
        f"clone   {clone[a:lo]}[{clone[lo:hi]}]{clone[hi:b]}",
        f"        {' ' * (lo - a)} {marker[lo - a:hi - a]}",
    ]
    head = f"{name}: {len(diffs)} base(s) differ, at {lo + 1}-{hi} of {len(parent)} bp"
    return f'<pre class="out-text">{escape(head)}\n\n{escape(chr(10).join(rows))}</pre>'


def art_overhang_matrix() -> str:
    return _svg(tiled_optimized_library().plot_overhangs())


def art_overhang_pairs() -> str:
    lib = tiled_optimized_library()
    pairs = lib.overhang_pairs()[
        ["end_a", "overhang_a", "end_b", "overhang_b", "same_reaction",
         "shared", "shared_flipped", "risk"]
    ]
    return df_preview_html(pairs.head(6), numeric={"shared", "shared_flipped"})


def art_mispriming() -> str:
    """Where the constant flanks of the bundled glucokinase pool anneal but should not."""
    lib = tiled_optimized_library()
    tab = lib.mispriming()[["handle", "region", "where", "position", "strand", "paired",
                            "aligned", "mismatches", "matched", "risk"]]
    return df_preview_html(tab.head(6),
                           numeric={"position", "paired", "aligned", "mismatches"})


def art_boundary_search() -> str:
    """What moving the boundaries buys, on the bundled glucokinase example.

    Counts only the comparisons that can misfire: a tile's two ends against each other, and
    each end against its own reverse complement. Two tiles never share a tube."""
    from library_designer.checks.overhangs import risk_of, self_risk, shared_bases
    from library_designer.layout.boundaries import layout_ends
    from library_designer.layout.tiled import compute_tiles
    from library_designer.regions import reverse_complement

    lib = tiled_optimized_library()
    params = lib.tiled_params
    searched = [(t.start, t.end) for t in lib.tiles]
    lines = []
    for label, windows in (("balanced split", compute_tiles(len(lib.reference), params)),
                           ("after the search", searched)):
        ends = layout_ends(lib.reference, params, windows)
        tally = {"collision": 0, "high": 0, "watch": 0}
        for i, (tile_a, a) in enumerate(ends):
            if self_risk(shared_bases(a, reverse_complement(a)), len(a)) == "collision":
                tally["collision"] += 1                      # palindromic
            for tile_b, b in ends[i + 1:]:
                if tile_a != tile_b:
                    continue                                 # separate tubes, never meet
                worst = max(shared_bases(a, b), shared_bases(a, reverse_complement(b)))
                r = risk_of(worst, len(a))
                if r in tally:
                    tally[r] += 1
        lines += [
            f"{label:>16}  {' '.join(s for _, s in ends)}",
            f"{'':>16}  tiles {[e - s for s, e in windows]}",
            f"{'':>16}  per reaction: {tally['collision']} unusable, "
            f"{tally['high']} at one mismatch, {tally['watch']} above target",
        ]
    padded = sorted({int(x) for x in lib.df["oligo_length"].dropna()})
    lines += ["", f"{'oligos, padded':>16}  {padded} bp, one length across the pool"]
    return f'<pre class="out-text">{escape(chr(10).join(lines))}</pre>'


def art_sequence_set() -> str:
    from library_designer import LibrarySpec, SequenceSet

    spec = LibrarySpec(name="my_designs", platform="pooled")
    lib = (SequenceSet.from_fasta(spec, EXAMPLES / "example_designs.faa")
           .generate().codon_optimize())
    df = lib.df[["name", "protein", "variable_dna"]].head(5)
    return df_preview_html(df, maxlen=34)


# --------------------------------------------------------------------------- #
def _safe(name: str, fn, *args):
    try:
        return fn(*args)
    except Exception as exc:                       # noqa: BLE001 - media are optional
        print(f"  (skipping {name}: {type(exc).__name__}: {exc})")
        return None


def render_all(tmp: Path) -> dict[str, str]:
    """Every guide artifact, keyed by the anchor it belongs to. Missing keys are fine."""
    out = {
        "spec-table": _safe("spec-table", art_spec_table),
        "variant-table": _safe("variant-table", art_variant_table),
        "summary-text": _safe("summary-text", art_summary_text),
        "qc-text": _safe("qc-text", art_qc_text),
        "codon-usage": _safe("codon-usage", art_codon_usage),
        "tiling": _safe("tiling", art_tiling),
        "overhang-matrix": _safe("overhang-matrix", art_overhang_matrix),
        "overhang-pairs": _safe("overhang-pairs", art_overhang_pairs),
        "mispriming": _safe("mispriming", art_mispriming),
        "boundary-search": _safe("boundary-search", art_boundary_search),
        "oligo-pool": _safe("oligo-pool", art_oligo_pool, tmp),
        "primer-order": _safe("primer-order", art_primer_order, tmp),
        "vector-table": _safe("vector-table", art_vector_table, tmp),
        "vector-features": _safe("vector-features", art_vector_features, tmp),
        "assembly-text": _safe("assembly-text", art_assembly_text, tmp),
        "clone-diff": _safe("clone-diff", art_clone_diff, tmp),
        "sequence-set": _safe("sequence-set", art_sequence_set),
    }
    return {k: v for k, v in out.items() if v}


if __name__ == "__main__":
    import tempfile

    with tempfile.TemporaryDirectory() as d:
        got = render_all(Path(d))
    print(f"{len(got)} artifacts: {', '.join(sorted(got))}")
