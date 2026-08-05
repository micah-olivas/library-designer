"""What the printed summary and QC report say.

``print(lib.summary())`` is the block a reviewer reads before sending an order, so the
lines it keeps and the ones it drops are behavior, not cosmetics. A line that carries no
information (an empty adaptor pair, an unset GC bound, a "no hit(s)" per pattern) pushes
the lines that matter off the screen.
"""
from __future__ import annotations

from pathlib import Path

import pytest

from library_designer import (
    CodonOptimizationParams,
    LibrarySpec,
    SequenceSet,
    SubstitutionScan,
    TiledAssemblyParams,
)

REPO = Path(__file__).resolve().parents[1]
PROTEIN = "AEGNTLISVDYEIFGKVQGVFFRKHTQAEGKKLGLVGWVQNTDRGTVQGQLQGPISKVRHMQEWLETRG"


def _spec(**kw) -> LibrarySpec:
    kw.setdefault("name", "demo")
    kw.setdefault("protein_sequence", PROTEIN)
    kw.setdefault("substitutions", ["A"])
    return LibrarySpec(**kw)


@pytest.fixture(scope="module")
def scan():
    spec = _spec(adaptor_5="ggtctccaagc", adaptor_3="ggtgggagacc",
                 platform="twist_oligo_pools")
    return SubstitutionScan(spec).generate().codon_optimize()


# --- the header ---------------------------------------------------------------

def test_header_says_all_optimized_when_nothing_failed(scan):
    head = str(scan.summary()).splitlines()[0]
    assert head.endswith("variants, all optimized")


def test_header_says_not_optimized_before_codon_optimize():
    lib = SubstitutionScan(_spec()).generate()
    head = str(lib.summary()).splitlines()[0]
    assert head.endswith("variants, not optimized yet")


# --- lines that only appear when they say something ---------------------------

def test_empty_adaptors_give_a_construct_line_and_no_adaptor_line():
    lib = SubstitutionScan(_spec()).generate().codon_optimize()
    text = str(lib.summary())
    assert "adaptors:" not in text
    assert re_line(text, "  construct:") is not None


def test_adaptors_are_printed_unquoted_with_the_construct_length(scan):
    # Labels are padded to one width and the parts separated by spaces rather than commas, so
    # a block of rows reads as a column.
    line = re_line(str(scan.summary()), "  adaptors:")
    assert line.split() == [
        "adaptors:", "5'", "ggtctccaagc", "(11", "bp)",
        "3'", "ggtgggagacc", "(11", "bp)", "construct", "229", "bp",
    ]


def test_one_construct_length_is_not_printed_as_a_range(scan):
    # Every member of a substitution scan is the same length, so "229-229 bp" is noise.
    assert "-229 bp" not in str(scan.summary())


def test_sequence_set_drops_the_single_members_bucket():
    lib = SequenceSet(_spec(), {"d1": PROTEIN, "d2": PROTEIN[:60]}).generate()
    assert "sublibraries:" not in str(lib.summary())


def test_optimization_line_keeps_only_what_was_set(scan):
    assert re_line(str(scan.summary()), "  optimization:").split(":", 1)[1].strip() == (
        "e_coli, use_best_codon"
    )


def test_optimization_line_reports_gc_and_iters_when_they_apply():
    spec = _spec(optimization=CodonOptimizationParams(
        method="match_codon_usage", gc_min=0.35, gc_max=0.65, max_random_iters=2000))
    line = re_line(str(SubstitutionScan(spec).generate().summary()), "  optimization:")
    assert line.split(":", 1)[1].strip() == (
        "e_coli, match_codon_usage, gc 0.35-0.65, 2000 iters"
    )


def test_starting_vector_is_named_by_file_not_by_path(tmp_path):
    plasmid = _genbank(tmp_path)
    spec = _spec(starting_vector=str(plasmid))
    lib = SubstitutionScan(spec).generate()
    assert re_line(str(lib.summary()), "  vector:").split(":", 1)[1].strip() == (
        "demo_backbone.gb"
    )


# --- the QC block -------------------------------------------------------------

def test_clean_enzyme_and_motif_checks_collapse_to_one_line(scan):
    text = str(scan.summary())
    assert "no hit(s)" not in text
    line = re_line(text, "    forbidden sequences:")
    assert line.split(":", 1)[1].strip() == "none (checked BsaI, 2 motifs)"


def test_qc_header_in_the_summary_does_not_repeat_the_variant_count(scan):
    text = str(scan.summary())
    assert re_line(text, "  QC report: ") == "  QC report: PASS"
    # On its own the report has no other line giving the count, so it keeps it.
    assert str(scan.check()).splitlines()[0] == (
        f"QC report: {len(scan)} variants, PASS"
    )


def test_failures_are_named_on_the_same_line():
    report = _report_with_failures(["I7F", "S8Y", "V9A", "D10F", "Y11M", "E12A"])
    line = re_line(str(report), "  translation:")
    assert line.endswith("(I7F, S8Y, V9A, D10F, and 2 more)")
    assert "    x " not in str(report)


# --- tiled layout replaces the (empty) adaptor pair ---------------------------

def test_tiled_summary_reports_tiles_and_oligo_lengths():
    spec = LibrarySpec.from_toml(REPO / "examples" / "gck_tiled.toml")
    lib = SubstitutionScan(spec).generate().codon_optimize().tile()
    s = lib.summary()
    text = str(s)
    assert s.tiles["n_tiles"] == len(lib.tiles)
    assert s.tiles["cds_len"] == len(lib.reference)
    line = re_line(text, "  tiles:")
    assert line.split(":", 1)[1].strip().startswith(
        f"{len(lib.tiles)} over a {len(lib.reference)} bp CDS"
    )
    assert f"primer set {lib.tiled_params.primer_set}" in line
    assert "adaptors:" not in text and "construct" not in text


def test_untiled_library_has_no_tiles_line():
    lib = SubstitutionScan(_spec()).generate()
    assert lib.summary().tiles is None
    assert "  tiles:" not in str(lib.summary())


def test_tiles_line_survives_a_library_tiled_with_explicit_params():
    # tile(params) with no [tiled] block on the spec: the summary must read the params off
    # the library, the same way QC does.
    lib = SubstitutionScan(_spec()).generate().codon_optimize()
    lib = lib.tile(TiledAssemblyParams(oligo_budget=200))
    assert lib.summary().tiles["primer_set"] == lib.tiled_params.primer_set


# --- helpers ------------------------------------------------------------------

def re_line(text: str, prefix: str) -> str | None:
    """The one line of ``text`` starting with ``prefix``, or None."""
    hits = [ln for ln in text.splitlines() if ln.startswith(prefix)]
    assert len(hits) <= 1, f"{prefix!r} matched {len(hits)} lines"
    return hits[0] if hits else None


def _report_with_failures(names: list[str]):
    from library_designer.checks import CheckReport

    return CheckReport(n_variants=100, translation_fail=names)


def _genbank(tmp_path: Path) -> Path:
    """A minimal circular plasmid carrying the scan's CDS as its sole CDS feature."""
    from Bio import SeqIO
    from Bio.Seq import Seq
    from Bio.SeqFeature import SeqFeature, SimpleLocation
    from Bio.SeqRecord import SeqRecord

    lib = SubstitutionScan(_spec()).generate().codon_optimize()
    flank = "TTGACAATTAATCATCGGCTCGTATAATG"
    seq = flank + lib.reference + flank
    rec = SeqRecord(Seq(seq), id="demo", name="demo", description="test backbone",
                    annotations={"molecule_type": "DNA", "topology": "circular"})
    rec.features = [
        SeqFeature(SimpleLocation(len(flank), len(flank) + len(lib.reference), strand=1),
                   type="CDS", qualifiers={"label": ["insert"]})
    ]
    path = tmp_path / "demo_backbone.gb"
    SeqIO.write(rec, str(path), "genbank")
    return path


# --- reading the report in code -----------------------------------------------

def test_the_report_echoes_as_the_readable_block_not_a_field_dump(scan):
    """A notebook echoing ``lib.check()`` used to print every field, empty ones included.
    The repr is the same block ``print()`` gives."""
    r = repr(scan.check())
    assert r == str(scan.check())
    assert r.startswith("QC report: ")
    assert "n_variants=" not in r and "motif_hits=" not in r
    # And it renders as the same text in a notebook.
    assert scan.check()._repr_html_().startswith("<pre")


def test_issues_is_empty_exactly_when_the_report_passes(scan):
    rep = scan.check()
    assert rep.passed and rep.issues == {}

    failing = _report_with_failures(["I7F", "S8Y"])
    assert not failing.passed
    assert failing.issues == {"translation_fail": ["I7F", "S8Y"]}


def test_issues_names_the_enzyme_or_motif_that_was_hit():
    """The two per-pattern checks get one key each, so a script can tell which enzyme
    or which motif was hit without re-deriving it."""
    spec = LibrarySpec(name="j", protein_sequence="MAAK", substitutions=["GGT"],
                       adaptor_3="CTCAAA", avoid_enzymes=["BsaI"], avoid_patterns=[],
                       seed=1)
    rep = SubstitutionScan(spec).generate().codon_optimize().check()
    assert rep.issues == {"enzyme_hits:BsaI": ["K4G"]}
    assert not rep.passed


def test_an_assembly_count_mismatch_shows_up_as_an_issue():
    from library_designer.checks import CheckReport

    rep = CheckReport(n_variants=10, assembly_checked=10, assembly_correct=8)
    assert not rep.passed
    (msg,) = rep.issues["assembly_incorrect"]
    assert msg == "2 of 10 members do not rebuild their intended variant"


def test_advisories_join_both_informational_lists_and_never_fail_the_report():
    from library_designer.checks import CheckReport

    rep = CheckReport(n_variants=1, reference_advisories=["kept a motif"],
                      overhang_advisories=["two ends share 3 bases"])
    assert rep.advisories == ["kept a motif", "two ends share 3 bases"]
    assert rep.passed and rep.issues == {}      # advisories are not failures


def test_check_fmt_returns_text_or_plain_data(scan):
    import json

    assert scan.check(fmt="text") == scan.check().text()

    d = scan.check(fmt="dict")
    assert isinstance(d, dict)
    assert d["passed"] is True and d["issues"] == {} and d["advisories"] == []
    assert d["n_variants"] == len(scan)
    assert json.loads(json.dumps(d)) == d       # straight into a run record

    with pytest.raises(ValueError, match="Unknown fmt"):
        scan.check(fmt="parsed")
