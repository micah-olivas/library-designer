"""End-to-end golden path for the main non-tiled substitution scan.

Loads ``examples/mbo038.toml`` and runs the full pipeline (protein -> generate ->
codon_optimize -> QC/export), which exercises the DNA Chisel codon-optimization path
that the native-CDS tiled/OMEGA fixtures skip. Also regression-tests the junction
screen: QC judges the *assembled* construct (adaptor + variable + adaptor) against an
assembled-WT baseline, so a restriction site spelled across an adaptor<->variable
junction is caught even when the variable region alone is clean.
"""
from __future__ import annotations

import re
from pathlib import Path

import pandas as pd
import pytest

from library_designer import LibrarySpec, SubstitutionScan
from library_designer.checks.motifs import count_enzyme_sites
from library_designer.checks.translation import translates_to
from library_designer.regions import assemble

REPO = Path(__file__).resolve().parents[1]
MBO038_TOML = REPO / "examples" / "mbo038.toml"


@pytest.fixture(scope="module")
def mbo038():
    spec = LibrarySpec.from_toml(MBO038_TOML)
    return SubstitutionScan(spec).generate().codon_optimize()


# --- codon-optimization round trip -------------------------------------------

def test_every_variant_round_trips_to_intended_protein(mbo038):
    # The main path codon-optimizes a protein (no native CDS) and must not fail.
    assert mbo038.failed == {}
    for protein, dna in zip(mbo038.df["protein"], mbo038.df["variable_dna"]):
        assert isinstance(dna, str)
        assert translates_to(dna, protein)


def test_single_wt_reference_invariant(mbo038):
    # On the *codon-optimized* reference (not a native-CDS fixture): every member is
    # matches the frozen reference except within its own stamped codon.
    ref = mbo038.reference
    saw_mutant = False
    for dna, idx in zip(mbo038.df["variable_dna"], mbo038.df["mut_index"]):
        assert len(dna) == len(ref)
        diffs = [k for k in range(len(ref)) if dna[k] != ref[k]]
        if pd.isna(idx):
            assert diffs == []                        # WT control is the reference
            continue
        codon = range(int(idx) * 3, int(idx) * 3 + 3)
        assert all(k in codon for k in diffs)         # differs only within its codon
        saw_mutant = True
    assert saw_mutant


# --- QC, including the assembled-construct junction screen (fix #1) -----------

def test_qc_passes_and_no_variant_exceeds_assembled_baseline(mbo038):
    spec = mbo038.spec
    a5, a3 = spec.adaptor_5.upper(), spec.adaptor_3.upper()
    baseline = {e: count_enzyme_sites(assemble(a5, mbo038.reference, a3), e)
                for e in spec.avoid_enzymes}
    for dna in mbo038.df["variable_dna"]:
        construct = assemble(a5, dna, a3)
        for e in spec.avoid_enzymes:
            # No variant manufactures a site (junction included) beyond the intended ones.
            assert count_enzyme_sites(construct, e) <= baseline[e]

    rep = mbo038.check()
    assert rep.passed
    assert not rep.optimization_failed and not rep.translation_fail
    assert not any(rep.enzyme_hits.values())
    assert not any(rep.motif_hits.values())
    assert not rep.length_exceeded


def test_qc_flags_restriction_site_at_adaptor_junction():
    """A stamped edge codon can spell a BsaI site with adjacent adaptor bases. The
    variable region alone is clean, so the old variable-only screen missed it; QC on
    the assembled construct must catch it (this is the fix #1 regression)."""
    # Pin GGT (Gly) at every position. adaptor_3 begins "CTC...", so the last codon of
    # the variable region + adaptor_3 spells GGT|CTC = GGTCTC (BsaI) only across the
    # junction, never inside the variable region.
    spec = LibrarySpec(
        name="junction",
        protein_sequence="MAAK",
        substitutions=["GGT"],       # literal codon -> placed verbatim, residue symbol "G"
        adaptor_3="CTCAAA",
        avoid_enzymes=["BsaI"],
        avoid_patterns=[],           # isolate the enzyme check
        seed=1,
    )
    lib = SubstitutionScan(spec).generate().codon_optimize()

    row = lib.df[lib.df["name"] == "K4G"].iloc[0]
    assert isinstance(row["variable_dna"], str)
    # Clean when screened in isolation ...
    assert count_enzyme_sites(row["variable_dna"], "BsaI") == 0
    # ... but the assembled construct carries a junction BsaI site.
    construct = assemble("", row["variable_dna"], spec.adaptor_3.upper())
    assert count_enzyme_sites(construct, "BsaI") == 1

    rep = lib.check()
    assert "K4G" in rep.enzyme_hits["BsaI"]
    assert not rep.passed


# --- exporters ----------------------------------------------------------------

def test_usortm_export_contract(mbo038, tmp_path):
    p = tmp_path / "variants.csv"
    mbo038.to_usortm(p)
    out = pd.read_csv(p)
    assert list(out.columns) == ["name", "sequence"]
    assert not out["name"].astype(str).str.contains(r"[/|>\s]", regex=True).any()

    a5, a3 = mbo038.spec.adaptor_5, mbo038.spec.adaptor_3
    for seq in out["sequence"]:
        assert seq.startswith(a5.lower()) and seq.endswith(a3.lower())
        variable = seq[len(a5):len(seq) - len(a3)]
        assert variable.isupper() and set(variable) <= set("ACGT")


def test_vendor_pooled_schema(mbo038, tmp_path):
    p = tmp_path / "order.csv"
    mbo038.to_vendor(p)   # platform = twist_oligo_pools -> pooled
    out = pd.read_csv(p)
    assert list(out.columns) == ["Pool Name", "Insert Length", "Insert Sequence"]
    assert (out["Pool Name"].astype(str) == mbo038.spec.name).all()
    assert (out["Insert Length"] == out["Insert Sequence"].str.len()).all()


# --- run directories / provenance ---------------------------------------------

_RUN_DIR = re.compile(r"^hAcyP1_\d{8}_\d{6}$")


def test_export_all_writes_a_dated_run_directory(mbo038, tmp_path):
    """Files land in their own dated directory, and the directory, the record inside it,
    and the library all name the same run."""
    import json

    mbo038.export_all(tmp_path / "out", plots=False)

    (run,) = list((tmp_path / "out").iterdir())
    assert _RUN_DIR.match(run.name)
    assert mbo038.output_dir == run
    assert {p.name for p in run.iterdir()} == {
        "variants.csv", "hAcyP1_full_library.csv", "hAcyP1_order.csv",
        "hAcyP1_design_specs.json", "oligos",
    }
    specs = json.loads((run / "hAcyP1_design_specs.json").read_text())
    assert specs["output_dir"] == str(run)
    # The directory stamp is the moment the sequences were built, not the moment they were
    # written, so the name and the record agree.
    assert specs["created"] == mbo038.created
    assert run.name.endswith(mbo038.created[:19].replace("-", "").replace(":", "").replace("T", "_"))


def test_re_exporting_refreshes_one_directory(mbo038, tmp_path):
    mbo038.export_all(tmp_path / "out", plots=False)
    mbo038.export_all(tmp_path / "out", plots=False)
    assert len(list((tmp_path / "out").iterdir())) == 1


def test_a_later_run_gets_its_own_directory(mbo038, tmp_path):
    """Two runs must not share a directory, which is what the stamp is for."""
    from datetime import timedelta

    mbo038.export_all(tmp_path / "out", plots=False)
    first = mbo038.output_dir
    mbo038._created_at += timedelta(hours=1)          # stand in for a later run
    mbo038.export_all(tmp_path / "out", plots=False)
    assert mbo038.output_dir != first
    assert sorted(p.name for p in (tmp_path / "out").iterdir()) == sorted(
        [first.name, mbo038.output_dir.name]
    )
    mbo038._created_at -= timedelta(hours=1)          # leave the module fixture as found


def test_timestamp_false_writes_straight_into_the_named_directory(mbo038, tmp_path):
    mbo038.export_all(tmp_path / "flat", plots=False, timestamp=False)
    assert (tmp_path / "flat" / "variants.csv").is_file()
    assert mbo038.output_dir == tmp_path / "flat"


def test_export_before_optimizing_leaves_no_directory_behind(tmp_path):
    """The guard exists so a library that isn't ready leaves no empty run directory."""
    spec = LibrarySpec(name="early", protein_sequence="MKAILVDE", substitutions=["A"])
    lib = SubstitutionScan(spec).generate()
    with pytest.raises(ValueError, match="not codon-optimized"):
        lib.export_all(tmp_path / "out")
    assert not (tmp_path / "out").exists()


def test_run_dir_is_shared_by_the_piecemeal_exporters(mbo038, tmp_path):
    out = mbo038.run_dir(tmp_path / "out")
    assert out.is_dir() and _RUN_DIR.match(out.name)
    mbo038.to_usortm(out / "variants.csv")
    assert out == mbo038.run_dir(tmp_path / "out")     # same run, same directory
    assert (out / "variants.csv").is_file()


# --- determinism --------------------------------------------------------------

def test_codon_optimization_is_deterministic_given_seed():
    # The reproducibility claim, same seed gives the same output.
    def run():
        spec = LibrarySpec(name="det", protein_sequence="MKAILVDE",
                           substitutions=["A"], seed=3)
        return SubstitutionScan(spec).generate().codon_optimize()

    a, b = run(), run()
    assert a.reference == b.reference
    assert list(a.df["variable_dna"]) == list(b.df["variable_dna"])


# --- design-specs hygiene -----------------------------------------------------

def test_drop_failed_also_clears_the_design_specs_record():
    """The design-specs JSON is the handoff record, so it must not keep listing
    failures for variants the library no longer holds."""
    spec = LibrarySpec(name="dropped", protein_sequence="MKAILVDE", substitutions=["A"])
    lib = SubstitutionScan(spec).generate().codon_optimize()
    lib.failed = {"K2A": "synthetic failure"}
    lib.design_specs["failed"] = dict(lib.failed)

    lib.drop_failed()
    assert lib.failed == {}
    assert "failed" not in lib.design_specs


# --- spec rendering -----------------------------------------------------------

def test_spec_html_shows_readable_placeholders():
    """An unset field renders as words, not as the stray punctuation left behind when
    the placeholder glyph was stripped."""
    spec = LibrarySpec(name="x", protein_sequence="MKV", avoid_enzymes=[])
    html = spec._repr_html_()
    for label in ("platform", "max_oligo_length"):
        cell = re.search(rf">{label}</th><td[^>]*>(.*?)</td>", html).group(1)
        assert cell == "not set"
    assert "no enzymes" in html


# --- reproducibility ----------------------------------------------------------

# A stochastic configuration: match_codon_usage samples, and the GC ceiling plus the
# avoid patterns force the constraint solver to search, so an unpinned RNG would show.
_DET_PROTEIN = "MKAILVGADEQTRWYFNSHCPMKAILVGADEQ"


def _det_spec(**kwargs):
    from library_designer import CodonOptimizationParams

    return LibrarySpec(
        name="det", protein_sequence=_DET_PROTEIN, substitutions=["A"],
        optimization=CodonOptimizationParams(method="match_codon_usage", gc_max=0.62),
        **kwargs,
    )


def test_codon_optimize_ignores_the_ambient_rng():
    """The spec's seed is applied on every call, so where the caller's RNG happens to be
    cannot change the result. Calling the wrapper directly used to skip seeding."""
    import numpy as np

    from library_designer.optimize.codon import codon_optimize

    spec = _det_spec()
    np.random.seed(1)
    first = codon_optimize(_DET_PROTEIN, spec)
    np.random.seed(999)
    second = codon_optimize(_DET_PROTEIN, spec)
    assert first == second


def test_optimization_leaves_the_callers_rng_untouched():
    """Seeding numpy globally is how DNA Chisel is pinned, but it must not leak: the
    caller's own draws used to silently inherit our seed."""
    import numpy as np

    np.random.seed(4321)
    expected = np.random.rand(3).tolist()

    np.random.seed(4321)
    SubstitutionScan(_det_spec()).generate().codon_optimize()
    assert np.random.rand(3).tolist() == expected


def test_seed_none_opts_out_and_follows_the_ambient_rng():
    import numpy as np

    spec = _det_spec(seed=None)
    np.random.seed(7)
    before = np.random.get_state()[1].copy()
    a = SubstitutionScan(spec).generate().codon_optimize().reference
    assert not np.array_equal(before, np.random.get_state()[1])   # stream consumed, not restored

    np.random.seed(7)
    b = SubstitutionScan(spec).generate().codon_optimize().reference
    assert a == b                                                # ambient seed governs it


def test_seed_is_recorded_in_the_design_specs():
    lib = SubstitutionScan(_det_spec(seed=11)).generate().codon_optimize()
    assert lib.design_specs["seed"] == 11
    assert lib.design_specs["reference_cds"] == lib.reference


def test_same_spec_reproduces_in_a_fresh_process():
    """Reproducibility has to survive a new interpreter, including hash randomization,
    since that is what re-running a design next month actually looks like."""
    import os
    import subprocess
    import sys

    script = (
        "from library_designer import LibrarySpec, SubstitutionScan, CodonOptimizationParams;"
        f"spec = LibrarySpec(name='det', protein_sequence={_DET_PROTEIN!r}, substitutions=['A'],"
        "optimization=CodonOptimizationParams(method='match_codon_usage', gc_max=0.62));"
        "print(SubstitutionScan(spec).generate().codon_optimize().reference)"
    )
    runs = []
    for hashseed in ("0", "12345"):
        proc = subprocess.run([sys.executable, "-c", script], capture_output=True, text=True,
                              env={**os.environ, "PYTHONHASHSEED": hashseed})
        assert proc.returncode == 0, proc.stderr
        runs.append(proc.stdout.strip())
    assert runs[0] == runs[1]
    assert runs[0] == SubstitutionScan(_det_spec()).generate().codon_optimize().reference


# --- finding a library's own runs ---------------------------------------------

def test_a_library_finds_its_own_runs(mbo038, tmp_path):
    """The library knows its name, so nothing about a run has to be typed except which one."""
    from datetime import timedelta

    out = tmp_path / "out"
    mbo038.export_all(out, plots=False)
    first = mbo038.output_dir
    mbo038._created_at += timedelta(days=1)
    mbo038.export_all(out, plots=False)
    second = mbo038.output_dir
    (out / "someone_elses_library_20260101_000000").mkdir()
    (out / "notes").mkdir()

    assert mbo038.runs(out) == [first, second]          # only this library's, oldest first
    assert mbo038.latest_run(out) == second
    assert mbo038.run(first.name.split("_")[-2], out) == first     # by date alone
    mbo038._created_at -= timedelta(days=1)


def test_naming_a_run_that_is_not_there_lists_the_ones_that_are(mbo038, tmp_path):
    out = tmp_path / "out"
    mbo038.export_all(out, plots=False)
    stamp = mbo038.output_dir.name.split("hAcyP1_")[-1]

    with pytest.raises(FileNotFoundError, match=f"available: \\['{stamp}'\\]"):
        mbo038.run("19990101_000000", out)
    assert mbo038.latest_run(tmp_path / "nowhere") is None
    assert mbo038.runs(tmp_path / "nowhere") == []


def test_a_run_record_is_read_back_without_rebuilding_its_filename(mbo038, tmp_path):
    mbo038.export_all(tmp_path / "out", plots=False)
    run = mbo038.latest_run(tmp_path / "out")

    record = mbo038.run_record(run)
    assert record["reference_cds"] == mbo038.reference
    assert mbo038.matches_run(run)

    with pytest.raises(FileNotFoundError, match="No design-specs record"):
        mbo038.run_record(tmp_path)


def test_matches_run_says_no_when_the_design_has_moved_on(mbo038, tmp_path):
    """The guard before adding files to a directory somebody already ordered from."""
    mbo038.export_all(tmp_path / "out", plots=False)
    run = mbo038.latest_run(tmp_path / "out")
    assert mbo038.matches_run(run)

    other = SubstitutionScan(
        LibrarySpec(name="hAcyP1", protein_sequence="MKAILVDEQTRW", substitutions=["A"])
    ).generate().codon_optimize()
    assert not other.matches_run(run)                   # same name, different design


def test_exporters_default_to_the_libraries_own_run_directory(tmp_path, monkeypatch):
    import sys

    sys.path.insert(0, "tests")
    from test_vectors import CLEAN_CDS, _plasmid, _standard_lib, A5_BACKBONE, A3_BACKBONE

    monkeypatch.chdir(tmp_path)
    lib = _standard_lib(_plasmid(tmp_path, CLEAN_CDS), A5_BACKBONE, A3_BACKBONE)
    lib.to_assembled_vectors()                          # no path given
    lib.to_vector_maps()

    run = lib.latest_run()
    assert run is not None
    assert (run / "assembled_vectors" / "parent_WT.gb").is_file()
    assert (run / "vector" / "destination.gb").is_file()


def test_the_handoff_carries_the_run_identity_for_downstream_tools(mbo038, tmp_path):
    """uSort-M reads variants.csv plus the record beside it. The CSV stays exactly
    name,sequence because uSort-M parses it strictly, so the run identity travels in the
    record, along with a digest tying the two together."""
    import hashlib
    import json

    mbo038.export_all(tmp_path / "out", plots=False)
    run = mbo038.output_dir
    record = json.loads((run / "hAcyP1_design_specs.json").read_text())

    assert mbo038.run_id == run.name                    # one token, and it names the folder
    assert record["run_id"] == mbo038.run_id
    assert record["created"] == mbo038.created

    handoff = record["handoff"]
    assert handoff["run_id"] == mbo038.run_id
    assert handoff["variants_csv"] == "variants.csv"
    assert handoff["n_variants"] == len(mbo038)
    csv = (run / "variants.csv").read_bytes()
    assert handoff["sha256"] == hashlib.sha256(csv).hexdigest()
    # The contract itself is untouched: two columns, nothing prepended.
    assert csv.split(b"\n")[0] == b"name,sequence"


def test_the_run_identity_changes_with_the_run_and_not_with_the_export(mbo038, tmp_path):
    from datetime import timedelta

    first = mbo038.run_id
    mbo038.export_all(tmp_path / "a", plots=False)
    mbo038.export_all(tmp_path / "b", plots=False)
    assert mbo038.run_id == first                       # re-exporting is the same run

    mbo038._created_at += timedelta(hours=1)
    assert mbo038.run_id != first                       # a later build is a different one
    mbo038._created_at -= timedelta(hours=1)


def test_a_stale_handoff_record_does_not_survive_re_optimizing():
    spec = LibrarySpec(name="h", protein_sequence="MKAILVDEQTRW", substitutions=["A"])
    lib = SubstitutionScan(spec).generate().codon_optimize()
    lib.design_specs["handoff"] = {"run_id": "from the last run"}
    lib.codon_optimize()
    assert "handoff" not in lib.design_specs


# --- the synonymous fallback at a mutated position ----------------------------
#
# A hand-built reference where one substitution's preferred codon spells a BsaI site
# and a rarer synonymous one does not. Reference codon 3 is GTC and codon 4 starts TC,
# so stamping the E. coli-preferred Ala codon GCG at codon 2 completes G|GTC|TC =
# GGTCTC, while GCC (the next Ala codon down) leaves it clean. The CDS is given
# verbatim, so no optimizer runs and the case is exact.

FALLBACK_CDS = "ATGAAAGTCTCAAAA"        # M K V S K
FALLBACK_PROTEIN = "MKVSK"


def _fallback_spec(fallback: bool = True):
    from library_designer.spec import CodonOptimizationParams

    return LibrarySpec(
        name="fallback",
        protein_sequence=FALLBACK_PROTEIN,
        cds=FALLBACK_CDS,
        substitutions=["A"],
        avoid_enzymes=["BsaI"],
        optimization=CodonOptimizationParams(synonymous_fallback=fallback),
        seed=0,
    )


def _stamped_codon(lib, name):
    row = lib.df[lib.df["name"] == name].iloc[0]
    i = int(row["mut_index"])
    return row["variable_dna"][i * 3:i * 3 + 3]


def test_the_preferred_codon_really_would_introduce_a_site():
    # Guard the premise: if E. coli's preferred Ala codon ever stops being GCG, or the
    # site stops being spelled, the two tests below would pass for the wrong reason.
    from library_designer.optimize.backbone import ranked_codons

    assert ranked_codons("e_coli")["A"][0] == "GCG"
    assert count_enzyme_sites(FALLBACK_CDS, "BsaI") == 0
    assert count_enzyme_sites("ATG" + "GCG" + FALLBACK_CDS[6:], "BsaI") == 1


def test_synonymous_fallback_steps_down_to_keep_the_variant():
    lib = SubstitutionScan(_fallback_spec(True)).generate().codon_optimize()
    assert lib.failed == {}                        # K2A is still makeable
    assert _stamped_codon(lib, "K2A") == "GCC"     # at a rarer codon than GCG
    assert count_enzyme_sites(lib.df[lib.df["name"] == "K2A"]["variable_dna"].iloc[0],
                              "BsaI") == 0
    # Every other position took the preferred codon, so the step-down is local.
    others = [n for n in lib.df["name"] if n not in ("K2A", "WT")]
    assert all(_stamped_codon(lib, n) == "GCG" for n in others)


def test_synonymous_fallback_off_flags_the_position_instead():
    on = SubstitutionScan(_fallback_spec(True)).generate().codon_optimize()
    off = SubstitutionScan(_fallback_spec(False)).generate().codon_optimize()

    assert set(off.failed) == {"K2A"}
    msg = off.failed["K2A"]
    assert "GCG" in msg and "restricted motif" in msg and "synonymous_fallback is off" in msg
    assert pd.isna(off.df[off.df["name"] == "K2A"]["variable_dna"].iloc[0])
    assert "K2A" in off.check().optimization_failed

    # The flag changes only the position it had to compromise on: same reference, and
    # every other member is identical to the fallback-on run.
    assert off.reference == on.reference
    kept = {n: s for n, s in zip(off.df["name"], off.df["variable_dna"]) if n != "K2A"}
    was = {n: s for n, s in zip(on.df["name"], on.df["variable_dna"]) if n != "K2A"}
    assert kept == was


def test_turning_the_fallback_off_is_recorded_in_the_design_specs():
    from library_designer.spec import optimization_line

    off = SubstitutionScan(_fallback_spec(False)).generate().codon_optimize()
    assert off.design_specs["spec"]["optimization"]["synonymous_fallback"] is False
    assert "no synonymous fallback" in optimization_line(
        off.design_specs["spec"]["optimization"])
    # On by default, so the line stays uncluttered for the usual case.
    on = SubstitutionScan(_fallback_spec(True)).generate().codon_optimize()
    assert "synonymous fallback" not in optimization_line(
        on.design_specs["spec"]["optimization"])
    # An older record predates the field and described the unconditional stepping.
    assert "synonymous fallback" not in optimization_line({"species": "e_coli",
                                                           "method": "use_best_codon"})


# --- naming the designed protein ----------------------------------------------

def test_designed_sequence_is_the_whole_protein_when_nothing_is_truncated():
    spec = LibrarySpec(name="n", protein_sequence="MKAILV", substitutions=["A"])
    assert spec.truncation == 0
    assert spec.designed_sequence == "MKAILV"
    assert spec.truncated_sequence == spec.designed_sequence   # alias still works
    assert spec.protein_description() == "protein_sequence"


def test_designed_sequence_drops_the_truncation_and_says_so():
    spec = LibrarySpec(name="n", protein_sequence="MKAILV", substitutions=["A"],
                       truncation=2)
    assert spec.designed_sequence == "AILV"
    assert spec.truncated_sequence == "AILV"
    assert spec.protein_description() == "the truncated protein_sequence (truncation=2)"


def test_a_cds_mismatch_only_blames_truncation_when_the_spec_truncates():
    from library_designer.optimize.backbone import build_reference

    spec = LibrarySpec(name="n", protein_sequence="MKAILV", substitutions=["A"],
                       cds="ATGAAAGCTATTTTGGTC")          # encodes MKAILV
    assert build_reference(spec) == "ATGAAAGCTATTTTGGTC"

    wrong = LibrarySpec(name="n", protein_sequence="MKAILW", substitutions=["A"],
                        cds="ATGAAAGCTATTTTGGTC")
    with pytest.raises(ValueError) as exc:
        build_reference(wrong)
    assert "does not translate to protein_sequence" in str(exc.value)
    assert "truncated" not in str(exc.value)

    trunc = LibrarySpec(name="n", protein_sequence="QQMKAILW", substitutions=["A"],
                        truncation=2, cds="ATGAAAGCTATTTTGGTC")
    with pytest.raises(ValueError) as exc:
        build_reference(trunc)
    assert "the truncated protein_sequence (truncation=2)" in str(exc.value)


# --- one FASTA per oligo ------------------------------------------------------

def test_export_all_writes_one_fasta_per_member(mbo038, tmp_path):
    from library_designer.io import _file_stem

    mbo038.export_all(tmp_path / "out", plots=False, oligo_fmt="fasta")
    d = mbo038.output_dir / "oligos"
    files = {p.name for p in d.iterdir()}
    assert len(files) == len(mbo038.df)                 # every member, WT control included
    assert files == {f"{_file_stem(str(n))}.fasta" for n in mbo038.df["name"]}

    # Each file is one record holding the assembled construct, adaptors and all.
    for name, dna, a5, a3 in zip(mbo038.df["name"], mbo038.df["variable_dna"],
                                 mbo038.df["adaptor_5"], mbo038.df["adaptor_3"]):
        head, seq, tail = (d / f"{_file_stem(str(name))}.fasta").read_text().split("\n")
        assert head == f">{name}"                       # the real name, '*' and all
        assert seq == assemble(a5, dna, a3).upper()
        assert tail == ""


def test_an_amber_stop_is_spelled_out_in_the_filename(mbo038, tmp_path):
    mbo038.export_all(tmp_path / "out", plots=False, oligo_fmt="fasta")
    d = mbo038.output_dir / "oligos"
    ambers = [str(n) for n in mbo038.df["name"] if str(n).endswith("*")]
    assert ambers                                        # the fixture scans TAG
    for name in ambers:
        # '*' globs in a shell and is illegal on Windows, so the file says "stop" while
        # the record inside keeps the variant's real name.
        assert not (d / f"{name}.fasta").exists()
        f = d / f"{name[:-1]}stop.fasta"
        assert f.read_text().startswith(f">{name}\n")


def test_skipping_the_fastas_leaves_no_directory(mbo038, tmp_path):
    mbo038.export_all(tmp_path / "out", plots=False, oligos=False)
    assert not (mbo038.output_dir / "oligos").exists()
    assert (mbo038.output_dir / "variants.csv").exists()   # the rest still went out


def test_two_names_that_share_a_filename_are_refused(tmp_path):
    """``K2*`` and a literal ``K2stop`` would both want the same K2stop file. Overwriting
    would drop a variant from an order, so the export refuses instead."""
    from library_designer.io import to_oligo_files

    spec = LibrarySpec(name="c", protein_sequence="MKAILV", substitutions=["A"])
    lib = SubstitutionScan(spec).generate().codon_optimize()
    lib.df.loc[0, "name"] = "K2*"
    lib.df.loc[1, "name"] = "K2stop"

    with pytest.raises(ValueError, match=r"both map to the same filename \(K2stop\)"):
        to_oligo_files(lib, tmp_path / "oligos")
    # Checked before anything is written, so a refused export leaves no directory behind.
    assert not (tmp_path / "oligos").exists()


def test_genbank_oligos_are_annotated(mbo038, tmp_path):
    """The default per-oligo format. Each file opens in a plasmid editor already labelled
    with the coding stretch, the mutated codon, and the Type IIS sites."""
    from Bio import SeqIO

    from library_designer.io import _file_stem

    mbo038.export_all(tmp_path / "out", plots=False)         # genbank by default
    d = mbo038.output_dir / "oligos"
    assert {p.suffix for p in d.iterdir()} == {".gb"}

    row = mbo038.df[mbo038.df["name"] == "I7F"].iloc[0]
    rec = SeqIO.read(d / f"{_file_stem('I7F')}.gb", "genbank")
    assert str(rec.seq) == assemble(row["adaptor_5"], row["variable_dna"], row["adaptor_3"]).upper()

    by_type: dict[str, list] = {}
    for f in rec.features:
        by_type.setdefault(f.type, []).append(f)

    # The coding stretch, placed where it actually sits in the oligo.
    (cds,) = by_type["CDS"]
    assert str(cds.extract(rec.seq)) == row["variable_dna"]
    assert cds.qualifiers["label"] == [mbo038.spec.name]

    # The mutated codon, inside the CDS and in frame with it.
    (var,) = by_type["variation"]
    assert var.qualifiers["label"][0].startswith("I7F")
    assert (int(var.location.start) - int(cds.location.start)) % 3 == 0
    assert int(var.location.start) == int(cds.location.start) + int(row["mut_index"]) * 3

    # Both BsaI sites, one per adaptor, on opposite strands.
    sites = by_type["protein_bind"]
    assert {f.qualifiers["label"][0] for f in sites} == {"BsaI"}
    assert sorted(f.location.strand for f in sites) == [-1, 1]


def test_a_wild_type_control_gets_no_mutated_codon_feature(mbo038, tmp_path):
    from Bio import SeqIO

    mbo038.to_oligo_files(tmp_path / "oligos", fmt="genbank")
    rec = SeqIO.read(tmp_path / "oligos" / "WT.gb", "genbank")
    assert [f.type for f in rec.features].count("variation") == 0
    assert [f.type for f in rec.features].count("CDS") == 1


def test_both_formats_writes_two_files_per_oligo(mbo038, tmp_path):
    n = len(mbo038.df)
    mbo038.to_oligo_files(tmp_path / "oligos", fmt="both")
    d = tmp_path / "oligos"
    assert len(list(d.glob("*.gb"))) == n and len(list(d.glob("*.fasta"))) == n

    with pytest.raises(ValueError, match="Unknown fmt"):
        mbo038.to_oligo_files(tmp_path / "other", fmt="snapgene")


# --- the single-WT-reference invariant, as a runtime check --------------------

def test_off_target_edits_is_clean_on_a_well_formed_library(mbo038):
    from library_designer.checks.report import off_target_edits

    assert off_target_edits(mbo038) == []
    assert mbo038.check().off_target_edits == []


def test_a_variant_that_strays_from_the_reference_is_caught(mbo038):
    """The check that makes a single-mutant library interpretable. A base changed anywhere
    but the member's own codon has to be reported, since a phenotype could no longer be
    pinned on the substitution."""
    from library_designer.checks.report import off_target_edits

    lib = SubstitutionScan(mbo038.spec).generate().codon_optimize()
    row = lib.df.index[lib.df["name"] == "I7F"][0]
    dna = lib.df.at[row, "variable_dna"]
    far = (int(lib.df.at[row, "mut_index"]) + 20) * 3          # a codon it must not touch
    swapped = "A" if dna[far] != "A" else "C"
    lib.df.at[row, "variable_dna"] = dna[:far] + swapped + dna[far + 1:]

    assert off_target_edits(lib) == ["I7F"]
    rep = lib.check()
    assert rep.off_target_edits == ["I7F"] and not rep.passed
    assert rep.issues["off_target_edits"] == ["I7F"]
    assert "unintended edits:" in rep.text()
    assert "1 outside their own codon (I7F)" in rep.text()


def test_a_silent_edit_inside_the_intended_codon_is_not_a_finding(mbo038):
    """A different synonymous codon at the member's own position is the stamp doing its job
    (stepping down to avoid a motif), not an unintended edit."""
    from library_designer.checks.report import off_target_edits

    lib = SubstitutionScan(mbo038.spec).generate().codon_optimize()
    row = lib.df.index[lib.df["name"] == "I7F"][0]
    dna, i = lib.df.at[row, "variable_dna"], int(lib.df.at[row, "mut_index"])
    lib.df.at[row, "variable_dna"] = dna[:i * 3] + "TTC" + dna[i * 3 + 3:]   # Phe, other codon
    assert off_target_edits(lib) == []


def test_a_wild_type_control_must_match_the_reference_exactly(mbo038):
    from library_designer.checks.report import off_target_edits

    lib = SubstitutionScan(mbo038.spec).generate().codon_optimize()
    row = lib.df.index[lib.df["name"] == "WT"][0]
    dna = lib.df.at[row, "variable_dna"]
    lib.df.at[row, "variable_dna"] = ("A" if dna[0] != "A" else "C") + dna[1:]
    assert off_target_edits(lib) == ["WT"]


def test_a_sequence_set_has_no_shared_reference_to_check_against():
    from library_designer import SequenceSet
    from library_designer.checks.report import off_target_edits

    spec = LibrarySpec(name="s", substitutions=[])
    lib = SequenceSet(spec, proteins={"a": "MKAILV", "b": "MKGWQP"}).generate().codon_optimize()
    assert lib.reference is None
    assert off_target_edits(lib) == []              # independent genes, nothing to compare
    assert lib.check().off_target_edits == []


# --- the codon map ------------------------------------------------------------

def test_codon_matrix_groups_rows_by_amino_acid_and_counts_members(mbo038):
    """Rows are every codon of every residue in play, grouped by amino acid and ordered
    most- to least-used, so a codon low in its band is a compromise. Cells count members."""
    import matplotlib
    matplotlib.use("Agg")

    from library_designer.optimize.backbone import ranked_codons
    from library_designer.viz import _codon_rows

    ranked = ranked_codons(mbo038.spec.optimization.species)
    rows = _codon_rows("e_coli", {"A", "F"})
    assert rows == [("A", c) for c in ranked["A"]] + [("F", c) for c in ranked["F"]]
    assert rows[0][1] == ranked["A"][0]          # the preferred codon leads its group

    fig = mbo038.plot_codon_matrix()
    (ax,) = [a for a in fig.axes if a.get_xlabel()]      # the map, not the colorbar
    labels = [t.get_text() for t in ax.get_yticklabels()]
    assert labels[0].split()[0] <= labels[-1].split()[0]  # amino acids in order down the axis
    assert "F TTT" in labels and "*" in " ".join(labels)   # the amber scan gets its own band

    # One cell per position carries the reference codon, so the busiest cell is a whole
    # sublibrary's worth of members and the quietest is a single substitution.
    (im,) = ax.get_images()
    import numpy as np
    data = np.array(im.get_array(), dtype=float)
    assert data.shape == (len(labels), len(mbo038.reference) // 3)
    assert np.nanmin(data) == 1                            # a substitution, carried by one
    # The reference codon at a position is carried by everyone who does not mutate it, so the
    # busiest cell is most of the library and every column sums to the whole library.
    assert np.nanmax(data) > len(mbo038.df) / 2
    assert set(np.nansum(data, axis=0).tolist()) == {float(len(mbo038.df))}


def test_the_codon_map_of_the_reference_alone_has_one_cell_per_position(mbo038):
    import matplotlib
    matplotlib.use("Agg")
    import numpy as np

    fig = mbo038.plot_codon_matrix(reference_only=True)
    (ax,) = [a for a in fig.axes if a.get_xlabel()]
    (im,) = ax.get_images()
    data = np.array(im.get_array(), dtype=float)
    filled = np.isfinite(data).sum(axis=0)
    assert set(filled.tolist()) == {1}                    # exactly one codon per position
    assert np.nanmax(data) == 1


def test_export_all_writes_both_qc_plots(mbo038, tmp_path):
    mbo038.export_all(tmp_path / "out", plots=True)
    pngs = sorted(p.name for p in mbo038.output_dir.glob("*.png"))
    assert pngs == ["hAcyP1_codon_matrix.png", "hAcyP1_codon_usage.png"]


def test_dpi_is_settable_on_every_plot_and_reaches_the_file(mbo038, tmp_path):
    """One number governs a figure inline and on disk. The default is raised above
    matplotlib's 100, since the codon map draws cells a couple of pixels wide."""
    import matplotlib
    matplotlib.use("Agg")

    from library_designer.viz import DEFAULT_DPI

    assert DEFAULT_DPI > 100
    for make in (mbo038.plot_codon_usage, mbo038.plot_codon_matrix):
        assert make().dpi == DEFAULT_DPI
        assert make(dpi=400).dpi == 400

    # A saved file inherits the figure's dpi, so the same call sharpens both.
    from PIL import Image

    mbo038.to_qc_plots(tmp_path / "low.png", dpi=100)
    mbo038.to_qc_plots(tmp_path / "high.png", dpi=300)
    low, high = Image.open(tmp_path / "low.png"), Image.open(tmp_path / "high.png")
    assert high.width > low.width * 2.5

    mbo038.export_all(tmp_path / "out", dpi=300)
    matrix = Image.open(mbo038.output_dir / "hAcyP1_codon_matrix.png")
    mbo038.export_all(tmp_path / "out2", dpi=100)
    smaller = Image.open(mbo038.output_dir / "hAcyP1_codon_matrix.png")
    assert matrix.width > smaller.width * 2.5
