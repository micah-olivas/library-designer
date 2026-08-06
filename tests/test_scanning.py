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
    # On the *codon-optimized* reference (not a native-CDS fixture): every member matches
    # the frozen reference except within its own stamped codon.
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
    assert spec.protein_description() == (
        "the truncated protein_sequence (truncation=2 from the N terminus)"
    )


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
    assert "the truncated protein_sequence (truncation=2 from the N terminus)" in str(exc.value)


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
    (ax,) = [a for a in fig.axes if a.get_images()]      # the map, not the colourbar
    # Ticks name the codon. The amino acid is one larger letter per group out in the margin,
    # so a residue is read once rather than repeated down every row of its block.
    labels = [t.get_text() for t in ax.get_yticklabels()]
    assert "TTT" in labels and "GCG" in labels
    assert all(len(t) == 3 for t in labels)              # codon only, no residue prefix

    letters = [a.get_text() for a in ax.texts]
    assert letters == sorted({aa for aa, _ in _codon_rows("e_coli", set(
        mbo038.df["mut_residue"].dropna().astype(str)) | set(mbo038.spec.designed_sequence))})
    assert "*" in letters and letters == sorted(letters)  # amber included, residues in order
    assert len(letters) < len(labels)                     # one letter per group, not per row

    # Grey banding is gone; the grouping is a black rule at each boundary instead.
    rules = [ln for ln in ax.get_lines() if ln.get_color() == "black"]
    assert len(rules) == len(letters) - 1                 # a rule between groups, not above the first
    assert not [p for p in ax.patches if p.get_facecolor()[:3] == (0.92, 0.92, 0.92)]

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
    (ax,) = [a for a in fig.axes if a.get_images()]
    (im,) = ax.get_images()
    data = np.array(im.get_array(), dtype=float)
    filled = np.isfinite(data).sum(axis=0)
    assert set(filled.tolist()) == {1}                    # exactly one codon per position
    assert np.nanmax(data) == 1


def test_export_all_writes_the_qc_plots_in_a_qc_subdir(mbo038, tmp_path):
    mbo038.export_all(tmp_path / "out", plots=True)
    qc = mbo038.output_dir / "qc"
    assert sorted(p.name for p in qc.iterdir()) == [
        "codon_matrix.pdf", "codon_usage.pdf", "gc_distribution.pdf",
    ]
    assert not list(mbo038.output_dir.glob("*.p*g"))   # nothing left loose in the run dir
    for plot in qc.iterdir():
        assert plot.read_bytes()[:4] == b"%PDF"


def test_exported_plots_name_their_font_in_the_pdf(mbo038, tmp_path):
    """Text goes out as embedded TrueType, not Type 3 outlines, so a label stays selectable
    and the face is named in the file. Skipped where the first choice is not installed, since
    the stack then legitimately falls back."""
    from matplotlib import font_manager

    from library_designer.viz import FONT_STACK

    wanted = FONT_STACK[0]
    if wanted not in {f.name for f in font_manager.fontManager.ttflist}:
        pytest.skip(f"{wanted} is not installed on this machine")

    mbo038.export_all(tmp_path / "out")
    for plot in (mbo038.output_dir / "qc").iterdir():
        assert wanted.encode() in plot.read_bytes(), f"{plot.name} is not typeset in {wanted}"


def test_export_all_skips_the_qc_subdir_when_plots_are_off(mbo038, tmp_path):
    mbo038.export_all(tmp_path / "out", plots=False)
    assert not (mbo038.output_dir / "qc").exists()


def test_dpi_is_settable_on_every_plot_and_reaches_the_file(mbo038, tmp_path):
    """One number governs a figure inline and in a raster file. The default is raised above
    matplotlib's 100, since the codon map draws cells a couple of pixels wide."""
    import matplotlib
    matplotlib.use("Agg")

    from library_designer.viz import DEFAULT_DPI

    assert DEFAULT_DPI > 100
    for make in (mbo038.plot_codon_usage, mbo038.plot_codon_matrix):
        assert make().dpi == DEFAULT_DPI
        assert make(dpi=400).dpi == 400

    # A saved image inherits the figure's dpi, so the same call sharpens both.
    from PIL import Image

    mbo038.to_qc_plots(tmp_path / "low.png", dpi=100)
    mbo038.to_qc_plots(tmp_path / "high.png", dpi=300)
    low, high = Image.open(tmp_path / "low.png"), Image.open(tmp_path / "high.png")
    assert high.width > low.width * 2.5

    # The exported QC plots are vector, so dpi cannot blur them and the file stands either way.
    for value in (100, 300):
        mbo038.export_all(tmp_path / f"out{value}", dpi=value)
        assert (mbo038.output_dir / "qc" / "codon_matrix.pdf").read_bytes()[:4] == b"%PDF"


# --- which terminus the truncation comes off ----------------------------------

_TRUNC_PROTEIN = "MKAILVDEQ"                        # 9 aa
_TRUNC_CDS = "ATGAAAGCTATTTTGGTTGATGAACAA"         # 27 bp, encodes it


def _trunc_spec(**kw):
    kw.setdefault("cds", _TRUNC_CDS)
    return LibrarySpec(name="t", protein_sequence=_TRUNC_PROTEIN, substitutions=["A"], **kw)


def test_truncation_defaults_to_the_n_terminus():
    spec = _trunc_spec(truncation=2)
    assert spec.terminus == "N"
    assert spec.designed_sequence == "AILVDEQ"      # MK removed
    assert spec.numbering_offset == 2


def test_a_c_terminal_truncation_takes_residues_off_the_other_end():
    spec = _trunc_spec(truncation=2, truncation_terminus="C")
    assert spec.terminus == "C"
    assert spec.designed_sequence == "MKAILVD"      # EQ removed
    # Nothing was removed before residue 1, so full-protein numbering is unchanged.
    assert spec.numbering_offset == 0


def test_the_terminus_is_named_in_messages_and_the_repr():
    n = _trunc_spec(truncation=2)
    c = _trunc_spec(truncation=2, truncation_terminus="C")
    assert "from the N terminus" in n.protein_description()
    assert "from the C terminus" in c.protein_description()
    assert "at the C terminus" in repr(c)
    assert _trunc_spec().protein_description() == "protein_sequence"   # silent when 0


def test_an_unknown_terminus_is_refused():
    with pytest.raises(ValueError, match="truncation_terminus must be 'N' or 'C'"):
        _trunc_spec(truncation=1, truncation_terminus="middle").designed_sequence


def test_a_full_length_cds_is_trimmed_to_the_designed_region():
    """The truncation is stated once and applies to the protein and to the DNA encoding it,
    so a full-length CDS is trimmed rather than refused."""
    from library_designer.optimize.backbone import build_reference

    n = build_reference(_trunc_spec(truncation=2))
    assert n == _TRUNC_CDS[6:]                       # first two codons dropped
    assert len(n) == 3 * len("AILVDEQ")

    c = build_reference(_trunc_spec(truncation=2, truncation_terminus="C"))
    assert c == _TRUNC_CDS[:-6]                      # last two codons dropped


def test_an_already_trimmed_cds_is_still_accepted():
    """Both forms work, so a spec that used to pass a pre-trimmed CDS keeps working."""
    from library_designer.optimize.backbone import build_reference

    spec = _trunc_spec(truncation=2, cds=_TRUNC_CDS[6:])
    assert build_reference(spec) == _TRUNC_CDS[6:]


def test_a_cds_matching_neither_form_is_still_refused():
    from library_designer.optimize.backbone import build_reference

    wrong = "ATG" + _TRUNC_CDS[3:]                   # right length, wrong residues
    spec = LibrarySpec(name="t", protein_sequence="MWAILVDEQ", substitutions=["A"],
                       truncation=2, cds=wrong)
    with pytest.raises(ValueError, match="does not translate to"):
        build_reference(spec)


def test_variant_names_follow_the_truncated_end():
    """Names stay on full-protein numbering. An N-terminal truncation shifts them; a
    C-terminal one cannot, since it removes residues past the end."""
    n = SubstitutionScan(_trunc_spec(truncation=2)).generate()
    c = SubstitutionScan(_trunc_spec(truncation=2, truncation_terminus="C")).generate()

    # An alanine scan skips the position that is already Ala, which is residue 3 in both, so
    # compare against the designed span minus that one.
    n_pos = sorted(int(p) for p in n.df["position"].dropna())
    c_pos = sorted(int(p) for p in c.df["position"].dropna())
    assert n_pos == [4, 5, 6, 7, 8, 9]               # designed span is residues 3..9
    assert c_pos == [1, 2, 4, 5, 6, 7]               # designed span is residues 1..7

    # And the residue letter in each name is the one at that full-protein position.
    for df in (n.df, c.df):
        for name, pos, wt in zip(df["name"], df["position"], df["wt_residue"]):
            if pd.isna(pos):
                continue
            assert _TRUNC_PROTEIN[int(pos) - 1] == wt == str(name)[0]


def test_compare_reference_trims_a_full_length_paste_to_the_designed_region():
    """The same rule as spec.cds: a truncation applies to the protein and to any DNA
    encoding it. Comparing untrimmed would read the two out of frame by the truncation and
    report almost every codon as different."""
    spec = LibrarySpec(name="c", protein_sequence="MKAILVDEQ", substitutions=["A"],
                       cds="ATGAAAGCTATTTTGGTTGATGAACAA", truncation=2)
    lib = SubstitutionScan(spec).generate().codon_optimize()
    assert len(lib.reference) == 21                       # MK dropped

    full = "ATGAAAGCTATTTTGGTTGATGAACAA"                  # what a tool given the whole protein returns
    rep = lib.compare_reference(full, label="ext")
    assert rep.protein_match                              # trimmed, not a mismatch
    assert "trimmed by 2 codon(s) at the N terminus" in rep.protein_note
    assert rep.codon_agreement == 1.0 and rep.n_codons == 7
    assert "note:" in str(rep) and "!" not in str(rep).split("\n")[1]

    # An already-trimmed paste still works and says nothing about trimming.
    same = lib.compare_reference(full[6:], label="ext")
    assert same.protein_match and same.protein_note == ""


def test_compare_reference_still_reports_a_genuinely_different_protein():
    spec = LibrarySpec(name="c", protein_sequence="MKAILVDEQ", substitutions=["A"],
                       cds="ATGAAAGCTATTTTGGTTGATGAACAA", truncation=2)
    lib = SubstitutionScan(spec).generate().codon_optimize()
    rep = lib.compare_reference("ATGTGGTGGTGGTGGTGGTGGTGG", label="ext")
    assert not rep.protein_match
    assert "the designed region" in rep.protein_note or "differing at" in rep.protein_note


def test_the_plots_number_codons_the_way_variant_names_do():
    """A truncated library's plots have to agree with its labels. The codon map spans the whole
    protein in full-protein numbering, with the truncated residues greyed, and the codon-usage
    trace starts at the first designed residue."""
    import matplotlib
    matplotlib.use("Agg")

    plain = SubstitutionScan(_trunc_spec()).generate().codon_optimize()
    n = SubstitutionScan(_trunc_spec(truncation=2)).generate().codon_optimize()
    c = SubstitutionScan(
        _trunc_spec(truncation=2, truncation_terminus="C")).generate().codon_optimize()

    def x_span(lib):
        fig = lib.plot_codon_matrix()
        (ax,) = [a for a in fig.axes if a.get_images()]
        lo, hi = ax.get_xlim()
        return round(lo + 0.5), round(hi - 0.5), ax.get_xlabel()

    # The protein is 9 aa, so every map covers 1..9 whichever end the truncation came off.
    assert x_span(plain) == (1, 9, "Codon position (CDS)")
    lo, hi, label = x_span(n)
    assert (lo, hi) == (1, 9) and "starts at 3" in label
    assert x_span(c) == (1, 9, "Codon position (CDS)")             # numbering unshifted

    # The data itself still covers the designed window only; the rest is a grey span.
    def designed_and_grey(lib):
        (ax,) = [a for a in lib.plot_codon_matrix().axes if a.get_images()]
        (im,) = ax.get_images()
        x0, x1, *_ = im.get_extent()
        grey = [(round(p.get_x() + 0.5), round(p.get_x() + p.get_width() - 0.5))
                for p in ax.patches]
        return (round(x0 + 0.5), round(x1 - 0.5)), grey

    assert designed_and_grey(plain) == ((1, 9), [])                # nothing truncated
    assert designed_and_grey(n) == ((3, 9), [(1, 2)])              # residues 1-2 dropped
    assert designed_and_grey(c) == ((1, 7), [(8, 9)])              # residues 8-9 dropped

    # And the first designed position is the first variant's own position.
    assert 3 == min(int(p) for p in n.df["position"].dropna()) - 1   # residue 3 is already Ala


def test_the_usage_overlay_is_trimmed_to_the_designed_region():
    """A full-length paste overlays in frame on a truncated reference rather than shifted."""
    import matplotlib
    matplotlib.use("Agg")

    from library_designer.compare import trim_to_designed

    spec = _trunc_spec(truncation=2)
    lib = SubstitutionScan(spec).generate().codon_optimize()
    trimmed, did = trim_to_designed(spec, _TRUNC_CDS)
    assert did and trimmed == _TRUNC_CDS[6:]

    fig = lib.plot_codon_usage(compare=_TRUNC_CDS, compare_label="ext")
    (line,) = [ln for ln in fig.axes[0].get_lines() if ln.get_label() == "ext"]
    xs = list(line.get_xdata())
    assert xs[0] == 3 and len(xs) == len(lib.reference) // 3       # in frame, full width

    # An already-designed-region sequence is left alone.
    assert trim_to_designed(spec, _TRUNC_CDS[6:]) == (_TRUNC_CDS[6:], False)


def test_the_codon_map_writes_the_wt_residue_over_each_column(mbo038):
    """The WT amino acid above each position, so a codon row can be read against what the
    wild type has there. The truncated residues are lettered too, in grey, since their columns
    are drawn."""
    import matplotlib
    matplotlib.use("Agg")
    from dnachisel import translate

    fig = mbo038.plot_codon_matrix()
    (ax,) = [a for a in fig.axes if a.get_images()]
    (track,) = ax.child_axes          # the secondary x-axis carrying the residues
    labels = track.get_xticklabels()
    assert "".join(t.get_text() for t in labels) == mbo038.spec.protein_sequence

    # The designed stretch comes from the reference the map counts, and is the part left in ink.
    inked = "".join(t.get_text() for t in labels if t.get_color() != "0.45")
    assert inked == translate(mbo038.reference)
    greyed = "".join(t.get_text() for t in labels if t.get_color() == "0.45")
    assert greyed == mbo038.spec.protein_sequence[:mbo038.spec.truncation]

    # Positioned on the same full-protein numbering as the map itself, from residue 1.
    assert [round(t) for t in track.get_xticks()][:3] == [1, 2, 3]
    # And the title still clears it.
    assert ax.get_title(loc="left").startswith("Codon map,")

    off = mbo038.plot_codon_matrix(wt_track=False)
    assert not [a for a in off.axes if a.get_images()][0].child_axes


def test_the_wt_track_is_dropped_when_the_columns_are_too_narrow_to_letter():
    """A long CDS cannot carry a legible letter per column, so the track is left off rather
    than drawn as a smear."""
    import matplotlib
    matplotlib.use("Agg")

    spec = LibrarySpec.from_toml(REPO / "examples" / "gck_tiled.toml")
    gck = SubstitutionScan(spec).generate().codon_optimize()
    assert len(gck.reference) // 3 > 400
    fig = gck.plot_codon_matrix()
    (ax,) = [a for a in fig.axes if a.get_images()]
    assert not ax.child_axes                                         # no track
    assert ax.get_title(loc="left").startswith("Codon map,")


# --- masking positions out of the scan ----------------------------------------
#
# Masking and truncating are different mechanisms. A masked residue is still encoded and still
# synthesized on every oligo; it just gets no variants. A truncated one leaves the designed
# region altogether and, with a vector, is supplied by the plasmid.

def test_masked_positions_get_no_variants_but_stay_in_the_construct():
    from dnachisel import translate

    plain = SubstitutionScan(_trunc_spec()).generate().codon_optimize()
    masked = SubstitutionScan(_trunc_spec(mask_positions=[1, 2])).generate().codon_optimize()

    # Same sequence, fewer members.
    assert masked.reference == plain.reference == _TRUNC_CDS
    assert translate(masked.reference) == _TRUNC_PROTEIN
    assert len(masked.df) < len(plain.df)

    pos = {int(p) for p in masked.df["position"].dropna()}
    assert 1 not in pos and 2 not in pos
    assert pos == {int(p) for p in plain.df["position"].dropna()} - {1, 2}
    # And no variant name refers to a masked position.
    assert not [n for n in masked.df["name"] if str(n)[1:-1] in ("1", "2")]


def test_masking_and_truncating_differ_in_what_is_encoded():
    """The distinction that matters: masking keeps the residues in the reference, truncating
    takes them out of it."""
    mask = SubstitutionScan(_trunc_spec(mask_positions=[1, 2])).generate().codon_optimize()
    trunc = SubstitutionScan(_trunc_spec(truncation=2)).generate().codon_optimize()

    assert len(mask.reference) == len(_TRUNC_CDS)
    assert len(trunc.reference) == len(_TRUNC_CDS) - 6
    # Both leave the same positions unscanned, so the member counts agree.
    assert len(mask.df) == len(trunc.df)
    assert ({int(p) for p in mask.df["position"].dropna()}
            == {int(p) for p in trunc.df["position"].dropna()})


def test_mask_positions_are_normalized_and_range_checked():
    spec = _trunc_spec(mask_positions=[3, 1, 3])
    assert spec.mask_positions == [1, 3] and spec.masked == {1, 3}      # sorted, deduped

    for bad in ([0], [len(_TRUNC_PROTEIN) + 1], [-2]):
        with pytest.raises(ValueError, match="fall outside protein_sequence"):
            _trunc_spec(mask_positions=bad).masked


def test_masking_every_position_is_refused():
    """Otherwise the library is nothing but the wild-type control, which is not a library."""
    with pytest.raises(ValueError, match="leaves no position to scan"):
        SubstitutionScan(
            _trunc_spec(mask_positions=list(range(1, len(_TRUNC_PROTEIN) + 1)))
        ).generate()


def test_masking_is_recorded_and_shown():
    spec = _trunc_spec(mask_positions=[1, 2])
    assert "masked:" in repr(spec) and "1, 2" in repr(spec)
    assert "encoded, not scanned" in repr(spec)
    assert _trunc_spec().mask_positions == [] and "masked:" not in repr(_trunc_spec())

    lib = SubstitutionScan(spec).generate().codon_optimize()
    assert lib.design_specs["spec"]["mask_positions"] == [1, 2]


def test_masking_uses_full_protein_numbering_even_when_truncated():
    """Both are full-protein numbers, so they compose without the caller converting."""
    lib = SubstitutionScan(
        _trunc_spec(truncation=2, mask_positions=[3, 4])).generate().codon_optimize()
    pos = {int(p) for p in lib.df["position"].dropna()}
    assert pos == {5, 6, 7, 8, 9}          # 1-2 truncated away, 3-4 masked


def test_the_spec_table_reports_the_mask():
    """The HTML table is the notebook's review surface, so it has to say what is masked, not
    only the text repr."""
    import re

    spec = _trunc_spec(mask_positions=[1, 2])
    html = spec._repr_html_()
    (row,) = re.findall(r">masked</th><td[^>]*>(.*?)</td>", html)
    text = re.sub("<[^>]+>", "", row)
    # The residue each masked position holds, not just the number, so the row says what is
    # being left out.
    assert text.startswith(f"{_TRUNC_PROTEIN[0]}1 {_TRUNC_PROTEIN[1]}2")
    assert "encoded, not scanned" in text

    # Absent entirely when nothing is masked, so an ordinary spec's table stays short.
    assert ">masked</th>" not in _trunc_spec()._repr_html_()


def test_a_long_mask_is_elided_in_the_table():
    spec = LibrarySpec(name="m", protein_sequence="MKAILVDEQTRWYFGH",
                       substitutions=["A"], mask_positions=list(range(1, 15)))
    import re

    html = spec._repr_html_()
    (row,) = re.findall(r">masked</th><td[^>]*>(.*?)</td>", html)
    assert row.count("<code>") == 12         # capped, so a long mask cannot flood the table
    assert "and 2 more" in row


# --- the GC window on the ordered molecule ------------------------------------

def test_gc_bounds_accepts_fractions_or_percentages():
    for given in ((0.35, 0.65), (35, 65), [35, 65]):
        assert _trunc_spec(gc_bounds=given).gc_bounds == (0.35, 0.65)
    assert _trunc_spec().gc_bounds is None                    # off by default

    for bad in ((0.65, 0.35), (-0.1, 0.5), (0.5, 2.0), (0.5, 0.5)):
        with pytest.raises(ValueError, match="gc_bounds must be"):
            _trunc_spec(gc_bounds=bad)


def test_the_gc_gate_judges_the_ordered_molecule_not_the_coding_region(mbo038):
    """A vendor's GC window refers to what it receives, so the gate is on the whole molecule.
    With adaptors on, that differs from the variable region's own GC."""
    from library_designer.checks.report import gc_fraction, ordered_molecules

    ordered = ordered_molecules(mbo038)
    assert set(ordered) == set(mbo038.df["name"].astype(str))
    row = mbo038.df.iloc[0]
    assert ordered[str(row["name"])] == assemble(
        row["adaptor_5"], row["variable_dna"], row["adaptor_3"])
    assert gc_fraction(ordered[str(row["name"])]) != gc_fraction(row["variable_dna"])


def test_the_gc_gate_is_off_until_bounds_are_set():
    """No vendor registry ships with the package, so the window is the caller's to state."""
    spec = LibrarySpec.from_toml(MBO038_TOML)
    assert spec.gc_bounds is None
    lib = SubstitutionScan(spec).generate().codon_optimize()
    assert lib.check().gc_out_of_range == []

    from dataclasses import replace
    tight = SubstitutionScan(replace(spec, gc_bounds=(0.50, 0.52))).generate().codon_optimize()
    rep = tight.check()
    assert rep.gc_out_of_range and not rep.passed
    assert "GC window" in rep.text() and "outside gc_bounds" in rep.text()
    assert rep.issues["gc_out_of_range"] == rep.gc_out_of_range


def test_every_member_inside_the_window_passes():
    from dataclasses import replace

    spec = replace(LibrarySpec.from_toml(MBO038_TOML), gc_bounds=(0.0, 1.0))
    lib = SubstitutionScan(spec).generate().codon_optimize()
    rep = lib.check()
    assert rep.gc_out_of_range == [] and rep.passed


def test_ordered_gc_is_written_beside_the_variable_region_gc(mbo038, tmp_path):
    import pandas as pd

    from library_designer.checks.report import gc_fraction, ordered_molecules

    mbo038.to_full_csv(tmp_path / "full.csv")
    out = pd.read_csv(tmp_path / "full.csv")
    assert {"gc_content", "ordered_gc"} <= set(out.columns)

    ordered = ordered_molecules(mbo038)
    for name, got in zip(out["name"], out["ordered_gc"]):
        assert got == pytest.approx(round(gc_fraction(ordered[str(name)]), 3))
    # The two columns are different numbers, so neither can be mistaken for the other.
    assert not out["gc_content"].equals(out["ordered_gc"])


def test_a_tiled_library_is_gated_on_its_oligo(tmp_path):
    """A tiled pool orders the assembled oligo, primers and sites included, so that is the
    molecule the window applies to rather than the bare CDS."""
    from dataclasses import replace

    from library_designer.checks.report import ordered_molecules

    spec = replace(LibrarySpec.from_toml(REPO / "examples" / "gck_tiled.toml"),
                   gc_bounds=(0.35, 0.65))
    lib = SubstitutionScan(spec).generate().codon_optimize().drop_failed().tile()
    ordered = ordered_molecules(lib)
    placed = {str(n): o for n, o in zip(lib.df["name"], lib.df["oligo"]) if isinstance(o, str)}
    assert ordered == placed                      # the oligo, not adaptor+CDS+adaptor
    assert "WT" not in ordered                    # the global WT row rides on no oligo


# --- the GC distribution figure -----------------------------------------------

def test_gc_table_reports_both_gc_numbers_per_member(mbo038):
    from library_designer.checks.report import gc_fraction, ordered_molecules

    t = mbo038.gc_table()
    assert list(t.columns) == ["name", "sublibrary", "ordered_gc", "variable_gc", "in_bounds"]
    assert len(t) == len(ordered_molecules(mbo038))
    row = t[t["name"] == "WT"].iloc[0]
    assert row["sublibrary"] == "WT"                     # the control is its own bucket
    assert row["ordered_gc"] == pytest.approx(
        gc_fraction(ordered_molecules(mbo038)["WT"]))
    # The ordered molecule carries the adaptors, so it is not the coding region's number.
    assert (t["ordered_gc"] != t["variable_gc"]).all()
    assert t["in_bounds"].isna().all()                   # no bounds set on this spec


def test_in_bounds_follows_gc_bounds():
    from dataclasses import replace

    spec = replace(LibrarySpec.from_toml(MBO038_TOML), gc_bounds=(0.50, 0.52))
    lib = SubstitutionScan(spec).generate().codon_optimize()
    t = lib.gc_table()
    assert set(t["in_bounds"]) <= {True, False}
    outside = set(t[~t["in_bounds"].astype(bool)]["name"])
    assert outside == set(lib.check().gc_out_of_range)   # the table and the gate agree


def test_the_gc_figure_has_two_panels_and_marks_the_window():
    import matplotlib
    matplotlib.use("Agg")
    from dataclasses import replace

    spec = replace(LibrarySpec.from_toml(MBO038_TOML), gc_bounds=(0.35, 0.65))
    lib = SubstitutionScan(spec).generate().codon_optimize()
    fig = lib.plot_gc_distribution()
    detail, window = fig.axes
    # Titles are set loc="left", and get_title() reads the centre one.
    assert "own scale" in detail.get_title(loc="left")
    assert "35% to 65%" in window.get_title(loc="left")

    # The right panel shows the whole window; the left one does not have to.
    lo, hi = window.get_xlim()
    assert lo <= 0.35 and hi >= 0.65
    assert (detail.get_xlim()[1] - detail.get_xlim()[0]) < (hi - lo)
    # Both panels mark the bounds, so the margin is readable on either.
    for ax in (detail, window):
        assert [round(ln.get_xdata()[0], 2) for ln in ax.get_lines()] == [0.35, 0.65]

    # Stacked by sublibrary, plus the coding-region outline.
    labels = [t.get_text() for t in window.get_legend().get_texts()]
    assert "coding region only" in labels
    assert sum(lab.startswith("to ") for lab in labels) == lib.df["mut_residue"].nunique()
    assert "WT control" in labels


def test_the_window_panel_says_so_when_no_bounds_are_set(mbo038):
    import matplotlib
    matplotlib.use("Agg")

    fig = mbo038.plot_gc_distribution()
    _detail, window = fig.axes
    assert window.get_title(loc="left") == "gc_bounds not set"
    assert not window.get_lines()                        # nothing to mark

    off = mbo038.plot_gc_distribution(show_variable_region=False)
    labels = [t.get_text() for t in off.axes[1].get_legend().get_texts()]
    assert "coding region only" not in labels


# --- homopolymer runs ---------------------------------------------------------
#
# Both prevented and checked. The optimizer avoids long runs in the coding region, the stamp
# will not introduce one, and QC screens the finished molecule, which is where a run spelled
# across an adaptor junction shows up: the optimizer never sees the flanks.

_RUN_PROTEIN = "MKKKKKKKKKAILVDEQTRW"      # lysine runs, which AAA/AAG can spell as a long A run


def _run_lib(limit=None, **kw):
    from library_designer.checks.report import longest_run  # noqa: F401

    if limit:
        kw["max_homopolymer"] = limit
    kw.setdefault("adaptor_5", "gcgtcggtctccaagc")
    kw.setdefault("adaptor_3", "ggtgagagaccgacgc")
    spec = LibrarySpec(name="h", protein_sequence=_RUN_PROTEIN, substitutions=["A"],
                       seed=3, **kw)
    return SubstitutionScan(spec).generate().codon_optimize()


def test_the_limit_is_a_run_length_and_yields_one_pattern_per_base():
    spec = LibrarySpec(name="h", protein_sequence="MKV", max_homopolymer=7)
    assert spec.homopolymer_patterns == ["A{8,}", "C{8,}", "G{8,}", "T{8,}"]
    assert LibrarySpec(name="h", protein_sequence="MKV").homopolymer_patterns == []

    with pytest.raises(ValueError, match="max_homopolymer must be at least 1"):
        LibrarySpec(name="h", protein_sequence="MKV", max_homopolymer=0)


def test_the_optimizer_breaks_up_a_run_that_the_protein_invites():
    """A lysine stretch reverse-translates to a long A run unless something stops it."""
    from library_designer.checks.report import longest_run

    assert longest_run(_run_lib().reference) > 7          # unconstrained, the run survives
    assert longest_run(_run_lib(limit=7).reference) <= 7  # constrained, it is broken up


def test_a_stamped_codon_cannot_introduce_a_run():
    from library_designer.checks.report import longest_run, ordered_molecules

    lib = _run_lib(limit=7)
    assert all(longest_run(dna) <= 7 for dna in lib.df["variable_dna"] if isinstance(dna, str))
    assert lib.check().homopolymer_hits == []
    assert max(longest_run(s) for s in ordered_molecules(lib).values()) <= 7


def test_a_run_spelled_across_the_adaptor_junction_is_caught():
    """The screen QC adds on top of the constraint. The optimizer only sees the coding region,
    so a run finished off by the adaptor is invisible to it."""
    from library_designer.checks.report import longest_run, ordered_molecules

    lib = _run_lib(limit=7, adaptor_5="ggtctccAAAAAAA")   # 7 A's meeting the CDS
    assert longest_run(lib.reference) <= 7                # the coding region is clean
    rep = lib.check()
    assert rep.homopolymer_hits and not rep.passed
    assert "with a run over max_homopolymer" in rep.text()

    # Every member whose first codon still starts with A crosses the limit; the one mutated to
    # GCG does not, which is why this is judged per member and not per adaptor.
    ordered = ordered_molecules(lib)
    for name in lib.df["name"].astype(str):
        crossed = longest_run(ordered[name]) > 7
        assert (name in rep.homopolymer_hits) is crossed
    assert set(lib.df["name"]) - set(rep.homopolymer_hits) == {"M1A"}


def test_the_gate_is_off_until_a_limit_is_set(mbo038):
    assert mbo038.spec.max_homopolymer is None
    assert mbo038.check().homopolymer_hits == []
