"""End-to-end golden path for the flagship non-tiled substitution scan.

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
    # The flagship path codon-optimizes a protein (no native CDS) and must not fail.
    assert mbo038.failed == {}
    for protein, dna in zip(mbo038.df["protein"], mbo038.df["variable_dna"]):
        assert isinstance(dna, str)
        assert translates_to(dna, protein)


def test_single_wt_reference_invariant(mbo038):
    # On the *codon-optimized* reference (not a native-CDS fixture): every member is
    # byte-identical to the frozen reference except within its own stamped codon.
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
        "hAcyP1_design_specs.json",
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
    """Two runs must not share a directory, which is the whole point of the stamp."""
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
    """The guard exists so a library that isn't ready cannot litter an empty run dir."""
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
    # The design-specs/reproducibility claim: same seed -> byte-identical output.
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
    cannot change the design. Calling the wrapper directly used to skip seeding."""
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
