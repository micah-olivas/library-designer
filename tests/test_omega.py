"""Tests for the OMEGA integration (integrations/omega.py).

OMEGA itself is a separately-installed, GPL-3.0 tool that library_designer never
imports. These tests exercise the arm's-length seam, the FASTA/primer inputs we
write, the command we build, and the outputs we parse, with ``subprocess.run``
stubbed so no real OMEGA (or its conda env) is required.
"""
from __future__ import annotations

import subprocess
from pathlib import Path

import pandas as pd
import pytest

from library_designer import LibrarySpec, OmegaParams, SubstitutionScan
from library_designer.integrations import omega

# A tiny native-CDS library: MKAILV, alanine/tyrosine scan. A native CDS skips codon
# optimization, so the fixture is instant and every member gets a full-length CDS.
PROTEIN = "MKAILV"
CDS = "ATGAAAGCGATTCTGGTT"   # M K A I L V, no BsaI site or SD-like motif


@pytest.fixture(scope="module")
def lib():
    spec = LibrarySpec(name="mini", protein_sequence=PROTEIN, cds=CDS,
                       substitutions=["A", "Y"])
    return SubstitutionScan(spec).generate().codon_optimize()


def _fake_run_factory():
    """A stub for subprocess.run that plays OMEGA: reads --output_dir out of the
    argv and drops the three expected CSVs there, then reports success."""
    def fake_run(cmd, cwd=None, capture_output=False, text=False):
        out = Path(cmd[cmd.index("--output_dir") + 1])
        out.mkdir(parents=True, exist_ok=True)
        pd.DataFrame({"name": ["mini_1"], "oligo_sequence": ["ACGTACGT"]}
                     ).to_csv(out / "oligo_order.csv", index=False)
        pd.DataFrame({"name": ["mini_1"], "oligo_sequence": ["ACGTACGT"], "fidelity": [0.95]}
                     ).to_csv(out / "optimization_results.csv", index=False)
        pd.DataFrame({"pool": [0], "fidelity": [0.95], "n_genes": [1], "random_seed": [42]}
                     ).to_csv(out / "pool_stats.csv", index=False)
        return subprocess.CompletedProcess(cmd, 0, stdout="OMEGA done", stderr="")
    return fake_run


# --- inputs we write ----------------------------------------------------------

def test_write_fasta_emits_full_cds(lib, tmp_path):
    path = tmp_path / "seqs.fasta"
    n = omega.write_fasta(lib, path)
    lines = path.read_text().splitlines()
    headers = [ln for ln in lines if ln.startswith(">")]
    seqs = [ln for ln in lines if not ln.startswith(">")]
    assert n == len(headers) == len(seqs) == len(lib)
    assert all(s.isupper() and len(s) == len(CDS) and len(s) % 3 == 0 for s in seqs)
    # coding region only, no adaptor bases leak in
    assert all(set(s) <= set("ACGT") for s in seqs)


def test_write_fasta_requires_optimized(tmp_path):
    spec = LibrarySpec(name="raw", protein_sequence=PROTEIN, cds=CDS, substitutions=["A"])
    raw = SubstitutionScan(spec).generate()          # not codon-optimized
    with pytest.raises(ValueError, match="codon_optimize"):
        omega.write_fasta(raw, tmp_path / "x.fasta")


def test_write_fasta_rejects_bad_names(lib, tmp_path):
    bad = lib.df.copy()
    bad.loc[0, "name"] = "has space"
    lib.df, saved = bad, lib.df
    try:
        with pytest.raises(ValueError, match="FASTA header"):
            omega.write_fasta(lib, tmp_path / "x.fasta")
    finally:
        lib.df = saved


def test_write_primers_omega_format(tmp_path):
    path = tmp_path / "primers.csv"
    n = omega.write_primers(path, source="subramanian2018", enzyme="BsaI")
    df = pd.read_csv(path)
    assert list(df.columns) == ["fwd_name", "fwd_sequence", "rev_name", "rev_sequence"]
    assert n == 164 // 2               # 164 screened primers -> 82 orthogonal pairs
    assert len(df) == n
    # capped by n_pairs
    omega.write_primers(path, n_pairs=5)
    assert len(pd.read_csv(path)) == 5


# --- command we build ---------------------------------------------------------

def test_build_command_carries_params(tmp_path):
    params = OmegaParams(njunctions=40, nopt_steps=500, add_primers=False)
    cmd = omega.build_command(
        "python", Path("/omega/code/omega.py"),
        tmp_path / "s.fasta", tmp_path / "p.csv", tmp_path / "out", params,
    )
    assert cmd[:3] == ["python", "/omega/code/omega.py", "genes"]
    assert cmd[cmd.index("--njunctions") + 1] == "40"
    assert cmd[cmd.index("--nopt_steps") + 1] == "500"
    assert cmd[cmd.index("--add_primers") + 1] == "false"
    assert cmd[cmd.index("--output_dir") + 1] == str(tmp_path / "out")
    assert "--opt_seeds" not in cmd                       # omitted when None


def test_build_command_opt_seeds(tmp_path):
    params = OmegaParams(njunctions=20, opt_seeds=[1, 2, 3])
    cmd = omega.build_command("python", Path("omega.py"),
                              tmp_path / "s", tmp_path / "p", tmp_path / "o", params)
    assert cmd[cmd.index("--opt_seeds") + 1] == "[1,2,3]"


# --- outputs we parse ---------------------------------------------------------

def test_parse_output_reads_three_csvs(tmp_path):
    _fake_run_factory()(["x", "--input_seqs", str(tmp_path / "s.fasta"),
                         "--output_dir", str(tmp_path / "out")], cwd=str(tmp_path))
    res = omega.parse_output(tmp_path / "out")
    assert len(res.oligos) == 1 and len(res.genes) == 1 and len(res.pools) == 1
    specs = res.design_specs()
    assert specs["n_pools"] == 1 and specs["pool_seeds"] == [42]


def test_parse_output_missing_raises(tmp_path):
    (tmp_path / "out").mkdir()
    with pytest.raises(FileNotFoundError, match="output is missing"):
        omega.parse_output(tmp_path / "out")


# --- resolution ---------------------------------------------------------------

def test_resolve_omega_needs_location(monkeypatch):
    monkeypatch.delenv("OMEGA_HOME", raising=False)
    with pytest.raises(ValueError, match="OMEGA location unknown"):
        omega._resolve_omega(None, None)


def test_resolve_omega_needs_cli(tmp_path, monkeypatch):
    monkeypatch.delenv("OMEGA_HOME", raising=False)
    with pytest.raises(FileNotFoundError, match="No OMEGA CLI"):
        omega._resolve_omega(tmp_path, None)   # dir exists but has no code/omega.py


# --- full path with subprocess stubbed ---------------------------------------

def test_assemble_with_omega_end_to_end(lib, tmp_path, monkeypatch):
    home = tmp_path / "omega"
    (home / "code").mkdir(parents=True)
    (home / "code" / "omega.py").write_text("# stub\n")
    monkeypatch.setattr(omega.subprocess, "run", _fake_run_factory())

    result = lib.assemble_with_omega(
        OmegaParams(njunctions=20),
        omega_home=home, work_dir=tmp_path / "work",
    )
    assert len(result.oligos) == 1
    assert (tmp_path / "work" / "input_seqs.fasta").is_file()
    assert (tmp_path / "work" / "primers.csv").is_file()
    # design specs recorded on the library
    assert lib.omega is result
    assert lib.design_specs["omega"]["params"]["njunctions"] == 20
    assert lib.design_specs["omega"]["pool_seeds"] == [42]
