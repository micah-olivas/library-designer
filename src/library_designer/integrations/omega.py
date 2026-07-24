"""OMEGA integration: assemble long genes from oligopools via Golden Gate.

OMEGA (RomeroLab, https://github.com/RomeroLab/omega) takes a set of codon-optimized
full-length genes, computationally fragments each into short Golden Gate oligos, and
picks fragmentation junctions whose Type IIS overhangs maximize predicted ligation
fidelity. It is the natural downstream companion to a variant generator: where
``.tile()`` rides one mutant tile into a vector holding the rest of the WT CDS,
OMEGA assembles *whole, distinct* genes from many oligos, so it fits libraries of
fully independent sequences (homologs, generative designs, deep multi-mutants).

**library_designer does not depend on, import, or bundle OMEGA.** OMEGA is GPL-3.0;
library_designer is MIT. To keep them separate works, this module runs a *separately
installed* OMEGA as an arm's-length subprocess: we write its inputs (a FASTA of
codon-optimized CDSs and a primer CSV), invoke its CLI, and read back its three
output CSVs. No OMEGA code is ever loaded into this process.

Point the runner at your OMEGA checkout with ``omega_home`` (or the ``OMEGA_HOME``
environment variable). If OMEGA lives in its own Python environment (it ships a
conda env, not a pip package), name that interpreter with ``omega_python`` (or
``OMEGA_PYTHON``); it defaults to ``"python"``.

Typical use::

    from library_designer.integrations.omega import OmegaParams

    lib = SubstitutionScan(spec).generate().codon_optimize()
    result = lib.assemble_with_omega(
        OmegaParams(njunctions=40),
        omega_home="~/repos/omega",
    )
    result.oligos.to_csv("oligo_order.csv", index=False)
"""
from __future__ import annotations

import os
import re
import subprocess
import tempfile
from dataclasses import asdict, dataclass, field
from pathlib import Path

import pandas as pd

from ..primers import load_primer_set

# Characters that break a FASTA header (``>``, whitespace) or the downstream
# uSort-M contract (``/``, ``|``). Kept identical to the export validation in io.py.
_BAD_NAME = re.compile(r"[/|>\s]")

# OMEGA writes these three files into its ``--output_dir``.
_OUTPUT_FILES = {
    "oligos": "oligo_order.csv",           # the oligos to order
    "genes": "optimization_results.csv",   # per-gene: sequence, oligo, primers, fidelity
    "pools": "pool_stats.csv",             # per-pool: fidelity, sites, primers, seed
}


@dataclass
class OmegaParams:
    """Design knobs for an OMEGA run, mirroring its ``genes`` CLI. Captured in the
    design specs so an assembly is reproducible. Only ``njunctions`` (the number of
    Golden Gate sites per subpool) is required; the rest match OMEGA's defaults.

    Runtime *location* (where OMEGA is installed, which interpreter runs it) is kept
    out of this record on purpose, it is passed to the runner separately so the
    design parameters stay portable across machines.
    """

    njunctions: int                                 # Golden Gate sites per subpool (the GG budget)
    upstream_bbsite: str = "AATG"                   # 5' backbone overhang (MoClo CDS start context)
    downstream_bbsite: str = "TTAG"                 # 3' backbone overhang (MoClo CDS stop context)
    enzyme: str = "BsaI"                            # Type IIS enzyme
    oligo_len: int = 300                            # oligo length (bp); fragments are padded up to this
    min_size: int = 40                              # minimum fragment size (bp)
    nopt_steps: int = 1000                          # simulated-annealing steps per run
    nopt_runs: int = 5                              # independent SA runs (best kept)
    njobs: int = 1                                  # parallel worker processes
    ligation_data: str = "T4_18h_37C"              # empirical fidelity dataset (see OMEGA)
    optimization: str = "simulated_annealing"       # or "greedy"
    add_primers: bool = True                        # append subpool amplification primers to oligos
    pad_oligos: bool = True                         # pad fragments to a uniform length
    opt_seeds: list[int] | None = None              # fixed per-run seeds (else OMEGA chooses)


@dataclass
class OmegaResult:
    """Parsed output of an OMEGA run: three DataFrames plus run metadata."""

    oligos: pd.DataFrame                # oligo_order.csv
    genes: pd.DataFrame                 # optimization_results.csv
    pools: pd.DataFrame                 # pool_stats.csv
    output_dir: Path
    command: list[str] = field(default_factory=list)
    stdout: str = ""

    def __repr__(self) -> str:
        return (
            f"<OmegaResult: {len(self.genes)} gene(s), {len(self.pools)} pool(s), "
            f"{len(self.oligos)} oligo(s) -> {self.output_dir}>"
        )

    def design_specs(self) -> dict:
        """A JSON-serializable record for ``Library.design_specs``: run sizes, the
        per-pool seeds OMEGA used (found case-insensitively), and the command."""
        specs: dict = {
            "n_genes": int(len(self.genes)),
            "n_pools": int(len(self.pools)),
            "n_oligos": int(len(self.oligos)),
            "command": list(self.command),
        }
        seed_col = next(
            (c for c in self.pools.columns if "seed" in c.lower()), None
        )
        if seed_col is not None:
            specs["pool_seeds"] = self.pools[seed_col].tolist()
        return specs


def _validate_names(names: pd.Series) -> None:
    bad = names[names.str.contains(_BAD_NAME)]
    if not bad.empty:
        raise ValueError(
            "Gene names contain characters that break a FASTA header or the uSort-M "
            "contract (/ | > whitespace): " + ", ".join(bad.head(10))
        )
    dups = names[names.duplicated()].unique()
    if len(dups):
        raise ValueError(f"Duplicate gene names would collide in the FASTA: {', '.join(dups[:10])}")


def write_fasta(library, path: str | Path) -> int:
    """Write each library member's full codon-optimized CDS as a FASTA record for
    OMEGA (``>name`` / uppercase DNA). Returns the number of records written.

    Only the coding region (``variable_dna``) is emitted, not the adaptors, OMEGA
    adds its own Type IIS sites, backbone overhangs, and amplification primers.
    Variants that failed optimization (NA ``variable_dna``) are skipped.
    """
    if "variable_dna" not in library.df.columns:
        raise ValueError("Library is not codon-optimized yet, call codon_optimize() first.")
    if getattr(library, "tiles", None) is not None:
        raise ValueError(
            "This library is already laid out for tiled assembly. OMEGA is an "
            "alternative whole-gene assembly route, use one or the other, not both."
        )

    df = library.df
    names = df["name"].astype(str)
    _validate_names(names)

    lines: list[str] = []
    for name, cds in zip(names, df["variable_dna"]):
        if not isinstance(cds, str):
            continue
        lines.append(f">{name}")
        lines.append(cds.upper())
    if not lines:
        raise ValueError("No optimized sequences to write (every variable_dna is NA).")

    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text("\n".join(lines) + "\n")
    return len(lines) // 2


def write_primers(path: str | Path, source: str = "subramanian2018",
                  enzyme: str = "BsaI", n_pairs: int | None = None) -> int:
    """Write subpool amplification primers in OMEGA's format
    (``fwd_name,fwd_sequence,rev_name,rev_sequence``) from a bundled or file primer
    set. Returns the number of pairs written.

    OMEGA recommends the Subramanian 2018 set, which library_designer already ships,
    so by default this reuses the same screened primers as tiled assembly. A flat
    orthogonal pool is chunked into consecutive pairs; a paired set is used as given.
    """
    ps = load_primer_set(source, enzyme)
    if ps.kind == "pool":
        prs = ps.primers
        rows = [
            {"fwd_name": prs[i][0], "fwd_sequence": prs[i][1],
             "rev_name": prs[i + 1][0], "rev_sequence": prs[i + 1][1]}
            for i in range(0, len(prs) - 1, 2)
        ]
    else:  # paired
        rows = [
            {"fwd_name": pid, "fwd_sequence": fwd, "rev_name": pid, "rev_sequence": rev}
            for pid, fwd, rev in ps.pairs
        ]
    if n_pairs is not None:
        rows = rows[:n_pairs]
    if not rows:
        raise ValueError(f"Primer set {source!r} yielded no usable pairs.")

    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(p, index=False)
    return len(rows)


def _resolve_omega(omega_home: str | Path | None,
                   omega_python: str | None) -> tuple[Path, Path, str]:
    """Locate the OMEGA CLI (``code/omega.py``) and the interpreter that runs it."""
    home = omega_home or os.environ.get("OMEGA_HOME")
    if not home:
        raise ValueError(
            "OMEGA location unknown. Pass omega_home=... or set OMEGA_HOME to your "
            "OMEGA checkout (https://github.com/RomeroLab/omega). library_designer "
            "does not bundle OMEGA."
        )
    home = Path(home).expanduser().resolve()
    script = home / "code" / "omega.py"
    if not script.is_file():
        raise FileNotFoundError(
            f"No OMEGA CLI at {script} (expected <omega_home>/code/omega.py). "
            f"Is omega_home={home} the root of an OMEGA checkout?"
        )
    python = omega_python or os.environ.get("OMEGA_PYTHON") or "python"
    return home, script, python


def build_command(python: str, script: Path, fasta: Path, primers: Path,
                  output_dir: Path, params: OmegaParams,
                  config: str | Path | None = None) -> list[str]:
    """Assemble the ``omega.py genes ...`` argv for a run. Paths are passed absolute
    so they resolve regardless of the working directory OMEGA runs in."""
    cmd = [python, str(script), "genes"]
    if config is not None:
        cmd += ["--config", str(config)]
    cmd += [
        "--input_seqs", str(fasta),
        "--primers", str(primers),
        "--output_dir", str(output_dir),
        "--njunctions", str(params.njunctions),
        "--upstream_bbsite", params.upstream_bbsite,
        "--downstream_bbsite", params.downstream_bbsite,
        "--enzyme", params.enzyme,
        "--oligo_len", str(params.oligo_len),
        "--min_size", str(params.min_size),
        "--nopt_steps", str(params.nopt_steps),
        "--nopt_runs", str(params.nopt_runs),
        "--njobs", str(params.njobs),
        "--ligation_data", params.ligation_data,
        "--optimization", params.optimization,
        "--add_primers", str(params.add_primers).lower(),
        "--pad_oligos", str(params.pad_oligos).lower(),
    ]
    if params.opt_seeds is not None:
        cmd += ["--opt_seeds", "[" + ",".join(str(s) for s in params.opt_seeds) + "]"]
    return cmd


def run_command(cmd: list[str], cwd: Path) -> subprocess.CompletedProcess:
    """Invoke OMEGA. Runs from ``cwd`` (the OMEGA checkout) so its bundled ligation
    data resolves. Raises with captured stderr on a non-zero exit."""
    proc = subprocess.run(cmd, cwd=str(cwd), capture_output=True, text=True)
    if proc.returncode != 0:
        raise RuntimeError(
            f"OMEGA exited with status {proc.returncode}.\n"
            f"command: {' '.join(cmd)}\n"
            f"--- stderr ---\n{proc.stderr}"
        )
    return proc


def parse_output(output_dir: str | Path) -> OmegaResult:
    """Read OMEGA's three output CSVs from ``output_dir`` into an ``OmegaResult``."""
    out = Path(output_dir)
    paths = {k: out / v for k, v in _OUTPUT_FILES.items()}
    missing = [str(p) for p in paths.values() if not p.is_file()]
    if missing:
        raise FileNotFoundError(
            "OMEGA finished but expected output is missing: " + ", ".join(missing)
        )
    return OmegaResult(
        oligos=pd.read_csv(paths["oligos"]),
        genes=pd.read_csv(paths["genes"]),
        pools=pd.read_csv(paths["pools"]),
        output_dir=out,
    )


def assemble(library, params: OmegaParams, *,
             omega_home: str | Path | None = None,
             omega_python: str | None = None,
             primer_source: str = "subramanian2018",
             work_dir: str | Path | None = None,
             config: str | Path | None = None) -> OmegaResult:
    """Run OMEGA end-to-end on a codon-optimized library and return the parsed result.

    Writes the FASTA and primer CSV into ``work_dir`` (a fresh temp dir if omitted),
    invokes the separately-installed OMEGA CLI, and reads its three output CSVs from
    ``work_dir/output``.
    """
    home, script, python = _resolve_omega(omega_home, omega_python)

    work = Path(work_dir) if work_dir is not None else Path(tempfile.mkdtemp(prefix="omega_"))
    work.mkdir(parents=True, exist_ok=True)
    fasta = work / "input_seqs.fasta"
    primers = work / "primers.csv"
    output_dir = work / "output"
    output_dir.mkdir(parents=True, exist_ok=True)

    write_fasta(library, fasta)
    write_primers(primers, source=primer_source, enzyme=params.enzyme)

    cmd = build_command(python, script, fasta, primers, output_dir, params, config=config)
    proc = run_command(cmd, cwd=home)

    result = parse_output(output_dir)
    result.command = cmd
    result.stdout = proc.stdout
    return result
