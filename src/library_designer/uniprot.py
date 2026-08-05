"""Look a protein up in UniProt instead of pasting it.

``LibrarySpec(uniprot="P07311")`` fills ``protein_sequence`` from the UniProt entry, and
records which entry it was, so the design specs say where the sequence came from rather
than leaving a bare string of residues in the file.

The download is cached on disk, so a spec resolves once and every later run reads the
cached FASTA. That keeps re-running a spec offline and repeatable, which matters because
UniProt entries do change: the fetched sequence is stored on the spec and in the design
specs, so the library stays reproducible even after the entry is revised. Pass
``refresh=True`` to fetch again.

Only the stdlib is used for the request, so nothing new is installed to make this work.
"""
from __future__ import annotations

import os
import re
import urllib.error
import urllib.request
from dataclasses import asdict, dataclass
from datetime import datetime, timezone
from pathlib import Path

BASE_URL = "https://rest.uniprot.org/uniprotkb"

# Accessions are 6 or 10 characters, optionally with an isoform suffix (P07311-2). Checking
# the shape keeps a pasted protein sequence or a typo from becoming an HTTP request.
_ACCESSION = re.compile(r"^[A-NR-Z][0-9](?:[A-Z][A-Z0-9]{2}[0-9]){1,2}(?:-\d+)?$|"
                        r"^[OPQ][0-9][A-Z0-9]{3}[0-9](?:-\d+)?$")


@dataclass
class UniProtEntry:
    """One UniProt record, reduced to what a design needs and what its provenance needs."""

    accession: str
    entry_name: str            # e.g. ACYP1_HUMAN
    protein_name: str          # e.g. Acylphosphatase-1
    organism: str
    gene: str | None
    sequence_version: int | None
    reviewed: bool             # a Swiss-Prot (sp) entry rather than TrEMBL (tr)
    sequence: str
    fetched: str               # when the FASTA was downloaded, ISO 8601

    def record(self) -> dict:
        """The provenance dict stored on the spec and written to the design specs. Holds
        everything but the sequence, which lives on the spec as ``protein_sequence``."""
        return {k: v for k, v in asdict(self).items() if k != "sequence"}

    def __str__(self) -> str:
        sv = f", SV {self.sequence_version}" if self.sequence_version else ""
        return (f"{self.accession} ({self.entry_name}{sv}): {self.protein_name}, "
                f"{self.organism}, {len(self.sequence)} aa")


def cache_dir() -> Path:
    """Where fetched FASTA files are kept. ``LIBRARY_DESIGNER_CACHE`` overrides it."""
    base = os.environ.get("LIBRARY_DESIGNER_CACHE")
    if base:
        return Path(base).expanduser() / "uniprot"
    xdg = os.environ.get("XDG_CACHE_HOME")
    root = Path(xdg).expanduser() if xdg else Path.home() / ".cache"
    return root / "library-designer" / "uniprot"


def _download(url: str, timeout: float) -> str:
    """The one function that touches the network, kept separate so tests can replace it."""
    with urllib.request.urlopen(url, timeout=timeout) as response:
        return response.read().decode("utf-8")


def parse_fasta(text: str) -> tuple[str, str]:
    """``(header, sequence)`` from a single-record FASTA, the sequence joined and uppercased."""
    lines = [ln.strip() for ln in text.strip().splitlines() if ln.strip()]
    if not lines or not lines[0].startswith(">"):
        raise ValueError("UniProt did not return a FASTA record.")
    return lines[0][1:], "".join(lines[1:]).upper()


def parse_entry(text: str, fetched: str) -> UniProtEntry:
    """Read a UniProt FASTA into a ``UniProtEntry``.

    The header reads ``sp|P07311|ACYP1_HUMAN Acylphosphatase-1 OS=Homo sapiens OX=9606
    GN=ACYP1 PE=1 SV=2``: database, accession, entry name, protein name, then ``KEY=value``
    fields that run to the next key."""
    header, sequence = parse_fasta(text)
    if not sequence or set(sequence) - set("ACDEFGHIKLMNPQRSTVWYXBZUO*"):
        bad = "".join(sorted(set(sequence) - set("ACDEFGHIKLMNPQRSTVWYXBZUO*")))
        raise ValueError(f"UniProt returned a sequence with non-residue characters: {bad!r}.")

    db = accession = entry_name = ""
    rest = header
    if header.count("|") >= 2:
        db, accession, rest = header.split("|", 2)
    entry_name, _, described = rest.partition(" ")
    fields = dict(re.findall(r"(\w{2})=(.*?)(?=\s+\w{2}=|$)", described))
    protein_name = described.split(" OS=")[0].strip() if " OS=" in described else described.strip()
    sv = fields.get("SV")
    return UniProtEntry(
        accession=accession or entry_name,
        entry_name=entry_name,
        protein_name=protein_name,
        organism=fields.get("OS", ""),
        gene=fields.get("GN"),
        sequence_version=int(sv) if sv and sv.isdigit() else None,
        reviewed=(db == "sp"),
        sequence=sequence,
        fetched=fetched,
    )


def fetch(accession: str, *, refresh: bool = False, timeout: float = 30.0,
          cache: str | Path | None = None) -> UniProtEntry:
    """Look ``accession`` up in UniProt and return its entry.

    Served from the on-disk cache when it is already there, so only the first call for an
    accession needs the network. ``refresh=True`` downloads again and replaces the cached
    copy. Raises with a readable reason when the accession is malformed, unknown, or the
    network is unreachable, since the fix differs in each case."""
    acc = str(accession).strip().upper()
    if not _ACCESSION.match(acc):
        raise ValueError(
            f"{accession!r} does not look like a UniProt accession (e.g. P07311, or "
            "P07311-2 for an isoform). Pass protein_sequence directly if you have the "
            "residues already."
        )
    directory = Path(cache).expanduser() if cache else cache_dir()
    path = directory / f"{acc}.fasta"
    if path.is_file() and not refresh:
        stamp = datetime.fromtimestamp(path.stat().st_mtime, timezone.utc)
        return parse_entry(path.read_text(), stamp.isoformat(timespec="seconds"))

    url = f"{BASE_URL}/{acc}.fasta"
    try:
        text = _download(url, timeout)
    except urllib.error.HTTPError as exc:
        if exc.code == 404:
            raise ValueError(f"No UniProt entry {acc!r} ({url} returned 404).") from exc
        raise ValueError(f"UniProt returned HTTP {exc.code} for {acc!r} ({url}).") from exc
    except urllib.error.URLError as exc:
        raise ValueError(
            f"Could not reach UniProt to look up {acc!r} ({exc.reason}). Set "
            "protein_sequence on the spec to work offline."
        ) from exc

    entry = parse_entry(text, datetime.now(timezone.utc).isoformat(timespec="seconds"))
    try:
        directory.mkdir(parents=True, exist_ok=True)
        path.write_text(text)
    except OSError:
        pass          # a read-only cache directory is no reason to fail the lookup
    return entry
