"""Sequence-set generator: a library of independent full-length proteins.

Where ``SubstitutionScan`` walks one wild-type protein and makes single mutants,
``SequenceSet`` takes a set of already-distinct proteins, orthologs, generative
designs, deep multi-mutants, and treats each as its own member. Every member is
codon-optimized on its own (see ``optimize/independent.py``); there is no shared
reference. This is the input a whole-gene assembler such as OMEGA expects, where
each member becomes a complete gene built from many oligos.

Typical use::

    from library_designer import LibrarySpec, SequenceSet

    spec = LibrarySpec(name="my_designs", platform="pooled")
    lib = SequenceSet(spec, proteins={
        "ortholog_ec": "MKAIL...",
        "design_01":   "MKGIL...",
    }).generate().codon_optimize()

    # or load them from a FASTA of protein sequences:
    lib = SequenceSet.from_fasta(spec, "designs.faa").generate().codon_optimize()
"""
from __future__ import annotations

import re
from pathlib import Path
from typing import Mapping

import pandas as pd

from ..library import Library
from ..spec import LibrarySpec

_AMINO_ACIDS = frozenset("ACDEFGHIKLMNPQRSTVWY")
# Characters that break a FASTA header (``>``, whitespace) or the uSort-M contract
# (``/``, ``|``); the same set io.py and integrations/omega.py reject in names.
_BAD_NAME = re.compile(r"[/|>\s]")


def _validate_protein(name: str, protein: str) -> str:
    seq = protein.strip().upper()
    if not seq:
        raise ValueError(f"Member {name!r} has an empty protein sequence.")
    bad = sorted(set(seq) - _AMINO_ACIDS)
    if bad:
        raise ValueError(
            f"Member {name!r} has non-amino-acid character(s) {bad} in its protein "
            "sequence. Give one-letter residues from the 20 canonical amino acids, "
            "with no stop ('*') or ambiguity ('X') codes."
        )
    return seq


def _read_fasta(path: str | Path) -> dict[str, str]:
    """Parse a FASTA of protein sequences into ``{name: sequence}``. The record name
    is the header up to the first whitespace."""
    name: str | None = None
    chunks: dict[str, list[str]] = {}
    order: list[str] = []
    for line in Path(path).read_text().splitlines():
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            name = line[1:].split()[0] if len(line) > 1 else ""
            if name in chunks:
                raise ValueError(f"Duplicate FASTA record name {name!r} in {path}.")
            chunks[name] = []
            order.append(name)
        elif name is None:
            raise ValueError(f"{path}: sequence data before the first '>' header.")
        else:
            chunks[name].append(line)
    if not order:
        raise ValueError(f"No FASTA records found in {path}.")
    return {n: "".join(chunks[n]) for n in order}


class SequenceSet:
    """Build a library from a set of full-length protein sequences.

    ``proteins`` maps each member's name to its one-letter protein sequence. Names
    must be unique and free of characters that break a FASTA header or the uSort-M
    handoff (``/ | >`` and whitespace). Sequences are the coding proteins only; the
    spec's adaptors (if any) flank every member, as for a scan.
    """

    def __init__(self, spec: LibrarySpec, proteins: Mapping[str, str]):
        self.spec = spec
        self.proteins = dict(proteins)

    @classmethod
    def from_fasta(cls, spec: LibrarySpec, path: str | Path) -> "SequenceSet":
        """Build from a FASTA file of protein sequences (one record per member)."""
        return cls(spec, _read_fasta(path))

    def generate(self) -> Library:
        spec = self.spec
        if not self.proteins:
            raise ValueError("SequenceSet needs at least one protein sequence.")

        names = [str(n) for n in self.proteins]
        bad = [n for n in names if _BAD_NAME.search(n)]
        if bad:
            raise ValueError(
                "Member names contain characters that break a FASTA header or the "
                "uSort-M contract (/ | > whitespace): " + ", ".join(bad[:10])
            )

        a5, a3 = spec.adaptor_5.upper(), spec.adaptor_3.upper()
        rows: list[dict] = []
        for name, protein in self.proteins.items():
            rows.append(
                {
                    "name": str(name),
                    "protein": _validate_protein(str(name), protein),
                    # Scan-specific columns kept as NA so every exporter, QC check, and
                    # summary that reads them stays happy across both library kinds.
                    "position": pd.NA,
                    "wt_residue": pd.NA,
                    "mut_residue": pd.NA,
                    "codon": pd.NA,
                    "mut_index": pd.NA,
                    "adaptor_5": a5,
                    "adaptor_3": a3,
                }
            )

        df = pd.DataFrame(rows)
        return Library(df, spec=spec, kind="sequence_set")
