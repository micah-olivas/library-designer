"""Single-substitution scan, the MBO-038 amber + missense workflow."""
from __future__ import annotations

import pandas as pd

from ..library import Library
from ..spec import LibrarySpec

_AMINO_ACIDS = frozenset("ACDEFGHIKLMNPQRSTVWY")
_DNA = frozenset("ACGT")


def _classify(entry: str) -> tuple[str, str]:
    """Resolve a substitution entry to ``(residue_symbol, fixed_codon)``.

    - amino-acid letter (``"A"``)  -> ``("A", "")``   : optimizer picks the codon
    - DNA codon (``"TAG"``)        -> ``("*", "TAG")`` : placed verbatim + protected

    The residue symbol is what appears in the variant name and protein sequence
    (``translate("TAG") == "*"``), so amber variants stay ``K7*``.
    """
    from dnachisel import translate

    e = entry.upper()
    if len(e) == 1 and e in _AMINO_ACIDS:
        return e, ""
    if len(e) == 3 and set(e) <= _DNA:
        return translate(e), e
    raise ValueError(
        f"Substitution {entry!r} is neither an amino acid (1 letter) "
        f"nor a DNA codon (3 of ACGT)."
    )


class SubstitutionScan:
    """Introduce each requested residue/codon at every position of the designed protein
    (``spec.designed_sequence``, the whole protein unless ``truncation`` is set).

    - Positions where the substitution matches the wild-type residue are skipped.
    - Positions in ``spec.mask_positions`` get no variants, though the residue is still
      encoded and synthesized on every oligo.
    - A codon substitution is placed verbatim and protected during optimization;
      an amino-acid substitution is codon-optimized freely.
    - The wild-type sequence is appended as a control (name ``"WT"``).

    Variant names use full-protein numbering (an N-terminal ``truncation`` added back),
    e.g. ``"K7A"`` or ``"K7*"``.
    """

    def __init__(self, spec: LibrarySpec):
        self.spec = spec

    def generate(self) -> Library:
        """Build the variant table, one row per substitution per position, plus a WT control.

        Positions where the substitution is already the wild-type residue are skipped, so an
        alanine scan of a protein carrying 12 alanines gives 12 fewer members than there are
        positions. A ``*`` in the input protein is skipped as a position too, and so is any
        position in ``spec.mask_positions``, which leaves the residue in the construct but out
        of the scan. Masking every scannable position raises rather than handing back a library
        of nothing but the wild-type control. Names carry full-protein numbering even when
        ``truncation`` is set, so ``K7A`` means residue 7 of
        the protein you supplied. Two substitutions that resolve to the same residue symbol
        at one position would collide, and that raises rather than silently dropping one.

        Sequences are not built here. The returned library has no ``variable_dna`` column
        until ``codon_optimize()`` runs.
        """
        spec = self.spec
        seq = spec.designed_sequence
        masked = spec.masked
        subs = [_classify(s) for s in spec.substitutions]
        a5, a3 = spec.adaptor_5.upper(), spec.adaptor_3.upper()  # explicit region columns
        rows: list[dict] = []

        for i, wt in enumerate(seq):
            if wt == "*":
                continue
            position = i + 1 + spec.numbering_offset   # 1-indexed on the full protein
            if position in masked:
                # Still encoded and still synthesized on every oligo, just not varied.
                continue
            for symbol, codon in subs:
                if symbol == wt:
                    continue
                rows.append(
                    {
                        "name": f"{wt}{position}{symbol}",
                        "protein": seq[:i] + symbol + seq[i + 1:],
                        "position": position,
                        "wt_residue": wt,
                        "mut_residue": symbol,
                        "codon": codon or pd.NA,   # literal codon to protect, if any
                        "mut_index": i,            # 0-based index in the designed protein
                        "adaptor_5": a5,
                        "adaptor_3": a3,
                    }
                )

        rows.append(
            {
                "name": "WT",
                "protein": seq,
                "position": pd.NA,
                "wt_residue": pd.NA,
                "mut_residue": pd.NA,
                "codon": pd.NA,
                "mut_index": pd.NA,
                "adaptor_5": a5,
                "adaptor_3": a3,
            }
        )

        if masked and len(rows) == 1:            # only the WT control was appended
            raise ValueError(
                f"mask_positions {sorted(masked)} leaves no position to scan. Mask fewer "
                "positions, or widen spec.substitutions."
            )
        df = pd.DataFrame(rows)
        dupes = df["name"][df["name"].duplicated()].unique()
        if len(dupes):
            raise ValueError(
                "Substitutions produced duplicate variant names "
                f"({', '.join(map(str, dupes[:5]))}). Two entries map to the same "
                "residue symbol at a position (e.g. two different stop codons)."
            )
        return Library(df, spec=spec)
