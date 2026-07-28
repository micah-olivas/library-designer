"""Cross-check an externally codon-optimized CDS against the library's frozen WT
reference. Built for pasting in the output of IDT's Codon Optimization Tool (the
lab's long-time default) and seeing where it agrees or diverges. In particular,
whether IDT's choice would introduce a restriction/ribosome-binding motif that our
design deliberately avoids (IDT's optimizer doesn't know about your Golden Gate
enzymes).

Both sequences are compared as coding regions only (no adaptors). Submit the
**truncated** protein to IDT, the one the reference encodes (``spec.truncated_sequence``);
the comparison warns if the pasted sequence encodes a different protein.
"""
from __future__ import annotations

import re
from dataclasses import dataclass, field


@dataclass
class ReferenceComparison:
    """Where an outside CDS agrees with the library's frozen reference and where it does
    not. Printing it gives the readable report.

    ``label`` names the outside sequence (``"IDT"`` unless you say otherwise) and
    ``n_codons`` is how many codons were compared, which is the shorter of the two
    sequences. ``codon_matches`` and the ``codon_agreement`` property say how many of
    those are identical, and ``codon_diffs`` lists the rest as ``(position, our codon,
    their codon)``, position 1-based.

    ``gc_ref`` and ``gc_other`` are the GC fraction of each sequence. ``mean_adapt_ref``
    and ``mean_adapt_other`` are the mean relative codon adaptiveness in the spec's
    species, where 1.0 would mean every codon is the most-used one for its residue.

    ``enzyme_sites`` and ``motif_hits`` count what the pasted sequence carries of the
    things the design avoids, the enzymes in ``spec.avoid_enzymes`` and the patterns in
    ``spec.avoid_patterns``. Both count the pasted sequence only, since the reference was
    optimized to have none. ``introduces_sites`` is True when any count is nonzero.

    ``protein_match`` says whether the pasted sequence translates to the reference
    protein, and ``protein_note`` says how it differs when it does not.
    """

    label: str
    n_codons: int
    codon_matches: int
    codon_diffs: list          # (position_1based, our_codon, their_codon)
    protein_match: bool
    gc_ref: float
    gc_other: float
    mean_adapt_ref: float      # mean relative codon adaptiveness (1.0 = all optimal)
    mean_adapt_other: float
    enzyme_sites: dict         # enzyme -> count in the pasted seq (enzymes we avoid)
    motif_hits: dict           # forbidden pattern -> count in the pasted seq
    protein_note: str = ""

    @property
    def codon_agreement(self) -> float:
        return self.codon_matches / self.n_codons if self.n_codons else 0.0

    @property
    def introduces_sites(self) -> bool:
        return any(self.enzyme_sites.values()) or any(self.motif_hits.values())

    def __str__(self) -> str:
        pct = f"{self.codon_agreement:.0%}"
        lines = [f"{self.label} vs reference:  "
                 + ("protein OK" if self.protein_match else "PROTEIN MISMATCH")]
        if self.protein_note:
            lines.append(f"  ! {self.protein_note}")
        lines.append(f"  codon agreement:  {self.codon_matches}/{self.n_codons} ({pct})")
        lines.append(f"  GC:               reference {self.gc_ref:.3f} | {self.label} {self.gc_other:.3f}")
        lines.append(f"  mean adaptiveness: reference {self.mean_adapt_ref:.3f} | {self.label} {self.mean_adapt_other:.3f}")
        flagged = False
        for enz, n in self.enzyme_sites.items():
            if n:
                lines.append(f"  ⚠ {enz} site(s) in {self.label} seq: {n}  (our design avoids this)")
                flagged = True
        for pat, n in self.motif_hits.items():
            if n:
                lines.append(f"  ⚠ motif /{pat}/ in {self.label} seq: {n}  (our design avoids this)")
                flagged = True
        if not flagged:
            lines.append(f"  ✓ no avoided restriction/RBS motifs in the {self.label} sequence")
        return "\n".join(lines)

    def _repr_html_(self) -> str:
        from html import escape as e
        return f"<pre style='margin:0;line-height:1.4'>{e(str(self))}</pre>"


def _clean_dna(seq: str) -> str:
    """Accept raw DNA or FASTA; strip headers/whitespace, uppercase, keep only ACGT."""
    body = "\n".join(ln for ln in seq.splitlines() if not ln.lstrip().startswith(">"))
    return "".join(c for c in body.upper() if c in "ACGT")


def _gc(seq: str) -> float:
    return (seq.count("G") + seq.count("C")) / len(seq) if seq else 0.0


def compare_reference(library, other_dna: str, label: str = "IDT") -> ReferenceComparison:
    """Cross-check an externally codon-optimized WT CDS against this library's frozen
    reference, returning a ``ReferenceComparison``. ``Library.compare_reference`` calls
    this.

    ``other_dna`` can be raw DNA or a FASTA record; headers, whitespace, and anything
    outside ACGT are stripped first. ``label`` is the name the report gives the outside
    sequence. Raises if the library has not been codon-optimized, since there is no
    reference to compare against, and raises if the pasted text holds no DNA.

    Codons are compared position by position over the shorter of the two sequences, so
    extra length on either side is ignored rather than shifting the frame. There is no
    alignment step, so an insertion or deletion makes every codon after it read as a
    difference.

    The protein is checked on its own. The pasted sequence is translated only when its
    length is a multiple of 3, and a mismatch is written into ``protein_note`` rather
    than raised, so the rest of the report is still there to read. The note says which
    kind of mismatch it is, a length that is not a multiple of 3, a protein of the wrong
    length (usually the full protein where the truncated one was wanted), or one that
    differs at some number of residues.
    """
    if library.reference is None:
        raise ValueError("Library is not codon-optimized yet, call codon_optimize() first.")
    from dnachisel import translate

    from .checks.motifs import count_enzyme_sites
    from .optimize.backbone import relative_adaptiveness

    spec = library.spec
    ref = library.reference
    other = _clean_dna(other_dna)
    if not other:
        raise ValueError("No DNA found in the pasted sequence.")
    w = relative_adaptiveness(spec.optimization.species)

    prot_ref = translate(ref)
    prot_other = translate(other) if len(other) % 3 == 0 else None
    protein_match = prot_other == prot_ref
    note = ""
    if not protein_match:
        if prot_other is None:
            note = f"{label} length {len(other)} nt is not a multiple of 3."
        elif len(prot_other) != len(prot_ref):
            note = (f"{label} encodes {len(prot_other)} aa but the reference is {len(prot_ref)} aa. "
                    f"Submit the truncated protein (truncation={spec.truncation}), not the full one.")
        else:
            d = sum(a != b for a, b in zip(prot_ref, prot_other))
            note = f"{label} encodes a protein differing at {d} residue(s) from the reference."

    n = min(len(ref), len(other)) // 3
    matches, diffs = 0, []
    for i in range(n):
        a, b = ref[i * 3:i * 3 + 3], other[i * 3:i * 3 + 3]
        if a == b:
            matches += 1
        else:
            diffs.append((i + 1, a, b))

    def _mean_adapt(seq: str) -> float:
        m = len(seq) // 3
        return sum(w.get(seq[i * 3:i * 3 + 3], 0.0) for i in range(m)) / m if m else 0.0

    return ReferenceComparison(
        label=label,
        n_codons=n,
        codon_matches=matches,
        codon_diffs=diffs,
        protein_match=protein_match,
        gc_ref=round(_gc(ref), 3),
        gc_other=round(_gc(other), 3),
        mean_adapt_ref=round(_mean_adapt(ref), 3),
        mean_adapt_other=round(_mean_adapt(other), 3),
        enzyme_sites={enz: count_enzyme_sites(other, enz) for enz in spec.avoid_enzymes},
        motif_hits={pat: len(re.findall(pat, other)) for pat in spec.avoid_patterns},
        protein_note=note,
    )
