"""Cross-check an externally codon-optimized CDS against the library's frozen WT
reference. Built for pasting in the output of IDT's Codon Optimization Tool (the
lab's long-time default) and seeing where it agrees or diverges — in particular
whether IDT's choice would introduce a restriction/ribosome-binding motif that our
design deliberately avoids (IDT's optimizer doesn't know about your Golden Gate
enzymes).

Both sequences are compared as coding regions only (no adaptors). Submit the
**truncated** protein to IDT — the one the reference encodes (``spec.truncated_sequence``);
the comparison warns if the pasted sequence encodes a different protein.
"""
from __future__ import annotations

import re
from dataclasses import dataclass, field


@dataclass
class ReferenceComparison:
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
            note = (f"{label} encodes {len(prot_other)} aa but the reference is {len(prot_ref)} aa "
                    f"— submit the truncated protein (truncation={spec.truncation}), not the full one.")
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
