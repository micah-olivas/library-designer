"""Cross-check an externally codon-optimized CDS against the library's frozen WT
reference. Built for pasting in the output of IDT's Codon Optimization Tool (the
lab's long-time default) and seeing where it agrees or diverges. In particular,
whether IDT's choice would introduce a restriction/ribosome-binding motif that our
design deliberately avoids (IDT's optimizer doesn't know about your Golden Gate
enzymes).

Both sequences are compared as coding regions only (no adaptors). Submit either the protein
the reference encodes (``spec.designed_sequence``) or, with a truncation set, the full-length
protein: a pasted full-length CDS is trimmed to the designed region before comparing, so the
two are not read out of frame. The comparison warns if the pasted sequence encodes a
different protein than either form.
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
    things the spec avoids, the enzymes in ``spec.avoid_enzymes`` and the patterns in
    ``spec.avoid_patterns``. Both count the pasted sequence only, since the reference was
    optimized to have none. ``introduces_sites`` is True when any count is nonzero.

    ``protein_match`` says whether the pasted sequence translates to the reference protein.
    ``protein_note`` says how it differs when it does not, and when it does match after being
    trimmed to the designed region, that it was.
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
        """Fraction of compared codons that match, 0.0 when there was nothing to compare.
        Codons are compared position by position with no alignment, so an indel makes
        everything downstream count as a difference."""
        return self.codon_matches / self.n_codons if self.n_codons else 0.0

    @property
    def introduces_sites(self) -> bool:
        """True when the pasted sequence carries any restricted enzyme site or avoid-motif
        that our reference was optimized to avoid."""
        return any(self.enzyme_sites.values()) or any(self.motif_hits.values())

    def __str__(self) -> str:
        pct = f"{self.codon_agreement:.0%}"
        lines = [f"{self.label} vs reference:  "
                 + ("protein OK" if self.protein_match else "PROTEIN MISMATCH")]
        if self.protein_note:
            lines.append(f"  {'note:' if self.protein_match else '!'} {self.protein_note}")
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


def trim_to_designed(spec, dna: str) -> tuple[str, bool]:
    """``(dna, trimmed)``: an outside CDS narrowed to the designed region when it is the
    full-length one and the spec truncates.

    One rule, used wherever an outside sequence meets a truncated reference (the comparison
    and the codon-usage overlay), so they cannot disagree about the frame. Off by default: a
    sequence that is already the designed region, or that encodes a different protein, comes
    back untouched, and the caller reports the mismatch as it would have anyway.
    """
    if not spec.truncation or len(dna) != 3 * len(spec.protein_sequence):
        return dna, False
    from dnachisel import translate

    if translate(dna) != spec.protein_sequence:
        return dna, False
    drop = spec.truncation * 3
    return (dna[drop:] if spec.terminus == "N" else dna[:-drop]), True


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
    length, or one that differs at some number of residues. A full-length CDS against a
    truncated reference is not a mismatch: it is trimmed to the designed region and the note
    records that it was.
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
    note = ""
    # A pasted full-length CDS against a truncated reference is trimmed the same way
    # ``build_reference`` trims ``spec.cds``: the truncation applies to the protein and to any
    # DNA encoding it. Otherwise the two would be compared out of frame by the truncation,
    # which reads as almost every codon differing.
    other, trimmed = trim_to_designed(spec, other)
    if trimmed:
        prot_other = translate(other)
        note = (f"{label} was the full-length CDS, trimmed by {spec.truncation} codon(s) at "
                f"the {spec.terminus} terminus to match the designed region.")
    protein_match = prot_other == prot_ref
    if not protein_match:
        if prot_other is None:
            note = f"{label} length {len(other)} nt is not a multiple of 3."
        elif len(prot_other) != len(prot_ref):
            # Only blame truncation when the spec truncates; otherwise the two are just
            # different proteins and telling the user to truncate would send them wrong.
            hint = (
                f"Submit either the designed region ({spec.protein_description()}) or the "
                "full-length protein, which is trimmed to match."
                if spec.truncation
                else "Submit the protein the reference encodes (spec.designed_sequence)."
            )
            note = (f"{label} encodes {len(prot_other)} aa but the reference is "
                    f"{len(prot_ref)} aa. {hint}")
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
