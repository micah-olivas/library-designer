"""Orthogonal primer sets for tiled assembly.

Each tile in a pooled tiled library carries a tile-specific primer pair so that its
sublibrary can be selectively amplified out of the shared pool. Those primers must be
*mutually orthogonal* (no cross-amplification), which is why we ship validated sets
rather than designing primers ad hoc:

- ``subramanian2018`` (default), 165 experimentally validated, mutually orthogonal
  20-mers (Subramanian, Russ & Ranganathan, *Synth. Biol.* 2018; PMC7445780). A
  flat *pool*: any two are orthogonal, so the tiler draws two per tile.
- ``gck800``, 7 synthetic forward/reverse pairs kept as a worked example for the
  glucokinase tiled design. A *paired* set. These are illustrative, not validated for
  orthogonality, so screen your own before ordering.

Sets can also be loaded from a CSV path: a pool has a ``sequence`` column (plus an
optional ``primer_id``); a paired set has ``forward`` and ``reverse`` columns (plus an
optional ``pair_id``). Primers are given 5'->3' as ordered.

Whatever the set, on load we drop any primer carrying the Golden Gate enzyme's
recognition site (on either strand), a published set is validated for orthogonality,
not for *your* enzyme, so a few primers routinely collide with it.
"""
from __future__ import annotations

import io
from dataclasses import dataclass, field
from importlib import resources
from pathlib import Path

import pandas as pd

from .checks.motifs import contains_enzyme_site


@dataclass
class PrimerSet:
    """A screened primer set. ``kind`` is ``"pool"`` (interchangeable orthogonal
    primers; the tiler pairs them up) or ``"paired"`` (explicit forward/reverse)."""

    name: str
    kind: str
    primers: list[tuple[str, str]] = field(default_factory=list)     # pool: (id, seq 5'->3')
    pairs: list[tuple[str, str, str]] = field(default_factory=list)  # paired: (id, fwd, rev), 5'->3'
    dropped: list[str] = field(default_factory=list)                 # ids removed for carrying the enzyme site

    @property
    def capacity(self) -> int:
        """Maximum number of tiles this set can flank."""
        return len(self.primers) // 2 if self.kind == "pool" else len(self.pairs)


def _read_set_csv(source: str) -> pd.DataFrame:
    """Load a primer-set CSV, resolving ``source`` as a file path first, then as the
    name of a set bundled under ``library_designer/data/primers/``."""
    p = Path(source)
    if p.exists():
        return pd.read_csv(p)
    res = resources.files("library_designer").joinpath(f"data/primers/{source}.csv")
    if not res.is_file():
        raise FileNotFoundError(
            f"Primer set {source!r} is neither a file nor a bundled set "
            "(known bundled sets: 'subramanian2018', 'gck800')."
        )
    return pd.read_csv(io.StringIO(res.read_text()))


def load_primer_set(source: str = "subramanian2018", enzyme: str = "BsaI") -> PrimerSet:
    """Load and screen a primer set. Primers carrying ``enzyme``'s recognition site
    (either strand) are dropped and listed in ``PrimerSet.dropped``."""
    df = _read_set_csv(source)
    cols = {c.lower().strip(): c for c in df.columns}

    if "sequence" in cols:  # flat orthogonal pool
        id_col = cols.get("primer_id")
        primers, dropped = [], []
        for i, row in df.iterrows():
            seq = str(row[cols["sequence"]]).strip().upper()
            pid = str(row[id_col]) if id_col else f"P{i}"
            if contains_enzyme_site(seq, enzyme):
                dropped.append(pid)
            else:
                primers.append((pid, seq))
        return PrimerSet(name=source, kind="pool", primers=primers, dropped=dropped)

    if "forward" in cols and "reverse" in cols:  # explicit fwd/rev pairs
        id_col = cols.get("pair_id")
        pairs, dropped = [], []
        for i, row in df.iterrows():
            fwd = str(row[cols["forward"]]).strip().upper()
            rev = str(row[cols["reverse"]]).strip().upper()
            pid = str(row[id_col]) if id_col else f"pair{i}"
            if contains_enzyme_site(fwd, enzyme) or contains_enzyme_site(rev, enzyme):
                dropped.append(pid)
            else:
                pairs.append((pid, fwd, rev))
        return PrimerSet(name=source, kind="paired", pairs=pairs, dropped=dropped)

    raise ValueError(
        f"Primer set {source!r} must have a 'sequence' column (pool) or "
        "'forward'+'reverse' columns (paired)."
    )
