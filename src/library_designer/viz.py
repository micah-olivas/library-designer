"""QC visualization: codon usage, the tile layout, and overhang homology.

The codon-usage figure plots relative codon adaptiveness (w = freq / best-synonymous-freq;
1.0 = optimal codon) along the CDS: the WT CDS as a line, and each variant's stamped
codon as a point coloured by sublibrary. Points below the threshold flag codons the
motif-avoidance had to compromise on. Amber (`*`) stamps sit low, since a stop codon is
pinned rather than optimized.

The overhang figure is the review surface for Golden Gate specificity, every fused overhang
against every other (see `checks/overhangs.py`).
"""
from __future__ import annotations

# Resolution every figure is built at, and the one a saved image inherits unless the caller
# says otherwise. Matplotlib defaults to 100, which is legible for a line plot but soft for a
# dense figure like the codon map, where a cell is a couple of pixels wide.
DEFAULT_DPI = 200

# Okabe-Ito colourblind-safe categorical palette (for sublibraries).
_PALETTE = ["#E69F00", "#56B4E9", "#009E73", "#CC79A7", "#0072B2", "#D55E00", "#F0E442", "#999999"]


def _new_figure(**kwargs):
    """A Figure created through pyplot and immediately closed.

    Importing pyplot registers the notebook inline PNG formatter, so the returned figure
    renders as an image without needing a ``%matplotlib inline`` in the notebook (a bare
    object-oriented Figure otherwise shows only as ``<Figure ...>`` text). Closing it
    removes it from pyplot's registry so it is not auto-shown a second time on top of the
    returned value. A closed figure still accepts axes/content and still supports
    ``savefig``. ``dpi`` sets how many pixels the figure renders to, inline and on disk, so
    raising it is what sharpens a dense figure rather than resizing it.
    """
    try:
        import matplotlib.pyplot as plt
    except ImportError as exc:
        raise ImportError(
            "Plotting requires matplotlib, a base dependency. "
            "Reinstall library-designer if it is missing."
        ) from exc
    fig = plt.figure(**kwargs)
    plt.close(fig)
    return fig


def codon_usage_figure(library, metric: str = "frequency", low_usage: float = 0.1,
                       compare: str | None = None, compare_label: str = "IDT",
                       dpi: int = DEFAULT_DPI):
    """Return a matplotlib Figure of per-position codon usage along the CDS.

    ``metric="frequency"`` (default) plots absolute host codon-usage frequency. The
    WT CDS varies position to position and each substitution's codon sits at its own
    frequency, so the plot is readable under any optimization method. Rare codons
    sit near the low-usage guide. ``metric="adaptiveness"`` plots relative
    adaptiveness w (freq / best-synonymous freq); useful for spotting motif-forced
    compromises, but under ``use_best_codon`` every optimized codon is 1.0.

    ``compare`` overlays an external coding sequence (e.g. IDT's WT optimization; raw
    DNA or FASTA) as a dashed line, so you can see where its codon choices diverge from
    the WT CDS.
    """
    from .optimize.backbone import codon_frequency, relative_adaptiveness

    if library.reference is None:
        if getattr(library, "kind", "scan") == "sequence_set":
            # The plot draws each variant's stamped codon against the shared WT backbone, and
            # a sequence set has no shared backbone: its members are independent genes.
            raise ValueError(
                "The codon-usage plot compares each variant's stamped codon against one "
                "shared reference, which a sequence set does not have (every member is its "
                "own gene). Plot a member's sequence yourself, or use a SubstitutionScan."
            )
        raise ValueError("Library is not codon-optimized yet, call codon_optimize() first.")

    spec = library.spec
    if metric == "frequency":
        table = codon_frequency(spec.optimization.species)
        ylabel = f"Codon usage frequency ({spec.optimization.species})"
        guide = low_usage
    elif metric == "adaptiveness":
        table = relative_adaptiveness(spec.optimization.species)
        ylabel = "Relative codon adaptiveness (w)"
        guide = 0.5
    else:
        raise ValueError(f"metric must be 'frequency' or 'adaptiveness', got {metric!r}")

    ref = library.reference
    n_codons = len(ref) // 3
    ref_y = [table.get(ref[i * 3:i * 3 + 3], 0.0) for i in range(n_codons)]

    fig = _new_figure(figsize=(12.5, 4), dpi=dpi)
    ax = fig.subplots()
    ax.plot(range(1, n_codons + 1), ref_y, color="0.45", lw=1.5, zorder=3, label="WT CDS")

    df = library.df
    muts = df[df["mut_index"].notna() & df["variable_dna"].notna()]
    subs = sorted(muts["mut_residue"].astype(str).unique())
    colours = {s: _PALETTE[i % len(_PALETTE)] for i, s in enumerate(subs)}
    for sub, grp in muts.groupby(muts["mut_residue"].astype(str)):
        xs = grp["mut_index"].astype(int) + 1
        ys = [table.get(v[int(i) * 3:int(i) * 3 + 3], 0.0)
              for v, i in zip(grp["variable_dna"], grp["mut_index"])]
        ax.scatter(xs, ys, s=16, alpha=0.65, color=colours[sub], edgecolors="none",
                   zorder=2, label=f"to {sub}")

    if compare:
        from .compare import _clean_dna
        other = _clean_dna(compare)
        m = min(len(other) // 3, n_codons)
        other_y = [table.get(other[i * 3:i * 3 + 3], 0.0) for i in range(m)]
        ax.plot(range(1, m + 1), other_y, color="#111111", lw=1.3, ls="--",
                alpha=0.9, zorder=4, label=compare_label)

    if guide:
        ax.axhline(guide, ls="--", lw=1, color="0.8", zorder=1)
    ax.set_ylim(-0.02, 1.05)
    ax.set_xlim(0, n_codons + 1)
    ax.set_xlabel("Codon position (CDS)")
    ax.set_ylabel(ylabel)
    ax.set_title(f"Codon-usage QC, {spec.name}  ({spec.optimization.species})", loc="left")
    # Legend outside the axes on the right, one entry per row. Inside the plot it covered
    # the low-usage codons in the lower-right corner, which are the ones worth reading.
    # tight_layout reserves the right margin for it, so a saved figure is not clipped.
    ax.legend(loc="center left", bbox_to_anchor=(1.02, 0.5), ncol=1, fontsize=8,
              frameon=False, borderaxespad=0)
    fig.tight_layout(rect=(0, 0, 0.86, 1))
    return fig


def tiling_figure(library, dpi: int = DEFAULT_DPI):
    """Return a matplotlib Figure of the tile layout for a tiled library.

    Each tile is a bar spanning the coding-sequence region it covers, labelled with
    its variant count and primer pair, above a histogram of where the mutations fall.
    Shows tile sizes, boundaries, and how variants distribute across sublibraries.
    """
    if library.tiles is None:
        raise ValueError("Library is not tiled yet, call tile() first.")

    ref = library.reference
    tiles = library.tiles
    df = library.df

    fig = _new_figure(figsize=(11, 0.45 * len(tiles) + 2.6), dpi=dpi)
    ax_top, ax = fig.subplots(
        2, 1, sharex=True, gridspec_kw={"height_ratios": [1, 2.6], "hspace": 0.12},
    )
    fig.subplots_adjust(left=0.14, right=0.98, top=0.9, bottom=0.13)

    muts = df[df["mut_index"].notna()]
    positions = muts["mut_index"].astype(int) * 3 + 1
    ax_top.hist(positions, bins=min(80, max(10, len(ref) // 45)), color="0.6")
    ax_top.set_ylabel("mutations", fontsize=8)
    ax_top.set_title(
        f"Tile layout, {library.spec.name}  ({len(tiles)} tiles, {len(ref)} nt)", loc="left"
    )

    for t in tiles:
        colour = _PALETTE[t.index % len(_PALETTE)]
        # Mutants only. The tile also carries a WT control, which has no position, so the
        # histogram above leaves it out too.
        n = int(df.loc[df["tile"] == t.index, "mut_index"].notna().sum()) \
            if "tile" in df.columns else 0
        ax.barh(t.index, t.end - t.start, left=t.start, height=0.72,
                color=colour, alpha=0.9, edgecolor="white")
        ax.text((t.start + t.end) / 2, t.index, f"tile {t.index}  ·  {n} variants",
                ha="center", va="center", fontsize=8, color="white", weight="bold")

    ax.set_yticks(range(len(tiles)))
    ax.set_yticklabels([f"{t.fwd_id} / {t.rev_id}" for t in tiles], fontsize=7)
    ax.invert_yaxis()
    ax.set_ylim(len(tiles) - 0.4, -0.6)
    ax.set_xlim(0, len(ref))
    ax.set_xlabel("CDS position (nt)")
    ax.set_ylabel("tile (fwd / rev primer)")
    return fig


# Risk tiers from checks/overhangs.py, palest to loudest.
_RISK_COLOURS = {"ok": "#EDF2F4", "watch": "#F7E3A1", "high": "#E8A33D", "collision": "#C1392B"}
# Two tiles never share a tube, so their cell is shown but not graded.
_CROSS_COLOUR = "#F7F8F9"


def overhang_figure(library, dpi: int = DEFAULT_DPI):
    """Return a matplotlib Figure of overhang homology, every fused overhang against every
    other.

    Each cell is how many of the four bases two overhangs share, taking the worse of the two
    orientations: as written, which is how the cut vector re-closes on itself, and against
    the reverse complement, which is how a fragment goes in backwards. The diagonal compares
    an overhang with its own reverse complement, so a hot diagonal cell is a palindrome.
    Cells are coloured by the risk tier the count falls in, and the ones inside a single
    reaction are boxed, since those are the pairs that meet in one tube.

    A design with no Golden Gate reaction has nothing to draw, so this raises. Call
    ``tile()`` first, or set a starting vector.
    """
    from .checks.overhangs import (
        MAX_SHARED,
        overhang_ends,
        risk_of,
        self_risk,
        self_shared,
        shared_bases,
    )
    from .regions import reverse_complement

    ends = overhang_ends(library)
    if not ends:
        raise ValueError(
            "No fused overhangs to plot. Tile the library with tile(), or set "
            "spec.starting_vector so there is a destination vector to cut."
        )

    n = len(ends)
    width = len(ends[0].seq)
    counts = [[0] * n for _ in range(n)]
    for i, a in enumerate(ends):
        for j, b in enumerate(ends):
            counts[i][j] = (
                self_shared(a) if i == j
                else max(shared_bases(a.seq, b.seq),
                         shared_bases(a.seq, reverse_complement(b.seq)))
            )

    side = max(4.2, 0.52 * n + 2.4)
    fig = _new_figure(figsize=(side, side * 0.92), dpi=dpi)
    ax = fig.add_subplot(111)
    fig.subplots_adjust(left=0.22, right=0.98, top=0.86, bottom=0.24)

    for i in range(n):
        for j in range(n):
            c = counts[i][j]
            # Only a cell that can actually misfire is graded: a tile's two ends against each
            # other, or the diagonal, an end against its own reverse complement (scored on its
            # own scale, see checks/overhangs.self_risk). Two different tiles are assembled in
            # separate tubes, so their cell carries the count and no verdict.
            together = ends[i].reaction == ends[j].reaction
            graded = together or i == j
            risk = (self_risk(c, width) if i == j else risk_of(c, width)) if graded else None
            ax.add_patch(_rect(
                j - 0.5, i - 0.5, 1, 1, linewidth=1.2, edgecolor="white",
                facecolor=_RISK_COLOURS[risk] if graded else _CROSS_COLOUR,
            ))
            ax.text(j, i, str(c), ha="center", va="center", fontsize=8,
                    color="white" if risk == "collision" else "#22303C",
                    alpha=1.0 if graded else 0.35,
                    weight="bold" if risk in ("high", "collision") else "normal")
            # Box the pairs that share a tube.
            if i != j and together:
                ax.add_patch(_rect(j - 0.5, i - 0.5, 1, 1, fill=False,
                                   edgecolor="#22303C", linewidth=1.6))

    labels = [f"{e.label}  {e.seq}" for e in ends]
    ax.set_xticks(range(n))
    ax.set_yticks(range(n))
    ax.set_xticklabels(labels, rotation=90, fontsize=7, family="monospace")
    ax.set_yticklabels(labels, fontsize=7, family="monospace")
    ax.set_xlim(-0.5, n - 0.5)
    ax.set_ylim(n - 0.5, -0.5)
    ax.set_aspect("equal")
    for s in ax.spines.values():
        s.set_visible(False)
    ax.tick_params(length=0)
    ax.set_title(
        f"Overhang homology, {library.spec.name}\nbases shared out of {width}, "
        f"worse orientation; boxed pairs share a reaction",
        loc="left", fontsize=9,
    )

    from matplotlib.patches import Patch

    tiers = [("ok", f"≤ {MAX_SHARED} shared"), ("watch", "above target"),
             ("high", f"{width - 1} shared, one mismatch"), ("collision", "identical")]
    handles = [Patch(facecolor=_RISK_COLOURS[k], label=v) for k, v in tiers]
    handles.append(Patch(facecolor=_CROSS_COLOUR, edgecolor="#DDE2E6",
                         label="separate reactions, not graded"))
    # On the figure rather than the axes, so the rotated tick labels above it do not push
    # it off the canvas as the overhang count grows.
    fig.legend(
        handles=handles,
        loc="lower center", bbox_to_anchor=(0.5, 0.005), ncol=3, frameon=False, fontsize=7.5,
    )
    return fig


def _rect(x, y, w, h, **kwargs):
    from matplotlib.patches import Rectangle

    return Rectangle((x, y), w, h, **kwargs)


# Codons grouped for the matrix figure: amino acid, then that residue's codons from most to
# least used in the host. A codon sitting low in its own band was a compromise.
def _codon_rows(species: str, symbols: set[str]) -> list[tuple[str, str]]:
    """``[(amino acid, codon), ...]`` for the residues in play, in row order top to bottom."""
    from .optimize.backbone import ranked_codons

    ranked = ranked_codons(species)
    rows: list[tuple[str, str]] = []
    for aa in sorted(symbols):
        for codon in ranked.get(aa, []):
            rows.append((aa, codon))
    return rows


def codon_matrix_figure(library, log: bool = True, reference_only: bool = False,
                        dpi: int = DEFAULT_DPI):
    """Return a matplotlib Figure of which codon sits at every position of the CDS.

    Codons run down the y axis grouped by the amino acid they encode, each group in a shaded
    band and ordered most- to least-used in the host, so a codon drawn low in its own band is
    one the design had to compromise on. Codon position along the CDS runs across the x axis.
    A cell counts how many library members carry that codon at that position, so the frozen
    reference reads as a continuous path and each stamped substitution shows up as a mark off
    it. ``reference_only=True`` counts the reference alone, which reduces the plot to that
    path.

    Counts span the whole library, so one codon carried by every member sits orders of
    magnitude above a substitution carried by one. ``log=True`` (the default) scales the
    colour logarithmically so both are visible at once; pass ``log=False`` for a linear scale,
    which suits a ``SequenceSet`` whose members are independent and spread more evenly.
    """
    from matplotlib.colors import LinearSegmentedColormap, LogNorm
    from matplotlib.patches import Rectangle
    import matplotlib.pyplot as plt

    if library.reference is None and reference_only:
        raise ValueError(
            "There is no shared reference to plot. A sequence set has no frozen WT, so pass "
            "reference_only=False to count across its members instead."
        )
    df = library.df
    if "variable_dna" not in df.columns:
        raise ValueError("Library is not codon-optimized yet, call codon_optimize() first.")

    seqs = [library.reference] if reference_only else [
        s for s in df["variable_dna"] if isinstance(s, str)
    ]
    if not seqs:
        raise ValueError("No sequences to plot (every variable_dna is NA).")
    n_codons = min(len(s) for s in seqs) // 3

    spec = library.spec
    from dnachisel import translate

    # Rows cover every residue the library actually encodes, including a pinned stop, so an
    # amber scan gets its own band rather than dropping off the figure.
    present = {translate(s[i * 3:i * 3 + 3]) for s in seqs for i in range(n_codons)}
    rows = _codon_rows(spec.optimization.species, present)
    if not rows:
        raise ValueError("No codons to plot for the residues in this library.")
    row_of = {codon: i for i, (_aa, codon) in enumerate(rows)}

    counts = [[0] * n_codons for _ in rows]
    for s in seqs:
        for j in range(n_codons):
            r = row_of.get(s[j * 3:j * 3 + 3])
            if r is not None:
                counts[r][j] += 1

    height = max(3.0, 0.16 * len(rows) + 1.4)
    width = max(7.0, min(22.0, 0.055 * n_codons + 3.2))
    fig = _new_figure(figsize=(width, height), dpi=dpi)
    ax = fig.subplots()

    # Alternating bands behind whole amino-acid groups, so the eye can find a residue's rows
    # without reading every tick label.
    start = 0
    for i, (aa, _codon) in enumerate(rows + [("", "")]):
        if i == len(rows) or aa != rows[start][0]:
            if (rows[start][0] and
                    sorted({a for a, _ in rows}).index(rows[start][0]) % 2 == 0):
                ax.add_patch(Rectangle((-0.5, start - 0.5), n_codons, i - start,
                                       facecolor="0.92", edgecolor="none", zorder=0))
            start = i

    biggest = max((v for row in counts for v in row), default=1)
    masked = [[(v if v else float("nan")) for v in row] for row in counts]
    norm = LogNorm(vmin=1, vmax=max(biggest, 2)) if log and biggest > 1 else None
    # A single-hue ramp, palest at one member and darkest at all of them, so the frozen
    # reference reads as the strong path and the one-off substitutions stay quiet. A diverging
    # or rainbow map inverts that and the eye lands on the substitution bands instead. The
    # ramp starts partway in, because Blues begins at white and a single-member cell drawn
    # white is a cell you cannot see at all.
    base = plt.get_cmap("Blues")
    cmap = LinearSegmentedColormap.from_list(
        "codonmap", [base(v) for v in [0.22 + 0.78 * i / 255 for i in range(256)]]
    )
    mesh = ax.imshow(masked, aspect="auto", interpolation="nearest", origin="upper",
                     cmap=cmap, norm=norm, zorder=2,
                     extent=(-0.5, n_codons - 0.5, len(rows) - 0.5, -0.5))

    ax.set_yticks(range(len(rows)))
    ax.set_yticklabels([f"{aa} {codon}" for aa, codon in rows], fontsize=6,
                       fontfamily="monospace")
    ax.set_xlabel("Codon position (CDS)")
    ax.set_ylabel("Codon, grouped by amino acid")
    what = "the reference" if reference_only else f"{len(seqs)} members"
    ax.set_title(f"Codon map, {spec.name}  ({what}, {spec.optimization.species})", loc="left")
    ax.set_xlim(-0.5, n_codons - 0.5)
    ax.set_ylim(len(rows) - 0.5, -0.5)
    bar = fig.colorbar(mesh, ax=ax, pad=0.01, fraction=0.025)
    bar.set_label("members carrying this codon here", fontsize=7)
    bar.ax.tick_params(labelsize=6)
    fig.tight_layout()
    return fig
