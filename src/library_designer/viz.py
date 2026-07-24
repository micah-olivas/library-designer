"""Codon-usage QC visualization.

Plots relative codon adaptiveness (w = freq / best-synonymous-freq; 1.0 = optimal
codon) along the CDS: the WT backbone as a line, and each variant's stamped codon
as a point coloured by sublibrary. Points below the threshold flag codons the
motif-avoidance had to compromise on. Amber (`*`) stamps are expected to sit low.
a stop codon is pinned, not optimized.
"""
from __future__ import annotations

# Okabe-Ito colourblind-safe categorical palette (for sublibraries).
_PALETTE = ["#E69F00", "#56B4E9", "#009E73", "#CC79A7", "#0072B2", "#D55E00", "#F0E442", "#999999"]


def _new_figure(**kwargs):
    """A Figure created through pyplot and immediately closed.

    Importing pyplot registers the notebook inline PNG formatter, so the returned figure
    renders as an image without needing a ``%matplotlib inline`` in the notebook (a bare
    object-oriented Figure otherwise shows only as ``<Figure ...>`` text). Closing it
    removes it from pyplot's registry so it is not auto-shown a second time on top of the
    returned value. A closed figure still accepts axes/content and still supports
    ``savefig``.
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
                       compare: str | None = None, compare_label: str = "IDT"):
    """Return a matplotlib Figure of per-position codon usage along the CDS.

    ``metric="frequency"`` (default) plots absolute host codon-usage frequency. The
    WT backbone is a weaving landscape and each substitution's codon sits at its own
    frequency, so the plot is informative under any optimization method. Rare codons
    sit near the low-usage guide. ``metric="adaptiveness"`` plots relative
    adaptiveness w (freq / best-synonymous freq); useful for spotting motif-forced
    compromises, but under ``use_best_codon`` every optimized codon is 1.0.

    ``compare`` overlays an external coding sequence (e.g. IDT's WT optimization; raw
    DNA or FASTA) as a dashed line, so you can see where its codon choices diverge from
    the WT backbone.
    """
    from .optimize.backbone import codon_frequency, relative_adaptiveness

    if library.reference is None:
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

    fig = _new_figure(figsize=(11, 4))
    ax = fig.subplots()
    ax.plot(range(1, n_codons + 1), ref_y, color="0.45", lw=1.5, zorder=3, label="WT backbone")

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
    ax.legend(loc="lower right", ncol=4, fontsize=8, framealpha=0.9)
    fig.tight_layout()
    return fig


def tiling_figure(library):
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

    fig = _new_figure(figsize=(11, 0.45 * len(tiles) + 2.6))
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
        n = int((df["tile"] == t.index).sum()) if "tile" in df.columns else 0
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
