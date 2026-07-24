# Notebooks

Both sets run with the project environment active:

```bash
uv run --with jupyterlab jupyter lab
```

## Tutorials (`tutorials/`)

Worked examples that run end to end on the bundled data in `examples/`. Start here to
see a full workflow and to check that the package works.

- `01-standard-scan.ipynb`. Single-substitution scan of a CDS that fits on one oligo:
  generate variants, codon-optimize onto a frozen WT reference, QC, and export for
  synthesis and for uSort-M.
- `02-tiled-assembly.ipynb`. Tiled assembly for a CDS longer than one oligo (glucokinase):
  codon-aligned tiles with orthogonal primer pairs, then Golden Gate into per-tile
  destination vectors.
- `03-omega-sequence-set.ipynb`. A library of independent full-length sequences
  (`examples/example_designs.faa`) built into whole genes from oligopools via
  [OMEGA](https://github.com/RomeroLab/omega). Needs a separately installed OMEGA.

## Templates (`templates/`)

Lean, copy-and-edit skeletons of the same three workflows. Fill in the input cell and
run top to bottom.

- `standard-scan.ipynb`
- `tiled-assembly.ipynb`
- `omega-sequence-set.ipynb`

Copy a template into `personal/` for your own design. That folder is git-ignored, so your
working notebooks stay out of version control.
