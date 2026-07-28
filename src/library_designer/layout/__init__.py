"""Where a design lands on real DNA. ``tiled`` cuts a long CDS into oligo-sized tile
windows and builds the synthesis oligos, ``destination`` builds the one plasmid a
standard library clones into, and ``vector_io`` reads the starting backbone both of
them clone against.
"""
from .tiled import TileInfo, assemble_oligo, assign_tile, compute_tiles, tile_library

__all__ = ["TileInfo", "assemble_oligo", "assign_tile", "compute_tiles", "tile_library"]
