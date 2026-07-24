"""Generator seam. A generator turns a spec into a Library of variants; new
library types (combinatorial, saturation/NNK, barcodes, ...) implement this
without touching the rest of the pipeline.
"""
from __future__ import annotations

from typing import Protocol, runtime_checkable

from ..library import Library


@runtime_checkable
class Generator(Protocol):
    def generate(self) -> Library: ...
