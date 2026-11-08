from __future__ import annotations

from dataclasses import dataclass
from typing import Any

from .segments import Segment, normalize


@dataclass
class ARGLineage:
    lineage_id: int
    population: str
    segments: tuple[Segment, ...]
    birth_time: float
    background: Any = None

    def __post_init__(self):
        self.segments = normalize(self.segments, float("inf"))

    @property
    def material_length(self) -> float:
        return sum(s.right - s.left for s in self.segments)
