from __future__ import annotations

from dataclasses import dataclass
from typing import Any

from .segments import Segment, interval_union, normalize


@dataclass
class ARGLineage:
    lineage_id: int
    population: str
    segments: tuple[Segment, ...]
    birth_time: float
    background: Any = None
    # Active intervals still have two or more ancestral carriers. Resolved
    # segments remain in ``segments`` for output but are excluded from hazards.
    active: tuple[tuple[float, float], ...] | None = None

    def __post_init__(self):
        self.segments = normalize(self.segments, float("inf"))
        if self.active is None:
            self.active = interval_union((s.left, s.right) for s in self.segments)
        else:
            self.active = interval_union(self.active)
        stored = interval_union((s.left, s.right) for s in self.segments)
        for left, right in self.active:
            if not any(a <= left and right <= b for a, b in stored):
                raise ValueError("active ancestral material must be stored in segments")

    @property
    def material_length(self) -> float:
        return sum(s.right - s.left for s in self.segments)

    @property
    def active_length(self) -> float:
        return sum(right - left for left, right in self.active)

    @property
    def is_active(self) -> bool:
        return bool(self.active)

    def with_active(self, intervals):
        return ARGLineage(self.lineage_id, self.population, self.segments, self.birth_time,
                          self.background, interval_union(intervals))
