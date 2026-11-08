from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class ARGEvent:
    index: int
    time: float
    event_type: str
    lineage_ids: tuple[int, ...]
    population_before: tuple[str, ...]
    population_after: tuple[str, ...]
    position: float | None = None
    rate: float | None = None
    material_length: float | None = None
    details: str = ""
