from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any


@dataclass
class SpatialModelRun:
    model_name: str
    chromosome_length_bp: float
    loci: list[dict[str, Any]]
    windows: list[dict[str, Any]]
    summary: dict[str, Any]
    marginal_q: tuple[float, float, float]
    feature_intervals: list[dict[str, Any]]
    run_dir: Path


def marginal_l1_distance(a: tuple[float, float, float], b: tuple[float, float, float]) -> float:
    return float(sum(abs(float(x) - float(y)) for x, y in zip(a, b)))
