from __future__ import annotations

from dataclasses import dataclass
from math import sqrt
from typing import Any, Callable, Iterable, Mapping

import numpy as np

from .analytic import TOPOLOGY_NAMES


TOPOLOGY_LABEL_TO_INDEX = {name: i for i, name in enumerate(TOPOLOGY_NAMES)}


@dataclass(frozen=True)
class QuartetInferenceResult:
    strategy: str
    inferred_topology_index: int
    inferred_topology: str
    support: tuple[float, float, float]
    total_weight: float


def topology_index(value: str | int) -> int:
    if isinstance(value, str):
        if value in TOPOLOGY_LABEL_TO_INDEX:
            return TOPOLOGY_LABEL_TO_INDEX[value]
        return int(value)
    return int(value)


def maximum_quartet_support(rows: Iterable[Mapping[str, Any]], weights: Iterable[float] | None = None) -> QuartetInferenceResult:
    rows = list(rows)
    if weights is None:
        weights = [1.0] * len(rows)
    support = np.zeros(3, dtype=float)
    total = 0.0
    for row, weight in zip(rows, weights):
        top = topology_index(row.get("topology_index", row.get("topology")))
        weight = float(weight)
        if weight < 0.0:
            raise ValueError("quartet weights must be nonnegative")
        support[top] += weight
        total += weight
    inferred = int(np.argmax(support)) if total > 0.0 else -1
    return QuartetInferenceResult(
        "custom",
        inferred,
        TOPOLOGY_NAMES[inferred] if inferred >= 0 else "",
        tuple(float(x) for x in support),
        float(total),
    )


def contribution_weights(
    rows: Iterable[Mapping[str, Any]],
    strategy: str,
    soft_weight: Callable[[Mapping[str, Any]], float] | None = None,
) -> list[float]:
    rows = list(rows)
    if strategy == "all_windows":
        return [1.0] * len(rows)
    if strategy == "oracle_filter":
        return [0.0 if bool(row.get("is_rearranged", False)) else 1.0 for row in rows]
    if strategy in {"block_collapse", "genealogy_block_collapse"}:
        counts: dict[int, int] = {}
        for row in rows:
            block_id = int(row["block_id"])
            counts[block_id] = counts.get(block_id, 0) + 1
        return [1.0 / counts[int(row["block_id"])] for row in rows]
    if strategy == "rearrangement_interval_collapse":
        counts: dict[tuple[str, str], int] = {}
        for row in rows:
            if bool(row.get("is_rearranged", False)):
                key = ("rearrangement", str(row.get("rearrangement_id", "")))
            else:
                key = ("background_block", str(row["block_id"]))
            counts[key] = counts.get(key, 0) + 1
        weights = []
        for row in rows:
            if bool(row.get("is_rearranged", False)):
                key = ("rearrangement", str(row.get("rearrangement_id", "")))
            else:
                key = ("background_block", str(row["block_id"]))
            weights.append(1.0 / counts[key])
        return weights
    if strategy == "soft_weight":
        if soft_weight is None:
            def soft_weight(row: Mapping[str, Any]) -> float:
                if "weight" in row:
                    return float(row["weight"])
                if "msrc_probability" in row:
                    return 1.0 - float(row["msrc_probability"])
                if "p_msrc" in row:
                    return 1.0 - float(row["p_msrc"])
                return 0.0 if bool(row.get("is_rearranged", False)) else 1.0
        return [float(soft_weight(row)) for row in rows]
    raise ValueError("strategy must be all_windows, oracle_filter, genealogy_block_collapse, rearrangement_interval_collapse, or soft_weight")


def infer_with_strategy(
    rows: Iterable[Mapping[str, Any]],
    strategy: str,
    soft_weight: Callable[[Mapping[str, Any]], float] | None = None,
) -> QuartetInferenceResult:
    rows = list(rows)
    weights = contribution_weights(rows, strategy, soft_weight)
    result = maximum_quartet_support(rows, weights)
    return QuartetInferenceResult(strategy, result.inferred_topology_index, result.inferred_topology, result.support, result.total_weight)


def support_fraction(result: QuartetInferenceResult, topology: int) -> float:
    if result.total_weight <= 0.0:
        return float("nan")
    return float(result.support[int(topology)] / result.total_weight)


def binomial_confidence_interval(successes: int, trials: int, z: float = 1.959963984540054) -> tuple[float, float]:
    if trials <= 0:
        return (float("nan"), float("nan"))
    p = successes / trials
    denom = 1.0 + z * z / trials
    center = (p + z * z / (2.0 * trials)) / denom
    half = z * sqrt((p * (1.0 - p) + z * z / (4.0 * trials)) / trials) / denom
    return (max(0.0, center - half), min(1.0, center + half))


def theoretical_flip_threshold(t1_msc: float, t2_msc: float, t1_msrc: float, t2_msrc: float) -> float | None:
    delta_msc = float(t1_msc) - float(t2_msc)
    beta = float(t2_msrc) - float(t1_msrc)
    denom = delta_msc + beta
    if delta_msc <= 0.0 or beta <= 0.0 or denom <= 0.0:
        return None
    return delta_msc / denom


def dominant_quartet_threshold(tau: float, beta: float) -> float:
    """Analytic failure threshold for the dominant-quartet benchmark grid."""
    tau = float(tau)
    beta = float(beta)
    if tau < 0.0:
        raise ValueError("tau must be nonnegative")
    if beta <= 0.0:
        raise ValueError("beta must be positive")
    delta_msc = 1.0 - float(np.exp(-tau))
    return float(delta_msc / (delta_msc + beta))


def msrc_probabilities_from_beta(beta: float) -> np.ndarray:
    """MSRC marginal with T2 exceeding T1 by beta."""
    beta = float(beta)
    if not (0.0 <= beta <= 1.0):
        raise ValueError("beta must be between 0 and 1")
    minor = (1.0 - beta) / 3.0
    return np.asarray([minor, minor + beta, minor], dtype=float)


def interpolate_first_crossing(xs: Iterable[float], ys: Iterable[float], target: float = 0.0) -> float | None:
    """Linearly interpolate the first in-grid crossing of y - target."""
    points = sorted((float(x), float(y) - float(target)) for x, y in zip(xs, ys))
    if not points:
        return None
    for x, y in points:
        if y == 0.0:
            return x
    for (x0, y0), (x1, y1) in zip(points, points[1:]):
        if y0 == 0.0:
            return x0
        if (y0 < 0.0 < y1) or (y0 > 0.0 > y1):
            return float(x0 + (0.0 - y0) * (x1 - x0) / (y1 - y0))
        if y1 == 0.0:
            return x1
    return None
