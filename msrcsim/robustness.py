from __future__ import annotations

from dataclasses import dataclass
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
    if strategy == "block_collapse":
        counts: dict[int, int] = {}
        for row in rows:
            block_id = int(row["block_id"])
            counts[block_id] = counts.get(block_id, 0) + 1
        return [1.0 / counts[int(row["block_id"])] for row in rows]
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
    raise ValueError("strategy must be all_windows, oracle_filter, block_collapse, or soft_weight")


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
