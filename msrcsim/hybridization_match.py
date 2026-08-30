from __future__ import annotations

from dataclasses import dataclass, asdict
from math import sqrt
from typing import Any

import numpy as np
from scipy.optimize import minimize

from .hybridization import topology_index
from .model_fitting import msc_probabilities, network_probabilities


@dataclass(frozen=True)
class HybridizationMatch:
    major_topology: int
    introgressed_topology: int
    gamma: float
    t_major: float
    t_introgressed: float
    target_q: tuple[float, float, float]
    fitted_q: tuple[float, float, float]
    l1_error: float
    l2_error: float
    success: bool
    message: str

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


def mixture_q(major_topology: int, introgressed_topology: int, gamma: float, t_major: float, t_introgressed: float) -> np.ndarray:
    if major_topology == introgressed_topology:
        return (1.0 - gamma) * msc_probabilities(major_topology, t_major) + gamma * msc_probabilities(introgressed_topology, t_introgressed)
    return network_probabilities(major_topology, introgressed_topology, gamma, t_major, t_introgressed)


def fit_hybridization_to_q(
    target_q: tuple[float, float, float] | list[float] | np.ndarray,
    major_topology: str | int,
    introgressed_topology: str | int,
    *,
    max_branch_length: float = 20.0,
) -> HybridizationMatch:
    target = np.asarray(target_q, dtype=float)
    if target.shape != (3,) or np.any(target < 0.0) or not np.isclose(target.sum(), 1.0, atol=1e-6):
        raise ValueError("target_q must contain three nonnegative values summing to one")
    major = topology_index(major_topology)
    intro = topology_index(introgressed_topology)

    def objective(theta: np.ndarray) -> float:
        q = mixture_q(major, intro, float(theta[0]), float(theta[1]), float(theta[2]))
        diff = q - target
        return float(np.dot(diff, diff))

    starts = []
    for gamma in (0.05, 0.2, 0.5, 0.8, 0.95):
        for t1 in (0.02, 0.25, 0.75, 1.5, 4.0):
            for t2 in (0.02, 0.25, 0.75, 1.5, 4.0):
                starts.append((gamma, t1, t2))
    best = None
    for start in starts:
        result = minimize(
            objective,
            np.asarray(start, dtype=float),
            method="L-BFGS-B",
            bounds=((0.0, 1.0), (0.0, max_branch_length), (0.0, max_branch_length)),
        )
        if best is None or result.fun < best.fun:
            best = result
    assert best is not None
    gamma, t_major, t_intro = [float(x) for x in best.x]
    fitted = mixture_q(major, intro, gamma, t_major, t_intro)
    diff = fitted - target
    return HybridizationMatch(
        major_topology=major,
        introgressed_topology=intro,
        gamma=gamma,
        t_major=t_major,
        t_introgressed=t_intro,
        target_q=(float(target[0]), float(target[1]), float(target[2])),
        fitted_q=(float(fitted[0]), float(fitted[1]), float(fitted[2])),
        l1_error=float(np.sum(np.abs(diff))),
        l2_error=float(sqrt(float(np.dot(diff, diff)))),
        success=bool(best.success),
        message=str(best.message),
    )
