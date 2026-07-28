from __future__ import annotations

from dataclasses import dataclass, asdict
from math import erfc, exp, inf, log, sqrt
from typing import Iterable

import numpy as np
from scipy.optimize import minimize

_EPS = 1e-15


def _safe_log(x: float) -> float:
    return log(max(float(x), _EPS))


def multinomial_log_likelihood(counts: Iterable[int], probabilities: Iterable[float]) -> float:
    return float(sum(int(n) * _safe_log(float(p)) for n, p in zip(counts, probabilities)))


def msc_probabilities(topology: int, branch_length: float) -> np.ndarray:
    if topology not in (0, 1, 2):
        raise ValueError("topology must be 0, 1, or 2")
    if branch_length < 0:
        raise ValueError("branch_length must be nonnegative")
    x = exp(-branch_length) / 3.0 if np.isfinite(branch_length) else 0.0
    q = np.full(3, x, dtype=float)
    q[topology] = 1.0 - 2.0 * x
    return q


@dataclass(frozen=True)
class MSCFit:
    topology: int
    branch_length: float
    q1: float
    q2: float
    q3: float
    log_likelihood: float
    aic: float
    distance: float


def fit_msc(counts: Iterable[int], topology: int) -> MSCFit:
    n = np.asarray(tuple(counts), dtype=float)
    if n.shape != (3,) or np.any(n < 0) or n.sum() <= 0:
        raise ValueError("counts must contain three nonnegative values with positive total")
    total = float(n.sum())
    minor = total - n[topology]
    x = min(1.0 / 3.0, max(0.0, minor / (2.0 * total)))
    branch = inf if x == 0.0 else max(0.0, -log(3.0 * x))
    q = np.full(3, x, dtype=float)
    q[topology] = 1.0 - 2.0 * x
    ll = multinomial_log_likelihood(n.astype(int), q)
    empirical = n / total
    return MSCFit(topology, branch, float(q[0]), float(q[1]), float(q[2]), ll, 2.0 - 2.0 * ll, float(np.linalg.norm(empirical - q)))


def fit_all_msc(counts: Iterable[int]) -> tuple[MSCFit, tuple[MSCFit, ...]]:
    fits = tuple(fit_msc(counts, i) for i in range(3))
    return max(fits, key=lambda f: f.log_likelihood), fits


def network_probabilities(parent_1: int, parent_2: int, gamma: float, t1: float, t2: float) -> np.ndarray:
    if parent_1 == parent_2:
        raise ValueError("parental topologies must differ")
    if not (0.0 <= gamma <= 1.0) or t1 < 0.0 or t2 < 0.0:
        raise ValueError("invalid network parameters")
    return (1.0 - gamma) * msc_probabilities(parent_1, t1) + gamma * msc_probabilities(parent_2, t2)


@dataclass(frozen=True)
class NetworkFit:
    parent_1: int
    parent_2: int
    gamma: float
    t1: float
    t2: float
    q1: float
    q2: float
    q3: float
    log_likelihood: float
    aic: float
    representation_error: float
    nondegenerate: bool
    well_interior: bool
    gamma_near_boundary: bool
    t1_near_zero: bool
    t2_near_zero: bool
    t1_near_upper_bound: bool
    t2_near_upper_bound: bool
    boundary_warning: bool


def fit_network_pair(counts: Iterable[int], parent_1: int, parent_2: int, max_branch: float = 20.0) -> NetworkFit:
    n = np.asarray(tuple(counts), dtype=float)
    if n.shape != (3,) or np.any(n < 0) or n.sum() <= 0:
        raise ValueError("counts must contain three nonnegative values with positive total")
    empirical = n / n.sum()

    def objective(theta: np.ndarray) -> float:
        gamma, t1, t2 = theta
        q = network_probabilities(parent_1, parent_2, float(gamma), float(t1), float(t2))
        return -multinomial_log_likelihood(n.astype(int), q)

    starts = [
        (0.25, 0.5, 0.5), (0.5, 0.5, 0.5), (0.75, 0.5, 0.5),
        (0.25, 2.0, 2.0), (0.5, 2.0, 2.0), (0.75, 2.0, 2.0),
        (0.5, 0.05, 5.0), (0.5, 5.0, 0.05),
    ]
    best = None
    for start in starts:
        result = minimize(objective, np.asarray(start), method="L-BFGS-B", bounds=((0.0, 1.0), (0.0, max_branch), (0.0, max_branch)))
        if best is None or result.fun < best.fun:
            best = result
    assert best is not None
    gamma, t1, t2 = map(float, best.x)
    q = network_probabilities(parent_1, parent_2, gamma, t1, t2)
    ll = multinomial_log_likelihood(n.astype(int), q)
    error = float(np.linalg.norm(empirical - q))
    formal_eps = 1e-4
    practical_gamma_eps = 0.02
    practical_t_eps = 0.02
    upper_margin = max(0.05, 0.01 * max_branch)
    nondegenerate = formal_eps < gamma < 1.0 - formal_eps and t1 > formal_eps and t2 > formal_eps
    gamma_near_boundary = gamma <= practical_gamma_eps or gamma >= 1.0 - practical_gamma_eps
    t1_near_zero = t1 <= practical_t_eps
    t2_near_zero = t2 <= practical_t_eps
    t1_near_upper_bound = t1 >= max_branch - upper_margin
    t2_near_upper_bound = t2 >= max_branch - upper_margin
    boundary_warning = gamma_near_boundary or t1_near_zero or t2_near_zero or t1_near_upper_bound or t2_near_upper_bound
    well_interior = nondegenerate and not boundary_warning
    return NetworkFit(
        parent_1, parent_2, gamma, t1, t2, float(q[0]), float(q[1]), float(q[2]),
        ll, 6.0 - 2.0 * ll, error, nondegenerate, well_interior,
        gamma_near_boundary, t1_near_zero, t2_near_zero,
        t1_near_upper_bound, t2_near_upper_bound, boundary_warning,
    )


def fit_all_networks(counts: Iterable[int]) -> tuple[NetworkFit, tuple[NetworkFit, ...]]:
    fits = tuple(fit_network_pair(counts, a, b) for a, b in ((0, 1), (0, 2), (1, 2)))
    return max(fits, key=lambda f: f.log_likelihood), fits


def off_arm_statistics(counts: Iterable[int], fitted_topology: int | None = None) -> dict[str, float | int]:
    n = np.asarray(tuple(counts), dtype=float)
    total = int(n.sum())
    if total <= 0:
        raise ValueError("positive counts are required")
    q = n / total
    best, _ = fit_all_msc(n.astype(int))
    topology = best.topology if fitted_topology is None else int(fitted_topology)
    minors = [i for i in range(3) if i != topology]
    difference = float(q[minors[0]] - q[minors[1]])
    variance = max(0.0, (q[minors[0]] + q[minors[1]] - difference * difference) / total)
    se = sqrt(variance)
    z = difference / se if se > 0 else (inf if difference > 0 else -inf if difference < 0 else 0.0)
    p = erfc(abs(z) / sqrt(2.0)) if np.isfinite(z) else 0.0
    arm_q = np.asarray((best.q1, best.q2, best.q3)) if topology == best.topology else np.asarray((fit_msc(n.astype(int), topology).q1, fit_msc(n.astype(int), topology).q2, fit_msc(n.astype(int), topology).q3))
    return {
        "nearest_msc_topology": topology,
        "nearest_msc_q1": float(arm_q[0]),
        "nearest_msc_q2": float(arm_q[1]),
        "nearest_msc_q3": float(arm_q[2]),
        "distance_to_nearest_msc_arm": float(np.linalg.norm(q - arm_q)),
        "off_arm_difference": difference,
        "off_arm_standard_error": se,
        "off_arm_z_score": z,
        "off_arm_p_value": p,
        "off_arm_ci95_low": difference - 1.959963984540054 * se,
        "off_arm_ci95_high": difference + 1.959963984540054 * se,
    }


def compare_models(counts: Iterable[int]) -> dict[str, object]:
    best_msc, _ = fit_all_msc(counts)
    best_net, _ = fit_all_networks(counts)
    result: dict[str, object] = {}
    result.update(off_arm_statistics(counts, best_msc.topology))
    result.update({f"best_msc_{k}": v for k, v in asdict(best_msc).items()})
    result.update({f"best_network_{k}": v for k, v in asdict(best_net).items()})
    delta_aic = best_net.aic - best_msc.aic
    result["delta_aic_network_vs_msc"] = delta_aic
    result["network_loglik_gain"] = best_net.log_likelihood - best_msc.log_likelihood

    network_representable = best_net.representation_error < 1e-6
    off_arm_supported = float(result["off_arm_p_value"]) < 0.05
    network_aic_preferred = delta_aic < 0.0
    network_strongly_preferred = delta_aic < -4.0
    result["network_representable"] = network_representable
    result["off_arm_supported"] = off_arm_supported
    result["network_aic_preferred"] = network_aic_preferred
    result["network_strongly_preferred"] = network_strongly_preferred

    if abs(float(result["off_arm_difference"])) < 1e-12:
        geometry = "single_tree_msc_arm"
    elif network_representable and best_net.well_interior:
        geometry = "network_interior_well_parameterized"
    elif network_representable and best_net.nondegenerate:
        geometry = "network_interior_boundary_warning"
    elif network_representable:
        geometry = "network_boundary"
    else:
        geometry = "poor_fit"

    if not off_arm_supported:
        evidence = "off_arm_not_supported"
    elif network_strongly_preferred:
        evidence = "network_strongly_preferred"
    elif network_aic_preferred:
        evidence = "network_weakly_preferred"
    else:
        evidence = "off_arm_supported_but_msc_aic_preferred"

    result["model_geometry_classification"] = geometry
    result["model_evidence_classification"] = evidence
    # Backward-compatible column; geometry and evidence columns should be preferred.
    result["model_classification"] = geometry
    return result
