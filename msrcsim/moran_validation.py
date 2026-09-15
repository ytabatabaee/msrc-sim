from __future__ import annotations

import numpy as np
from scipy.linalg import expm

from .moran import moran_rates


def moran_generator(total_chromosomes: int, selection: float = 0.0) -> np.ndarray:
    """Return the exact small-population Moran CTMC generator on K=0..M."""
    total = int(total_chromosomes)
    if total <= 0:
        raise ValueError("total_chromosomes must be positive")
    q = np.zeros((total + 1, total + 1), dtype=float)
    for k in range(1, total):
        rates = moran_rates(k, total, selection)
        q[k, k + 1] = rates.q_plus
        q[k, k - 1] = rates.q_minus
        q[k, k] = -rates.total
    return q


def moran_distribution_after_t(
    total_chromosomes: int,
    initial_count: int,
    elapsed_time: float,
    selection: float = 0.0,
) -> np.ndarray:
    """Exact finite-time Moran count distribution for validation only."""
    total = int(total_chromosomes)
    initial = int(initial_count)
    if not 0 <= initial <= total:
        raise ValueError(f"initial_count must be between 0 and {total}")
    if elapsed_time < 0:
        raise ValueError("elapsed_time must be non-negative")
    dist0 = np.zeros(total + 1, dtype=float)
    dist0[initial] = 1.0
    return dist0 @ expm(moran_generator(total, selection) * float(elapsed_time))


def moran_event_count_sample(
    total_chromosomes: int,
    initial_count: int,
    elapsed_time: float,
    rng: np.random.Generator,
    selection: float = 0.0,
) -> int:
    """Sample terminal K from the production Gillespie dynamics on one branch."""
    total = int(total_chromosomes)
    k = int(initial_count)
    t = 0.0
    while 0 < k < total and t < elapsed_time:
        rates = moran_rates(k, total, selection)
        if rates.total <= 0.0:
            break
        t += float(rng.exponential(1.0 / rates.total))
        if t >= elapsed_time:
            break
        k += 1 if rng.random() < rates.q_plus / rates.total else -1
    return k
