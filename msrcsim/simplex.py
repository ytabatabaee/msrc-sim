from __future__ import annotations

from math import sqrt


def barycentric_to_cartesian(q1: float, q2: float, q3: float) -> tuple[float, float]:
    """Map quartet probabilities to an equilateral simplex.

    Vertices: T1=(0,0), T2=(1,0), T3=(1/2,sqrt(3)/2).
    """
    total = q1 + q2 + q3
    if abs(total - 1.0) > 1e-8:
        raise ValueError(f"Quartet probabilities must sum to one, got {total}")
    return q2 + 0.5 * q3, (sqrt(3.0) / 2.0) * q3
