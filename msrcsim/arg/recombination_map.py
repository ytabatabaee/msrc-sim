from __future__ import annotations

from dataclasses import dataclass
import numpy as np
from .segments import interval_union


@dataclass(frozen=True)
class PiecewiseConstantRecombinationMap:
    """Per-base, per-generation map over [0, length)."""

    length: float
    breakpoints: tuple[float, ...] = ()
    rates: tuple[float, ...] = (0.0,)

    def __post_init__(self):
        if self.length <= 0 or len(self.rates) != len(self.breakpoints) + 1:
            raise ValueError("invalid recombination map")
        if tuple(sorted(self.breakpoints)) != self.breakpoints or any(x <= 0 or x >= self.length for x in self.breakpoints):
            raise ValueError("breakpoints must be sorted and internal to the chromosome")
        if any(r < 0 for r in self.rates):
            raise ValueError("recombination rates must be nonnegative")

    @classmethod
    def constant(cls, length: float, rate: float):
        return cls(float(length), (), (float(rate),))

    def _pieces(self):
        points = (0.0,) + self.breakpoints + (self.length,)
        return zip(points[:-1], points[1:], self.rates)

    def integrated_rate(self, segments) -> float:
        total = 0.0
        for left, right in segments:
            for a, b, rate in self._pieces():
                total += max(0.0, min(right, b) - max(left, a)) * rate
        return total

    def breakpoint_hazard(self, segments) -> float:
        """Hazard for ancestry-changing cuts.

        A cut is possible anywhere between the leftmost and rightmost active
        ancestral bases of a lineage. Thus gaps between disjoint segments are
        included, while bases outside the material span are excluded.
        """
        intervals = interval_union(segments)
        if not intervals:
            return 0.0
        left, right = intervals[0][0], intervals[-1][1]
        return self.integrated_rate(((left, right),))

    def sample_position(self, segments, rng: np.random.Generator) -> float:
        intervals = interval_union(segments)
        if not intervals:
            raise ValueError("cannot sample from empty ancestral material")
        segments = ((intervals[0][0], intervals[-1][1]),)
        pieces = []
        total = 0.0
        for left, right in segments:
            for a, b, rate in self._pieces():
                lo, hi = max(left, a), min(right, b)
                weight = max(0.0, hi - lo) * rate
                if weight:
                    pieces.append((lo, hi, weight))
                    total += weight
        if total <= 0:
            raise ValueError("cannot sample from a zero-hazard recombination map")
        draw = rng.random() * total
        for lo, hi, weight in pieces:
            if draw <= weight:
                return float(lo + (draw / weight) * (hi - lo))
            draw -= weight
        return float(pieces[-1][1])
