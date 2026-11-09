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

    @staticmethod
    def eligible_crossover_region(active):
        """Return the shared ancestry-changing crossover region.

        Hudson's continuous-genome construction permits a cut anywhere in
        the span between the leftmost and rightmost *unresolved* ancestral
        material carried by a lineage.  Thus gaps between active segments are
        eligible, while stored material outside that span (including material
        whose local genealogy is already complete) is not.  The endpoints are
        half-open and have zero probability under the continuous map.
        """
        intervals = interval_union(active)
        if not intervals:
            return ()
        return ((intervals[0][0], intervals[-1][1]),)

    def crossover_hazard_and_region(self, active):
        region = self.eligible_crossover_region(active)
        return self.integrated_rate(region), region

    def breakpoint_hazard(self, segments) -> float:
        """Hazard for ancestry-changing cuts over ``segments`` as active material."""
        hazard, _ = self.crossover_hazard_and_region(segments)
        return hazard

    def sample_position(self, segments, rng: np.random.Generator) -> float:
        """Sample from exactly the map-weighted region used by the hazard."""
        total, eligible = self.crossover_hazard_and_region(segments)
        if not eligible:
            raise ValueError("cannot sample from empty ancestral material")
        pieces = []
        for left, right in eligible:
            for a, b, rate in self._pieces():
                lo, hi = max(left, a), min(right, b)
                weight = max(0.0, hi - lo) * rate
                if weight:
                    pieces.append((lo, hi, weight))
        if total <= 0:
            raise ValueError("cannot sample from a zero-hazard recombination map")
        draw = rng.random() * total
        for lo, hi, weight in pieces:
            if draw < weight:
                return float(lo + (draw / weight) * (hi - lo))
            draw -= weight
        return float(pieces[-1][1])
