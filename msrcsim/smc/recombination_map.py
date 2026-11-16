from dataclasses import dataclass
import numpy as np


@dataclass(frozen=True)
class PiecewiseConstantMap:
    """Rates per base per generation on [0, length); boundaries are half-open."""
    length: float
    breakpoints: tuple[float, ...]
    rates: tuple[float, ...]

    def __post_init__(self):
        if not np.isfinite(self.length) or self.length <= 0 or len(self.rates) != len(self.breakpoints) + 1:
            raise ValueError("invalid recombination map dimensions")
        if (tuple(sorted(self.breakpoints)) != self.breakpoints or
                len(set(self.breakpoints)) != len(self.breakpoints) or
                any(not np.isfinite(x) or x <= 0 or x >= self.length for x in self.breakpoints)):
            raise ValueError("breakpoints must be strictly inside chromosome and sorted")
        if any(x < 0 or not np.isfinite(x) for x in self.rates):
            raise ValueError("recombination rates must be finite and nonnegative")

    @classmethod
    def constant(cls, length, rate):
        return cls(float(length), (), (float(rate),))

    def rate_at(self, x):
        return self.rates[int(np.searchsorted(self.breakpoints, x, side="right"))]

    def next_boundary(self, x):
        i = int(np.searchsorted(self.breakpoints, x, side="right"))
        return self.breakpoints[i] if i < len(self.breakpoints) else self.length

    def advance_hazard(self, x, hazard):
        """Invert integrated map rate; correctly crosses any number of intervals."""
        if not np.isfinite(x) or x < 0 or x > self.length:
            raise ValueError("starting position must lie within chromosome")
        if not np.isfinite(hazard) or hazard < 0:
            raise ValueError("integrated hazard must be finite and nonnegative")
        if hazard == 0:
            return float(x)
        while x < self.length:
            end = self.next_boundary(x)
            rate = self.rate_at(x)
            if rate > 0:
                available = (end - x) * rate
                # Generalized inverse convention: at an exact cumulative
                # hazard boundary, return the left boundary of any following
                # zero-rate plateau.
                if hazard <= available:
                    return x + hazard / rate
                hazard -= available
            x = end
        return self.length

    def to_json(self):
        return {"length": self.length, "breakpoints": list(self.breakpoints), "rates_per_bp_per_generation": list(self.rates)}
