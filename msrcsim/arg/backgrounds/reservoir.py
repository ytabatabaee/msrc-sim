"""Version 1 background reservoir for the shared Hudson ARG engine.

Moved without changing rates or boundary conventions from spatial_reservoir.
The old import path remains an alias. See docs/msrc_spatial_recommended_model.md.
"""
from dataclasses import dataclass
from typing import ClassVar
import math
from .msrc_two_state import MSRCTwoBackground


@dataclass(frozen=True)
class ReservoirBackground(MSRCTwoBackground):
    name: str = "msrc_reservoir_v1"
    version: ClassVar[str] = "reservoir_v1"
    coalescence_boundary_policy: ClassVar[str] = "reject_absent_occupied_background"
    marker_exchange_policy: ClassVar[str] = "reservoir_independent_of_stored_marker"

    def __post_init__(self):
        super().__post_init__()
        origin = self.structured_origin_age
        if origin is None or not math.isfinite(origin) or origin <= 0:
            raise ValueError("reservoir model requires an explicit positive finite origin age")
        if not math.isfinite(self.background_anchor):
            raise ValueError("marker coordinate must be finite")

    def frequency(self, population, time):
        p = float(self.history.frequency_at(population, time))
        if not math.isfinite(p) or not 0 <= p <= 1:
            raise ValueError("background frequency must be finite and in [0,1]")
        return p

    def coalescence_rate(self, population, state_a, state_b, effective_size, time):
        if effective_size <= 0:
            raise ValueError("effective population size must be positive")
        if not self._structured(time):
            return 1 / (2 * effective_size)
        p = self.frequency(population, time)
        for state in (state_a, state_b):
            if state not in (0, 1):
                raise ValueError("background must be 0 or 1")
            ph = p if state else 1 - p
            if ph == 0:
                raise ValueError("lineage occupies an absent background; boundary bridge unspecified")
        if state_a != state_b:
            return 0.0
        ph = p if state_a else 1 - p
        return 1 / (2 * effective_size * ph)

    def recombination_hazard(self, population, state, active, segments, time, rec_map):
        # Validate even a lineage with no same-population partner.
        self.coalescence_rate(population, state, state,
                              1.0, time)
        return super().recombination_hazard(population, state, active, segments, time, rec_map)

    def _region(self, active, segments):
        # Crucial versioned change: never condition exchange on stored marker
        # ancestry. The reservoir remains available after marker separation.
        neutral = self._map_region(active)
        if not neutral:
            return (), ()
        left, right = neutral[0][0], neutral[-1][1]
        return neutral, ((min(left, self.background_anchor), max(right, self.background_anchor)),)
