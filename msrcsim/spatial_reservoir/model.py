"""A spatial extension with a provable one-site MSRC generator.

The unobserved marker is a background reservoir, not a retained marker tree.
Crossovers assign the distal parent a draw from that reservoir even when no
stored segment contains the marker. See docs/msrc_spatial_recommended_model.md.
"""
from dataclasses import dataclass
import math
import numpy as np

from ..arg.backgrounds.msrc_two_state import MSRCTwoBackground
from ..arg.engine import ARGConfig, simulate_arg


@dataclass(frozen=True)
class ReservoirBackground(MSRCTwoBackground):
    name: str = "msrc_reservoir_v1"

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


def simulate_reservoir_arg(tree, rec_map, history, marker, suppressed_region,
                           kappa=1.0, seed=1, sampled_backgrounds=None):
    """Simulate one chromosome, holding one supplied history fixed throughout.

    Nonzero-rate history jumps to an absent occupied background raise instead
    of silently regularizing or rejecting and resampling the trajectory.
    """
    if not 0 <= marker <= rec_map.length:
        raise ValueError("marker must lie inside or at an endpoint of the chromosome")
    a, b = suppressed_region
    if not 0 <= a < b <= rec_map.length:
        raise ValueError("suppressed region must lie on the chromosome")
    if any(branch.effective_population_size <= 0 for branch in tree.branches.values()):
        raise ValueError("effective population sizes must be positive")
    points = (0.0,) + tuple(rec_map.breakpoints) + (rec_map.length,)
    if (not all(math.isfinite(x) for x in points)
            or any(b <= a for a, b in zip(points, points[1:]))
            or not all(math.isfinite(r) and r >= 0 for r in rec_map.rates)):
        raise ValueError("map requires finite ordered coordinates and nonnegative finite rates")
    rng = np.random.default_rng(seed)
    if sampled_backgrounds is None:
        sampled_backgrounds = {}
        for taxon in tree.taxa:
            p = float(history.terminal_frequency(taxon))
            if not np.isfinite(p) or not 0 <= p <= 1:
                raise ValueError("invalid terminal background frequency")
            sampled_backgrounds[taxon] = int(rng.random() < p)
    if set(sampled_backgrounds) != set(tree.taxa):
        raise ValueError("one background assignment is required for every taxon")
    if any(h not in (0, 1) for h in sampled_backgrounds.values()):
        raise ValueError("backgrounds must be 0 or 1")
    model = ReservoirBackground(history, marker, suppressed_region, kappa,
                                sampled_backgrounds)
    result = simulate_arg(tree, ARGConfig(rec_map.length, rec_map, model, True,
                                          history, sampled_backgrounds), rng)
    return result
