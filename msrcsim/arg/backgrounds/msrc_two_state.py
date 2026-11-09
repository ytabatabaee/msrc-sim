from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Mapping

from ..segments import interval_union


@dataclass(frozen=True)
class MSRCTwoBackground:
    """Two-state structured background model for the ARG.

    The model intentionally exposes the physical recombination measure rather
    than adding an independent background-switching clock.
    """

    history: Any
    background_anchor: float
    inversion_interval: tuple[float, float]
    cross_background_suppression: float = 1.0
    initial_backgrounds: Mapping[str, int] | None = None

    name: str = "msrc_two_state"

    def __post_init__(self):
        if not 0.0 <= self.cross_background_suppression <= 1.0:
            raise ValueError("cross_background_suppression must be in [0, 1]")
        left, right = self.inversion_interval
        if left < 0 or right <= left:
            raise ValueError("inversion_interval must be non-empty and ordered")

    @property
    def structured_origin_age(self):
        return getattr(self.history, "structured_origin_age", None)

    def initial_state(self, population: str, segment: Any, taxon: str | None = None) -> int:
        if taxon is not None and self.initial_backgrounds is not None and taxon in self.initial_backgrounds:
            state = int(self.initial_backgrounds[taxon])
            if state not in (0, 1):
                raise ValueError("background states must be 0 or 1")
            return state
        if taxon is not None:
            return int(self.history.terminal_frequency(taxon) >= 0.5)
        return 0

    def frequency(self, population: str, time: float) -> float:
        branch = population
        return min(1.0, max(0.0, float(self.history.frequency_at(branch, time))))

    def _structured(self, time: float) -> bool:
        origin = self.structured_origin_age
        return origin is None or time < float(origin) - 1e-10

    def coalescence_rate(self, population, state_a, state_b, effective_size, time):
        if not self._structured(time):
            return 1.0 / (2.0 * effective_size)
        if int(state_a) != int(state_b):
            return 0.0
        p = self.frequency(population, time)
        ph = p if int(state_a) == 1 else 1.0 - p
        return 1.0 / (2.0 * effective_size * max(ph, 1.0 / (2.0 * effective_size)))

    def next_boundary(self, population: str, time: float) -> float:
        values = [self.history.next_frequency_boundary(population, time)]
        origin = self.structured_origin_age
        if origin is not None and origin > time + 1e-10:
            values.append(float(origin))
        return min(values)

    def _kappa(self, position: float) -> float:
        left, right = self.inversion_interval
        return self.cross_background_suppression if left <= position < right else 1.0

    def _stored_anchor(self, segments) -> bool:
        return any(s.left <= self.background_anchor < s.right for s in segments)

    def _pieces(self, rec_map, left, right):
        cuts = [left, right]
        cuts.extend(x for x in self.inversion_interval if left < x < right)
        cuts = sorted(set(cuts))
        for lo, hi in zip(cuts, cuts[1:]):
            midpoint = (lo + hi) / 2.0
            for a, b, rate in rec_map._pieces():
                x, y = max(lo, a), min(hi, b)
                if y > x:
                    yield x, y, rate, self._kappa(midpoint)

    def _region(self, active, segments):
        neutral = self._map_region(active)
        if not neutral or not self._stored_anchor(segments):
            return neutral, neutral
        left, right = neutral[0][0], neutral[-1][1]
        marker_region = ((min(left, self.background_anchor), max(right, self.background_anchor)),)
        return neutral, marker_region

    @staticmethod
    def _map_region(active):
        intervals = interval_union(active)
        if not intervals:
            return ()
        return ((intervals[0][0], intervals[-1][1]),)

    def recombination_hazard(self, population, state, active, segments, time, rec_map):
        if not self._structured(time):
            return rec_map.breakpoint_hazard(active)
        _neutral, region = self._region(active, segments)
        p = self.frequency(population, time)
        same = p if int(state) == 1 else 1.0 - p
        cross = p if int(state) == 0 else 1.0 - p
        total = 0.0
        for left, right in region:
            for lo, hi, rate, kappa in self._pieces(rec_map, left, right):
                total += (hi - lo) * rate * (same + cross * kappa)
        return total

    def sample_recombination(self, population, state, active, segments, time, rec_map, rng):
        if not self._structured(time):
            position = rec_map.sample_position(active, rng)
            return {"position": position, "parent_background": int(state), "kind": "recombination",
                    "neutral_region": True, "background_changed": False}
        neutral, region = self._region(active, segments)
        if not region:
            raise ValueError("cannot recombine without stored active material")
        p = self.frequency(population, time)
        same = p if int(state) == 1 else 1.0 - p
        cross = p if int(state) == 0 else 1.0 - p
        pieces = []
        total = 0.0
        for left, right in region:
            for lo, hi, rate, kappa in self._pieces(rec_map, left, right):
                if rate <= 0:
                    continue
                weight = rate * (same + cross * kappa) * (hi - lo)
                pieces.append((lo, hi, rate, weight)); total += weight
        if total <= 0:
            raise ValueError("cannot sample from a zero-hazard MSRC recombination measure")
        draw = rng.random() * total
        for lo, hi, rate, weight in pieces:
            if draw < weight:
                position = lo + (draw / weight) * (hi - lo)
                kappa = self._kappa(position)
                same_weight = same
                cross_weight = cross * kappa
                parent = int(state) if rng.random() * (same_weight + cross_weight) < same_weight else 1 - int(state)
                neutral_event = any(a < position < b for a, b in neutral)
                return {"position": float(position), "parent_background": parent,
                        "kind": "recombination" if neutral_event else "background_transfer",
                        "neutral_region": neutral_event, "background_changed": parent != int(state)}
            draw -= weight
        raise RuntimeError("failed to sample MSRC recombination event")

    def recombination_states(self, population, state, left, right, time, position=None, parent_background=None):
        if parent_background is None or position is None or not self._structured(time):
            return state, state
        g = int(parent_background); h = int(state)
        if g == h:
            return h, h
        if position <= self.background_anchor:
            return g, h
        return h, g

    def transfer_state(self, state, parent_background):
        return int(parent_background)
