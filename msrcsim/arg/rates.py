"""Bounded memoization of the existing built-in ARG rate functions.

No rate formula, summation order, event order, or RNG operation is changed.
Custom adapters/maps keep their original call behavior. Shared histories must
be fixed functions of population and evolutionary time, as the model requires.
"""
from .backgrounds.neutral import NeutralBackground
from .backgrounds.msrc_two_state import MSRCTwoBackground
from .backgrounds.reservoir import ReservoirBackground
from .recombination_map import PiecewiseConstantRecombinationMap


class RateCache:
    def __init__(self, background, rec_map, capacity=8192):
        self.background = background
        self.rec_map = rec_map
        self.capacity = capacity
        self.hazards = {}
        self.epoch()

    @classmethod
    def supported(cls, background, rec_map):
        return (type(background) in (NeutralBackground, MSRCTwoBackground, ReservoirBackground)
                and type(rec_map) is PiecewiseConstantRecombinationMap)

    def epoch(self):
        self.pairs = {}
        self.frequencies = {}
        self.boundaries = {}

    def recombination(self, lineage, time):
        bg = self.background
        if type(bg) is MSRCTwoBackground:
            # Legacy marker-storage semantics keep the original calculation.
            return bg.recombination_hazard(lineage.population, lineage.background,
                                            lineage.active, lineage.segments, time, self.rec_map)
        if type(bg) is ReservoirBackground and bg._structured(time):
            if lineage.population not in self.frequencies:
                self.frequencies[lineage.population] = bg.frequency(lineage.population, time)
            p = self.frequencies[lineage.population]
        else:
            p = None
        hull = (lineage.active[0][0], lineage.active[-1][1]) if lineage.active else None
        key = (lineage.population, lineage.background, p, hull)
        if key not in self.hazards:
            if len(self.hazards) >= self.capacity:
                self.hazards.clear()
            if type(bg) is NeutralBackground:
                value, _ = self.rec_map.crossover_hazard_and_region(lineage.active)
            else:
                value = bg.recombination_hazard(lineage.population, lineage.background,
                                                lineage.active, lineage.segments, time, self.rec_map)
            self.hazards[key] = value
        return self.hazards[key]

    def coalescence(self, a, b, effective_size, time):
        key = (a.population, a.background, b.background, effective_size)
        if key not in self.pairs:
            self.pairs[key] = self.background.coalescence_rate(
                a.population, a.background, b.background, effective_size, time)
        return self.pairs[key]

    def boundary(self, population, time):
        if population not in self.boundaries:
            self.boundaries[population] = self.background.next_boundary(population, time)
        return self.boundaries[population]
