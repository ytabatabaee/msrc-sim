from __future__ import annotations

from dataclasses import dataclass
from typing import Iterable

import numpy as np


@dataclass(frozen=True)
class AncestryTract:
    tract_id: int
    start_bp: float
    end_bp: float
    ancestry_state: int

    @property
    def length_bp(self) -> float:
        return self.end_bp - self.start_bp

    @property
    def ancestry_label(self) -> str:
        return ancestry_label(self.ancestry_state)


def ancestry_label(state: int) -> str:
    if int(state) == 0:
        return "major"
    if int(state) == 1:
        return "introgressed"
    raise ValueError("ancestry state must be 0 or 1")


def simulate_ancestry_tracts(
    length_bp: float,
    gamma: float,
    generations_since_pulse: float,
    recombination_rate_per_bp_per_generation: float,
    rng: np.random.Generator,
) -> list[AncestryTract]:
    """Simulate the pulse-hybridization ancestry mosaic along a chromosome.

    Hybridization sets the foreign ancestry fraction gamma. Recombination in
    the h generations after that pulse fragments the ancestry into tracts; it
    is not the origin of introgressed ancestry.
    """
    length = float(length_bp)
    if length <= 0.0:
        raise ValueError("chromosome length must be positive")
    gamma = float(gamma)
    h = float(generations_since_pulse)
    r = float(recombination_rate_per_bp_per_generation)
    if not (0.0 <= gamma <= 1.0):
        raise ValueError("gamma must be between 0 and 1")
    if h < 0.0:
        raise ValueError("generations_since_pulse must be nonnegative")
    if r < 0.0:
        raise ValueError("recombination rate must be nonnegative")

    if gamma == 0.0:
        return [AncestryTract(0, 0.0, length, 0)]
    if gamma == 1.0:
        return [AncestryTract(0, 0.0, length, 1)]

    initial_state = int(rng.random() < gamma)
    if h == 0.0 or r == 0.0:
        return [AncestryTract(0, 0.0, length, initial_state)]

    kappa = h * r
    alpha = kappa * gamma
    beta = kappa * (1.0 - gamma)
    tracts: list[AncestryTract] = []
    state = initial_state
    start = 0.0
    tract_id = 0
    while start < length:
        rate = alpha if state == 0 else beta
        tract_length = float(rng.exponential(1.0 / rate))
        end = min(length, start + tract_length)
        tracts.append(AncestryTract(tract_id, start, end, state))
        start = end
        state = 1 - state
        tract_id += 1
    return tracts


def introgressed_fraction_by_length(tracts: Iterable[AncestryTract]) -> float:
    total = 0.0
    introgressed = 0.0
    for tract in tracts:
        length = tract.length_bp
        total += length
        if tract.ancestry_state == 1:
            introgressed += length
    return introgressed / total if total > 0.0 else float("nan")
