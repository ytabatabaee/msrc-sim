from __future__ import annotations

from typing import Protocol

import numpy as np

from .rearrangement import Rearrangement
from .species_tree import SpeciesTree


class PopulationFrequencyHistory(Protocol):
    """Frequency-history interface consumed by the backward MSRC simulator."""

    def frequency_at(self, branch_id: str, age: float) -> float: ...

    def next_frequency_boundary(self, branch_id: str, age: float) -> float: ...

    def terminal_frequency(self, taxon: str) -> float: ...


def population_process_model(config: dict | None) -> str:
    section = (config or {}).get("population_process", {}) or {}
    model = str(section.get("model", "wright_fisher")).lower().replace("-", "_")
    if model not in {"wright_fisher", "moran"}:
        raise ValueError("population_process.model must be 'wright_fisher' or 'moran'")
    return model


def resolved_population_process(config: dict | None) -> dict:
    model = population_process_model(config)
    if model == "moran":
        return {
            "model": "moran",
            "time_scale": "generation_equivalent",
            "chromosome_copies": "M=2Ne",
            "attempted_replacement_rate_per_generation": "M/2",
            "selection": "native_continuous_time_moran",
            "split_initialization": "binomial_from_parent_terminal_frequency",
        }
    return {"model": "wright_fisher"}


def simulate_population_history(
    tree: SpeciesTree,
    rearrangement: Rearrangement,
    rng: np.random.Generator,
    config: dict | None = None,
):
    model = population_process_model(config)
    if model == "wright_fisher":
        from .wright_fisher import simulate_frequency_history

        return simulate_frequency_history(tree, rearrangement, rng)
    from .moran import simulate_frequency_history

    return simulate_frequency_history(tree, rearrangement, rng)
