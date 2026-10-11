from __future__ import annotations

from typing import Protocol

import numpy as np

from .rearrangement import Rearrangement
from .species_tree import SpeciesTree
from .speciation_split import ORIENTATION_RULE, resolve_speciation_split


class PopulationFrequencyHistory(Protocol):
    """Frequency-history interface consumed by the backward MSRC simulator."""

    def frequency_at(self, branch_id: str, age: float) -> float: ...

    def next_frequency_boundary(self, branch_id: str, age: float) -> float: ...

    def terminal_frequency(self, taxon: str) -> float: ...


def population_process_model(config: dict | None) -> str:
    section = (config or {}).get("population_process", {}) or {}
    model = str(section.get("model", "wright_fisher")).lower().replace("-", "_")
    if model not in {"wright_fisher", "moran", "wf_diffusion"}:
        raise ValueError("population_process.model must be 'wright_fisher', 'moran', or 'wf_diffusion'")
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
    if model == "wf_diffusion":
        section = (config or {}).get("population_process", {}) or {}
        diffusion = section.get("diffusion", {}) or {}
        return {
            "model": "wf_diffusion",
            "time_scale": "generations; branch-local coalescent units t/(2Ne)",
            "chromosome_copies": "continuous diffusion (reference scale M=2Ne)",
            "selection": "unsupported; neutral only",
            "split_initialization": str(diffusion.get("split_initialization", "deterministic_parent_frequency")),
            "tolerance": float(diffusion.get("tolerance", 0.001)),
            "max_step_cu": float(diffusion.get("max_step_cu", 0.01)),
            "boundary_treatment": "absorbing boundary-layer snapping",
        }
    return {"model": "wright_fisher"}


def resolved_speciation_split(config: dict | None) -> dict | None:
    section = (config or {}).get("speciation_split")
    if section is None:
        return None
    return {
        "enabled": True,
        "model": str(section.get("model", "reciprocal_log_odds")),
        "eta": float(section.get("eta", 0.0)),
        "orientation_rule": ORIENTATION_RULE,
        "by_split": dict(section.get("by_split", {}) or {}),
    }


def simulate_population_history(
    tree: SpeciesTree,
    rearrangement: Rearrangement,
    rng: np.random.Generator,
    config: dict | None = None,
):
    model = population_process_model(config)
    split_policy = resolve_speciation_split(config, tree)
    if model == "wright_fisher":
        from .wright_fisher import simulate_frequency_history

        return simulate_frequency_history(tree, rearrangement, rng, split_policy)
    if model == "wf_diffusion":
        from .wf_diffusion import simulate_frequency_history

        return simulate_frequency_history(tree, rearrangement, rng, config, split_policy)
    from .moran import simulate_frequency_history

    return simulate_frequency_history(tree, rearrangement, rng, split_policy)
