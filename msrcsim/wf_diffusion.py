"""Neutral Wright--Fisher diffusion histories.

The numerical transition is a bounded beta transition.  For a branch-local
coalescent-time increment ``du``, it has mean ``p`` and variance
``p*(1-p)*du``.  A boundary layer of width ``tolerance`` is snapped to an
absorbing boundary layer with a mean-preserving projection onto the nearest
boundary-layer endpoint.  This avoids the systematic boundary and fixation
errors caused by clipping Gaussian Euler steps.  The method is intended as a
practical, experimental approximation;
the recorded trajectory is piecewise constant for the genealogy engine.
"""

from __future__ import annotations

from bisect import bisect_right
from dataclasses import dataclass, field
from typing import Any, Mapping

import numpy as np

from .rearrangement import Rearrangement
from .species_tree import PopulationBranch, SpeciesTree
from .speciation_split import (
    SpeciationSplitPolicy,
    SplitInitializationEvent,
    daughter_target_frequency,
    mass_conserving_daughter_frequencies,
)
from .wright_fisher import FrequencyRecord


@dataclass
class DiffusionFrequencyHistory:
    records: list[FrequencyRecord]
    by_branch: dict[str, list[FrequencyRecord]]
    split_events: list[Any] = field(default_factory=list)
    structured_origin_age: float | None = None
    structured_origin_branch: str | None = None
    population_process: str = "wf_diffusion"
    numerical_method: str = "bounded beta transition with absorbing boundary layer"
    _sorted_by_branch: dict[str, list[FrequencyRecord]] = field(default_factory=dict, init=False, repr=False)
    _ages_by_branch: dict[str, list[float]] = field(default_factory=dict, init=False, repr=False)

    def _ascending(self, branch_id: str):
        if branch_id not in self._sorted_by_branch:
            rows = sorted(self.by_branch[branch_id], key=lambda r: r.absolute_age)
            self._sorted_by_branch[branch_id] = rows
            self._ages_by_branch[branch_id] = [r.absolute_age for r in rows]
        return self._sorted_by_branch[branch_id]

    def frequency_at(self, branch_id: str, age: float) -> float:
        rows = self._ascending(branch_id)
        if not rows:
            return 0.0
        ages = self._ages_by_branch[branch_id]
        return float(rows[max(0, bisect_right(ages, float(age) + 1e-12) - 1)].frequency_A1)

    def next_frequency_boundary(self, branch_id: str, age: float) -> float:
        self._ascending(branch_id)
        ages = self._ages_by_branch[branch_id]
        i = bisect_right(ages, float(age) + 1e-10)
        return ages[i] if i < len(ages) else float("inf")

    def terminal_frequency(self, taxon: str) -> float:
        return self.frequency_at(taxon, 0.0)


# Public alias matching the historical WF/Moran naming convention.
FrequencyHistory = DiffusionFrequencyHistory


def _validate_options(config: Mapping[str, Any] | None) -> tuple[float, float, str]:
    section = (config or {}).get("population_process", {}) or {}
    opts = section.get("diffusion", {}) or {}
    tolerance = float(opts.get("tolerance", 0.001))
    max_step_cu = float(opts.get("max_step_cu", 0.01))
    split_mode = str(opts.get("split_initialization", "deterministic_parent_frequency"))
    if not (0.0 < tolerance < 0.25):
        raise ValueError("population_process.diffusion.tolerance must be in (0, 0.25)")
    if not (0.0 < max_step_cu <= 1.0):
        raise ValueError("population_process.diffusion.max_step_cu must be in (0, 1]")
    if split_mode not in {"deterministic_parent_frequency", "finite_copy_sampling"}:
        raise ValueError("diffusion.split_initialization must be deterministic_parent_frequency or finite_copy_sampling")
    return tolerance, max_step_cu, split_mode


def _step(p: float, du: float, rng: np.random.Generator, tolerance: float) -> float:
    """One bounded, mean/variance matched neutral diffusion transition."""
    p = min(1.0, max(0.0, float(p)))
    # Project the unresolved boundary layer onto {0, tolerance}.  The
    # probability of retaining tolerance is p/tolerance, so this projection
    # has conditional mean p and cannot turn a rare variant into immediate
    # fixation.  The upper boundary is symmetric.  At tolerance itself the
    # beta transition is used, avoiding a zero-length artificial holding step.
    if p < tolerance:
        return tolerance if rng.random() < p / tolerance else 0.0
    if p > 1.0 - tolerance:
        return 1.0 - tolerance if rng.random() < (1.0 - p) / tolerance else 1.0
    # Beta(a,b) has mean p and variance p(1-p)/(a+b+1).
    concentration = max(1e-9, 1.0 / max(float(du), 1e-15) - 1.0)
    draw = float(rng.beta(p * concentration, (1.0 - p) * concentration))
    if draw < tolerance:
        return tolerance if rng.random() < draw / tolerance else 0.0
    if draw > 1.0 - tolerance:
        return 1.0 - tolerance if rng.random() < (1.0 - draw) / tolerance else 1.0
    return draw


def _split_frequency(tree: SpeciesTree, policy: SpeciationSplitPolicy, parent_id: str,
                    child_id: str, parent_p: float) -> tuple[float, SplitInitializationEvent]:
    rule = policy.rule_for(parent_id)
    children = next(n for n in tree.nodes() if n.name == parent_id).children
    child_names = [n.name for n in children]
    if len(child_names) == 2:
        offset = rule.eta if child_id == rule.positive_child else -rule.eta
        side = "positive" if child_id == rule.positive_child else "negative"
    else:
        offset, side = 0.0, "neutral"
    if policy.model == "mass_conserving" and len(child_names) == 2:
        plus, minus = mass_conserving_daughter_frequencies(parent_p, rule.eta)
        target = plus if side == "positive" else minus
    else:
        target = daughter_target_frequency(parent_p, offset)
    event = SplitInitializationEvent(
        parent_id, child_id, side, rule.positive_child,
        "lost" if parent_p == 0 else ("fixed" if parent_p == 1 else "segregating"),
        None, None, parent_p, rule.eta, offset,
        float(np.exp(offset)) if policy.model == "reciprocal_log_odds" else None,
        target, 2 * tree.branches[child_id].effective_population_size, None, policy.model,
    )
    return float(target), event


def _record(rearrangement: Rearrangement, branch: PopulationBranch, index: int,
            age: float, p: float, *, origin=False, start=False, end=False) -> FrequencyRecord:
    s = branch.selection_coefficient if branch.selection_coefficient != 0.0 else rearrangement.selection_coefficient
    status = "not_present" if p <= 0.0 and not origin else ("newly_originated" if origin else ("lost" if p <= 0 else ("fixed" if p >= 1 else "segregating")))
    total = 2 * branch.effective_population_size
    return FrequencyRecord(rearrangement.rearrangement_id, branch.branch_id, branch.parent_branch_id,
                           index, float(age), p * total, (1.0 - p) * total, total, p, 1.0 - p,
                           status, origin, start, end, float(s), branch.effective_population_size)


def simulate_frequency_history(tree: SpeciesTree, rearrangement: Rearrangement,
                              rng: np.random.Generator, config: Mapping[str, Any] | None = None,
                              split_policy: SpeciationSplitPolicy | None = None) -> DiffusionFrequencyHistory:
    tolerance, max_step_cu, split_mode = _validate_options(config)
    if rearrangement.selection_coefficient != 0.0 or any(b.selection_coefficient != 0.0 for b in tree.branches.values()):
        raise ValueError("wf_diffusion currently supports neutral rearrangements only (selection must be zero)")
    if rearrangement.origin_branch not in tree.branches:
        raise ValueError("Unknown origin branch")
    origin_branch = tree.branches[rearrangement.origin_branch]
    origin_age = origin_branch.older_age - rearrangement.origin_time_from_branch_start
    if not (origin_branch.younger_age <= origin_age <= origin_branch.older_age):
        raise ValueError("Origin time falls outside origin branch")

    records: list[FrequencyRecord] = []
    by = {bid: [] for bid in tree.branches}
    terminal: dict[str, float] = {}
    split_events: list[Any] = []
    order = sorted(tree.branches.values(), key=lambda b: b.older_age, reverse=True)
    for branch in order:
        p = None if branch.branch_id == rearrangement.origin_branch else terminal.get(branch.parent_branch_id)
        if p is not None and split_policy is not None:
            p, event = _split_frequency(tree, split_policy, branch.parent_branch_id, branch.branch_id, p)
            split_events.append(event)
        if p is not None and split_mode == "finite_copy_sampling":
            p = float(rng.binomial(2 * branch.effective_population_size, p) / (2 * branch.effective_population_size))
        age = float(branch.older_age)
        idx = 0
        rows: list[FrequencyRecord] = []
        is_origin_branch = branch.branch_id == rearrangement.origin_branch
        if is_origin_branch:
            p = None
        if is_origin_branch and origin_age >= age - 1e-12:
            p = min(1.0, rearrangement.initial_copy_count / (2 * branch.effective_population_size))
            rows.append(_record(rearrangement, branch, idx, age, p, origin=True, start=True))
        else:
            rows.append(_record(rearrangement, branch, idx, age, 0.0 if p is None else p, start=True))
        while age > branch.younger_age + 1e-12:
            if is_origin_branch and p is None and age <= origin_age + 1e-12:
                if origin_age < branch.older_age - 1e-12:
                    rows.append(_record(rearrangement, branch, idx, origin_age + 1e-10, 0.0))
                p = min(1.0, rearrangement.initial_copy_count / (2 * branch.effective_population_size))
                idx += 1
                rows.append(_record(rearrangement, branch, idx, origin_age, p, origin=True,
                                    end=origin_age <= branch.younger_age + 1e-12))
                age = origin_age
                continue
            next_age = max(branch.younger_age, age - max_step_cu * 2.0 * branch.effective_population_size)
            # The origin is a hard event boundary.  A pending origin may not
            # be crossed by an integration step.
            if is_origin_branch and p is None and next_age < origin_age:
                next_age = origin_age
            if p is not None and next_age < age:
                du = (age - next_age) / (2.0 * branch.effective_population_size)
                p = _step(p, du, rng, tolerance)
            age = next_age
            if is_origin_branch and p is None and abs(age - origin_age) <= 1e-12:
                if origin_age < branch.older_age - 1e-12:
                    rows.append(_record(rearrangement, branch, idx, origin_age + 1e-10, 0.0))
                p = min(1.0, rearrangement.initial_copy_count / (2 * branch.effective_population_size))
                idx += 1
                rows.append(_record(rearrangement, branch, idx, age, p, origin=True,
                                    end=age <= branch.younger_age + 1e-12))
            else:
                idx += 1
                rows.append(_record(rearrangement, branch, idx, age, 0.0 if p is None else p, end=age <= branch.younger_age + 1e-12))
        if rows:
            last = rows[-1]
            if abs(last.absolute_age - branch.younger_age) <= 1e-12 and not last.is_branch_end:
                rows[-1] = _record(rearrangement, branch, last.forward_generation, last.absolute_age,
                                   last.frequency_A1, origin=last.is_origin, end=True)
        records.extend(rows)
        by[branch.branch_id].extend(rows)
        terminal[branch.branch_id] = 0.0 if p is None else float(p)
    history = DiffusionFrequencyHistory(records, by, split_events, origin_age, rearrangement.origin_branch)
    history.diffusion_tolerance = tolerance
    history.diffusion_max_step_cu = max_step_cu
    history.split_initialization = split_mode
    return history
