from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .rearrangement import Rearrangement
from .species_tree import PopulationBranch, SpeciesTree
from .wright_fisher import FrequencyHistory, FrequencyRecord


TIME_CONVENTION = (
    "Continuous-time Moran histories use M=2Ne chromosome copies and an "
    "attempted replacement-event rate M/2 per generation-equivalent time. "
    "Under neutrality this gives Var[dp]/dt=p(1-p)/M=p(1-p)/(2Ne)."
)


@dataclass(frozen=True)
class MoranRates:
    q_plus: float
    q_minus: float

    @property
    def total(self) -> float:
        return self.q_plus + self.q_minus


def moran_rates(k: int, total: int, selection: float = 0.0) -> MoranRates:
    """State-changing Moran CTMC rates under the generation-equivalent scale.

    With ``M=total=2Ne`` and ``p=K/M``, neutral rates are
    ``q_plus=q_minus=(M/2) p(1-p)``. With native genic Moran selection,
    A1 reproductive fitness is ``1+s`` and A0 fitness is 1, giving
    ``p_sel = p(1+s)/(1+s*p)`` and rates
    ``q_plus=(M/2) p_sel (1-p)`` and
    ``q_minus=(M/2) (1-p_sel) p``.
    """
    if 1.0 + float(selection) <= 0.0:
        raise ValueError("Moran selection requires 1 + s_M > 0")
    if k <= 0 or k >= total:
        return MoranRates(0.0, 0.0)
    p = k / total
    if float(selection) == 0.0:
        p_sel = p
    else:
        p_sel = p * (1.0 + float(selection)) / (1.0 + float(selection) * p)
    scale = total / 2.0
    return MoranRates(scale * p_sel * (1.0 - p), scale * (1.0 - p_sel) * p)


def _status(k: int, total: int, ever: bool) -> str:
    if not ever:
        return "not_present"
    if k == 0:
        return "lost"
    if k == total:
        return "fixed"
    return "segregating"


def _record(
    rearrangement: Rearrangement,
    branch: PopulationBranch,
    event_index: int,
    age: float,
    k: int | None,
    ever: bool,
    *,
    is_origin: bool = False,
    is_branch_start: bool = False,
    is_branch_end: bool = False,
) -> FrequencyRecord:
    total = 2 * branch.effective_population_size
    kk = 0 if k is None else int(k)
    status = "newly_originated" if is_origin else _status(kk, total, ever)
    s = branch.selection_coefficient if branch.selection_coefficient != 0.0 else rearrangement.selection_coefficient
    return FrequencyRecord(
        rearrangement.rearrangement_id,
        branch.branch_id,
        branch.parent_branch_id,
        int(event_index),
        float(age),
        kk,
        total - kk,
        total,
        kk / total,
        1.0 - kk / total,
        status,
        bool(is_origin),
        bool(is_branch_start),
        bool(is_branch_end),
        float(s),
        branch.effective_population_size,
    )


def _simulate_branch_events(
    branch: PopulationBranch,
    rearrangement: Rearrangement,
    rng: np.random.Generator,
    start_count: int | None,
    origin_age: float | None,
) -> tuple[list[FrequencyRecord], int | None, int]:
    records: list[FrequencyRecord] = []
    total = 2 * branch.effective_population_size
    s = branch.selection_coefficient if branch.selection_coefficient != 0.0 else rearrangement.selection_coefficient
    if 1.0 + float(s) <= 0.0:
        raise ValueError("Moran selection requires 1 + s_M > 0")

    age = float(branch.older_age)
    younger = float(branch.younger_age)
    k = start_count
    ever = k is not None and k > 0
    event_index = 0
    records.append(_record(rearrangement, branch, event_index, age, k, ever, is_branch_start=True))

    if origin_age is not None and k is None:
        if origin_age < age - 1e-12:
            event_index += 1
            records.append(_record(rearrangement, branch, event_index, origin_age, None, False))
        age = float(origin_age)
        k = min(int(rearrangement.initial_copy_count), total)
        ever = True
        # The origin record marks the first interval younger than the origin.
        # Continuous-time queries just older than origin still see absence via
        # the preceding not-present record.
        event_index += 1
        records.append(_record(rearrangement, branch, event_index, np.nextafter(age, younger), k, ever, is_origin=True))

    while k is not None and younger < age and 0 < k < total:
        rates = moran_rates(k, total, s)
        if rates.total <= 0.0:
            break
        dt = float(rng.exponential(1.0 / rates.total))
        next_age = age - dt
        if next_age <= younger:
            break
        old_k = k
        if rng.random() < rates.q_plus / rates.total:
            k += 1
        else:
            k -= 1
        event_index += 1
        # A record at an event age describes the interval immediately older
        # than that event; the changed count is carried by the next younger
        # event record or the branch-end record.
        records.append(_record(rearrangement, branch, event_index, next_age, old_k, ever))
        age = next_age

    event_index += 1
    records.append(_record(rearrangement, branch, event_index, younger, k, ever, is_branch_end=True))
    return sorted(records, key=lambda r: r.absolute_age, reverse=True), k, max(0, event_index - 1)


def simulate_frequency_history(tree: SpeciesTree, rearrangement: Rearrangement, rng: np.random.Generator):
    """Simulate a continuous-time Moran structural-frequency history.

    The rearrangement is introduced once using the same origin-branch,
    origin-time, and initial-copy-count logic as Wright-Fisher. Daughter
    branches are initialized by binomial sampling from the parental terminal
    frequency, matching the existing species-split rule.
    """
    if rearrangement.origin_branch not in tree.branches:
        raise ValueError("Unknown origin branch")
    origin_branch = tree.branches[rearrangement.origin_branch]
    origin_age = origin_branch.older_age - rearrangement.origin_time_from_branch_start
    if not (origin_branch.younger_age <= origin_age <= origin_branch.older_age):
        raise ValueError("Origin time falls outside origin branch")

    records: list[FrequencyRecord] = []
    by: dict[str, list[FrequencyRecord]] = {branch_id: [] for branch_id in tree.branches}
    end_counts: dict[str, tuple[int, int] | None] = {}
    event_count = 0

    order = sorted(tree.branches.values(), key=lambda b: b.older_age, reverse=True)
    for branch in order:
        total = 2 * branch.effective_population_size
        if branch.branch_id == rearrangement.origin_branch:
            start_count = None
            branch_origin_age = origin_age
        else:
            parent = end_counts.get(branch.parent_branch_id)
            start_count = None if parent is None else int(rng.binomial(total, parent[0] / parent[1]))
            branch_origin_age = None
        branch_records, end_count, n_events = _simulate_branch_events(
            branch, rearrangement, rng, start_count, branch_origin_age
        )
        records.extend(branch_records)
        by[branch.branch_id].extend(branch_records)
        end_counts[branch.branch_id] = None if end_count is None else (int(end_count), total)
        event_count += n_events

    history = FrequencyHistory(records, by)
    history.population_process = "moran"
    history.moran_time_convention = TIME_CONVENTION
    history.num_population_events = event_count
    return history
