import numpy as np

from msrcsim.rearrangement import Rearrangement
from msrcsim.species_tree import SpeciesTree
from msrcsim.structured_coalescent import simulate_genealogy
from msrcsim.wright_fisher import FrequencyHistory, FrequencyRecord


def _constant_history(tree, M=20, origin_age=1.0, p=0.5):
    records = []
    for age in (origin_age, 0.0):
        records.append(
            FrequencyRecord(
                "x", "ROOT", None, int(age * 10), age, M * p, M * (1 - p), M,
                p, 1 - p, "segregating", False, age == origin_age,
                age == 0.0, 0.0, M // 2,
            )
        )
    by = {branch_id: [] for branch_id in tree.branches}
    by["ROOT"] = records
    return FrequencyHistory(records, by, [], origin_age, "ROOT")


def test_structured_origin_is_recorded_from_explicit_rearrangement_age():
    tree = SpeciesTree("(1:0,2:0)ROOT;", 10, 4.0)
    rearrangement = Rearrangement("x", "inversion", "ROOT", 2, 10)
    origin_age = tree.branches["ROOT"].older_age - rearrangement.origin_time_from_branch_start
    rng = np.random.default_rng(4)
    from msrcsim.moran import simulate_frequency_history

    history = simulate_frequency_history(tree, rearrangement, rng)
    assert history.structured_origin_age == origin_age
    assert history.structured_origin_branch == "ROOT"


def test_wright_fisher_history_carries_the_same_origin_boundary():
    tree = SpeciesTree("(1:0,2:0)ROOT;", 10, 4.0)
    rearrangement = Rearrangement("x", "inversion", "ROOT", 2, 10)
    origin_age = tree.branches["ROOT"].older_age - rearrangement.origin_time_from_branch_start
    from msrcsim.wright_fisher import simulate_frequency_history

    history = simulate_frequency_history(tree, rearrangement, np.random.default_rng(5))
    assert history.structured_origin_age == origin_age
    assert history.structured_origin_branch == "ROOT"


def test_opposite_backgrounds_coalesce_after_boundary_without_switching():
    tree = SpeciesTree("(1:0,2:0)ROOT;", 10, 1.0)
    history = _constant_history(tree, origin_age=1.0)
    rng = np.random.default_rng(8)
    result = simulate_genealogy(
        0, tree, history, {"1": 0, "2": 1}, 0.5, 1.0, rng, record_events=True
    )
    assert len(result.coalescences) == 1
    assert result.coalescences[0].time >= 1.0
    assert result.coalescences[0].arrangement == 0
    assert all(event.event_type != "switch" for event in result.events)
    coal_event = next(event for event in result.events if event.event_type == "coalescence")
    assert coal_event.frequency_A1 is None
    assert coal_event.event_rate == 1 / 20


def test_finite_structured_interval_has_uniform_ordinary_ancestor_fallback():
    M = 20
    ne = M // 2
    T = 1.0
    tree = SpeciesTree("(1:0,2:0,3:0,4:0)ROOT;", ne, T)
    history = _constant_history(tree, M=M, origin_age=T)
    sampled = {"1": 1, "2": 0, "3": 1, "4": 0}
    rng = np.random.default_rng(19)
    counts = np.zeros(3, dtype=int)
    reps = 12000
    for locus in range(reps):
        result = simulate_genealogy(locus, tree, history, sampled, 0.0, 1.0, rng, record_events=False)
        counts[result.topology_index] += 1

    lam = 2 / M
    a0 = 1.0 - np.exp(-2 * lam * T)
    expected = np.array([(1 - a0) / 3, (1 + 2 * a0) / 3, (1 - a0) / 3])
    observed = counts / reps
    assert np.all(np.abs(observed - expected) < 0.025)
    assert observed[1] < 0.55
