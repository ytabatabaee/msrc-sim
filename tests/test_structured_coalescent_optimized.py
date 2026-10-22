import numpy as np

from msrcsim.species_tree import SpeciesTree
from msrcsim.structured_coalescent import (
    simulate_genealogy,
    simulate_genealogy_reference,
)
from msrcsim.wright_fisher import FrequencyHistory, FrequencyRecord


def _fixture():
    tree = SpeciesTree("(A:0,B:0,C:0,D:0)ROOT;", 10, 4.0)
    records = [FrequencyRecord("x", "ROOT", None, int(age), age, 10, 10, 20,
                               0.5, 0.5, "segregating", False, age == 4,
                               age == 0, 0.0, 10) for age in (4., 3., 2., 1., 0.)]
    history = FrequencyHistory(records, {key: records for key in tree.branches},
                               [], 4.0, "ROOT")
    return tree, history, {"A": 1, "B": 0, "C": 1, "D": 0}


def test_optimized_process_returns_complete_genealogy():
    tree, history, sampled = _fixture()
    result = simulate_genealogy(0, tree, history, sampled, 0.02, 1.0,
                                np.random.default_rng(12), True)
    assert len(result.coalescences) == 3
    assert result.newick.endswith(";")
    assert set(result.newick) >= {"A", "B", "C", "D"}
    assert all(event.event_type in {"switch", "coalescence", "species_boundary"}
               for event in result.events)


def test_optimized_and_reference_have_same_event_invariants():
    tree, history, sampled = _fixture()
    for function in (simulate_genealogy, simulate_genealogy_reference):
        result = function(0, tree, history, sampled, 0.02, 1.0,
                          np.random.default_rng(42), True)
        assert len(result.coalescences) == len(tree.taxa) - 1
        assert len(result.coalescence_times) == len(tree.taxa) - 1
        assert all(x >= 0 for x in result.coalescence_times)


def test_optimized_topology_distribution_matches_reference_within_mc_error():
    tree, history, sampled = _fixture()
    counts = []
    for function in (simulate_genealogy, simulate_genealogy_reference):
        observed = np.zeros(3, dtype=int)
        for seed in range(1500):
            result = function(seed, tree, history, sampled, 0.02, 1.0,
                              np.random.default_rng(seed), False)
            observed[result.topology_index] += 1
        counts.append(observed)
    difference = np.abs(counts[0] / 1500 - counts[1] / 1500)
    assert np.all(difference < 0.06)
