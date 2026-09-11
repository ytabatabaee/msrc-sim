from __future__ import annotations

from math import comb

import numpy as np

from msrcsim.pattern_probabilities import (
    PATTERNS,
    equal_frequency_pattern_class_probabilities,
    exact_pattern_probabilities,
    monte_carlo_pattern_probabilities,
    pattern_class_summary,
    run_pattern_probability_analysis,
    wf_distribution_after_t,
    wf_transition_matrix,
)
from msrcsim.rearrangement import Rearrangement
from msrcsim.species_tree import SpeciesTree


def _small_tree(ne: int = 4) -> SpeciesTree:
    return SpeciesTree(
        "((1:3,2:3)A:2,(3:3,4:3)B:2)ROOT;",
        ne,
        5,
    )


def _small_rearrangement() -> Rearrangement:
    return Rearrangement("x", "inversion", "ROOT", 3, 4, 0.0)


def _small_config(seed: int = 11) -> dict:
    return {
        "mode": "mechanistic",
        "seed": seed,
        "num_loci": 10,
        "species_tree": {
            "newick": "((1:3,2:3)A:2,(3:3,4:3)B:2)ROOT;",
            "root_extension": 5,
            "default_effective_population_size": 4,
        },
        "rearrangement": {
            "id": "x",
            "type": "inversion",
            "origin_branch": "ROOT",
            "origin_time_from_branch_start": 3,
            "initial_copy_count": 4,
            "selection": {"model": "genic", "coefficient": 0.0},
        },
        "recombination": {
            "baseline_rate": 0.01,
            "effective_cross_arrangement_fraction": 0.05,
        },
        "sampling": {"samples_per_species": 1},
    }


def test_wf_transition_rows_sum_to_one_and_absorbing_states():
    matrix = wf_transition_matrix(4)
    assert np.allclose(matrix.sum(axis=1), 1.0)
    assert matrix[0, 0] == 1.0
    assert matrix[0, 1:].sum() == 0.0
    assert matrix[-1, -1] == 1.0
    assert matrix[-1, :-1].sum() == 0.0


def test_one_generation_transition_matches_binomial_probabilities():
    ne = 3
    initial_count = 2
    total = 2 * ne
    p = initial_count / total
    dist = wf_distribution_after_t(ne, initial_count, 1)
    expected = np.array(
        [comb(total, j) * p**j * (1.0 - p) ** (total - j) for j in range(total + 1)]
    )
    assert np.allclose(dist, expected)


def test_pattern_probabilities_sum_to_one_and_weights_match_2_2():
    result = exact_pattern_probabilities(_small_tree(), _small_rearrangement())
    assert set(result.pattern_probabilities) == set(PATTERNS)
    assert np.isclose(sum(result.pattern_probabilities.values()), 1.0)
    summary = pattern_class_summary(result.pattern_probabilities)
    assert np.isclose(summary["w1"] + summary["w2"] + summary["w3"], summary["P_2_2"])
    assert result.summary["required_speciation_nodes"] == ["ROOT", "A", "B"]


def test_equal_frequency_closed_form_pattern_classes():
    p = 0.35
    pattern_probs = {}
    for pattern in PATTERNS:
        ones = pattern.count("1")
        pattern_probs[pattern] = p**ones * (1.0 - p) ** (4 - ones)
    observed = pattern_class_summary(pattern_probs)
    expected = equal_frequency_pattern_class_probabilities(p)
    for key in ["P_2_2", "P_3_1", "P_4_0"]:
        assert np.isclose(observed[key], expected[key])


def test_symmetric_daughter_branches_have_symmetric_patterns():
    result = exact_pattern_probabilities(_small_tree(), _small_rearrangement())
    probs = result.pattern_probabilities
    assert np.isclose(probs["1000"], probs["0100"])
    assert np.isclose(probs["0010"], probs["0001"])
    assert np.isclose(probs["1100"], probs["0011"])
    assert np.isclose(probs["0101"], probs["1010"])
    assert np.isclose(probs["0110"], probs["1001"])


def test_same_seed_reproduces_monte_carlo_comparison():
    tree = _small_tree()
    rearrangement = _small_rearrangement()
    a = monte_carlo_pattern_probabilities(tree, rearrangement, replicates=200, seed=5)
    b = monte_carlo_pattern_probabilities(tree, rearrangement, replicates=200, seed=5)
    assert a == b


def test_theory_and_simulation_agree_within_monte_carlo_tolerance(tmp_path):
    result = run_pattern_probability_analysis(
        _small_config(),
        tmp_path,
        monte_carlo_replicates=20000,
        seed=19,
    )
    assert (tmp_path / "theoretical_pattern_probabilities.csv").exists()
    assert (tmp_path / "theoretical_summary.json").exists()
    assert (tmp_path / "theory_vs_simulation.csv").exists()
    assert (tmp_path / "theory_vs_simulation_pattern_probabilities.pdf").exists()
    assert result["max_abs_error"] < 0.02
