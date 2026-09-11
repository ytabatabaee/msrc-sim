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
    assert result.summary["taxon_order"] == ["1", "2", "3", "4"]


def test_joint_and_conditional_persistent_probabilities_are_consistent():
    result = exact_pattern_probabilities(_small_tree(), _small_rearrangement())
    persistent = result.summary["persistent_at_all_required_speciation_events"]
    for pattern in PATTERNS:
        assert (
            result.joint_persistent_pattern_probabilities[pattern]
            <= result.pattern_probabilities[pattern] + 1e-12
        )
    assert np.isclose(sum(result.joint_persistent_pattern_probabilities.values()), persistent)
    assert persistent > 0.0
    assert np.isclose(sum(result.conditional_persistent_pattern_probabilities.values()), 1.0)
    assert np.isclose(
        result.summary["P_2_2_given_persistent"]
        + result.summary["P_3_1_given_persistent"]
        + result.summary["P_4_0_given_persistent"],
        1.0,
    )
    assert np.isclose(
        result.summary["w1_given_persistent"]
        + result.summary["w2_given_persistent"]
        + result.summary["w3_given_persistent"],
        result.summary["P_2_2_given_persistent"],
    )
    assert np.isclose(
        result.summary["discordant_2_2_given_persistent"],
        result.summary["w2_given_persistent"] + result.summary["w3_given_persistent"],
    )


def test_zero_persistent_probability_has_null_conditionals():
    result = exact_pattern_probabilities(
        _small_tree(),
        Rearrangement("x", "inversion", "ROOT", 3, 8, 0.0),
    )
    assert result.summary["persistent_at_all_required_speciation_events"] == 0.0
    assert all(value is None for value in result.conditional_persistent_pattern_probabilities.values())
    assert result.summary["P_2_2_given_persistent"] is None
    assert result.summary["conditional_persistent_diagnostic"]


def test_invalid_initial_copy_count_raises_error():
    with np.testing.assert_raises_regex(ValueError, "initial_copy_count"):
        exact_pattern_probabilities(
            _small_tree(),
            Rearrangement("x", "inversion", "ROOT", 3, 9, 0.0),
        )


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
    for key in [
        "P_2_2_given_persistent",
        "P_3_1_given_persistent",
        "P_4_0_given_persistent",
        "w1_given_persistent",
        "w2_given_persistent",
        "w3_given_persistent",
        "discordant_2_2_given_persistent",
    ]:
        theory = result["theory"].summary[key]
        simulation = result["simulation"]["summary"][key]
        assert abs(theory - simulation) < 0.04
