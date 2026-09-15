import numpy as np

from msrcsim.robustness import dominant_quartet_threshold
from msrcsim.wf_moran_robustness import (
    block_normalized_epsilon,
    correction_epsilon_fields,
    correction_recovered,
    dense_window_epsilon,
    four_taxon_tree,
    matched_discordant_history,
    oracle_filter_effective_epsilon,
    simulate_forward_attempt,
    simulate_quartet_summary,
    soft_weight_effective_epsilon,
    threshold_from_beta,
)


def test_matched_0101_history_favors_t2_under_strong_suppression():
    tree = four_taxon_tree(ne=20, tau=0.6)
    history, sampled = matched_discordant_history(tree, "0101")
    result = simulate_quartet_summary(
        tree,
        history,
        sampled,
        num_loci=1200,
        seed=7,
        recombination_rate=0.02,
        effective_fraction=0.001,
    )
    assert result.q2 > result.q1
    assert result.q2 > result.q3
    assert result.unique_favored_dominant


def test_process_labels_share_matched_history_downstream_path():
    tree = four_taxon_tree(ne=20, tau=0.6)
    history, sampled = matched_discordant_history(tree, "1010")
    wf = simulate_quartet_summary(tree, history, sampled, 600, 10, 0.02, 0.001)
    moran = simulate_quartet_summary(tree, history, sampled, 600, 10, 0.02, 0.001)
    assert (wf.n1, wf.n2, wf.n3) == (moran.n1, moran.n2, moran.n3)
    assert wf.q2 > wf.q1


def test_threshold_calculation_uses_process_specific_beta():
    wf_beta = 0.30
    moran_beta = 0.18
    tau = 0.6
    assert threshold_from_beta(tau, wf_beta) == dominant_quartet_threshold(tau, wf_beta)
    assert threshold_from_beta(tau, wf_beta) != threshold_from_beta(tau, moran_beta)


def test_oracle_filter_effective_epsilon_is_zero_after_filtering():
    assert oracle_filter_effective_epsilon(0.55) == 0.0


def test_soft_weight_effective_epsilon_formula():
    assert np.isclose(soft_weight_effective_epsilon(0.55, w0=1.0, w1=0.5), 0.3793103448275862)
    assert np.isclose(soft_weight_effective_epsilon(0.55, w0=1.0, w1=0.1), 0.10891089108910891)


def test_block_normalized_epsilon_is_block_fraction():
    pi = 0.25
    assert block_normalized_epsilon(pi) == pi
    assert dense_window_epsilon(pi, mu0=1.0, mu1=3.0) == 0.5


def test_correction_reporting_keeps_epsilon_fields_distinct():
    naive = correction_epsilon_fields("all_windows", 0.55, 0.312)
    oracle = correction_epsilon_fields("oracle_filter", 0.55, 0.312)
    block = correction_epsilon_fields("genealogy_block_collapse", 0.55, 0.312, block_epsilon=0.25)
    soft = correction_epsilon_fields("soft_weight", 0.55, 0.312, soft_ratio=0.5)
    assert naive["raw_epsilon"] == 0.55
    assert naive["effective_epsilon"] == 0.55
    assert oracle["raw_epsilon"] == 0.55
    assert oracle["effective_epsilon"] == 0.0
    assert block["raw_epsilon"] == 0.55
    assert block["block_epsilon"] == 0.25
    assert block["effective_epsilon"] == 0.25
    assert np.isclose(soft["weighted_epsilon"], 0.3793103448275862)
    assert np.isclose(soft["effective_epsilon"], 0.3793103448275862)


def test_correction_utilities_accept_process_agnostic_rows():
    rows = []
    for block_id in range(6):
        rows.append({"population_process": "moran", "block_id": block_id, "topology_index": 0, "is_rearranged": False})
    for window in range(24):
        rows.append({"population_process": "moran", "block_id": 99, "topology_index": 1, "is_rearranged": True, "rearrangement_id": "inv"})
    assert not correction_recovered(rows, "all_windows")
    assert correction_recovered(rows, "genealogy_block_collapse")
    assert correction_recovered(rows, "oracle_filter")
    assert correction_recovered(rows, "soft_weight", 0.1)


def test_forward_attempt_runs_for_both_population_processes():
    for process in ("wright_fisher", "moran"):
        tree, history, sampled, pattern = simulate_forward_attempt(process, ne=12, tau=0.3, rng=np.random.default_rng(4))
        assert tree.taxa == ("1", "2", "3", "4")
        assert set(history.by_branch) == set(tree.branches)
        assert set(sampled) == set(tree.taxa)
        assert len(pattern) == 4
