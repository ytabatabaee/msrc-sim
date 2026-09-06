import argparse

import numpy as np

from msrcsim.robustness import dominant_quartet_threshold, interpolate_first_crossing, msrc_probabilities_from_beta
from msrcsim.validation_grid import (
    CORRECTION_SUMMARY_FIELDS,
    CONSISTENCY_SUMMARY_FIELDS,
    KAPPA_REPLICATE_FIELDS,
    KAPPA_SUMMARY_FIELDS,
    THRESHOLD_REGRESSION_FIELDS,
    THRESHOLD_REPLICATE_FIELDS,
    THRESHOLD_SUMMARY_FIELDS,
    fit_margin_regression,
    run_consistency_benchmark,
    run_kappa_grid,
    run_threshold_grid,
    summarize_threshold_regression,
)


def _args(**overrides):
    values = dict(
        seed=11,
        chrom="chr1",
        chrom_length=10000.0,
        windows=120,
        block_windows=4,
        kappa=0.25,
        tau_values="0.1",
        beta_values="0.5",
        rearrangement_fractions="0.0,0.1,0.2,0.3,0.4,0.5,0.6",
        threshold_replicates=6,
        threshold_bootstrap_replicates=50,
        kappa_values="1.0,0.5,0.25",
        kappa_rearrangement_fractions="0.5",
        kappa_replicates=80,
        kappa_tau=0.5,
        kappa_beta=0.5,
        correction_cells="",
        correction_rearrangement_fractions="0.0,0.25,0.5",
        correction_replicates=4,
        soft_probability_mode="oracle",
        soft_sensitivity=1.0,
        soft_specificity=1.0,
        soft_noise_sd=0.0,
        consistency_tau=0.25,
        consistency_beta=0.5,
        consistency_n_blocks="10,25,50,100,250",
        consistency_replicates=160,
        consistency_soft_probability_mode="oracle",
        consistency_soft_sensitivity=1.0,
        consistency_soft_specificity=1.0,
        consistency_soft_noise_sd=0.0,
    )
    values.update(overrides)
    return argparse.Namespace(**values)


def test_dominant_quartet_threshold_formula_and_beta_probabilities():
    tau = 0.25
    beta = 0.5
    expected = (1.0 - np.exp(-tau)) / (1.0 - np.exp(-tau) + beta)
    assert dominant_quartet_threshold(tau, beta) == expected
    q = msrc_probabilities_from_beta(beta)
    assert np.isclose(q.sum(), 1.0)
    assert np.isclose(q[1] - q[0], beta)


def test_empirical_crossing_interpolation_synthetic_examples():
    assert np.isclose(interpolate_first_crossing([0.0, 0.1, 0.2], [0.3, 0.1, -0.1]), 0.15)
    assert np.isclose(interpolate_first_crossing([0.0, 0.1, 0.2], [0.2, 0.0, -0.2]), 0.1)


def test_no_extrapolated_threshold_is_reported():
    assert interpolate_first_crossing([0.0, 0.1, 0.2], [0.3, 0.2, 0.1]) is None
    assert interpolate_first_crossing([0.0, 0.1, 0.2], [-0.3, -0.2, -0.1]) is None


def test_regression_threshold_is_correct_on_synthetic_linear_data():
    fit = fit_margin_regression([0.0, 0.25, 0.5, 0.75], [0.4, 0.2, 0.0, -0.2])
    assert np.isclose(fit["intercept"], 0.4)
    assert np.isclose(fit["slope"], -0.8)
    assert np.isclose(fit["threshold"], 0.5)
    assert np.isclose(fit["r2"], 1.0)


def test_regression_threshold_not_reported_outside_sampled_grid():
    rows = []
    for replicate_id in range(3):
        for fraction, margin in [(0.0, 0.5), (0.25, 0.4), (0.5, 0.3)]:
            rows.append({
                "tau": 0.25,
                "beta": 0.5,
                "rearrangement_fraction": fraction,
                "support_margin_t1_minus_t2": margin,
                "replicate_id": replicate_id,
            })
    summary = summarize_threshold_regression(rows, bootstrap_replicates=10, seed=5)
    assert summary[0]["regression_threshold"] == ""


def test_bootstrap_regression_ci_is_reproducible_by_seed():
    rows, _ = run_threshold_grid(_args(threshold_replicates=10))
    a = summarize_threshold_regression(rows, bootstrap_replicates=40, seed=99)
    b = summarize_threshold_regression(rows, bootstrap_replicates=40, seed=99)
    assert a == b


def test_output_schemas_are_stable():
    assert THRESHOLD_REPLICATE_FIELDS == [
        "seed", "replicate_id", "tau", "beta", "rearrangement_fraction",
        "theory_threshold", "support_t1", "support_t2", "support_t3",
        "support_margin_t1_minus_t2", "inferred_topology_index", "recovered_true_t1",
        "num_windows", "num_genealogy_blocks", "kappa", "mode",
    ]
    assert THRESHOLD_SUMMARY_FIELDS == [
        "tau", "beta", "theory_threshold", "empirical_support_threshold",
        "empirical_recovery50_threshold", "abs_error", "n_replicates",
    ]
    assert "observed_ratio" in KAPPA_SUMMARY_FIELDS
    assert "seed" in KAPPA_REPLICATE_FIELDS
    assert THRESHOLD_REGRESSION_FIELDS[:4] == ["tau", "beta", "theory_threshold", "regression_threshold"]
    assert "ci_low" in CORRECTION_SUMMARY_FIELDS and "ci_high" in CORRECTION_SUMMARY_FIELDS
    assert "n_blocks" in CONSISTENCY_SUMMARY_FIELDS


def test_same_seed_reproduces_identical_threshold_grid_rows():
    rows_a, summary_a = run_threshold_grid(_args())
    rows_b, summary_b = run_threshold_grid(_args())
    assert rows_a == rows_b
    assert summary_a == summary_b


def test_kappa_one_observed_breakpoint_ratio_near_one():
    _, summary = run_kappa_grid(_args(kappa_values="1.0", kappa_replicates=160))
    ratio = float(summary[0]["observed_ratio"])
    assert abs(ratio - 1.0) < 0.25


def test_observed_ratio_decreases_monotonically_with_kappa_in_aggregate():
    _, summary = run_kappa_grid(_args(kappa_values="1.0,0.5,0.25,0.1", kappa_replicates=180))
    by_kappa = {float(row["kappa"]): float(row["observed_ratio"]) for row in summary}
    assert by_kappa[1.0] > by_kappa[0.5] > by_kappa[0.25] > by_kappa[0.1]


def test_pure_msc_recovery_increases_with_n_blocks_in_aggregate():
    rows = run_consistency_benchmark(_args(consistency_n_blocks="10,25,100,250"))
    pure = [row for row in rows if row["regime"] == "pure_msc" and row["strategy"] == "all_windows"]
    assert pure[-1]["recovery_probability"] >= pure[0]["recovery_probability"]


def test_below_threshold_recovery_approaches_one():
    rows = run_consistency_benchmark(_args(consistency_n_blocks="50,250,1000", consistency_replicates=220))
    below = [row for row in rows if row["regime"] == "below_threshold" and row["strategy"] == "all_windows"]
    assert below[-1]["recovery_probability"] > 0.9


def test_above_threshold_naive_recovery_decreases_toward_zero():
    rows = run_consistency_benchmark(_args(consistency_n_blocks="25,250,1000", consistency_replicates=220))
    above = [row for row in rows if row["regime"] == "above_threshold" and row["strategy"] == "all_windows"]
    assert above[-1]["recovery_probability"] < above[0]["recovery_probability"]
    assert above[-1]["recovery_probability"] < 0.2


def test_corrected_methods_recover_when_weighted_contamination_below_threshold():
    rows = run_consistency_benchmark(_args(consistency_n_blocks="1000", consistency_replicates=220))
    above = {row["strategy"]: row["recovery_probability"] for row in rows if row["regime"] == "above_threshold"}
    assert above["all_windows"] < 0.2
    assert above["oracle_filter"] > 0.95
    assert above["rearrangement_interval_collapse"] > 0.95
    assert above["soft_weight"] > 0.95
