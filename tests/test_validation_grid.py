import argparse

import numpy as np

from msrcsim.robustness import dominant_quartet_threshold, interpolate_first_crossing, msrc_probabilities_from_beta
from msrcsim.validation_grid import (
    CORRECTION_SUMMARY_FIELDS,
    KAPPA_REPLICATE_FIELDS,
    KAPPA_SUMMARY_FIELDS,
    THRESHOLD_REPLICATE_FIELDS,
    THRESHOLD_SUMMARY_FIELDS,
    run_kappa_grid,
    run_threshold_grid,
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
    assert "ci_low" in CORRECTION_SUMMARY_FIELDS and "ci_high" in CORRECTION_SUMMARY_FIELDS


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
