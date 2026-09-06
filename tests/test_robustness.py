from copy import deepcopy

import numpy as np

from msrcsim.model_fitting import msc_probabilities
from msrcsim.robustness import binomial_confidence_interval, contribution_weights, infer_with_strategy
from msrcsim.robustness_cli import _assign_background, _linkage_diagnostics, _make_window_skeleton, _rows_for_fraction


def _constructed_rows():
    rows = []
    for block in range(5):
        rows.append({"block_id": block, "topology_index": 0, "topology": "12|34", "is_rearranged": False, "rearrangement_id": "", "msrc_probability": 0.0})
    for window in range(20):
        rows.append({"block_id": 5, "topology_index": 1, "topology": "13|24", "is_rearranged": True, "rearrangement_id": "inv", "msrc_probability": 1.0})
    return rows


def test_block_collapse_gives_each_block_total_weight_one():
    rows = _constructed_rows()
    weights = contribution_weights(rows, "block_collapse")
    totals = {}
    for row, weight in zip(rows, weights):
        totals[row["block_id"]] = totals.get(row["block_id"], 0.0) + weight
    assert set(round(v, 12) for v in totals.values()) == {1.0}


def test_filtering_and_weighting_never_rewrite_topology():
    rows = _constructed_rows()
    original = deepcopy(rows)
    infer_with_strategy(rows, "oracle_filter")
    infer_with_strategy(rows, "soft_weight")
    infer_with_strategy(rows, "block_collapse")
    assert rows == original


def test_constructed_example_naive_flips_but_corrections_recover_t1():
    rows = _constructed_rows()
    assert infer_with_strategy(rows, "all_windows").inferred_topology_index == 1
    assert infer_with_strategy(rows, "oracle_filter").inferred_topology_index == 0
    assert infer_with_strategy(rows, "block_collapse").inferred_topology_index == 0
    assert infer_with_strategy(rows, "soft_weight").inferred_topology_index == 0


def test_genealogy_block_and_rearrangement_interval_collapse_can_differ():
    rows = []
    for block in range(2):
        rows.append({"block_id": block, "topology_index": 0, "is_rearranged": False, "rearrangement_id": "", "msrc_probability": 0.0})
    for block in range(2, 8):
        rows.append({"block_id": block, "topology_index": 1, "is_rearranged": True, "rearrangement_id": "inv", "msrc_probability": 1.0})
    assert contribution_weights(rows, "genealogy_block_collapse") != contribution_weights(rows, "rearrangement_interval_collapse")
    assert infer_with_strategy(rows, "genealogy_block_collapse").inferred_topology_index == 1
    assert infer_with_strategy(rows, "rearrangement_interval_collapse").inferred_topology_index == 0


def test_rearrangement_interval_weights_sum_to_one_per_interval():
    rows = []
    for i in range(3):
        rows.append({"block_id": i, "topology_index": 1, "is_rearranged": True, "rearrangement_id": "inv_a"})
    for i in range(3, 8):
        rows.append({"block_id": i, "topology_index": 1, "is_rearranged": True, "rearrangement_id": "inv_b"})
    weights = contribution_weights(rows, "rearrangement_interval_collapse")
    totals = {}
    for row, weight in zip(rows, weights):
        totals[row["rearrangement_id"]] = totals.get(row["rearrangement_id"], 0.0) + weight
    assert {round(value, 12) for value in totals.values()} == {1.0}


def test_oracle_soft_weights_reproduce_oracle_filtering():
    rows = _constructed_rows()
    assert contribution_weights(rows, "soft_weight") == contribution_weights(rows, "oracle_filter")


def test_nonbinary_soft_probabilities_differ_from_oracle_filtering():
    rows = _constructed_rows()
    for row in rows:
        row["msrc_probability"] = 0.8 if row["is_rearranged"] else 0.2
        row["weight"] = 1.0 - row["msrc_probability"]
    assert contribution_weights(rows, "soft_weight") != contribution_weights(rows, "oracle_filter")


def test_ci_bounds_are_constrained_to_unit_interval():
    for successes, trials in [(0, 10), (10, 10), (47, 100)]:
        low, high = binomial_confidence_interval(successes, trials)
        assert 0.0 <= low <= high <= 1.0


def test_paired_mode_preserves_background_realization_where_intended():
    rng = np.random.default_rng(101)
    q_msc = msc_probabilities(0, 1.2)
    q_msrc = np.asarray([0.025, 0.95, 0.025])
    baseline = _assign_background(_make_window_skeleton("chr1", 1000.0, 20), q_msc=q_msc, block_windows=2, rng=rng)
    low = _rows_for_fraction(
        0.1, replicate_id=0, chrom="chr1", length=1000.0, windows=20,
        block_windows=2, msrc_block_windows=1, rng=rng, q_msc=q_msc, q_msrc=q_msrc,
        soft_probability_mode="oracle", soft_sensitivity=1.0, soft_specificity=1.0,
        soft_noise_sd=0.0, baseline_rows=baseline,
    )
    high = _rows_for_fraction(
        0.3, replicate_id=0, chrom="chr1", length=1000.0, windows=20,
        block_windows=2, msrc_block_windows=1, rng=rng, q_msc=q_msc, q_msrc=q_msrc,
        soft_probability_mode="oracle", soft_sensitivity=1.0, soft_specificity=1.0,
        soft_noise_sd=0.0, baseline_rows=baseline,
    )
    high_by_window = {row["window_id"]: row for row in high}
    for row in low:
        other = high_by_window[row["window_id"]]
        if not row["is_rearranged"] and not other["is_rearranged"]:
            assert row["topology_index"] == other["topology_index"]


def test_rows_for_fraction_reproducible_by_seed():
    q_msc = msc_probabilities(0, 1.2)
    q_msrc = np.asarray([0.025, 0.95, 0.025])
    kwargs = dict(
        fraction=0.25, replicate_id=0, chrom="chr1", length=1000.0, windows=20,
        block_windows=2, msrc_block_windows=1, q_msc=q_msc, q_msrc=q_msrc,
        soft_probability_mode="noisy", soft_sensitivity=0.8, soft_specificity=0.9,
        soft_noise_sd=0.02,
    )
    a = _rows_for_fraction(rng=np.random.default_rng(55), **kwargs)
    b = _rows_for_fraction(rng=np.random.default_rng(55), **kwargs)
    assert a == b


def test_replicate_block_coordinates_differ_across_seeds():
    q_msc = msc_probabilities(0, 1.2)
    q_msrc = np.asarray([0.025, 0.95, 0.025])
    kwargs = dict(
        fraction=0.5, replicate_id=0, chrom="chr1", length=1000.0, windows=50,
        block_windows=5, msrc_block_windows=1, q_msc=q_msc, q_msrc=q_msrc,
        soft_probability_mode="oracle", soft_sensitivity=1.0, soft_specificity=1.0,
        soft_noise_sd=0.0, kappa=0.5,
    )
    a = _rows_for_fraction(rng=np.random.default_rng(1), **kwargs)
    b = _rows_for_fraction(rng=np.random.default_rng(2), **kwargs)
    coords_a = {(row["block_id"], round(row["block_start"], 6), round(row["block_end"], 6)) for row in a}
    coords_b = {(row["block_id"], round(row["block_start"], 6), round(row["block_end"], 6)) for row in b}
    assert coords_a != coords_b


def test_benchmark_density_ratio_tracks_kappa_statistically():
    q_msc = msc_probabilities(0, 1.2)
    q_msrc = np.asarray([0.025, 0.95, 0.025])
    ratios = []
    for seed in range(160):
        rows = _rows_for_fraction(
            0.5, replicate_id=seed, chrom="chr1", length=10000.0, windows=200,
            block_windows=4, msrc_block_windows=1, rng=np.random.default_rng(seed),
            q_msc=q_msc, q_msrc=q_msrc, soft_probability_mode="oracle",
            soft_sensitivity=1.0, soft_specificity=1.0, soft_noise_sd=0.0,
            kappa=0.35,
        )
        diag = _linkage_diagnostics(rows, fraction=0.5, replicate_id=seed)
        ratio = float(diag[0]["observed_inside_outside_rate_ratio"])
        if np.isfinite(ratio):
            ratios.append(ratio)
    assert abs(float(np.mean(ratios)) - 0.35) < 0.15


def test_benchmark_kappa_near_zero_usually_one_rearranged_block():
    q_msc = msc_probabilities(0, 1.2)
    q_msrc = np.asarray([0.025, 0.95, 0.025])
    one_block = 0
    for seed in range(80):
        rows = _rows_for_fraction(
            0.5, replicate_id=seed, chrom="chr1", length=10000.0, windows=200,
            block_windows=5, msrc_block_windows=1, rng=np.random.default_rng(seed),
            q_msc=q_msc, q_msrc=q_msrc, soft_probability_mode="oracle",
            soft_sensitivity=1.0, soft_specificity=1.0, soft_noise_sd=0.0,
            kappa=0.001,
        )
        inside_blocks = {row["block_id"] for row in rows if row["is_rearranged"]}
        one_block += len(inside_blocks) == 1
    assert one_block / 80 > 0.9
