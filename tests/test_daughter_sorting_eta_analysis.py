from __future__ import annotations

from collections import Counter

import pytest

from msrcsim.statistics import wilson_interval
from scripts.run_daughter_sorting_eta import _summarize_eta, theoretical_targets
from scripts.revise_daughter_sorting_eta import recompute_summary, validate_patterns


def test_eta_summary_bookkeeping_and_wilson_intervals():
    row = _summarize_eta(
        eta=2.0,
        n_histories=10,
        pattern_counts=Counter({"1010": 3, "0101": 1, "0011": 2, "0000": 4}),
        persistence_n=5,
        strict_n=6,
        discordant_n=4,
        persistent_discordant_n=4,
        d_without_p_count=0,
        diagnostic={
            "mean_realized_daughter_difference": 0.4,
            "opposite_enrichment_p": 0.7,
            "n_polymorphic_split_events": 10,
        },
    )
    assert row["persistence_p"] == 0.5
    assert row["strict_2_2_p"] == 0.6
    assert row["discordant_2_2_p"] == 0.4
    assert row["p_discordant_given_persistent"] == 0.8
    assert row["p_discordant_given_2_2"] == pytest.approx(4 / 6)
    assert row["pattern_1010_p"] == 0.3
    assert row["d_without_p_count"] == 0
    assert (row["conditional_ci_low"], row["conditional_ci_high"]) == wilson_interval(4, 5)


def test_theoretical_targets_use_expected_reciprocal_direction():
    rows = theoretical_targets([0.0, 2.0])
    zero = next(row for row in rows if row["eta"] == 0.0 and row["parent_frequency_A1"] == 0.5)
    positive = next(row for row in rows if row["eta"] == 2.0 and row["parent_frequency_A1"] == 0.5)
    assert zero["target_frequency_positive"] == zero["target_frequency_negative"] == 0.5
    assert positive["target_frequency_positive"] > 0.5
    assert positive["target_frequency_negative"] < 0.5


def test_revised_pattern_validation_preserves_leading_zeroes():
    rows = [
        {"terminal_pattern": "0101", "persistent": "True", "strict_2_2": "True", "discordant_2_2": "True", "process": "wright_fisher", "split_mode": "descendant", "eta": "0.0"},
        {"terminal_pattern": "0011", "persistent": "True", "strict_2_2": "True", "discordant_2_2": "False", "process": "wright_fisher", "split_mode": "descendant", "eta": "0.0"},
    ]
    validate_patterns(rows)
    summary = [{
        "process": "wright_fisher",
        "split_mode": "descendant",
        "eta": "0.0",
        "persistence_p": "1.0",
        "strict_2_2_p": "1.0",
        "discordant_2_2_p": "0.5",
        "p_discordant_given_persistent": "0.5",
        "p_discordant_given_2_2": "0.5",
        "d_without_p_count": "0",
        **{f"pattern_{format(i, '04b')}_p": ("0.5" if format(i, "04b") in {"0101", "0011"} else "0.0") for i in range(16)},
    }]
    check = recompute_summary(rows, summary)
    assert check["max_abs_difference"] == 0.0

    with pytest.raises(ValueError):
        validate_patterns([{"terminal_pattern": "101"}])
