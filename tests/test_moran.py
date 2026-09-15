from __future__ import annotations

import numpy as np

from msrcsim.moran import moran_rates, simulate_frequency_history
from msrcsim.moran_validation import moran_distribution_after_t, moran_event_count_sample
from msrcsim.population_process import simulate_population_history
from msrcsim.rearrangement import Rearrangement
from msrcsim.species_tree import SpeciesTree


def test_moran_transition_rates_neutral_and_selection():
    rates = moran_rates(5, 20, 0.0)
    assert rates.q_plus == rates.q_minus == 20 / 2 * 0.25 * 0.75
    selected = moran_rates(5, 20, 0.1)
    assert selected.q_plus > selected.q_minus


def test_moran_absorbing_boundaries():
    assert moran_rates(0, 20).total == 0.0
    assert moran_rates(20, 20).total == 0.0


def test_continuous_frequency_boundaries_known_times():
    tree = SpeciesTree("(A:10,B:10)ROOT;", 10, 0)
    rr = Rearrangement("x", "inversion", "ROOT", 1, 4, 0.0)
    root = tree.branches["ROOT"]
    from msrcsim.wright_fisher import FrequencyHistory, FrequencyRecord

    def rec(age, k):
        return FrequencyRecord("x", "ROOT", None, int(age), age, k, 20 - k, 20, k / 20, 1 - k / 20, "segregating", False, False, False, 0.0, 10)

    history = FrequencyHistory([], {"ROOT": [rec(10.0, 4), rec(8.7, 4), rec(5.2, 5), rec(1.1, 6), rec(0.0, 7)]})
    assert history.frequency_at("ROOT", 8.700001) == 4 / 20
    assert history.frequency_at("ROOT", 8.699999) == 5 / 20
    assert history.frequency_at("ROOT", 5.200001) == 5 / 20
    assert history.frequency_at("ROOT", 5.199999) == 6 / 20
    assert history.frequency_at("ROOT", 1.100001) == 6 / 20
    assert history.frequency_at("ROOT", 1.099999) == 7 / 20
    assert history.next_frequency_boundary("ROOT", 5.0) == 5.2
    assert root.branch_id == "ROOT"


def test_moran_history_reproducibility_and_one_copy_changes():
    tree = SpeciesTree("((A:4,B:4)AB:4,(C:4,D:4)CD:4)ROOT;", 12, 6)
    rr = Rearrangement("x", "inversion", "ROOT", 2, 8, 0.0)
    cfg = {"population_process": {"model": "moran"}}
    h1 = simulate_population_history(tree, rr, np.random.default_rng(5), cfg)
    h2 = simulate_population_history(tree, rr, np.random.default_rng(5), cfg)
    rows1 = [(r.branch_id, r.absolute_age, r.copy_count_A1) for r in h1.records]
    rows2 = [(r.branch_id, r.absolute_age, r.copy_count_A1) for r in h2.records]
    assert rows1 == rows2
    for branch_id, recs in h1.by_branch.items():
        ages = [r.absolute_age for r in recs]
        assert ages == sorted(ages, reverse=True)
        counts = [r.copy_count_A1 for r in recs if r.status != "not_present"]
        for a, b in zip(counts, counts[1:]):
            assert abs(a - b) in {0, 1}


def test_neutral_moran_fixation_probability_and_martingale_mc():
    rng = np.random.default_rng(11)
    total = 10
    initial = 3
    reps = 2500
    terminal = np.array([moran_event_count_sample(total, initial, 200.0, rng) for _ in range(reps)])
    assert abs(np.mean(terminal == total) - initial / total) < 0.04
    rng = np.random.default_rng(12)
    short = np.array([moran_event_count_sample(total, initial, 0.8, rng) / total for _ in range(reps)])
    assert abs(float(short.mean()) - initial / total) < 0.025


def test_moran_mc_matches_exact_ctmc_terminal_persistence_fix_loss():
    total = 8
    initial = 3
    elapsed = 1.5
    exact = moran_distribution_after_t(total, initial, elapsed)
    rng = np.random.default_rng(13)
    reps = 12000
    samples = np.array([moran_event_count_sample(total, initial, elapsed, rng) for _ in range(reps)])
    empirical = np.bincount(samples, minlength=total + 1) / reps
    assert np.max(np.abs(empirical - exact)) < 0.025
    for statistic in (
        lambda x: x[0],
        lambda x: x[-1],
        lambda x: x[1:-1].sum(),
    ):
        assert abs(float(statistic(empirical)) - float(statistic(exact))) < 0.025


def test_moran_arbitrary_tree_more_than_four_taxa_and_split_rule():
    tree = SpeciesTree("((A:3,B:3)AB:3,((C:3,D:3)CD:2,(E:3,F:3)EF:2)Y:1)ROOT;", 10, 4)
    rr = Rearrangement("x", "inversion", "ROOT", 1, 10, 0.0)
    history = simulate_frequency_history(tree, rr, np.random.default_rng(19))
    assert len(tree.taxa) == 6
    assert set(history.by_branch) == set(tree.branches)
    for child in ("AB", "Y"):
        start = max(history.by_branch[child], key=lambda r: r.absolute_age)
        assert 0 <= start.copy_count_A1 <= start.population_chromosomes
