from __future__ import annotations

import numpy as np
import pytest

from msrcsim.population_process import population_process_model, simulate_population_history
from msrcsim.rearrangement import Rearrangement
from msrcsim.species_tree import SpeciesTree
from msrcsim.wf_diffusion import _step


def _history(ne=20, root_extension=20, seed=1, **diffusion):
    tree = SpeciesTree("(A:10,B:10)ROOT;", ne, root_extension)
    rr = Rearrangement("x", "inversion", "ROOT", 5, ne, 0.0)
    config = {"population_process": {"model": "wf_diffusion", "diffusion": diffusion}}
    return tree, simulate_population_history(tree, rr, np.random.default_rng(seed), config)


def test_diffusion_is_a_third_process_and_rejects_selection():
    assert population_process_model({"population_process": {"model": "wf_diffusion"}}) == "wf_diffusion"
    tree = SpeciesTree("(A:2,B:2)ROOT;", 10, 4)
    rr = Rearrangement("x", "inversion", "ROOT", 1, 1, 0.1)
    with pytest.raises(ValueError, match="neutral"):
        simulate_population_history(tree, rr, np.random.default_rng(2), {"population_process": {"model": "wf_diffusion"}})


def test_diffusion_boundaries_are_absorbing_and_records_are_continuous():
    tree = SpeciesTree("(A:5,B:5)ROOT;", 20, 10)
    for initial in (0, 40):
        rr = Rearrangement("x", "inversion", "ROOT", 1, initial, 0.0)
        h = simulate_population_history(tree, rr, np.random.default_rng(3), {"population_process": {"model": "wf_diffusion"}})
        assert all(r.frequency_A1 == initial / 40 for r in h.by_branch["ROOT"] if r.absolute_age <= 9)
    _, h = _history(seed=4)
    assert any(isinstance(r.copy_count_A1, float) for r in h.records)
    assert all(0.0 <= r.frequency_A1 <= 1.0 for r in h.records)


def test_diffusion_reproducibility_and_genealogy_interface():
    tree, h1 = _history(seed=8, max_step_cu=0.02)
    _, h2 = _history(seed=8, max_step_cu=0.02)
    assert [(r.absolute_age, r.frequency_A1) for r in h1.records] == [(r.absolute_age, r.frequency_A1) for r in h2.records]
    assert h1.next_frequency_boundary("ROOT", 2.0) > 2.0
    assert h1.terminal_frequency("A") == h1.frequency_at("A", 0.0)


def test_deterministic_diffusion_split_preserves_parent_frequency_at_split():
    tree = SpeciesTree("((A:4,B:4)AB:4,C:8)ROOT;", 20, 10)
    rr = Rearrangement("x", "inversion", "ROOT", 2, 10, 0.0)
    cfg = {"population_process": {"model": "wf_diffusion"}}
    h = simulate_population_history(tree, rr, np.random.default_rng(7), cfg)
    parent = h.frequency_at("ROOT", tree.branches["ROOT"].younger_age)
    assert h.frequency_at("AB", tree.branches["AB"].older_age) == pytest.approx(parent)
    assert h.frequency_at("C", tree.branches["C"].older_age) == pytest.approx(parent)


def test_finite_copy_split_is_explicit_binomial_sampling():
    tree = SpeciesTree("((A:4,B:4)AB:4,C:8)ROOT;", 20, 10)
    rr = Rearrangement("x", "inversion", "ROOT", 2, 10, 0.0)
    cfg = {"population_process": {"model": "wf_diffusion", "diffusion": {
        "split_initialization": "finite_copy_sampling",
    }}}
    h = simulate_population_history(tree, rr, np.random.default_rng(7), cfg)
    child_start = h.frequency_at("AB", tree.branches["AB"].older_age)
    assert child_start * (2 * tree.branches["AB"].effective_population_size) == pytest.approx(
        round(child_start * (2 * tree.branches["AB"].effective_population_size))
    )


def test_boundary_layer_does_not_create_immediate_fixation():
    rng = np.random.default_rng(101)
    draws = np.array([_step(0.0005, 0.01, rng, 0.001) for _ in range(20000)])
    assert draws.mean() == pytest.approx(0.0005, abs=0.00015)
    assert np.mean(draws == 1.0) < 0.01
    assert np.all((draws >= 0.0) & (draws <= 1.0))


def test_transition_mean_and_short_time_variance():
    rng = np.random.default_rng(102)
    p, du, n = 0.3, 0.001, 60000
    draws = np.array([_step(p, du, rng, 1e-5) for _ in range(n)])
    assert draws.mean() == pytest.approx(p, abs=0.0004)
    assert draws.var() == pytest.approx(p * (1 - p) * du, abs=0.00003)


def test_diffusion_martingale_and_heterozygosity_over_one_coalescent_unit():
    rng = np.random.default_rng(103)
    p0, du, steps, n = 0.3, 0.01, 100, 2500
    terminal = []
    heterozygosity = []
    for _ in range(n):
        p = p0
        for _ in range(steps):
            p = _step(p, du, rng, 1e-4)
        terminal.append(p)
        heterozygosity.append(2 * p * (1 - p))
    assert np.mean(terminal) == pytest.approx(p0, abs=0.025)
    assert np.mean(heterozygosity) == pytest.approx(2 * p0 * (1 - p0) * np.exp(-1), abs=0.012)


@pytest.mark.parametrize("origin_time", [0, 1, 9, 10])
def test_origin_is_recorded_at_exact_branch_age(origin_time):
    tree = SpeciesTree("(A:10,B:10)ROOT;", 20, 10)
    rr = Rearrangement("x", "inversion", "ROOT", origin_time, 4, 0.0)
    h = simulate_population_history(
        tree, rr, np.random.default_rng(202),
        {"population_process": {"model": "wf_diffusion", "diffusion": {"max_step_cu": 0.05}}},
    )
    origin_age = tree.branches["ROOT"].older_age - origin_time
    origin_rows = [r for r in h.by_branch["ROOT"] if r.is_origin]
    assert len(origin_rows) == 1
    assert origin_rows[0].absolute_age == pytest.approx(origin_age)
    if origin_time > 0:
        assert h.frequency_at("ROOT", origin_age + 1e-9) == 0.0
    assert h.frequency_at("ROOT", origin_age) == pytest.approx(4 / 40)
